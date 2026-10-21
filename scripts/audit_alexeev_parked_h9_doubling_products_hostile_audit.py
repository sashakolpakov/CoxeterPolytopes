#!/usr/bin/env python3
"""Independent hostile audit of the parked H9 doubling-product theorem.

This script does not import or call the product certifier or its scouting
helper.  It proves the all-parameter facet formulae from exact reflection
normal forms and global wall stabilizers, independently rebuilds the six GMP
cones and Vinberg gates, and performs a bounded comparison with pinned
CoxIter snapshots.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import reduce
import gzip
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
import tarfile
from typing import Any, Iterable, Sequence
import zipfile

import networkx as nx
import sympy as sp

from scripts.certify_k3_s3_h10_redoubling_family import _parse_coxiter_text
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import DiagramEnumerator
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "alexeev_parked_rank10_h9_examples.json"
PRODUCTS = ROOT / "data" / "enumerations" / "alexeev_parked_h9_doubling_products.json"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_parked_h9_doubling_products_hostile_audit.json"
)
COXITER_ZIP = Path("/private/tmp/coxiter-current.zip")
DOUBLING_TAR = Path(
    "/private/tmp/coxiter-current/CoxIter-master/graphs/testing_index2.tar.gz"
)

SCHEMA = "coxeter.alexeev_parked_h9_doubling_products.hostile_audit.v1"
SOURCE_SHA256 = "e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430"
PRODUCTS_SHA256 = "10e6323895bacd2d18f897f33ea5f9d4bb1ee4eedad061d545aa5302b447c694"
PRODUCTS_PAYLOAD_SHA256 = "3aebc89f3e91e1dd2f8f95569bc0d5ca0d01b81a5b76b95820e837b63ae9c0ea"
COXITER_ZIP_SHA256 = "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
DOUBLING_TAR_SHA256 = "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"

# Zero-based indices in the immutable seed ordering.
PATH_FAMILIES: dict[str, dict[str, Any]] = {
    "A14_8_1_C02": {
        "pair": (14, 15),
        "cube": (6,),
        "formula": "4*N+13",
        "slope": 4,
        "intercept": 13,
    },
    "A14_10_1_C02": {
        "pair": (5, 14),
        "cube": (0, 2),
        "formula": "3*N+14",
        "slope": 3,
        "intercept": 14,
    },
}
CUBE_SEED = "A14_10_1_C01"
CUBE_WALLS = (2, 4, 10)
RECTANGLE_SEED = "A16_14_1_C01"
RECTANGLE_PAIRS = ((0, 13), (7, 14))


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(matrix[row, column]) for column in range(matrix.cols)] for row in range(matrix.rows)]


def _pair(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> sp.Rational:
    return sp.Rational((left.T * gram * right)[0])


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    return sp.simplify(
        sp.eye(gram.rows)
        - sp.Rational(2, _pair(gram, root, root)) * root * (root.T * gram)
    )


def _primitive(values: Iterable[sp.Rational | int]) -> tuple[int, ...]:
    entries = [sp.Rational(value) for value in values]
    denominator = int(sp.ilcm(*(entry.q for entry in entries)))
    integers = [int(entry * denominator) for entry in entries]
    divisor = reduce(math.gcd, (abs(value) for value in integers if value))
    return tuple(value // divisor for value in integers)


def _unoriented_key(vector: sp.Matrix) -> tuple[int, ...]:
    key = _primitive(vector)
    first = next(value for value in key if value)
    return key if first > 0 else tuple(-value for value in key)


def _token(gram: sp.Matrix, left: int, right: int) -> str:
    if left == right:
        return "1"
    ratio = sp.factor(
        gram[left, right] ** 2 / (gram[left, left] * gram[right, right])
    )
    known = {
        sp.Rational(0): "2",
        sp.Rational(1, 4): "3",
        sp.Rational(1, 2): "4",
        sp.Rational(3, 4): "6",
        sp.Rational(1): "parallel",
    }
    if ratio in known:
        return known[ratio]
    if ratio > 1:
        return f"dotted:{ratio}"
    raise ArithmeticError(f"unexpected normalized product {ratio}")


def _relations(gram: sp.Matrix) -> tuple[tuple[str, ...], ...]:
    return tuple(
        tuple(_token(gram, left, right) for right in range(gram.rows))
        for left in range(gram.rows)
    )


def _histogram(relations: Sequence[Sequence[str]]) -> dict[str, int]:
    return dict(
        sorted(
            Counter(
                relations[left][right]
                for left in range(len(relations))
                for right in range(left)
            ).items()
        )
    )


def _doubling_wall_gate(
    gram: sp.Matrix, roots: Sequence[sp.Matrix], index: int
) -> bool:
    facet_gram = sp.Matrix(
        [[_pair(gram, left, right) for right in roots] for left in roots]
    )
    relations = _relations(facet_gram)
    return all(
        token in {"1", "2", "4", "6", "parallel"}
        or token.startswith("dotted:")
        for token in relations[index]
    )


def _graph(relations: Sequence[Sequence[str]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(relations)))
    for left in range(len(relations)):
        for right in range(left):
            if relations[left][right] != "2":
                graph.add_edge(left, right, relation=relations[left][right])
    return graph


EDGE_MATCH = nx.algorithms.isomorphism.categorical_edge_match("relation", None)


def _isomorphic(
    left: Sequence[Sequence[str]], right: Sequence[Sequence[str]]
) -> bool:
    return len(left) == len(right) and nx.is_isomorphic(
        _graph(left), _graph(right), edge_match=EDGE_MATCH
    )


def _load_seeds() -> dict[str, dict[str, Any]]:
    if _sha(SOURCE) != SOURCE_SHA256:
        raise ArithmeticError("the immutable parked-H9 seed archive changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    result: dict[str, dict[str, Any]] = {}
    for case in source["cases"]:
        section = case["section"]
        gram = _matrix(section["saturated_fixed_lattice"]["gram"])
        interior = sp.Matrix(
            [
                sp.Rational(value)
                for value in section["strict_fixed_interior"][
                    "fixed_lattice_coordinates"
                ]
            ]
        )
        active = tuple(
            row
            for row in section["restricted_mirrors"]["records"]
            if row["hyperbolic_status"] == "active_facet"
        )
        raw = tuple(
            sp.Matrix(row["primitive_integral_same_direction_coordinates"])
            for row in active
        )
        roots = tuple(
            root if _pair(gram, interior, root) > 0 else -root for root in raw
        )
        if not (
            gram.shape == (10, 10)
            and exact_inertia(gram) == (1, 9, 0)
            and all(_pair(gram, interior, root) > 0 for root in roots)
        ):
            raise ArithmeticError(f"{case['id']}: seed orientation changed")
        result[case["id"]] = {
            "gram": gram,
            "roots": roots,
            "interior": interior,
        }
    return result


def _matrix_key(matrix: sp.Matrix) -> tuple[sp.Rational, ...]:
    return tuple(sp.Rational(value) for value in matrix)


def _cube_elements(
    gram: sp.Matrix, roots: Sequence[sp.Matrix], indices: Sequence[int]
) -> tuple[tuple[int, sp.Matrix], ...]:
    reflections = [_reflection(gram, roots[index]) for index in indices]
    if any(
        left * right != right * left
        for position, left in enumerate(reflections)
        for right in reflections[position + 1 :]
    ):
        raise ArithmeticError("a claimed orthogonal cube does not commute")
    rows = []
    for mask in range(1 << len(indices)):
        action = sp.eye(gram.rows)
        for bit, reflection in enumerate(reflections):
            if mask & (1 << bit):
                action *= reflection
        rows.append((mask, action))
    if len({_matrix_key(action) for _, action in rows}) != len(rows):
        raise ArithmeticError("a claimed C2 cube action has a kernel")
    return tuple(rows)


def _unipotent_power(identity: sp.Matrix, nilpotent: sp.Matrix, exponent: int) -> sp.Matrix:
    return sp.simplify(
        identity
        + exponent * nilpotent
        + sp.Rational(exponent * (exponent - 1), 2) * nilpotent**2
    )


def _path_factor(
    gram: sp.Matrix, roots: Sequence[sp.Matrix], pair: tuple[int, int]
) -> dict[str, Any]:
    a, b = (roots[index] for index in pair)
    sa, sb = (_reflection(gram, root) for root in (a, b))
    delta = a + b
    translation = sa * sb
    nilpotent = translation - sp.eye(gram.rows)
    if not (
        _pair(gram, a, a) == _pair(gram, b, b) == -2
        and _pair(gram, a, b) == 2
        and _pair(gram, delta, delta) == 0
        and reduce(math.gcd, (abs(int(value)) for value in delta if value)) == 1
        and sp.Matrix.hstack(a, delta).rank() == 2
        and nilpotent.rank() == 2
        and (nilpotent**2).rank() == 1
        and nilpotent**3 == sp.zeros(gram.rows)
        and translation * a == a + 2 * delta
        and translation.det() == 1
        and sa.det() == sb.det() == -1
    ):
        raise ArithmeticError("a claimed infinite-dihedral factor failed exact replay")
    # N^3=0 proves the displayed binomial expression for every integer z.
    for exponent in range(-3, 4):
        if _unipotent_power(sp.eye(gram.rows), nilpotent, exponent) != translation**exponent:
            raise ArithmeticError("the exact parabolic power law failed")
    return {
        "a": a,
        "b": b,
        "sa": sa,
        "sb": sb,
        "delta": delta,
        "translation": translation,
        "nilpotent": nilpotent,
        "public": {
            "walls_one_based": [pair[0] + 1, pair[1] + 1],
            "generator_gram": [[-2, 2], [2, -2]],
            "primitive_null_sum": [str(value) for value in delta],
            "translation_definition": "T=s_a*s_b",
            "translation_power_identity": "T^z=I+z*N+z*(z-1)/2*N^2 for every integer z",
            "nilpotent_ranks_N_N2_N3": [2, 1, 0],
            "translation_on_a": "T^z(a)=a+2*z*(a+b)",
            "all_translation_powers_distinct": True,
            "normal_form": "T^z or T^z*s_a, z in Z",
            "normal_form_injective": True,
            "factor": "D_infinity",
        },
    }


def _translation_solutions(
    gram: sp.Matrix,
    factor: dict[str, Any],
    base: sp.Matrix,
    target: sp.Matrix,
) -> str | tuple[int, ...]:
    """Solve T^z base = target over all z in Z, without a finite search."""

    a = factor["a"]
    delta = factor["delta"]
    delta_pair = _pair(gram, delta, base)
    a_pair = _pair(gram, a, base)
    candidates: list[sp.Rational] = []
    if delta_pair:
        candidates = [
            sp.factor(
                (_pair(gram, a, base) - _pair(gram, a, target))
                / (2 * delta_pair)
            )
        ]
    elif a_pair:
        coordinate = next(index for index, value in enumerate(delta) if value)
        candidates = [
            sp.factor(-(target[coordinate] - base[coordinate]) / (a_pair * delta[coordinate]))
        ]
    else:
        return "all" if base == target else ()
    answers = []
    for candidate in candidates:
        if candidate.q != 1:
            continue
        exponent = int(candidate)
        power = _unipotent_power(
            sp.eye(gram.rows), factor["nilpotent"], exponent
        )
        if power * base == target:
            answers.append(exponent)
    return tuple(sorted(set(answers)))


def _expected_path_stabilizers(
    orthogonal_a: bool, orthogonal_b: bool
) -> tuple[tuple[int, str | tuple[int, ...]], ...]:
    if orthogonal_a and orthogonal_b:
        return ((0, "all"), (1, "all"))
    if orthogonal_a:
        return ((0, (0,)), (1, (0,)))
    if orthogonal_b:
        return ((0, (0,)), (1, (-1,)))
    return ((0, (0,)),)


def _path_cube_global_orbits(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    pair: tuple[int, int],
    cube: tuple[int, ...],
    factor: dict[str, Any],
) -> dict[str, Any]:
    selected = set((*pair, *cube))
    nonselected = tuple(index for index in range(len(roots)) if index not in selected)
    cubes = _cube_elements(gram, roots, cube)
    identity = sp.eye(gram.rows)
    actual: set[tuple[int, int, int, int, int, str | tuple[int, ...]]] = set()
    for target_index in nonselected:
        for source_index in nonselected:
            for epsilon in (0, 1):
                for mask, cube_action in cubes:
                    base = (
                        (factor["sa"] if epsilon else identity)
                        * cube_action
                        * roots[source_index]
                    )
                    for sign in (1, -1):
                        solutions = _translation_solutions(
                            gram, factor, base, sign * roots[target_index]
                        )
                        if solutions == "all" or solutions:
                            actual.add(
                                (
                                    target_index,
                                    source_index,
                                    epsilon,
                                    mask,
                                    sign,
                                    solutions,
                                )
                            )
    expected: set[tuple[int, int, int, int, int, str | tuple[int, ...]]] = set()
    wall_rows = []
    for index in nonselected:
        orthogonal_a = _pair(gram, roots[index], factor["a"]) == 0
        orthogonal_b = _pair(gram, roots[index], factor["b"]) == 0
        orthogonal_cube_bits = tuple(
            bit
            for bit, wall in enumerate(cube)
            if _pair(gram, roots[index], roots[wall]) == 0
        )
        allowed_masks = tuple(
            mask
            for mask in range(1 << len(cube))
            if all(not (mask & (1 << bit)) for bit in range(len(cube)) if bit not in orthogonal_cube_bits)
        )
        for epsilon, solutions in _expected_path_stabilizers(
            orthogonal_a, orthogonal_b
        ):
            for mask in allowed_masks:
                expected.add((index, index, epsilon, mask, 1, solutions))
        nonorthogonal_cube_count = len(cube) - len(orthogonal_cube_bits)
        wall_rows.append(
            {
                "wall_one_based": index + 1,
                "norm": str(_pair(gram, roots[index], roots[index])),
                "pairings_with_path_walls": [
                    str(_pair(gram, roots[index], factor["a"])),
                    str(_pair(gram, roots[index], factor["b"])),
                ],
                "orthogonal_to_path_a": orthogonal_a,
                "orthogonal_to_path_b": orthogonal_b,
                "nonorthogonal_cube_direction_count": nonorthogonal_cube_count,
                "cube_orbit_multiple": 1 << nonorthogonal_cube_count,
            }
        )
    if actual != expected:
        missing = expected - actual
        extra = actual - expected
        raise ArithmeticError(
            f"global path/cube wall orbit census failed: missing={len(missing)}, extra={len(extra)}"
        )
    endpoint_collisions = 0
    for target_index in pair:
        for source_index in nonselected:
            for epsilon in (0, 1):
                for _mask, cube_action in cubes:
                    base = (
                        (factor["sa"] if epsilon else identity)
                        * cube_action
                        * roots[source_index]
                    )
                    for sign in (1, -1):
                        solutions = _translation_solutions(
                            gram, factor, base, sign * roots[target_index]
                        )
                        if solutions == "all" or solutions:
                            endpoint_collisions += 1
    if endpoint_collisions:
        raise ArithmeticError("a path endpoint orbit collides with a nonselected wall")
    return {
        "nonselected_wall_count": len(nonselected),
        "all_integer_translation_exponents_solved_symbolically": True,
        "cross_seed_wall_orbit_identification_count": 0,
        "opposite_orientation_identification_count": 0,
        "path_endpoint_to_nonselected_orbit_collision_count": 0,
        "stabilizers_exactly_generated_by_orthogonal_selected_reflections": True,
        "wall_component_data": wall_rows,
    }


def _path_family_proof(
    identifier: str, seed: dict[str, Any], spec: dict[str, Any]
) -> dict[str, Any]:
    gram: sp.Matrix = seed["gram"]
    roots: tuple[sp.Matrix, ...] = seed["roots"]
    pair: tuple[int, int] = spec["pair"]
    cube: tuple[int, ...] = spec["cube"]
    factor = _path_factor(gram, roots, pair)
    cube_elements = _cube_elements(gram, roots, cube)
    if not all(
        _doubling_wall_gate(gram, roots, index) for index in (*pair, *cube)
    ):
        raise ArithmeticError(f"{identifier}: a selected seed wall is not a doubling wall")
    if not all(
        factor["translation"] * action == action * factor["translation"]
        and factor["sa"] * action == action * factor["sa"]
        and action * factor["a"] == factor["a"]
        and action * factor["b"] == factor["b"]
        for _, action in cube_elements
    ):
        raise ArithmeticError(f"{identifier}: path and cube factors do not commute")
    orbits = _path_cube_global_orbits(gram, roots, pair, cube, factor)
    sum_multiples = sum(
        row["cube_orbit_multiple"] for row in orbits["wall_component_data"]
    )
    coefficient_a = sum(
        row["cube_orbit_multiple"]
        for row in orbits["wall_component_data"]
        if row["orthogonal_to_path_a"]
    )
    coefficient_b = sum(
        row["cube_orbit_multiple"]
        for row in orbits["wall_component_data"]
        if row["orthogonal_to_path_b"]
    )
    if coefficient_a != coefficient_b:
        raise ArithmeticError(f"{identifier}: an unrecorded path phase remains")
    slope = sum_multiples - coefficient_a
    intercept = 2 + coefficient_a
    if (slope, intercept) != (spec["slope"], spec["intercept"]):
        raise ArithmeticError(f"{identifier}: all-parameter facet formula changed")
    return {
        "seed": identifier,
        "path_factor": factor["public"],
        "cube_walls_one_based": [index + 1 for index in cube],
        "cube_factor": f"C2^{len(cube)}",
        "cube_action_order": len(cube_elements),
        "path_cube_direct_product_injective": True,
        "selected_cube_mirrors_are_internal": True,
        "all_selected_seed_walls_pass_exact_even_angle_doubling_gate": True,
        "selected_path_endpoint_facet_count": 2,
        "global_wall_orbit_audit": orbits,
        "component_formula": {
            "A": sum_multiples,
            "B_a": coefficient_a,
            "B_b": coefficient_b,
            "unsimplified": (
                f"2+{sum_multiples}*N-{coefficient_a}*floor(N/2)"
                f"-{coefficient_b}*floor((N-1)/2)"
            ),
            "identity_used": "floor(N/2)+floor((N-1)/2)=N-1",
            "simplified": spec["formula"],
            "all_N_at_least_one": True,
        },
        "copy_index_and_volume_multiple": f"{len(cube_elements)}*N",
        "convex_Coxeter_induction": (
            "iterate exact doubling along the infinite-dihedral gallery; "
            "the commuting cube mirrors remain orthogonal doubling mirrors"
        ),
    }


def _finite_cube_proof(seed: dict[str, Any]) -> dict[str, Any]:
    gram: sp.Matrix = seed["gram"]
    roots: tuple[sp.Matrix, ...] = seed["roots"]
    full_elements = _cube_elements(gram, roots, CUBE_WALLS)
    full_reflections = [_reflection(gram, roots[index]) for index in CUBE_WALLS]
    if not (
        len(full_elements) == 8
        and all(reflection**2 == sp.eye(gram.rows) for reflection in full_reflections)
        and all(
            _pair(gram, roots[left], roots[right]) == 0
            for left, right in itertools.combinations(CUBE_WALLS, 2)
        )
        and all(
            _doubling_wall_gate(gram, roots, index) for index in CUBE_WALLS
        )
    ):
        raise ArithmeticError("the finite C2^3 factor failed")

    rows = []
    relations_by_subset: dict[tuple[int, ...], tuple[tuple[str, ...], ...]] = {}
    for size in range(1, len(CUBE_WALLS) + 1):
        for selected in itertools.combinations(CUBE_WALLS, size):
            selected_set = set(selected)
            elements = _cube_elements(gram, roots, selected)
            nonselected = tuple(
                index for index in range(len(roots)) if index not in selected_set
            )
            actual_maps = set()
            expected_maps = set()
            for target_index in nonselected:
                for source_index in nonselected:
                    for mask, action in elements:
                        for sign in (1, -1):
                            if action * roots[source_index] == sign * roots[target_index]:
                                actual_maps.add(
                                    (target_index, source_index, mask, sign)
                                )
                orthogonal_bits = tuple(
                    bit
                    for bit, wall in enumerate(selected)
                    if _pair(gram, roots[target_index], roots[wall]) == 0
                )
                for mask in range(1 << len(selected)):
                    if all(
                        not (mask & (1 << bit))
                        for bit in range(len(selected))
                        if bit not in orthogonal_bits
                    ):
                        expected_maps.add(
                            (target_index, target_index, mask, 1)
                        )
            if actual_maps != expected_maps:
                raise ArithmeticError("a finite cube wall orbit has an extra collision")
            boundary: dict[tuple[int, ...], sp.Matrix] = {}
            for index in nonselected:
                for _mask, action in elements:
                    vector = action * roots[index]
                    boundary.setdefault(_unoriented_key(vector), vector)
            boundary_roots = tuple(boundary[key] for key in sorted(boundary))
            expected_count = sum(
                1
                << sum(
                    _pair(gram, roots[index], roots[wall]) != 0
                    for wall in selected
                )
                for index in nonselected
            )
            if len(boundary_roots) != expected_count:
                raise ArithmeticError("finite cube component formula failed")
            facet_gram = sp.Matrix(
                [
                    [_pair(gram, left, right) for right in boundary_roots]
                    for left in boundary_roots
                ]
            )
            relations = _relations(facet_gram)
            relations_by_subset[selected] = relations
            rows.append(
                {
                    "selected_walls_one_based": [index + 1 for index in selected],
                    "copy_count": len(elements),
                    "facet_count": len(boundary_roots),
                    "component_formula": (
                        "sum over unselected seed walls of "
                        "2^(number of selected nonorthogonal generators)"
                    ),
                    "relation_histogram": _histogram(relations),
                }
            )
    classes = []
    for selected, relations in relations_by_subset.items():
        match = next(
            (row for row in classes if _isomorphic(relations, row["relations"])),
            None,
        )
        if match is None:
            classes.append(
                {
                    "relations": relations,
                    "copy_count": 1 << len(selected),
                    "facet_count": len(relations),
                    "subsets": [[index + 1 for index in selected]],
                }
            )
        else:
            match["subsets"].append([index + 1 for index in selected])
    summary = [
        (
            row["copy_count"],
            row["facet_count"],
            len(row["subsets"]),
        )
        for row in classes
    ]
    if summary != [(2, 14, 3), (4, 15, 3), (8, 16, 1)]:
        raise ArithmeticError("the hostile finite-cube isomorphism census changed")
    return {
        "seed": CUBE_SEED,
        "walls_one_based": [index + 1 for index in CUBE_WALLS],
        "generator_gram": [
            [
                str(_pair(gram, roots[left], roots[right]))
                for right in CUBE_WALLS
            ]
            for left in CUBE_WALLS
        ],
        "group": "C2^3",
        "group_action_order": len(full_elements),
        "group_action_injective": True,
        "all_selected_seed_walls_pass_exact_even_angle_doubling_gate": True,
        "complete_nonempty_subset_count": len(rows),
        "global_wall_orbit_census_has_no_extra_identifications": True,
        "subset_rows": rows,
        "labelled_isomorphism_types": [
            {key: value for key, value in row.items() if key != "relations"}
            for row in classes
        ],
    }


def _symbolic_power(nilpotent: sp.Matrix, variable: sp.Symbol) -> sp.Matrix:
    return sp.simplify(
        sp.eye(nilpotent.rows)
        + variable * nilpotent
        + variable * (variable - 1) * nilpotent**2 / 2
    )


def _rectangle_proof(seed: dict[str, Any]) -> dict[str, Any]:
    gram: sp.Matrix = seed["gram"]
    roots: tuple[sp.Matrix, ...] = seed["roots"]
    first, second = RECTANGLE_PAIRS
    x = _path_factor(gram, roots, first)
    y = _path_factor(gram, roots, second)
    generators = (*first, *second)
    if not (
        x["delta"] == y["delta"]
        and all(
            _pair(gram, roots[left], roots[right]) == 0
            for left in first
            for right in second
        )
        and x["translation"] * y["translation"]
        == y["translation"] * x["translation"]
        and x["nilpotent"] * y["nilpotent"] == sp.zeros(gram.rows)
        and y["nilpotent"] * x["nilpotent"] == sp.zeros(gram.rows)
        and sp.Matrix.hstack(x["a"], y["a"], x["delta"]).rank() == 3
        and all(
            _doubling_wall_gate(gram, roots, index) for index in generators
        )
    ):
        raise ArithmeticError("the two affine-A1 directions do not form the claimed rectangle")
    fixed = tuple(
        index
        for index, root in enumerate(roots)
        if index not in generators
        and all(_pair(gram, root, roots[wall]) == 0 for wall in generators)
    )
    corners = tuple(
        index
        for index in range(len(roots))
        if index not in (*generators, *fixed)
    )
    if fixed != (2, 3, 4, 5, 6, 9, 10) or corners != (1, 8, 11, 12):
        raise ArithmeticError("the seven-fixed/four-corner decomposition changed")

    k, ell = sp.symbols("k ell", integer=True)
    symbolic_x = _symbolic_power(x["nilpotent"], k)
    symbolic_y = _symbolic_power(y["nilpotent"], ell)
    corner_specs = (
        ("AC", 1, (-4 * k, 4 * k + 2, -4 * ell, 4 * ell + 2)),
        ("AD", 8, (-4 * k, 4 * k + 2, 2 - 4 * ell, 4 * ell)),
        ("BC", 11, (4 - 8 * k, 8 * k, -8 * ell, 8 * ell + 4)),
        ("BD", 12, (4 - 8 * k, 8 * k, 4 - 8 * ell, 8 * ell)),
    )
    corner_rows = []
    for parity, index, expected_pairings in corner_specs:
        vector = sp.simplify(symbolic_x * symbolic_y * roots[index])
        actual_pairings = tuple(
            sp.factor((roots[wall].T * gram * vector)[0]) for wall in generators
        )
        if not all(
            sp.expand(actual - expected) == 0
            for actual, expected in zip(actual_pairings, expected_pairings, strict=True)
        ):
            raise ArithmeticError(f"rectangle corner {parity} orbit formula changed")
        norm = _pair(gram, roots[index], roots[index])
        delta_pair = _pair(gram, x["delta"], roots[index])
        fixed_x = x["sa"] if parity[0] == "A" else x["sb"]
        other_x = x["sb"] if parity[0] == "A" else x["sa"]
        fixed_y = y["sa"] if parity[1] == "C" else y["sb"]
        other_y = y["sb"] if parity[1] == "C" else y["sa"]
        expected_other_x = (
            _unipotent_power(sp.eye(gram.rows), x["nilpotent"], -1)
            if parity[0] == "A"
            else x["translation"]
        )
        expected_other_y = (
            _unipotent_power(sp.eye(gram.rows), y["nilpotent"], -1)
            if parity[1] == "C"
            else y["translation"]
        )
        if not (
            delta_pair > 0
            and fixed_x * roots[index] == roots[index]
            and fixed_y * roots[index] == roots[index]
            and other_x * roots[index] == expected_other_x * roots[index]
            and other_y * roots[index] == expected_other_y * roots[index]
        ):
            raise ArithmeticError("a corner orbit has the wrong orientation")
        corner_rows.append(
            {
                "parity_class": parity,
                "seed_wall_one_based": index + 1,
                "norm": str(norm),
                "positive_pairing_with_common_delta": str(delta_pair),
                "pairings_with_a_b_c_d": [str(value) for value in actual_pairings],
                "grid_coordinates": {
                    "AC": "i=2*k+1, j=2*ell+1",
                    "AD": "i=2*k+1, j=2*ell",
                    "BC": "i=2*k, j=2*ell+1",
                    "BD": "i=2*k, j=2*ell",
                }[parity],
                "full_group_orbit_equals_displayed_translation_orbit": True,
            }
        )

    # These exact congruences distinguish all four infinite orbits.  The
    # displayed linear pairings then recover k and ell inside each orbit.
    if not (
        _pair(gram, roots[1], roots[1])
        == _pair(gram, roots[8], roots[8])
        == -4
        and _pair(gram, roots[11], roots[11])
        == _pair(gram, roots[12], roots[12])
        == -8
        and all(
            _unoriented_key(roots[left]) != _unoriented_key(roots[right])
            for left, right in itertools.combinations((*fixed, *corners), 2)
        )
    ):
        raise ArithmeticError("rectangle orbit separation failed")

    # Four selected endpoint mirrors are mutually distinct for every box:
    # each direction has a+2z*delta, and the two direction spans meet only in
    # the null line.  Fixed walls pair zero with a,b,c,d, while corners pair
    # positively with delta, so neither can collide with an endpoint.
    endpoints_separate = (
        sp.Matrix.hstack(x["a"], x["delta"]).rank() == 2
        and sp.Matrix.hstack(y["a"], y["delta"]).rank() == 2
        and sp.Matrix.hstack(x["a"], y["a"], x["delta"]).rank() == 3
        and all(
            _pair(gram, roots[index], x["a"]) == 0
            and _pair(gram, roots[index], x["b"]) == 0
            and _pair(gram, roots[index], y["a"]) == 0
            and _pair(gram, roots[index], y["b"]) == 0
            for index in fixed
        )
    )
    if not endpoints_separate:
        raise ArithmeticError("rectangle endpoints have not been separated")
    return {
        "seed": RECTANGLE_SEED,
        "x_factor": x["public"],
        "y_factor": y["public"],
        "common_primitive_null_direction": [str(value) for value in x["delta"]],
        "cross_pair_gram_zero": True,
        "translations_commute": True,
        "all_four_seed_walls_pass_exact_even_angle_doubling_gate": True,
        "mixed_nilpotent_products_zero": True,
        "direct_product": {
            "group": "D_infinity x D_infinity",
            "injective": True,
            "proof": (
                "the y factor fixes a,b; an x normal form fixing a,b is the identity; "
                "then the y normal form is the identity"
            ),
        },
        "fixed_walls_one_based": [index + 1 for index in fixed],
        "fixed_wall_count": len(fixed),
        "corner_orbits": corner_rows,
        "corner_orbit_separation": {
            "x_parity_from_norm": "-4 means A/odd; -8 means B/even",
            "y_parity_from_c_pairing": (
                "mod 4 for A-orbits and mod 8 for B-orbits"
            ),
            "translation_coordinates_recovered_linearly": True,
            "opposite_sign_collision_excluded_by_positive_delta_pairing": True,
            "all_integer_k_ell_covered_symbolically": True,
        },
        "endpoint_orbits": {
            "count": 4,
            "distinct_for_all_integer_translation_exponents": True,
            "distinct_from_fixed_and_corner_orbits": True,
        },
        "all_parameter_component_count": {
            "x_A_components": "m-floor(m/2)",
            "x_B_components": "m-floor((m-1)/2)",
            "x_sum": "m+1",
            "y_C_components": "n-floor(n/2)",
            "y_D_components": "n-floor((n-1)/2)",
            "y_sum": "n+1",
            "corner_grid": "(m+1)*(n+1)",
            "facet_decomposition": "4 endpoints + 7 fixed + (m+1)*(n+1)",
            "facet_formula": "m*n+m+n+12",
            "all_m_n_at_least_one": True,
        },
        "copy_index_and_volume_multiple": "m*n",
        "convex_Coxeter_induction": (
            "iterate exact doubling in either commuting direction; the other pair "
            "is fixed because all cross products vanish"
        ),
    }


DIRECT_EXPECTED: dict[str, dict[str, int]] = {
    "A14_10_1_C01_D2_W3": {
        "facets": 14,
        "finite": 94,
        "ideal": 8,
        "edges": 515,
        "affine": 17,
    },
    "A14_10_1_C01_D4_W3_5": {
        "facets": 15,
        "finite": 121,
        "ideal": 10,
        "edges": 661,
        "affine": 22,
    },
    "A14_10_1_C01_D8_W3_5_11": {
        "facets": 16,
        "finite": 148,
        "ideal": 14,
        "edges": 828,
        "affine": 31,
    },
    "A14_8_1_C02_P2_X_C2": {
        "facets": 21,
        "finite": 334,
        "ideal": 60,
        "edges": 2407,
        "affine": 137,
    },
    "A14_10_1_C02_P2_X_C2SQ": {
        "facets": 20,
        "finite": 312,
        "ideal": 49,
        "edges": 2096,
        "affine": 109,
    },
    "A16_14_1_C01_RECT_2X2": {
        "facets": 20,
        "finite": 210,
        "ideal": 24,
        "edges": 1294,
        "affine": 58,
    },
}

DIRECT_BUILD: dict[str, dict[str, Any]] = {
    "A14_10_1_C01_D2_W3": {"kind": "cube", "selected": (2,)},
    "A14_10_1_C01_D4_W3_5": {"kind": "cube", "selected": (2, 4)},
    "A14_10_1_C01_D8_W3_5_11": {
        "kind": "cube",
        "selected": (2, 4, 10),
    },
    "A14_8_1_C02_P2_X_C2": {
        "kind": "path_cube",
        "pair": (14, 15),
        "cube": (6,),
        "N": 2,
    },
    "A14_10_1_C02_P2_X_C2SQ": {
        "kind": "path_cube",
        "pair": (5, 14),
        "cube": (0, 2),
        "N": 2,
    },
    "A16_14_1_C01_RECT_2X2": {
        "kind": "rectangle",
        "first": (0, 13),
        "second": (7, 14),
        "m": 2,
        "n": 2,
    },
}


def _path_representatives(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    pair: tuple[int, int],
    copies: int,
) -> tuple[sp.Matrix, ...]:
    factor = _path_factor(gram, roots, pair)
    rows = []
    for position in range(copies):
        exponent, odd = divmod(position, 2)
        action = _unipotent_power(
            sp.eye(gram.rows), factor["nilpotent"], exponent
        )
        if odd:
            action *= factor["sa"]
        rows.append(action)
    if len({_matrix_key(action) for action in rows}) != copies:
        raise ArithmeticError("finite path representatives collide")
    return tuple(rows)


def _product_actions(*factors: Sequence[sp.Matrix]) -> tuple[sp.Matrix, ...]:
    rows = [sp.eye(factors[0][0].rows)]
    for factor in factors:
        rows = [left * right for left in rows for right in factor]
    if len({_matrix_key(action) for action in rows}) != len(rows):
        raise ArithmeticError("a finite product action has a kernel")
    return tuple(rows)


def _direct_boundary_keys(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    transforms: Sequence[sp.Matrix],
    generator_indices: Sequence[int],
) -> set[tuple[int, ...]]:
    transform_keys = {_matrix_key(transform) for transform in transforms}
    generators = {
        index: _reflection(gram, roots[index]) for index in generator_indices
    }
    internal = set()
    for transform in transforms:
        for index, generator in generators.items():
            if _matrix_key(transform * generator) in transform_keys:
                internal.add(_unoriented_key(transform * roots[index]))
    boundary = set()
    for transform in transforms:
        for root in roots:
            key = _unoriented_key(transform * root)
            if key not in internal:
                boundary.add(key)
    return boundary


def _reconstruct_direct_boundary(
    identifier: str,
    gram: sp.Matrix,
    seed_roots: Sequence[sp.Matrix],
    archived_roots: Sequence[sp.Matrix],
) -> int:
    spec = DIRECT_BUILD[identifier]
    if spec["kind"] == "cube":
        selected = spec["selected"]
        transforms = tuple(action for _mask, action in _cube_elements(gram, seed_roots, selected))
        generators = selected
    elif spec["kind"] == "path_cube":
        pair = spec["pair"]
        cube = spec["cube"]
        path = _path_representatives(gram, seed_roots, pair, spec["N"])
        cube_actions = tuple(
            action for _mask, action in _cube_elements(gram, seed_roots, cube)
        )
        transforms = _product_actions(path, cube_actions)
        generators = (*pair, *cube)
    elif spec["kind"] == "rectangle":
        first = _path_representatives(
            gram, seed_roots, spec["first"], spec["m"]
        )
        second = _path_representatives(
            gram, seed_roots, spec["second"], spec["n"]
        )
        transforms = _product_actions(first, second)
        generators = (*spec["first"], *spec["second"])
    else:
        raise ArithmeticError("unknown hostile direct construction kind")
    expected = _direct_boundary_keys(gram, seed_roots, transforms, generators)
    archived = {_unoriented_key(root) for root in archived_roots}
    if expected != archived:
        raise ArithmeticError(f"{identifier}: independent boundary roots differ")
    return len(transforms)


def _parse_inequality(path: Path) -> sp.Matrix:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    begin = lines.index("begin")
    rows, columns, kind = lines[begin + 1].split()
    if kind != "rational":
        raise ArithmeticError("a direct cdd input is not rational")
    payload = [
        [sp.Rational(value) for value in lines[begin + 2 + index].split()]
        for index in range(int(rows))
    ]
    if any(len(row) != int(columns) or row[0] != 0 for row in payload):
        raise ArithmeticError("unexpected inhomogeneous direct cdd input")
    return sp.Matrix([row[1:] for row in payload])


def _cdd_rays(path: Path, executable: str) -> tuple[tuple[int, ...], ...]:
    completed = subprocess.run(
        [executable, "--rep"],
        input=path.read_bytes(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    lines = completed.stdout.decode("utf-8", errors="strict").splitlines()
    marker = lines.index("V-representation")
    begin = lines.index("begin", marker)
    count, columns, kind = lines[begin + 1].split()
    if kind != "rational":
        raise ArithmeticError("fresh cdd output is not rational")
    result = []
    for index in range(int(count)):
        row = [sp.Rational(value) for value in lines[begin + 2 + index].split()]
        if len(row) != int(columns) or row[0] != 0:
            raise ArithmeticError("fresh homogeneous cdd output contains a non-ray")
        result.append(_primitive(row[1:]))
    return tuple(result)


def _ray_archive(path: Path) -> tuple[dict[str, Any], tuple[list[Any], ...]]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle if line.strip()]
    return rows[0], tuple(rows[1:])


def _direct_geometry(
    record: dict[str, Any], seeds: dict[str, dict[str, Any]], executable: str
) -> dict[str, Any]:
    identifier = record["id"]
    expected = DIRECT_EXPECTED[identifier]
    seed = seeds[record["seed"]]
    gram: sp.Matrix = seed["gram"]
    roots = tuple(
        sp.Matrix([sp.Rational(value) for value in row])
        for row in record["roots_in_seed_fixed_lattice_coordinates"]
    )
    reconstructed_copy_count = _reconstruct_direct_boundary(
        identifier, gram, seed["roots"], roots
    )
    facet_gram = sp.Matrix(
        [[_pair(gram, left, right) for right in roots] for left in roots]
    )
    archived_facet_gram = _matrix(record["facet_gram"])
    inequality_path = ROOT / record["artifacts"]["inequality_path"]
    ray_path = ROOT / record["artifacts"]["ray_archive_path"]
    if not (
        _sha(inequality_path) == record["artifacts"]["inequality_sha256"]
        and _sha(ray_path) == record["artifacts"]["ray_archive_sha256"]
    ):
        raise ArithmeticError(f"{identifier}: a binary artifact hash changed")
    inequality = _parse_inequality(inequality_path)
    if not (
        len(roots) == expected["facets"]
        and record["facet_count"] == expected["facets"]
        and reconstructed_copy_count == record["copy_count_and_volume_multiple"]
        and facet_gram == archived_facet_gram
        and inequality == sp.Matrix.vstack(*(root.T * gram for root in roots))
        and inequality.rank() == 10
        and all(reduce(math.gcd, (abs(int(value)) for value in root if value)) == 1 for root in roots)
        and _histogram(_relations(facet_gram)) == record["relation_histogram"]
    ):
        raise ArithmeticError(f"{identifier}: direct lattice/cone reconstruction failed")
    reflections = tuple(_reflection(gram, root) for root in roots)
    if not all(
        all(value.is_Integer for value in reflection)
        and abs(reflection.det()) == 1
        and reflection.T * gram * reflection == gram
        for reflection in reflections
    ):
        raise ArithmeticError(f"{identifier}: a direct facet reflection is nonintegral")

    header, ray_records = _ray_archive(ray_path)
    archived_rays = tuple(tuple(map(int, row[0])) for row in ray_records)
    fresh_rays = _cdd_rays(inequality_path, executable)
    if not (
        header["case_id"] == identifier
        and header["ray_count"] == len(ray_records)
        and len(set(archived_rays)) == len(archived_rays)
        and set(fresh_rays) == set(archived_rays)
    ):
        raise ArithmeticError(f"{identifier}: fresh GMP rays differ from the archive")

    vectors = tuple(sp.Matrix(row[0]) for row in ray_records)
    masks = []
    kinds = []
    for vector, row in zip(vectors, ray_records, strict=True):
        active_positions = row[1]
        square_payload = row[2]
        kind = row[3]
        slacks = inequality * vector
        actual_active = [index + 1 for index, value in enumerate(slacks) if value == 0]
        square = _pair(gram, vector, vector)
        expected_kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        active_matrix = sp.Matrix.vstack(
            *(inequality.row(index - 1) for index in active_positions)
        )
        if not (
            all(value >= 0 for value in slacks)
            and actual_active == active_positions
            and str(square) == str(square_payload)
            and kind == expected_kind
            and active_matrix.rank() == 9
        ):
            raise ArithmeticError(f"{identifier}: ray incidence/type replay failed")
        masks.append(sum(1 << (index - 1) for index in active_positions))
        kinds.append(kind)

    strict_interior = sum(vectors, sp.zeros(gram.rows, 1))
    interior_slacks = inequality * strict_interior
    time_pairings = [_pair(gram, strict_interior, vector) for vector in vectors]
    if not (
        _pair(gram, strict_interior, strict_interior) > 0
        and all(value > 0 for value in interior_slacks)
        and all(value > 0 for value in time_pairings)
    ):
        raise ArithmeticError(f"{identifier}: hostile time-sheet reconstruction failed")

    witness_squares = []
    for target in range(len(roots)):
        witness = sum(
            (
                vectors[index]
                for index, active in enumerate(row[1] for row in ray_records)
                if target + 1 in active
            ),
            sp.zeros(gram.rows, 1),
        )
        slacks = inequality * witness
        square = _pair(gram, witness, witness)
        if not (
            slacks[target] == 0
            and square > 0
            and all(slacks[index] > 0 for index in range(len(roots)) if index != target)
        ):
            raise ArithmeticError(f"{identifier}: an independent strict facet witness failed")
        witness_squares.append(square)

    positive_gram = -facet_gram
    enumerator = DiagramEnumerator(positive_gram, dimension=9)
    finite_masks = {mask for mask, kind in zip(masks, kinds, strict=True) if kind == "finite"}
    ideal_masks = {mask for mask, kind in zip(masks, kinds, strict=True) if kind == "ideal"}
    finite_vertices = enumerator.finite_subdiagrams(9)
    finite_edges = enumerator.finite_subdiagrams(8)
    if set(finite_vertices) != finite_masks:
        raise ArithmeticError(f"{identifier}: elliptic vertices differ from GMP")
    edge_endpoints: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for right in range(len(masks)):
        for left in range(right):
            intersection = masks[left] & masks[right]
            if intersection.bit_count() >= 8:
                edge_endpoints[intersection].append((left, right))
    if not (
        set(edge_endpoints) == set(finite_edges)
        and all(mask.bit_count() == 8 for mask in edge_endpoints)
        and all(len(pairs) == 1 for pairs in edge_endpoints.values())
    ):
        raise ArithmeticError(f"{identifier}: literal elliptic edge gate failed")

    affine = enumerator.connected_affine_subdiagrams()
    neighborhoods = {}
    for component_mask in affine:
        component = tuple(
            index for index in range(len(roots)) if component_mask & (1 << index)
        )
        zero = tuple(
            facet
            for facet in range(len(roots))
            if facet not in component
            and all(positive_gram[facet, node] == 0 for node in component)
        )
        neighborhood = component_mask | sum(1 << index for index in zero)
        if neighborhood not in ideal_masks:
            raise ArithmeticError(f"{identifier}: affine zero-neighborhood is not a cusp")
        neighborhoods[component_mask] = neighborhood
    for cusp in ideal_masks:
        components = enumerator._components(cusp, all_nonorthogonal=True)
        if not (
            sum(component.bit_count() - 1 for component in components) == 8
            and all(component in affine for component in components)
            and all(neighborhoods[component] == cusp for component in components)
        ):
            raise ArithmeticError(f"{identifier}: cusp is not a full parabolic neighborhood")

    counts = Counter(kinds)
    archived_exact = record["exact_cone"]["exact_double_description"]
    archived_vinberg = record["literal_vinberg"]
    if not (
        counts == {"finite": expected["finite"], "ideal": expected["ideal"]}
        and len(finite_edges) == expected["edges"]
        and len(affine) == expected["affine"]
        and archived_exact["finite_vertex_count"] == expected["finite"]
        and archived_exact["ideal_vertex_count"] == expected["ideal"]
        and archived_exact["hyperideal_vertex_count"] == 0
        and archived_vinberg["elliptic_edge_count"] == expected["edges"]
        and archived_vinberg["connected_affine_component_count"] == expected["affine"]
        and archived_vinberg["bad_edge_count"]
        == archived_vinberg["bad_parabolic_neighborhood_count"]
        == 0
    ):
        raise ArithmeticError(f"{identifier}: direct archived summary differs from replay")
    return {
        "id": identifier,
        "seed": record["seed"],
        "copy_count_and_volume_multiple": record["copy_count_and_volume_multiple"],
        "independent_group_boundary_reconstruction": True,
        "facet_count": len(roots),
        "integral_reflection_count": len(reflections),
        "independent_strict_timelike_facet_witness_count": len(witness_squares),
        "fresh_GMP_extreme_ray_count": len(fresh_rays),
        "finite_vertex_count": len(finite_masks),
        "ideal_vertex_count": len(ideal_masks),
        "hyperideal_vertex_count": 0,
        "elliptic_edge_count": len(finite_edges),
        "connected_affine_component_count": len(affine),
        "all_exact_checks_pass": True,
        "_relations": _relations(facet_gram),
    }


def _corpus_relations(matrix: Sequence[Sequence[str]]) -> tuple[tuple[str, ...], ...]:
    replacements = {
        "0": "parallel",
        "infinity": "parallel",
        "∞": "parallel",
        "1": "dotted:unspecified",
    }
    return tuple(
        tuple(
            "1" if left == right else replacements.get(str(value), str(value))
            for right, value in enumerate(row)
        )
        for left, row in enumerate(matrix)
    )


def _coxiter_comparison(
    targets: dict[str, tuple[tuple[str, ...], ...]]
) -> dict[str, Any]:
    if not (
        _sha(COXITER_ZIP) == COXITER_ZIP_SHA256
        and _sha(DOUBLING_TAR) == DOUBLING_TAR_SHA256
    ):
        raise ArithmeticError("a pinned CoxIter snapshot changed")
    target_facets = {len(relations) for relations in targets.values()}
    ordinary_histogram: Counter[int] = Counter()
    doubling_histogram: Counter[int] = Counter()
    records: list[dict[str, Any]] = []
    with zipfile.ZipFile(COXITER_ZIP) as archive:
        ordinary_names = sorted(name for name in archive.namelist() if name.endswith(".coxiter"))
        for name in ordinary_names:
            raw = archive.read(name)
            facets, dimension, matrix = _parse_coxiter_text(
                raw.decode("utf-8", errors="replace")
            )
            if dimension == 9:
                ordinary_histogram[facets] += 1
                if facets in target_facets:
                    records.append(
                        {
                            "lane": "ordinary",
                            "path": name,
                            "facets": facets,
                            "relations": _corpus_relations(matrix),
                        }
                    )
    doubling_count = 0
    with tarfile.open(DOUBLING_TAR, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            doubling_count += 1
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"cannot read {member.name}")
            raw = handle.read()
            facets, dimension, matrix = _parse_coxiter_text(
                raw.decode("utf-8", errors="replace")
            )
            if dimension == 9:
                doubling_histogram[facets] += 1
                if facets in target_facets:
                    records.append(
                        {
                            "lane": "index_two_doubling",
                            "path": member.name,
                            "facets": facets,
                            "relations": _corpus_relations(matrix),
                        }
                    )
    if not (
        len(ordinary_names) == 1214
        and doubling_count == 2959
        and dict(sorted(ordinary_histogram.items())) == {10: 38, 11: 10, 12: 13}
        and dict(sorted(doubling_histogram.items()))
        == {12: 8, 13: 8, 15: 8, 19: 8, 27: 8, 43: 8}
    ):
        raise ArithmeticError("the pinned H9 header inventory changed")
    comparisons = []
    total_matches = 0
    for identifier, target in targets.items():
        same = [record for record in records if record["facets"] == len(target)]
        matches = [
            {"lane": record["lane"], "path": record["path"]}
            for record in same
            if _isomorphic(target, record["relations"])
        ]
        histogram_collisions = [
            record["path"]
            for record in same
            if _histogram(target) == _histogram(record["relations"])
        ]
        total_matches += len(matches)
        comparisons.append(
            {
                "id": identifier,
                "facet_count": len(target),
                "same_header_record_count": len(same),
                "relation_histogram_collision_paths": histogram_collisions,
                "exact_encoded_labelled_graph_matches": matches,
            }
        )
    if total_matches:
        raise ArithmeticError("a direct product representative matches pinned CoxIter")
    return {
        "scope": (
            "bounded exact encoded-labelled-graph comparison; unspecified dotted "
            "weights are not promoted to exact distances; no worldwide priority claim"
        ),
        "ordinary_archive_sha256": COXITER_ZIP_SHA256,
        "doubling_archive_sha256": DOUBLING_TAR_SHA256,
        "ordinary_file_count": len(ordinary_names),
        "doubling_file_count": doubling_count,
        "ordinary_H9_header_histogram": {
            str(key): value for key, value in sorted(ordinary_histogram.items())
        },
        "doubling_H9_header_histogram": {
            str(key): value for key, value in sorted(doubling_histogram.items())
        },
        "comparisons": comparisons,
        "total_exact_encoded_labelled_graph_match_count": 0,
    }


def audit(executable: str | None = None) -> dict[str, Any]:
    if _sha(PRODUCTS) != PRODUCTS_SHA256:
        raise ArithmeticError("the pinned product certificate changed")
    products = json.loads(PRODUCTS.read_text(encoding="utf-8"))
    if not (
        products["certificate_payload_sha256"] == PRODUCTS_PAYLOAD_SHA256
        and products["schema"] == "coxeter.alexeev_parked_h9_doubling_products.v1"
        and {record["id"] for record in products["direct_exact_audits"]}
        == set(DIRECT_EXPECTED)
    ):
        raise ArithmeticError("the product-certificate denominator changed")
    seeds = _load_seeds()
    finite_cube = _finite_cube_proof(seeds[CUBE_SEED])
    path_families = {
        identifier: _path_family_proof(identifier, seeds[identifier], spec)
        for identifier, spec in PATH_FAMILIES.items()
    }
    rectangle = _rectangle_proof(seeds[RECTANGLE_SEED])
    for identifier, proof in path_families.items():
        archived = products["infinite_families"][identifier]
        if not (
            proof["component_formula"]["simplified"] == archived["facet_formula"]
            and proof["copy_index_and_volume_multiple"]
            == archived["index_and_volume_multiple"]
        ):
            raise ArithmeticError(f"{identifier}: hostile formula differs from archive")
    archived_rectangle = products["infinite_families"][RECTANGLE_SEED]
    if not (
        rectangle["all_parameter_component_count"]["facet_formula"]
        == archived_rectangle["facet_formula"]
        and rectangle["copy_index_and_volume_multiple"]
        == archived_rectangle["index_and_volume_multiple"]
    ):
        raise ArithmeticError("the hostile rectangle formula differs from archive")

    cdd = executable or shutil.which("cddexec_gmp") or "cddexec_gmp"
    direct_internal = [
        _direct_geometry(record, seeds, cdd)
        for record in products["direct_exact_audits"]
    ]
    targets = {record["id"]: record["_relations"] for record in direct_internal}
    coxiter = _coxiter_comparison(targets)
    direct = [
        {key: value for key, value in record.items() if not key.startswith("_")}
        for record in direct_internal
    ]
    return {
        "schema": SCHEMA,
        "pinned_inputs": {
            "seed_archive_path": str(SOURCE.relative_to(ROOT)),
            "seed_archive_sha256": SOURCE_SHA256,
            "product_certificate_path": str(PRODUCTS.relative_to(ROOT)),
            "product_certificate_sha256": PRODUCTS_SHA256,
            "product_certificate_internal_payload_sha256": PRODUCTS_PAYLOAD_SHA256,
        },
        "independence": {
            "imports_product_certifier": False,
            "imports_product_scout": False,
            "calls_product_build": False,
            "all_integer_translation_exponents_solved_symbolically": True,
            "fresh_GMP_double_description_for_six_rows": True,
            "independent_ray_incidence_time_facet_and_Vinberg_replay": True,
        },
        "finite_C2_cubed": finite_cube,
        "infinite_path_cube_families": path_families,
        "infinite_rectangle_family": rectangle,
        "direct_exact_replays": direct,
        "coxiter_bounded_comparison": coxiter,
        "theorem_status_split": {
            "six_explicit_descendants": {
                "status": "independently_reconstructed_and_exactly_replayed",
                "count": len(direct),
                "gates": [
                    "independent group-boundary root reconstruction",
                    "fresh GMP extreme-ray equality",
                    "independent strict timelike witness for every facet",
                    "one-time-cone and zero-hyperideal-ray check",
                    "literal Vinberg elliptic-edge and full-cusp-neighborhood replay",
                ],
            },
            "finite_C2_cubed_census": {
                "status": "complete_nonempty_subset_census",
                "marked_subset_count": finite_cube["complete_nonempty_subset_count"],
                "labelled_type_count": len(finite_cube["labelled_isomorphism_types"]),
                "gates": [
                    "faithful commuting C2^3 reflection action",
                    "exact global wall stabilizers for every subset",
                    "component facet formula",
                    "exact labelled-graph isomorphism quotient",
                ],
            },
            "three_all_parameter_families": {
                "status": "universal_closed_orbit_proofs",
                "A14_8_1_C02": {
                    "facet_formula": "4*N+13",
                    "index": "2*N",
                    "universal_gates": [
                        "N^3=0 parabolic power identity for every integer exponent",
                        "faithful D_infinity x C2 action",
                        "all-z exact wall-stabilizer and collision census",
                        "component identity using floor(N/2)+floor((N-1)/2)=N-1",
                        "iterated convex Coxeter doubling",
                    ],
                },
                "A14_10_1_C02": {
                    "facet_formula": "3*N+14",
                    "index": "4*N",
                    "universal_gates": [
                        "N^3=0 parabolic power identity for every integer exponent",
                        "faithful D_infinity x C2^2 action",
                        "all-z exact wall-stabilizer and collision census",
                        "component identity using floor(N/2)+floor((N-1)/2)=N-1",
                        "iterated convex Coxeter doubling",
                    ],
                },
                "A16_14_1_C01": {
                    "facet_formula": "m*n+m+n+12",
                    "index": "m*n",
                    "universal_gates": [
                        "two N^3=0 parabolic power identities",
                        "faithful D_infinity x D_infinity action",
                        "four symbolic parity-corner orbit formulae for all integer k,ell",
                        "linear recovery of every grid coordinate and exclusion of extra identifications",
                        "4 endpoint + 7 fixed + (m+1)(n+1) component count",
                        "iterated convex Coxeter doubling in two commuting directions",
                    ],
                },
            },
        },
        "scope": {
            "all_parameter_group_and_facet_formulas_proved": True,
            "finite_volume_noncompact_for_all_parameters": True,
            "finite_volume_reason": (
                "exact convex Coxeter doubling induction gives a finite union of "
                "the displayed number of seed chambers"
            ),
            "noncompactness_reason": (
                "each constructed reflection subgroup has finite index in its "
                "nonuniform seed group"
            ),
            "mechanism_is_classical_Allcock_doubling": True,
            "worldwide_novelty_claimed": False,
            "surface_automorphism_or_nef_claimed": False,
        },
        "decision": {
            "formula_audit": "all_parameter_closed_orbit_proofs_pass",
            "direct_geometry": "six_exact_finite_volume_noncompact_H9_replays_pass",
            "bounded_lineage": "zero_pinned_CoxIter_encoded_labelled_graph_matches",
            "priority_status": "new_to_project_worldwide_priority_open",
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--cdd", default=None)
    args = parser.parse_args()
    payload = audit(args.cdd)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if OUTPUT.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("hostile product audit archive changed")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
    print(json.dumps(payload["decision"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
