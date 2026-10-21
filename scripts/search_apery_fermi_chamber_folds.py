"""Parse and screen Shimada's 80-wall Apéry--Fermi K3 chamber.

The source is the GAP-format ``AperyFermiCompData.txt`` archived at
Zenodo record 14842583.  This program deliberately uses only the weighted
wall arrangement and the finite chamber stabilizer.  It does *not* require a
chosen wall permutation to have any additional geometric interpretation.

For every subgroup of the full 32-element weighted chamber stabilizer, the
program:

* derives the wall permutation from the published 19 by 19 lattice matrix;
* forms exact orbit-sum wall normals in the fixed Lorentzian space;
* removes a wall only with an exact nonnegative-span certificate;
* proves every retained wall is a genuine hyperbolic facet using a strict
  timelike point on that wall; and
* checks every normalized wall product against the Coxeter angle list.

SciPy is used only to discover small supports for redundancy certificates.
Every accepted support is then solved and verified over the rationals.  Thus
no mathematical conclusion rests on floating-point output.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, deque
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import re
from typing import Any, Iterable

import networkx as nx
import numpy as np
from scipy.optimize import linprog
import sympy as sp

from coxeter_polytopes.verify.vinberg_volume import (
    verify_vinberg_volume_criterion,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _coxeter_label


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data"
    / "external"
    / "apery_fermi_shimada_2025"
    / "AperyFermiCompData.txt"
)
STOCK_OUTPUT = ROOT / "data" / "known_examples" / "apery_fermi_chamber_stock.json"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "apery_fermi_chamber_subgroup_folds.json"
)

Permutation = tuple[int, ...]
Subgroup = frozenset[Permutation]


def _extract_assignment(text: str, name: str) -> str:
    """Return one top-level bracketed GAP assignment as text."""

    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    if text[start] != "[":
        raise ValueError(f"assignment {name} is not a bracketed list")

    depth = 0
    quote: str | None = None
    for end in range(start, len(text)):
        character = text[end]
        if quote is not None:
            if character == quote and text[end - 1] != "\\":
                quote = None
            continue
        if character in {'"', "'"}:
            quote = character
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1]
    raise ValueError(f"unterminated assignment {name}")


def _rational(text: str) -> sp.Rational:
    return sp.Rational(text.strip())


def _rational_text(value: sp.Rational) -> str:
    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


def _matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [_rational_text(sp.Rational(matrix[row, column])) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _parse_source(path: Path = SOURCE) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    gram = sp.Matrix(ast.literal_eval(_extract_assignment(text, "GramSX")))
    ample = sp.Matrix(1, 19, ast.literal_eval(_extract_assignment(text, "a32")))
    orthogonal_stabilizer = tuple(
        sp.Matrix(matrix)
        for matrix in ast.literal_eval(_extract_assignment(text, "OSXD0"))
    )
    geometric_stabilizer = tuple(
        sp.Matrix(matrix)
        for matrix in ast.literal_eval(_extract_assignment(text, "AutXD0"))
    )

    walls_text = _extract_assignment(text, "WallsD0")
    vector_texts = re.findall(r"\bvect\s*:=\s*\[([^\]]+)\]", walls_text)
    innout = re.findall(r'\binnout\s*:=\s*"([^"]+)"', walls_text)
    orbit_numbers = [
        int(value)
        for value in re.findall(r"\borbnumb\s*:=\s*(\d+)", walls_text)
    ]
    if not (len(vector_texts) == len(innout) == len(orbit_numbers) == 80):
        raise ValueError("the WallsD0 record count is not 80")
    wall_vectors = tuple(
        tuple(_rational(entry) for entry in vector.split(","))
        for vector in vector_texts
    )
    if any(len(vector) != 19 for vector in wall_vectors):
        raise ValueError("a wall vector does not have 19 coordinates")

    return {
        "gram": gram,
        "ample": ample,
        # The census uses the full lattice/chamber stabilizer.  The geometric
        # subgroup is retained as metadata, never as a search gate.
        "stabilizer": orthogonal_stabilizer,
        "orthogonal_stabilizer": orthogonal_stabilizer,
        "geometric_stabilizer": geometric_stabilizer,
        "wall_vectors": wall_vectors,
        "wall_types": tuple(innout),
        "wall_orbit_numbers": tuple(orbit_numbers),
    }


def _compose(left: Permutation, right: Permutation) -> Permutation:
    """Return ``left`` after ``right``."""

    return tuple(left[right[index]] for index in range(len(left)))


def _identity(size: int) -> Permutation:
    return tuple(range(size))


def _permutation_order(permutation: Permutation) -> int:
    identity = _identity(len(permutation))
    power = identity
    for exponent in range(1, 10_000):
        power = _compose(permutation, power)
        if power == identity:
            return exponent
    raise ArithmeticError("permutation-order guard exceeded")


def _generated_subgroup(generators: Iterable[Permutation]) -> Subgroup:
    generators = tuple(generators)
    if not generators:
        raise ValueError("at least one sized permutation is required")
    identity = _identity(len(generators[0]))
    subgroup = {identity, *generators}
    queue = deque(subgroup)
    while queue:
        left = queue.popleft()
        for right in tuple(subgroup):
            for product in (_compose(left, right), _compose(right, left)):
                if product not in subgroup:
                    subgroup.add(product)
                    queue.append(product)
    return frozenset(subgroup)


def _all_subgroups(group: tuple[Permutation, ...]) -> tuple[Subgroup, ...]:
    identity = _identity(len(group[0]))
    subgroups: set[Subgroup] = {frozenset({identity})}
    queue = deque(subgroups)
    while queue:
        subgroup = queue.popleft()
        for element in group:
            if element in subgroup:
                continue
            generated = _generated_subgroup((*subgroup, element))
            if generated not in subgroups:
                subgroups.add(generated)
                queue.append(generated)
    return tuple(
        sorted(subgroups, key=lambda subgroup: (len(subgroup), tuple(sorted(subgroup))))
    )


def _orbits(subgroup: Subgroup, size: int) -> tuple[tuple[int, ...], ...]:
    unseen = set(range(size))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        orbit = tuple(sorted({element[start] for element in subgroup}))
        unseen.difference_update(orbit)
        orbits.append(orbit)
    return tuple(orbits)


def _derive_wall_permutations(
    source: dict[str, Any], matrices: tuple[sp.Matrix, ...] | None = None
) -> tuple[Permutation, ...]:
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    if matrices is None:
        matrices = source["stabilizer"]
    wall_vectors: tuple[tuple[sp.Rational, ...], ...] = source["wall_vectors"]
    walls = tuple(sp.Matrix(1, gram.rows, vector) for vector in wall_vectors)
    index = {vector: position for position, vector in enumerate(wall_vectors)}
    if len(index) != len(walls):
        raise ValueError("the published wall vectors are not distinct")

    permutations: list[Permutation] = []
    for matrix in matrices:
        if matrix * gram * matrix.T != gram:
            raise ArithmeticError("a published chamber stabilizer is not an isometry")
        if ample * matrix != ample:
            raise ArithmeticError("a published chamber stabilizer does not fix a32")
        images = tuple(tuple(wall * matrix) for wall in walls)
        if any(image not in index for image in images):
            raise ArithmeticError("a chamber stabilizer does not permute WallsD0")
        permutations.append(tuple(index[image] for image in images))

    group = tuple(permutations)
    if len(set(group)) != len(matrices):
        raise ArithmeticError("the published chamber-stabilizer wall action is not faithful")
    if _generated_subgroup(group) != frozenset(group):
        raise ArithmeticError("the published chamber-stabilizer list is not closed")
    return group


def _weighted_wall_graph_automorphisms(
    gram: sp.Matrix, wall_vectors: tuple[tuple[sp.Rational, ...], ...]
) -> tuple[Permutation, ...]:
    """Enumerate every abstract symmetry of the exact weighted wall Gram."""

    walls = tuple(sp.Matrix(1, gram.rows, vector) for vector in wall_vectors)
    graph = nx.Graph()
    for index, wall in enumerate(walls):
        graph.add_node(index, norm=_rational_text((wall * gram * wall.T)[0]))
    for left in range(len(walls)):
        for right in range(left):
            product = sp.Rational((walls[left] * gram * walls[right].T)[0])
            if product:
                graph.add_edge(left, right, product=_rational_text(product))
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=nx.algorithms.isomorphism.categorical_node_match("norm", None),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("product", None),
    )
    return tuple(
        sorted(
            {
                tuple(int(mapping[index]) for index in range(len(walls)))
                for mapping in matcher.isomorphisms_iter()
            }
        )
    )


def _fixed_rank(matrices: Iterable[sp.Matrix], ambient_rank: int) -> int:
    constraints = [matrix.T - sp.eye(ambient_rank) for matrix in matrices]
    stacked = sp.Matrix.vstack(*constraints)
    return ambient_rank - stacked.rank()


def _exact_nonnegative_dependence(
    coordinates: sp.Matrix, target: int
) -> tuple[tuple[int, sp.Rational], ...] | None:
    """Find and verify a rational conic dependence, or return ``None``.

    Rows of ``coordinates`` are exact vectors in an injective coordinate
    system.  HiGHS discovers a basic support; SymPy then proves the result.
    """

    others = [index for index in range(coordinates.rows) if index != target]
    matrix_float = np.array(coordinates[others, :].T.tolist(), dtype=float)
    target_float = np.array(coordinates[target, :].tolist(), dtype=float).reshape(-1)
    solution = linprog(
        np.zeros(len(others)),
        A_eq=matrix_float,
        b_eq=target_float,
        bounds=(0, None),
        method="highs",
    )
    if not solution.success:
        return None

    support = [
        others[position]
        for position, coefficient in enumerate(solution.x)
        if coefficient > 1e-8
    ]
    if not support:
        raise ArithmeticError("floating conic support is unexpectedly empty")
    exact_matrix = coordinates[support, :].T
    exact_target = coordinates[target, :].T
    solution_set = sp.linsolve((exact_matrix, exact_target))
    if solution_set is sp.EmptySet:
        raise ArithmeticError("floating conic support has no exact solution")
    coefficients = tuple(next(iter(solution_set)))
    if any(not coefficient.is_Rational for coefficient in coefficients):
        raise ArithmeticError("floating support did not yield a unique rational solution")
    rational_coefficients = tuple(sp.Rational(value) for value in coefficients)
    if any(coefficient <= 0 for coefficient in rational_coefficients):
        raise ArithmeticError("recovered conic dependence is not strictly positive")
    if exact_matrix * sp.Matrix(rational_coefficients) != exact_target:
        raise ArithmeticError("recovered conic dependence failed exact replay")
    return tuple(zip(support, rational_coefficients, strict=True))


def _vinberg_payload(facet_gram: sp.Matrix, dimension: int) -> dict[str, Any]:
    positive_wall_gram = [
        [
            Fraction(
                int(-sp.Rational(facet_gram[left, right]).p),
                int(sp.Rational(facet_gram[left, right]).q),
            )
            for right in range(facet_gram.cols)
        ]
        for left in range(facet_gram.rows)
    ]
    verification = verify_vinberg_volume_criterion(
        positive_wall_gram,
        dimension,
        simple_roots_complete=True,
    )
    endpoint_histogram = Counter(
        audit.endpoint_count for audit in verification.edge_audits
    )
    return {
        "decision": verification.decision,
        "finite_volume": verification.finite_volume,
        "compact": verification.compact,
        "elliptic_edge_count": verification.elliptic_edge_count,
        "finite_vertex_count": verification.finite_vertex_count,
        "ideal_vertex_count": verification.ideal_vertex_count,
        "bad_edge_count": len(verification.finite_volume_bad_edges),
        "edge_endpoint_count_histogram": {
            str(key): value for key, value in sorted(endpoint_histogram.items())
        },
    }


def _screen_subgroup(
    source: dict[str, Any],
    group: tuple[Permutation, ...],
    subgroup: Subgroup,
    subgroup_index: int,
) -> dict[str, Any]:
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    walls = tuple(sp.Matrix(1, gram.rows, vector) for vector in source["wall_vectors"])
    permutation_index = {permutation: index for index, permutation in enumerate(group)}
    matrix_elements = tuple(
        source["stabilizer"][permutation_index[element]] for element in subgroup
    )
    wall_orbits = _orbits(subgroup, len(walls))
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, gram.rows))
        for orbit in wall_orbits
    )
    stacked_normals = sp.Matrix.vstack(*orbit_normals)
    orbit_rank = stacked_normals.rank()
    fixed_rank = _fixed_rank(matrix_elements, gram.rows)
    if orbit_rank != fixed_rank:
        raise ArithmeticError("orbit normals do not span the full fixed space")
    _, pivot_columns = stacked_normals.rref()
    coordinates = stacked_normals[:, list(pivot_columns)]
    if coordinates.cols != fixed_rank:
        raise ArithmeticError("the fixed-space coordinate projection is not injective")

    orbit_gram = sp.Matrix(
        [
            [(left * gram * right.T)[0] for right in orbit_normals]
            for left in orbit_normals
        ]
    )
    inertia = exact_inertia(orbit_gram)
    if inertia != (1, fixed_rank - 1, len(wall_orbits) - fixed_rank):
        raise ArithmeticError("the fixed wall span is not Lorentzian")

    redundant: list[dict[str, Any]] = []
    facets: list[int] = []
    facet_witnesses: list[dict[str, Any]] = []
    for target, normal in enumerate(orbit_normals):
        dependence = _exact_nonnegative_dependence(coordinates, target)
        if dependence is not None:
            redundant.append(
                {
                    "orbit_index": target + 1,
                    "terms": [
                        {
                            "orbit_index": support + 1,
                            "coefficient": _rational_text(coefficient),
                        }
                        for support, coefficient in dependence
                    ],
                }
            )
            continue

        norm = sp.Rational(orbit_gram[target, target])
        if norm >= 0:
            raise ArithmeticError("a nonspacelike extreme inequality remains unresolved")
        ample_pairing = sp.Rational((ample * gram * normal.T)[0])
        witness = ample - ample_pairing / norm * normal
        witness_square = sp.Rational((witness * gram * witness.T)[0])
        pairings = [sp.Rational((witness * gram * other.T)[0]) for other in orbit_normals]
        if not (
            witness_square > 0
            and pairings[target] == 0
            and all(pairing > 0 for index, pairing in enumerate(pairings) if index != target)
        ):
            raise ArithmeticError("an extreme inequality lacks a strict timelike facet witness")
        facets.append(target)
        facet_witnesses.append(
            {
                "orbit_index": target + 1,
                "witness_square": _rational_text(witness_square),
                "minimum_strict_pairing": _rational_text(
                    min(pairing for index, pairing in enumerate(pairings) if index != target)
                ),
            }
        )

    if len(facets) + len(redundant) != len(wall_orbits):
        raise ArithmeticError("the exact wall reduction is not exhaustive")

    labels: Counter[str] = Counter()
    bad_pairs: list[dict[str, Any]] = []
    for left_position, left in enumerate(facets):
        for right in facets[left_position + 1 :]:
            left_norm = sp.Rational(orbit_gram[left, left])
            right_norm = sp.Rational(orbit_gram[right, right])
            product = sp.Rational(orbit_gram[left, right])
            ratio = sp.factor(product**2 / (left_norm * right_norm))
            try:
                label = _coxeter_label(left_norm, right_norm, product)
            except ValueError as error:
                bad_pairs.append(
                    {
                        "facet_orbit_indices": [left + 1, right + 1],
                        "normalized_squared_product": _rational_text(sp.Rational(ratio)),
                        "reason": str(error),
                    }
                )
            else:
                labels[label] += 1

    dimension = fixed_rank - 1
    facet_gram = orbit_gram.extract(facets, facets)
    all_coxeter = not bad_pairs
    vinberg = None
    if all_coxeter and dimension >= 2:
        vinberg = _vinberg_payload(facet_gram, dimension)

    return {
        "id": f"subgroup_{subgroup_index:03d}",
        "group_element_indices": sorted(
            permutation_index[element] + 1 for element in subgroup
        ),
        "order": len(subgroup),
        "element_order_histogram": {
            str(order): count
            for order, count in sorted(
                Counter(_permutation_order(element) for element in subgroup).items()
            )
        },
        "wall_orbits": [[index + 1 for index in orbit] for orbit in wall_orbits],
        "wall_orbit_count": len(wall_orbits),
        "fixed_rank": fixed_rank,
        "hyperbolic_dimension": dimension,
        "orbit_gram_inertia": list(inertia),
        "facet_count": len(facets),
        "facet_orbit_indices": [index + 1 for index in facets],
        "redundant_orbit_count": len(redundant),
        "redundancy_certificates": redundant,
        "facet_witnesses": facet_witnesses,
        "coxeter_edge_counts": dict(sorted(labels.items())),
        "non_coxeter_pair_count": len(bad_pairs),
        "non_coxeter_pairs": bad_pairs,
        "all_facet_pairs_are_coxeter": all_coxeter,
        "facet_gram": _matrix_payload(facet_gram),
        "vinberg": vinberg,
        "checks": {
            "subgroup_action_is_exact": True,
            "orbit_normals_span_fixed_space": True,
            "fixed_space_is_lorentzian": True,
            "wall_reduction_has_exact_certificate_for_every_orbit": True,
            "every_retained_wall_has_strict_timelike_relative_interior": True,
        },
    }


@lru_cache(maxsize=4)
def build_artifacts(path: Path = SOURCE) -> tuple[dict[str, Any], dict[str, Any]]:
    source = _parse_source(path)
    gram: sp.Matrix = source["gram"]
    ample: sp.Matrix = source["ample"]
    wall_rows = tuple(sp.Matrix(1, gram.rows, row) for row in source["wall_vectors"])
    group = _derive_wall_permutations(source)
    geometric_group = _derive_wall_permutations(
        source, source["geometric_stabilizer"]
    )
    abstract_group = _weighted_wall_graph_automorphisms(
        gram, source["wall_vectors"]
    )
    if set(abstract_group) != set(group):
        raise ArithmeticError(
            "published OSXD0 is not the full weighted-wall symmetry group"
        )
    if not set(geometric_group) < set(group):
        raise ArithmeticError("AutXD0 is not the expected proper wall subgroup")
    subgroups = _all_subgroups(group)
    screens = [
        _screen_subgroup(source, group, subgroup, index)
        for index, subgroup in enumerate(subgroups)
    ]

    wall_norms = Counter(
        sp.Rational((wall * gram * wall.T)[0]) for wall in wall_rows
    )
    stock = {
        "schema": "apery_fermi_exact_chamber_stock_v1",
        "source": {
            "title": "The automorphism group of an Ap'ery-Fermi K3 surface: Computational data",
            "author": "Ichiro Shimada",
            "doi": "10.5281/zenodo.14842583",
            "record_url": "https://zenodo.org/records/14842583",
            "source_file": "AperyFermiCompData.txt",
            "source_archive_md5": "51403629e27a00607e4bed40ceeb1f39",
            "license": "CC-BY-4.0",
        },
        "lattice": {
            "description": "U + E8(-1) + E8(-1) + <-12>",
            "rank": gram.rows,
            "determinant": int(gram.det()),
            "inertia": list(exact_inertia(gram)),
            "gram": _matrix_payload(gram),
            "ample_class_a32": [_rational_text(value) for value in ample],
            "ample_square": _rational_text(sp.Rational((ample * gram * ample.T)[0])),
        },
        "walls": [
            {
                "index": index + 1,
                "vector": [_rational_text(value) for value in vector],
                "square": _rational_text(sp.Rational((wall_rows[index] * gram * wall_rows[index].T)[0])),
                "kind": source["wall_types"][index],
                "published_orbit_number": source["wall_orbit_numbers"][index],
            }
            for index, vector in enumerate(source["wall_vectors"])
        ],
        "wall_summary": {
            "count": len(wall_rows),
            "kind_histogram": dict(sorted(Counter(source["wall_types"]).items())),
            "square_histogram": {
                _rational_text(norm): count
                for norm, count in sorted(wall_norms.items())
            },
        },
        "weighted_wall_symmetry_group": {
            "order": len(group),
            "element_order_histogram": {
                str(order): count
                for order, count in sorted(
                    Counter(_permutation_order(element) for element in group).items()
                )
            },
            "wall_permutations": [
                [image + 1 for image in permutation] for permutation in group
            ],
            "all_matrices_fix_a32": True,
            "all_matrices_preserve_gram": True,
            "all_matrices_permute_walls": True,
            "equals_full_abstract_weighted_wall_graph_automorphism_group": True,
            "published_source_assignment": "OSXD0",
        },
        "geometric_chamber_subgroup": {
            "order": len(geometric_group),
            "index_in_weighted_wall_symmetry_group": len(group) // len(geometric_group),
            "published_source_assignment": "AutXD0",
            "used_as_a_search_gate": False,
        },
    }

    above_nine = [screen for screen in screens if screen["hyperbolic_dimension"] > 9]
    above_nine_nontrivial = [screen for screen in above_nine if screen["order"] > 1]
    coxeter_screens = [screen for screen in screens if screen["all_facet_pairs_are_coxeter"]]
    genuine_coxeter_screens = [
        screen for screen in coxeter_screens if screen["hyperbolic_dimension"] >= 2
    ]
    degenerate_one_dimensional_screens = [
        screen for screen in coxeter_screens if screen["hyperbolic_dimension"] == 1
    ]
    enumeration = {
        "schema": "apery_fermi_chamber_subgroup_fold_census_v1",
        "method": {
            "symmetry_scope": (
                "all subgroups of the full abstract weighted wall-graph symmetry "
                "group, equal to published OSXD0"
            ),
            "surface_automorphism_required": False,
            "normal_restriction": "sum wall normals over each subgroup orbit",
            "wall_reduction": "exact positive rational dependences or strict timelike facet witnesses",
            "coxeter_gate": "exact normalized squared products",
            "vinberg_gate": "run on every Coxeter survivor of dimension at least two",
        },
        "group_order": len(group),
        "subgroup_count": len(subgroups),
        "subgroup_order_histogram": {
            str(order): count
            for order, count in sorted(Counter(len(subgroup) for subgroup in subgroups).items())
        },
        "subgroups": screens,
        "summary": {
            "dimension_histogram": {
                str(dimension): count
                for dimension, count in sorted(
                    Counter(screen["hyperbolic_dimension"] for screen in screens).items(),
                    reverse=True,
                )
            },
            "nontrivial_dimension_above_nine_count": len(above_nine_nontrivial),
            "nontrivial_dimension_above_nine_coxeter_count": sum(
                screen["all_facet_pairs_are_coxeter"] for screen in above_nine_nontrivial
            ),
            "dimension_above_nine_nontrivial_bad_pair_counts": sorted(
                screen["non_coxeter_pair_count"] for screen in above_nine_nontrivial
            ),
            "coxeter_subgroup_count": len(coxeter_screens),
            "genuine_h2_through_h8_coxeter_count": len(genuine_coxeter_screens),
            "genuine_h2_through_h8_vinberg_finite_volume_count": sum(
                screen["vinberg"] is not None
                and screen["vinberg"]["finite_volume"] is True
                for screen in genuine_coxeter_screens
            ),
            "degenerate_h1_coxeter_compatible_count": len(
                degenerate_one_dimensional_screens
            ),
            "vinberg_was_not_run_in_h1": True,
            "coxeter_dimensions": sorted(
                Counter(screen["hyperbolic_dimension"] for screen in coxeter_screens).elements(),
                reverse=True,
            ),
            "all_exact_reduction_checks_pass": all(
                all(screen["checks"].values()) for screen in screens
            ),
        },
    }
    return stock, enumeration


def _write_or_check(path: Path, payload: dict[str, Any], check: bool) -> None:
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if check:
        if not path.is_file() or path.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rendered, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--stock-output", type=Path, default=STOCK_OUTPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    stock, enumeration = build_artifacts(arguments.source)
    _write_or_check(arguments.stock_output, stock, arguments.check_archive)
    _write_or_check(arguments.output, enumeration, arguments.check_archive)
    print(json.dumps(enumeration["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
