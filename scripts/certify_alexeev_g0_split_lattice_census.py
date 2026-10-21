#!/usr/bin/env python3
"""Certify the split-root obstruction for Alexeev's ``r+a=22`` chambers.

The input is only the pinned ancillary file of complete mixed-reflection
wall Gram matrices.  For every long wall ``v_i`` (square -4, divisibility
two), an exact Coble-nef interpretation would require two exchanged K3
roots ``c_i^+`` and ``c_i^-``.  Put ``d_i=c_i^+-c_i^-`` and
``H=-Gram(d_i)/2``.  If ``n_ij=(v_i,v_j)/2``, nonnegative intersections of
the lifted curves force

    H_ii = 2,
    H_ij in {-n_ij, -n_ij+2, ..., n_ij}.

Moreover H must be positive semidefinite.  This script solves that finite
completion problem exactly for all nine raw chambers, checks the saturated
radical parity condition, and computes the even lattice generated after
adjoining ``(v_i+d_i)/2``.

This is a decorated-carrier obstruction.  It does not assert that an
undecorated arithmetic Coxeter chamber fails to exist.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path("/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt")
SOURCE_SHA256 = "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
OUTPUT = ROOT / "data/enumerations/alexeev_g0_split_lattice_census.json"
SCHEMA = "coxeter.alexeev_g0_split_lattice_census.v1"

KEYS = (
    (11, 11, 1),
    (12, 10, 1),
    (13, 9, 1),
    (14, 8, 0),
    (14, 8, 1),
    (15, 7, 1),
    (16, 6, 1),
    (18, 4, 0),
    (20, 2, 1),
)

EXPECTED_WALL_AND_LONG_COUNTS = {
    (11, 11, 1): (12, 10),
    (12, 10, 1): (14, 11),
    (13, 9, 1): (19, 10),
    (14, 8, 0): (17, 10),
    (14, 8, 1): (24, 10),
    (15, 7, 1): (60, 10),
    (16, 6, 1): (118, 25),
    (18, 4, 0): (48, 24),
    (20, 2, 1): (50, 5),
}

FOREST_INFEASIBLE = {
    (11, 11, 1): ((9, 1, 0), -1),
    (12, 10, 1): ((10, 1, 0), -8),
    (13, 9, 1): ((9, 1, 0), -16),
    (14, 8, 0): ((9, 1, 0), -16),
    (14, 8, 1): ((9, 1, 0), -64),
}

MINIMAL_BASES = {
    (16, 6, 1): (0, 1, 2, 3, 4, 5),
    (18, 4, 0): (0, 1, 2, 5),
    (20, 2, 1): (0,),
}

ROOT_LATTICE_NAMES = {
    (15, 7, 1): "E6",
    (16, 6, 1): "A5 + A1",
    (18, 4, 0): "A2 + A2",
    (20, 2, 1): "A1",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _gf2_rank(rows: Iterable[Iterable[int]]) -> int:
    packed = [
        sum((int(value) & 1) << column for column, value in enumerate(row))
        for row in rows
    ]
    rank = 0
    width = max((value.bit_length() for value in packed), default=0)
    for column in range(width):
        pivot = next(
            (
                index
                for index in range(rank, len(packed))
                if (packed[index] >> column) & 1
            ),
            None,
        )
        if pivot is None:
            continue
        packed[rank], packed[pivot] = packed[pivot], packed[rank]
        for index in range(len(packed)):
            if index != rank and ((packed[index] >> column) & 1):
                packed[index] ^= packed[rank]
        rank += 1
    return rank


def _ambient_reconstruction(
    wall_gram: sp.Matrix, rank: int
) -> tuple[sp.Matrix, sp.Matrix]:
    ambient = wall_gram[:rank, :rank]
    walls = wall_gram[:, :rank] * ambient.inv()
    if not (
        ambient.det() != 0
        and walls[:rank, :] == sp.eye(rank)
        and walls * ambient * walls.T == wall_gram
        and all(value.q == 1 for value in walls)
    ):
        raise ArithmeticError("the leading wall block is not an integral basis")
    return ambient, walls


def _allowed_h(pairing: int) -> tuple[int, ...]:
    if pairing < 0 or pairing % 2:
        raise ArithmeticError("long-wall products must be nonnegative and even")
    n = pairing // 2
    return tuple(range(-n, n + 1, 2))


def _long_data(
    source: str, key: tuple[int, int, int]
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, list[int]]:
    wall_gram = _extract_matrix(source, key)
    rank = key[0]
    ambient, walls = _ambient_reconstruction(wall_gram, rank)
    long_indices = [
        index for index in range(wall_gram.rows) if wall_gram[index, index] == -4
    ]
    long_gram = wall_gram.extract(long_indices, long_indices)
    long_vectors = walls[long_indices, :]
    divisibilities = []
    for row in range(long_vectors.rows):
        pairings = long_vectors.row(row) * ambient
        divisibilities.append(
            int(sp.gcd_list([abs(int(value)) for value in pairings]))
        )
    if not (
        wall_gram.rows == EXPECTED_WALL_AND_LONG_COUNTS[key][0]
        and len(long_indices) == EXPECTED_WALL_AND_LONG_COUNTS[key][1]
        and set(divisibilities) == {2}
        and {
            int(long_gram[i, j])
            for i in range(long_gram.rows)
            for j in range(i)
        }
        <= {0, 2, 4}
    ):
        raise ArithmeticError(f"raw long-wall data changed for {key}")
    return ambient, long_vectors, long_gram, long_indices


def _pair_histogram(long_gram: sp.Matrix) -> dict[str, int]:
    return {
        str(value): count
        for value, count in sorted(
            Counter(
                int(long_gram[i, j])
                for i in range(long_gram.rows)
                for j in range(i)
            ).items()
        )
    }


def _forest_components_and_edges(long_gram: sp.Matrix) -> tuple[int, list[tuple[int, int]]]:
    size = long_gram.rows
    edges = [
        (left, right)
        for left in range(size)
        for right in range(left)
        if long_gram[left, right] == 2
    ]
    parent = list(range(size))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    for left, right in edges:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            raise ArithmeticError("the unit-intersection graph is not a forest")
        parent[root_left] = root_right
    components = len({find(index) for index in range(size)})
    return components, edges


def _forest_no_go(long_gram: sp.Matrix) -> dict[str, Any]:
    components, edges = _forest_components_and_edges(long_gram)
    size = long_gram.rows
    representative = 2 * sp.eye(size)
    for left, right in edges:
        representative[left, right] = representative[right, left] = -1
    inertia = exact_inertia(representative)
    return {
        "status": "no_positive_semidefinite_completion",
        "reason": (
            "the n=1 graph is a forest, so every lift signing is diagonally "
            "switching-congruent to the displayed indefinite representative"
        ),
        "unit_intersection_graph": {
            "vertices": size,
            "edges": len(edges),
            "components": components,
            "cycle_rank": len(edges) - size + components,
            "forest": True,
        },
        "representative_H": [[int(value) for value in row] for row in representative.tolist()],
        "representative_inertia_positive_negative_zero": list(inertia),
        "representative_determinant": int(representative.det()),
        "minimum_anti_invariant_rank": None,
    }


def _spanning_forest(
    size: int, edges: list[tuple[int, int]]
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    parent = list(range(size))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    tree: list[tuple[int, int]] = []
    chords: list[tuple[int, int]] = []
    for left, right in sorted(edges):
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            chords.append((left, right))
        else:
            parent[root_left] = root_right
            tree.append((left, right))
    return tree, chords


def _type_vi_completion(long_gram: sp.Matrix) -> tuple[sp.Matrix, dict[str, Any]]:
    size = long_gram.rows
    edges = [
        (left, right)
        for left in range(size)
        for right in range(left)
        if long_gram[left, right] == 2
    ]
    tree, chords = _spanning_forest(size, edges)
    distribution: Counter[tuple[int, int, int]] = Counter()
    positive_semidefinite: list[tuple[int, sp.Matrix]] = []
    for mask in range(1 << len(chords)):
        completion = 2 * sp.eye(size)
        for left, right in edges:
            completion[left, right] = completion[right, left] = -1
        for bit, (left, right) in enumerate(chords):
            if (mask >> bit) & 1:
                completion[left, right] = completion[right, left] = 1
        inertia = exact_inertia(completion)
        distribution[inertia] += 1
        if inertia[1] == 0:
            positive_semidefinite.append((mask, completion))
    if not (
        len(chords) == 6
        and len(positive_semidefinite) == 1
        and exact_inertia(positive_semidefinite[0][1]) == (6, 0, 4)
    ):
        raise ArithmeticError("the Type-VI switching census changed")
    mask, completion = positive_semidefinite[0]
    return completion, {
        "method": "all switching classes after fixing a spanning tree",
        "unit_intersection_edges": len(edges),
        "spanning_tree_edges": [list(edge) for edge in tree],
        "cycle_edges": [list(edge) for edge in chords],
        "switching_classes_tested": 1 << len(chords),
        "inertia_distribution": {
            f"{positive},{negative},{zero}": count
            for (positive, negative, zero), count in sorted(distribution.items())
        },
        "positive_semidefinite_switching_classes": 1,
        "canonical_cycle_mask": mask,
    }


def _minimal_root_completion(
    long_gram: sp.Matrix, base_indices: tuple[int, ...]
) -> tuple[sp.Matrix, dict[str, Any]]:
    """Exhaust rank-minimal completions using a forced positive basis block."""

    size = long_gram.rows
    base_size = len(base_indices)
    base = 2 * sp.eye(base_size)
    base_edges: list[tuple[int, int]] = []
    for left in range(base_size):
        for right in range(left):
            allowed = _allowed_h(
                int(long_gram[base_indices[left], base_indices[right]])
            )
            if len(allowed) == 1:
                base[left, right] = base[right, left] = allowed[0]
            elif allowed == (-1, 1):
                base[left, right] = base[right, left] = -1
                base_edges.append((left, right))
            else:
                raise ArithmeticError("the advertised lower-bound block has a long choice")

    # Every signing on a forest is obtained from this all-minus representative
    # by a diagonal +/-1 congruence.  Thus the positive-definite block, and
    # hence the rank lower bound, does not depend on an unexamined lift choice.
    parent = list(range(base_size))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    for left, right in base_edges:
        root_left, root_right = find(left), find(right)
        if root_left == root_right:
            raise ArithmeticError("the rank-witness signing support is not a forest")
        parent[root_left] = root_right
    base_components = len({find(index) for index in range(base_size)})
    if len(base_edges) != base_size - base_components:
        raise ArithmeticError("the rank-witness forest count changed")
    if exact_inertia(base) != (base_size, 0, 0):
        raise ArithmeticError("the advertised lower-bound block is not positive definite")

    inverse = base.inv()
    coordinates: dict[int, tuple[sp.Rational, ...]] = {
        index: tuple(sp.eye(base_size).row(position))
        for position, index in enumerate(base_indices)
    }
    candidate_counts: dict[int, int] = {}
    for index in range(size):
        if index in coordinates:
            continue
        pairing_choices = [
            _allowed_h(int(long_gram[index, base_index]))
            for base_index in base_indices
        ]
        candidates = []
        for pairings in itertools.product(*pairing_choices):
            column = sp.Matrix(pairings)
            if (column.T * inverse * column)[0] == 2:
                candidates.append(tuple(inverse * column))
        candidates = sorted(set(candidates))
        if not (
            len(candidates) == 2
            and candidates[1] == tuple(-value for value in candidates[0])
            and all(value.q == 1 for value in candidates[0])
        ):
            raise ArithmeticError(f"unexpected root candidates at long wall {index}")
        candidate_counts[index] = len(candidates)
        coordinates[index] = candidates[0]

    coordinate_matrix = sp.Matrix(
        [[coordinates[index][column] for column in range(base_size)] for index in range(size)]
    )
    completion = coordinate_matrix * base * coordinate_matrix.T
    for left in range(size):
        for right in range(left):
            if int(completion[left, right]) not in _allowed_h(
                int(long_gram[left, right])
            ):
                raise ArithmeticError("canonical root coordinates are incompatible")
            allowed_absolute = {
                abs(value) for value in _allowed_h(int(long_gram[left, right]))
            }
            if abs(int(completion[left, right])) not in allowed_absolute:
                raise ArithmeticError("an independent sheet swap can violate a pairing")
    if exact_inertia(completion) != (base_size, 0, size - base_size):
        raise ArithmeticError("constructed completion has the wrong rank")
    return completion, {
        "method": (
            "a forced positive-definite principal block gives the lower bound; "
            "all allowed pairing columns of norm two in its dual are enumerated"
        ),
        "lower_bound_basis_long_indices_zero_based": list(base_indices),
        "lower_bound_basis_gram": [[int(value) for value in row] for row in base.tolist()],
        "lower_bound_basis_determinant": int(base.det()),
        "lower_bound_sign_support": {
            "edges": [list(edge) for edge in base_edges],
            "components": base_components,
            "forest": True,
            "all_signings_switch_to_the_displayed_basis": True,
        },
        "candidate_count_for_each_nonbasis_wall": {
            str(index): count for index, count in sorted(candidate_counts.items())
        },
        "all_candidate_sets_are_one_opposite_pair": True,
        "opposite_choice_is_exactly_the_independent_sheet_swap": True,
        "all_opposite_choices_remain_pairwise_compatible": True,
        "exhaustion_reason": (
            "at the lower-bound rank the positive principal block spans the whole "
            "difference space; pairings with that basis determine every other vector, "
            "and the exhaustive norm-two search returns only one opposite pair"
        ),
        "rank_minimal_switching_orbits": 1,
    }


def _integral_root_basis(
    completion: sp.Matrix,
) -> tuple[tuple[int, ...], sp.Matrix, sp.Matrix]:
    rank = completion.rank()
    for indices in itertools.combinations(range(completion.rows), rank):
        basis = completion.extract(indices, indices)
        if basis.det() == 0:
            continue
        coordinates = completion[:, indices] * basis.inv()
        if all(value.q == 1 for value in coordinates):
            return indices, basis, coordinates
    raise ArithmeticError("the completed difference roots have no integral root basis")


def _completion_arithmetic(
    *,
    key: tuple[int, int, int],
    ambient: sp.Matrix,
    long_vectors: sp.Matrix,
    long_gram: sp.Matrix,
    completion: sp.Matrix,
) -> dict[str, Any]:
    indices, root_gram, coordinates = _integral_root_basis(completion)
    root_rank = root_gram.rows
    radical_relations = []
    parity_pass = True
    for index in range(completion.rows):
        if index in indices:
            continue
        relation = [0] * completion.rows
        relation[index] = 1
        for position, basis_index in enumerate(indices):
            relation[basis_index] -= int(coordinates[index, position])
        invariant_sum = sp.Matrix(1, len(relation), relation) * long_vectors
        parity = [int(value) % 2 for value in invariant_sum]
        parity_pass = parity_pass and not any(parity)
        radical_relations.append(
            {
                "nonbasis_long_index": index,
                "relation": relation,
                "invariant_sum_mod_2": parity,
            }
        )

    lift_intersections = Counter()
    integral_lift_gram = True
    nonnegative_lift_gram = True
    for left in range(completion.rows):
        for right in range(left):
            n = int(long_gram[left, right]) // 2
            h = int(completion[left, right])
            a = (n - h) // 2
            b = (n + h) // 2
            if not (n - h) % 2 == 0 or not (n + h) % 2 == 0:
                integral_lift_gram = False
            if min(a, b) < 0:
                nonnegative_lift_gram = False
            lift_intersections[(a, b)] += 1

    glue_rows = [
        [int(value) for value in long_vectors.row(index)]
        + [int(value) for value in coordinates.row(index)]
        for index in range(completion.rows)
    ]
    glue_rank = _gf2_rank(glue_rows)
    determinant_m = abs(int(ambient.det()))
    determinant_h = abs(int(root_gram.det()))
    determinant_a = (2**root_rank) * determinant_h
    determinant_k = sp.Rational(
        determinant_m * determinant_a, 2 ** (2 * glue_rank)
    )
    if determinant_k.q != 1:
        raise ArithmeticError("the split-component overlattice has nonintegral determinant")

    integral_pairings = True
    for left in range(completion.rows):
        for right in range(completion.rows):
            numerator = int(long_gram[left, right]) - 2 * int(completion[left, right])
            if numerator % 4:
                integral_pairings = False

    complex_budget = 20 - key[0]
    total_rank = key[0] + root_rank
    checks = {
        "H_positive_semidefinite": exact_inertia(completion)[1] == 0,
        "H_has_claimed_rank": completion.rank() == root_rank,
        "root_coordinates_integral": all(value.q == 1 for value in coordinates),
        "saturated_radical_parity": parity_pass,
        "lift_intersections_integral": integral_lift_gram,
        "lift_intersections_nonnegative": nonnegative_lift_gram,
        "split_component_overlattice_integral": integral_pairings,
        "ambient_determinant_is_2_power_a": determinant_m == 2 ** key[1],
        "glue_rank_equals_difference_rank": glue_rank == root_rank,
        "complex_rank_budget_fails": root_rank > complex_budget,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"completion arithmetic failed for {key}: {checks}")

    return {
        "status": "positive_semidefinite_completion",
        "minimum_anti_invariant_rank": root_rank,
        "complex_anti_invariant_algebraic_rank_budget": complex_budget,
        "complex_budget_exceeded": True,
        "difference_root_lattice_H": ROOT_LATTICE_NAMES[key],
        "difference_lattice_Gram": f"{ROOT_LATTICE_NAMES[key]}(-2)",
        "canonical_H": [[int(value) for value in row] for row in completion.tolist()],
        "canonical_H_inertia_positive_negative_zero": list(exact_inertia(completion)),
        "integral_root_basis_long_indices_zero_based": list(indices),
        "integral_root_basis_H_gram": [
            [int(value) for value in row] for row in root_gram.tolist()
        ],
        "integral_root_coordinates": [
            [int(value) for value in row] for row in coordinates.tolist()
        ],
        "root_lattice_determinant": determinant_h,
        "radical": {
            "rank": completion.rows - root_rank,
            "integral_basis_relations": radical_relations,
            "every_relation_has_invariant_sum_in_2M": parity_pass,
        },
        "canonical_lift_intersection_pair_histogram_a_b": {
            f"{a},{b}": count
            for (a, b), count in sorted(lift_intersections.items())
        },
        "split_component_overlattice": {
            "definition": "K=<M,A,(v_i+d_i)/2>",
            "rank": total_rank,
            "signature_positive_negative": [1, total_rank - 1],
            "absolute_determinant_M": determinant_m,
            "absolute_determinant_A": determinant_a,
            "F2_glue_rank": glue_rank,
            "index_over_M_plus_A": 2**glue_rank,
            "absolute_determinant_K": int(determinant_k),
            "even_and_integral": integral_pairings,
        },
        "checks": checks,
    }


def _row(source: str, key: tuple[int, int, int]) -> dict[str, Any]:
    ambient, long_vectors, long_gram, long_indices = _long_data(source, key)
    common: dict[str, Any] = {
        "key": list(key),
        "boundary_count_from_fixed_locus": 1 + (key[0] - key[1]) // 2,
        "hyperbolic_dimension": key[0] - 1,
        "complete_mixed_reflection_wall_count": EXPECTED_WALL_AND_LONG_COUNTS[key][0],
        "long_wall_count": len(long_indices),
        "long_wall_indices_zero_based": long_indices,
        "long_wall_pairing_histogram": _pair_histogram(long_gram),
        "all_long_walls_have_divisibility_two": True,
        "decorated_requirement": (
            "every square-minus-four facet is a downstairs square-minus-two curve "
            "and therefore splits into two exchanged K3 roots"
        ),
    }
    if key in FOREST_INFEASIBLE:
        obstruction = _forest_no_go(long_gram)
        expected_inertia, expected_determinant = FOREST_INFEASIBLE[key]
        if not (
            tuple(obstruction["representative_inertia_positive_negative_zero"])
            == expected_inertia
            and obstruction["representative_determinant"] == expected_determinant
        ):
            raise ArithmeticError(f"forest obstruction changed for {key}")
        common.update(obstruction)
        common["complex_anti_invariant_algebraic_rank_budget"] = 20 - key[0]
        common["complex_carrier_excluded"] = True
        common["all_characteristic_split_system_excluded"] = True
        return common

    if key == (15, 7, 1):
        completion, enumeration = _type_vi_completion(long_gram)
    else:
        completion, enumeration = _minimal_root_completion(
            long_gram, MINIMAL_BASES[key]
        )
    arithmetic = _completion_arithmetic(
        key=key,
        ambient=ambient,
        long_vectors=long_vectors,
        long_gram=long_gram,
        completion=completion,
    )
    common.update(arithmetic)
    common["rank_minimality_certificate"] = enumeration
    common["complex_carrier_excluded"] = True
    common["all_characteristic_split_system_excluded"] = key == (16, 6, 1)
    return common


def build(source_path: Path = SOURCE) -> dict[str, Any]:
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    digest = _sha256(source_path)
    if digest != SOURCE_SHA256:
        raise ArithmeticError("Alexeev ancillary source hash changed")
    source = source_path.read_text(encoding="utf-8")
    rows = [_row(source, key) for key in KEYS]
    by_key = {tuple(row["key"]): row for row in rows}

    mii = by_key[(18, 4, 0)]["split_component_overlattice"]
    rank16 = by_key[(16, 6, 1)]["split_component_overlattice"]
    global_checks = {
        "all_nine_rows_present": len(rows) == 9,
        "all_nine_complex_exact_carriers_excluded": all(
            row["complex_carrier_excluded"] for row in rows
        ),
        "five_rows_have_no_numerical_split_system": sum(
            row["status"] == "no_positive_semidefinite_completion" for row in rows
        )
        == 5,
        "type_VI_minimum_rank_six": by_key[(15, 7, 1)][
            "minimum_anti_invariant_rank"
        ]
        == 6,
        "type_MII_minimum_rank_four": by_key[(18, 4, 0)][
            "minimum_anti_invariant_rank"
        ]
        == 4,
        "type_MII_rank22_determinant9": mii["rank"] == 22
        and mii["absolute_determinant_K"] == 9,
        "rank16_rank22_determinant12": rank16["rank"] == 22
        and rank16["absolute_determinant_K"] == 12,
    }
    if not all(global_checks.values()):
        raise ArithmeticError(f"global split census checks failed: {global_checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "POLYTOPE_TO_SURFACE_necessary_split_root_criterion",
        "raw_source": {
            "path": str(source_path),
            "sha256": digest,
            "paper": (
                "Valery Alexeev, Reflective hyperbolic 2-elementary lattices, "
                "K3 surfaces and hyperkahler manifolds, arXiv:2209.09110v4"
            ),
            "input_fields": (
                "only the nine complete wall Gram matrices with r+a=22; "
                "no Coble type labels or surface curve lists"
            ),
        },
        "decorated_polytope_class": {
            "definition": (
                "the nine pinned complete mixed-reflection chambers of even "
                "2-elementary lattices (r,a,delta) with r+a=22, equipped with "
                "the canonical Coble quotient norm dictionary: long (-4), "
                "divisibility-two facets are required to descend from (-2)-curves"
            ),
            "not_assumed": [
                "a pre-existing Coble surface",
                "a Kondo type or template",
                "effectivity of any wall",
                "that short (-2) walls split",
            ],
            "short_wall_warning": (
                "short square-minus-two walls encode boundary halves or downstairs "
                "square-minus-one curves; stable/tangent alternatives are not put "
                "through the forced split-root test"
            ),
        },
        "necessary_criterion": {
            "formula": (
                "for n_ij=(v_i.v_j)/2 and a_ij,b_ij>=0 with a_ij+b_ij=n_ij, "
                "H_ii=2 and H_ij=n_ij-2a_ij must be positive semidefinite"
            ),
            "saturated_radical_condition": (
                "for every integral z in ker(H), sum z_i v_i must lie in 2M"
            ),
            "complex_budget": "rank(H) <= 20-r",
        },
        "rows": rows,
        "theorems": {
            "class_wide_complex_no_go": (
                "none of the nine complete mixed-reflection chambers can be the "
                "exact full nef chamber of a complex strong-Coble quotient under "
                "the stated norm decoration"
            ),
            "five_absolute_numerical_no_gos": (
                "the rows (11,11,1), (12,10,1), (13,9,1), (14,8,0), "
                "and (14,8,1) admit no nonnegative split-intersection Gram at all"
            ),
            "rank16_all_characteristic_no_go": (
                "the (16,6,1) split lattice has rank 22 and determinant 12; "
                "every finite-index even rank-22 overclosure has determinant 12 "
                "or 3, neither the p^(2 sigma) determinant of a supersingular K3"
            ),
            "type_MII_host_forced": (
                "the (18,4,0) split lattice is rank 22, even, signature (1,21), "
                "and determinant 9; its only possible proper even overclosure would "
                "be even unimodular of signature (1,21), impossible modulo 8.  Thus "
                "any K3 carrier is supersingular in characteristic 3 with Artin invariant 1"
            ),
            "type_VI_calibration": (
                "the (15,7,1) split differences span E6(-2) of rank 6, so the "
                "generated algebraic lattice has rank 21; a carrier cannot be complex "
                "and needs at least one additional positive-characteristic Picard class"
            ),
        },
        "scope": {
            "proved": (
                "a non-tautological necessary carrier criterion and exact no-go/host "
                "consequences on an independently defined nine-element raw-P class"
            ),
            "not_proved": [
                "nonexistence of any of the nine abstract arithmetic polytopes",
                "nonexistence of carriers of another algebraic species",
                "sufficiency of a passing split lattice for a K3 involution",
                "effectivity, ramification-frame uniqueness, or exact-nef return",
                "a contradiction with the distinct 22/30-facet complex Type-I/II nef chambers",
            ],
            "type_I_II_noncontradiction": (
                "the complex Type-I/II Coble nef chambers are different subgroup/full-nef "
                "wall systems (for example H10/22 and H10/30), not Alexeev's complete "
                "12-wall mixed-reflection chamber for (11,11,1)"
            ),
        },
        "checks": global_checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build(arguments.source)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "all_checks_pass": result["all_checks_pass"],
                "complex_no_go_rows": len(result["rows"]),
                "minimum_ranks": {
                    "15,7,1": result["rows"][5]["minimum_anti_invariant_rank"],
                    "16,6,1": result["rows"][6]["minimum_anti_invariant_rank"],
                    "18,4,0": result["rows"][7]["minimum_anti_invariant_rank"],
                    "20,2,1": result["rows"][8]["minimum_anti_invariant_rank"],
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
