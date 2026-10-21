#!/usr/bin/env python3
"""Freeze and diagnose the exact bound-8 OG10 golden Vinberg prefix.

The output is a bounded-search certificate.  Its hyperideal ray proves that
the recorded 171-wall outer cone is not finite-volume; roots at later heights
may still cut that ray, so this is not a nonreflectivity theorem.
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from typing import Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, PHI, ZERO
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.search.og10_golden_vinberg import (
    D,
    og10_field_gram,
    run_og10_golden_vinberg_search,
)
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray


IDENTITY_BOUND = 8
HYPERIDEAL_WITNESS = (0, 1, 2, 3, 4, 5, 6, 7, 11, 21)


def _vector(vector: Iterable[Golden]) -> list[dict[str, str]]:
    return [entry.to_json() for entry in vector]


def _rank(rows: Sequence[Sequence[Golden]]) -> int:
    matrix = [[entry for entry in row] for row in rows]
    if not matrix:
        return 0
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [entry / pivot_value for entry in matrix[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not matrix[row][column]:
                continue
            multiplier = matrix[row][column]
            matrix[row] = [
                matrix[row][index] - multiplier * matrix[pivot_row][index]
                for index in range(column_count)
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def _spanning_root_indices(roots: Sequence[Sequence[Golden]]) -> tuple[int, ...]:
    selected: list[int] = []
    current_rank = 0
    for index, root in enumerate(roots):
        candidate_rank = _rank([roots[value] for value in (*selected, index)])
        if candidate_rank > current_rank:
            selected.append(index)
            current_rank = candidate_rank
        if current_rank == 11:
            break
    return tuple(selected)


def _relation(inner: Golden) -> str:
    if inner == 0:
        return "m2"
    if inner == -1:
        return "m3"
    if inner == -PHI:
        return "m5"
    if (inner + 2).sign(0) <= 0:
        return "ultraparallel"
    raise ArithmeticError(f"unrecognized acute norm-two product {inner}")


def _digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def build() -> dict[str, object]:
    search = run_og10_golden_vinberg_search(IDENTITY_BOUND)
    gram = og10_field_gram()
    roots_json = [_vector(root) for root in search.roots]
    relation_counts = Counter(
        _relation(_pair(left, gram, right))
        for index, left in enumerate(search.roots)
        for right in search.roots[index + 1 :]
    )
    spanning = _spanning_root_indices(search.roots)
    witness = verify_hyperideal_extreme_ray(
        search.roots,
        gram,
        HYPERIDEAL_WITNESS,
        dimension=10,
    )

    level_records = [
        {
            "control_coordinate": level.control_coordinate.to_json(),
            "height": level.height.to_json(),
            "target": level.target.to_json(),
            "root_count": len(level.roots),
            "visited_by_depth": list(level.visited_by_depth),
            "parity_pruned": level.parity_pruned,
            "schur_pruned_by_depth": list(level.schur_pruned_by_depth),
            "wall_pruned_by_depth": list(level.wall_pruned_by_depth),
            "leaf_norm_mismatches": level.leaf_norm_mismatches,
        }
        for level in search.levels
    ]
    checks = {
        "control_candidates_complete_through_bound": (
            [str(value) for value in search.candidate_coordinates]
            == [
                "-1 + 1*phi",
                "1",
                "1*phi",
                "1 + 1*phi",
                "2*phi",
                "2 + 1*phi",
                "1 + 2*phi",
                "2 + 2*phi",
                "1 + 3*phi",
                "2 + 3*phi",
                "1 + 4*phi",
                "3 + 3*phi",
            ]
        ),
        "height_zero_simple_system_A1_E8_A1": len(search.cone_roots) == 10,
        "level_root_counts": [len(level.roots) for level in search.levels]
        == [1, 2, 2, 2, 3, 0, 16, 31, 16, 65, 3, 20],
        "facet_count_171": len(search.roots) == 171,
        "roots_distinct": len(set(search.roots)) == len(search.roots),
        "all_integral_norm_two": all(
            all(entry.is_integral for entry in root) and _pair(root, gram, root) == 2
            for root in search.roots
        ),
        "root_span_rank_11": len(spanning) == 11,
        "all_pairs_are_Coxeter_or_ultraparallel": sum(relation_counts.values())
        == len(search.roots) * (len(search.roots) - 1) // 2,
        "hyperideal_witness_active": witness.active_facets == HYPERIDEAL_WITNESS,
        "hyperideal_witness_spacelike": witness.norm.sign(0) > 0,
        "hyperideal_witness_satisfies_all_walls": all(
            value.sign(0) <= 0 for value in witness.wall_pairings
        ),
    }
    return {
        "schema_version": "0.1.0",
        "certificate_kind": "bounded_exact_vinberg_prefix_with_hyperideal_witness",
        "candidate_id": "og10_golden_h10_carrier",
        "field_lattice": "[[-2*phi,1],[1,2]] + E8 + <2>",
        "hyperbolic_dimension": 10,
        "control_vector": _vector(search.control_vector),
        "control_norm": (-2 * D).to_json(),
        "identity_control_bound": IDENTITY_BOUND,
        "candidate_control_coordinates": [
            value.to_json() for value in search.candidate_coordinates
        ],
        "height_zero": {
            "root_system": "A1+E8+A1",
            "full_root_count": 244,
            "simple_root_count": len(search.cone_roots),
        },
        "levels": level_records,
        "visited_node_count": sum(
            sum(level.visited_by_depth) for level in search.levels
        ),
        "facet_count": len(search.roots),
        "spanning_root_indices": list(spanning),
        "pair_relation_counts": dict(sorted(relation_counts.items())),
        "roots_sha256": _digest(roots_json),
        "simple_roots": roots_json,
        "hyperideal_witness": {
            "active_facets": list(witness.active_facets),
            "ray": _vector(witness.ray),
            "norm": witness.norm.to_json(),
            "wall_pairings_sha256": _digest(
                [value.to_json() for value in witness.wall_pairings]
            ),
        },
        "checks": checks,
        "verified": all(checks.values()),
        "conclusion": (
            "the_exact_bound_8_outer_cone_is_not_finite_volume;_later_"
            "norm_two_roots_or_other_OG10_wall_types_may_cut_the_witness"
        ),
    }


def main() -> None:
    print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
