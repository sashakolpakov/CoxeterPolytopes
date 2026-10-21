#!/usr/bin/env python3
"""Build arrays-only free/satellite fixtures for PB-Atlas numerical Stage B1."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/fixtures/pb_atlas_numerical_stage_b1_programs.json"
SCHEMA = "coxeter.pb_atlas.numerical_cluster_stage_b1_fixtures.v1"


def _hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _pair(left: Sequence[int], right: Sequence[int], kind: str, e: int | None) -> int:
    model_rank = 1 if kind == "P2" else 2
    if kind == "P2":
        answer = left[0] * right[0]
    else:
        assert e is not None
        answer = -e * left[0] * right[0] + left[0] * right[1] + left[1] * right[0]
    return answer - sum(
        left[index] * right[index] for index in range(model_rank, len(left))
    )


def _program(
    *,
    endpoint_kind: str,
    hirzebruch_e: int | None,
    parents_chronological: Sequence[Sequence[int]],
    extra_rows_chronological: Sequence[Sequence[int]],
    boundary_extra_positions_and_weights: Sequence[tuple[int, int]],
) -> dict[str, Any]:
    """Construct a self-contained Stage-A-shaped endpoint program.

    Extra rows are supplied in chronological total-exceptional coordinates;
    the frozen Stage-A packet stores total exceptionals in blowdown order, so
    their exceptional tails are reversed here.
    """

    center_count = len(parents_chronological)
    model_rank = 1 if endpoint_kind == "P2" else 2
    if endpoint_kind == "P2":
        e = None
        model_gram = [[1]]
        canonical_model = [-3]
        nef_rays = [[1]]
    else:
        if hirzebruch_e is None or hirzebruch_e < 0:
            raise ValueError("F_e fixture requires e >= 0")
        e = int(hirzebruch_e)
        model_gram = [[-e, 1], [1, 0]]
        canonical_model = [-2, -(e + 2)]
        nef_rays = [[0, 1], [1, e]]

    proximity = [[0] * center_count for _ in range(center_count)]
    for child, parent_list in enumerate(parents_chronological):
        for parent in parent_list:
            proximity[child][int(parent)] = 1
    strict_chronological = [
        [int(row == column) - proximity[column][row] for column in range(center_count)]
        for row in range(center_count)
    ]
    strict_blowdown = [
        [
            strict_chronological[center_count - 1 - row][center_count - 1 - column]
            for column in range(center_count)
        ]
        for row in range(center_count)
    ]
    contracted_rows = [[0] * model_rank + row for row in strict_blowdown]
    extra_rows = [
        [*map(int, row[:model_rank]), *reversed([int(value) for value in row[model_rank:]])]
        for row in extra_rows_chronological
    ]
    coordinates = [*contracted_rows, *extra_rows]
    contracted = list(range(center_count))
    strict_gram = [
        [_pair(coordinates[left], coordinates[right], endpoint_kind, e) for right in contracted]
        for left in contracted
    ]
    size = model_rank + center_count
    basis_gram = [[0] * size for _ in range(size)]
    for left, row in enumerate(model_gram):
        for right, value in enumerate(row):
            basis_gram[left][right] = value
    for index in range(center_count):
        basis_gram[model_rank + index][model_rank + index] = -1
    standard = [[int(row == column) for row in range(size)] for column in range(size)]
    boundary = []
    for extra_position, weight in boundary_extra_positions_and_weights:
        wall_index = center_count + int(extra_position)
        boundary.append(
            {
                "wall_index_zero_based": wall_index,
                "anti_bicanonical_weight": int(weight),
                "endpoint_blowup_coordinates": coordinates[wall_index],
                "endpoint_image_coordinates": coordinates[wall_index][:model_rank],
            }
        )
    trace = [
        {
            "step": step + 1,
            "source_wall_index_zero_based": step,
            "current_total_exceptional_in_input_basis": standard[model_rank + step],
            "current_square": -1,
            "current_canonical_pairing": -1,
        }
        for step in range(center_count)
    ]
    seed = {
        "endpoint_kind": endpoint_kind,
        "hirzebruch_e": e,
        "parents": [list(map(int, row)) for row in parents_chronological],
        "coordinates": coordinates,
        "boundary": list(boundary_extra_positions_and_weights),
    }
    return {
        "endpoint_kind": endpoint_kind,
        "hirzebruch_e": e,
        "hirzebruch_minimal": None if endpoint_kind == "P2" else e != 1,
        "contracted_wall_indices_zero_based": contracted,
        "sequential_contraction_order_zero_based": contracted,
        "strict_contracted_wall_gram": strict_gram,
        "model_basis_in_input_picard_coordinates": standard[:model_rank],
        "total_exceptional_basis_in_input_picard_coordinates": standard[model_rank:],
        "endpoint_blowup_basis_gram": basis_gram,
        "canonical_coordinates_in_endpoint_blowup_basis": [
            *canonical_model,
            *([1] * center_count),
        ],
        "nef_rays_in_model_basis": nef_rays,
        "all_wall_coordinates_in_endpoint_blowup_basis": coordinates,
        "all_wall_endpoint_image_coordinates": [row[:model_rank] for row in coordinates],
        "boundary": boundary,
        "contraction_trace": trace,
        "permutation_invariant_profile_sha256": _hash(seed),
        "stage_B_payload_complete": True,
    }


def build_fixtures() -> dict[str, Any]:
    free = _program(
        endpoint_kind="P2",
        hirzebruch_e=None,
        parents_chronological=[[], [], [0], [2]],
        extra_rows_chronological=[
            [3, -1, -1, -1, -1],  # primitive -K, doubled to -2K
            [2, -2, -1, -1, -1],  # proximity-consistent divisor
            [2, -1, 0, -2, -2],   # hostile proximity violation
        ],
        boundary_extra_positions_and_weights=[(0, 2)],
    )
    # The first two proper centres commute.  This packet exchanges them and
    # transports every divisor coefficient; its raw witness order differs,
    # while its unlabelled Stage-B1 profile must agree.
    free_commuted = _program(
        endpoint_kind="P2",
        hirzebruch_e=None,
        parents_chronological=[[], [], [1], [2]],
        extra_rows_chronological=[
            [3, -1, -1, -1, -1],
            [2, -1, -2, -1, -1],
            [2, 0, -1, -2, -2],
        ],
        boundary_extra_positions_and_weights=[(0, 2)],
    )
    satellite = _program(
        endpoint_kind="P2",
        hirzebruch_e=None,
        parents_chronological=[[], [0], [0, 1]],
        extra_rows_chronological=[
            [3, -1, -1, -1],
            [3, -2, -1, -1],  # equalities 2=1+1 and 1=1
            [3, -1, -2, -2],  # violates both proximity inequalities
        ],
        boundary_extra_positions_and_weights=[(0, 2)],
    )
    ruled = _program(
        endpoint_kind="F_e",
        hirzebruch_e=2,
        parents_chronological=[[], [0]],
        extra_rows_chronological=[
            [2, 4, -1, -1],
            [1, 2, -1, 0],
        ],
        boundary_extra_positions_and_weights=[(0, 2)],
    )
    invalid_three_parent = _program(
        endpoint_kind="P2",
        hirzebruch_e=None,
        parents_chronological=[[], [0], [0, 1], [0, 1, 2]],
        extra_rows_chronological=[[3, -1, -1, -1, -1]],
        boundary_extra_positions_and_weights=[(0, 2)],
    )
    invalid_disjoint_satellite = _program(
        endpoint_kind="P2",
        hirzebruch_e=None,
        parents_chronological=[[], [], [0], [1], [2, 3]],
        extra_rows_chronological=[[3, -1, -1, -1, -1, -1]],
        boundary_extra_positions_and_weights=[(0, 2)],
    )
    return {
        "schema": SCHEMA,
        "scope": "arrays-only Stage-A-shaped endpoint programs; no named surface/type fixtures",
        "valid_cases": [
            {
                "fixture_id": "free_forest_two_commuting_roots",
                "program": free,
                "expected_center_kind_histogram": {"free": 2, "proper": 2},
            },
            {
                "fixture_id": "same_free_forest_commuted_order",
                "program": free_commuted,
                "expected_center_kind_histogram": {"free": 2, "proper": 2},
            },
            {
                "fixture_id": "satellite_center_two_parents",
                "program": satellite,
                "expected_center_kind_histogram": {"free": 1, "proper": 1, "satellite": 1},
            },
            {
                "fixture_id": "ruled_free_center",
                "program": ruled,
                "expected_center_kind_histogram": {"free": 1, "proper": 1},
            },
        ],
        "invalid_cases": [
            {
                "fixture_id": "hostile_center_with_three_proximity_parents",
                "program": invalid_three_parent,
                "expected_error_substring": "more than two older centres",
            },
            {
                "fixture_id": "hostile_satellite_with_disjoint_exceptional_parents",
                "program": invalid_disjoint_satellite,
                "expected_error_substring": "do not intersect immediately before the blowup",
            }
        ],
        "mutation_pairs": [
            {
                "left_fixture_id": "free_forest_two_commuting_roots",
                "right_fixture_id": "same_free_forest_commuted_order",
                "expected_equal_invariant_profile": True,
            }
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    rendered = json.dumps(build_fixtures(), indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
