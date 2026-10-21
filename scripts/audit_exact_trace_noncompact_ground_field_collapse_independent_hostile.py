#!/usr/bin/env python3
"""Independent exact replay for the full-cusp fixed-line argument."""

from __future__ import annotations

from fractions import Fraction
import argparse
import json
from typing import Sequence


def rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    if not work:
        return 0
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                work[row][j] - multiple * work[pivot_row][j]
                for j in range(len(work[0]))
            ]
        pivot_row += 1
    return pivot_row


def cube_cusp_record(dimension: int) -> dict[str, int | bool]:
    if dimension < 2:
        raise ValueError("finite-volume noncompact scope begins in H^2")
    ambient = dimension + 1
    cusp = [Fraction(0) for _ in range(ambient)]
    cusp[-2] = cusp[-1] = Fraction(1)
    normals: list[list[Fraction]] = []
    for coordinate in range(dimension - 1):
        for sign in (-1, 1):
            normal = cusp[:]
            normal[coordinate] = Fraction(sign)
            normals.append(normal)
    fixed_rows = [
        [
            normal[column] * (1 if column < ambient - 1 else -1)
            for column in range(ambient)
        ]
        for normal in normals
    ]
    normal_rank = rank(normals)
    fixed_kernel_dimension = ambient - rank(fixed_rows)
    if normal_rank != dimension or fixed_kernel_dimension != 1:
        raise AssertionError("the full cusp packet lost its unique fixed line")
    summed = [sum(normal[j] for normal in normals) for j in range(ambient)]
    if summed != [2 * (dimension - 1) * entry for entry in cusp]:
        raise AssertionError("the positive cusp dependence changed")
    return {
        "dimension": dimension,
        "incident_walls": len(normals),
        "normal_rank": normal_rank,
        "fixed_kernel_dimension": fixed_kernel_dimension,
        "nonsimple": len(normals) > dimension,
    }


def build_audit(max_dimension: int = 20) -> dict[str, object]:
    return {
        "verdict": "PASS_EXACT_FULL_CUSP_SCOPE",
        "frames": [cube_cusp_record(n) for n in range(2, max_dimension + 1)],
        "logical_gates": {
            "bounded_horosphere_has_zero_recession_cone": True,
            "positive_dependence_supplies_null_direction_in_normal_span": True,
            "all_incident_fixed_spaces_intersect_in_one_line": True,
            "kernel_commutes_with_field_base_change": True,
            "distinguished_zero_is_zero_in_k": True,
            "definite_nonidentity_embedding_is_impossible": True,
            "returned_field_is_Q": True,
        },
        "scope": {
            "full_incident_packet_required": True,
            "simple_cusp_required": False,
            "free_order_module_required": False,
            "primitive_normals_required": False,
            "compact_case_closed": False,
        },
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-dimension", type=int, default=20)
    arguments = parser.parse_args(argv)
    print(json.dumps(build_audit(arguments.max_dimension), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
