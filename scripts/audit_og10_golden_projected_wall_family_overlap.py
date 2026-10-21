#!/usr/bin/env python3
"""Audit overlap between the frozen scalar-unit and mixed-E8 wall families.

The published family counts are exact, but they are not disjoint.  If
``w=phi^n*v`` and ``v`` is an ordinary E8 vector, then in the mixed notation
``w=a+phi*b`` the vector ``b`` is ``F_n*v``.  Hence the shell-3 condition
``b^2<=6`` overlaps the scalar family exactly for ``n=1,2`` and the indicated
E8 shells.  This script computes the exact overlap in signature-admissible
wall records and the resulting geometric union count.
"""

from __future__ import annotations

import json

from scripts.search_og10_golden_projected_wall import (
    PHI_PAIR,
    _add,
    _coefficient_records,
    _is_signature_admissible,
    _mul,
    _power,
    _scale,
)


ZERO_E8_WALLS = 66
SCALAR_UNIT_WALLS = 2_706_240
MIXED_THROUGH_B6_WALLS = 1_472_640
E8_SHELL_COUNTS = {2: 240, 4: 2_160, 6: 6_720, 8: 17_520}


def fibonacci(index: int) -> int:
    left, right = 0, 1
    for _ in range(index):
        left, right = right, left + right
    return left


def audit() -> dict[str, object]:
    coefficients = _coefficient_records(30)
    overlap_rows: list[dict[str, int]] = []
    overlap_wall_count = 0
    overlap_eligible_direction_count = 0

    for power in range(1, 9):
        scalar = _power(PHI_PAIR, power)
        b_multiplier = fibonacci(power)
        for vector_square, vector_count in E8_SHELL_COUNTS.items():
            b_square = b_multiplier * b_multiplier * vector_square
            if b_square > 6:
                continue
            twice_w_square = _scale(_mul(scalar, scalar), 2 * vector_square)
            admissible_coefficients = sum(
                _is_signature_admissible(_add(record[2], twice_w_square))
                for record in coefficients
            )
            # Signature admissibility implies 0<(w')^2<2, so precisely the
            # rows with a nonzero coefficient count occur in shell 3.
            if admissible_coefficients:
                wall_count = vector_count * admissible_coefficients
                overlap_wall_count += wall_count
                overlap_eligible_direction_count += vector_count
                overlap_rows.append(
                    {
                        "unit_power": power,
                        "E8_vector_square": vector_square,
                        "mixed_b_square": b_square,
                        "E8_vector_count": vector_count,
                        "admissible_coefficient_count": admissible_coefficients,
                        "overlap_wall_count": wall_count,
                    }
                )

    family_indexed_total = (
        ZERO_E8_WALLS + SCALAR_UNIT_WALLS + MIXED_THROUGH_B6_WALLS
    )
    return {
        "kind": "og10_golden_projected_wall_family_overlap_audit",
        "coefficient_bound": 30,
        "mixed_b_square_bound": 6,
        "family_indexed_wall_evaluation_count": family_indexed_total,
        "scalar_mixed_overlap_rows": overlap_rows,
        "scalar_mixed_overlap_eligible_direction_count": (
            overlap_eligible_direction_count
        ),
        "scalar_mixed_overlap_wall_count": overlap_wall_count,
        "distinct_geometric_wall_evaluation_count": (
            family_indexed_total - overlap_wall_count
        ),
        "zero_E8_is_disjoint_from_nonzero_E8_families": True,
        "checksums_are_not_used_for_deduplication": True,
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
