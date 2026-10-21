"""Exact bookkeeping for the global complete K3-root return theorem.

This does not recompute the Nikulin--Vinberg or Wang classifications.  It
checks their published rank counts and the exact intersection of Wang's
18-element no-correction list with the project's nineteen rank-11--19 rows.
"""

from __future__ import annotations

import json


RANK_COUNTS = {
    3: 26,
    4: 14,
    5: 9,
    6: 10,
    7: 9,
    8: 12,
    9: 10,
    10: 9,
    11: 4,
    12: 4,
    13: 3,
    14: 3,
    15: 1,
    16: 1,
    17: 1,
    18: 1,
    19: 1,
}


NO_CORRECTION = {
    "U+E8+E7",
    "U+E8+D6",
    "U+E8+D4+A1",
    "U+E8+D4",
    "U+D8+D4",
    "U+E8+4A1",
    "U+E8+3A1",
    "U+D8+3A1",
    "U+E8+A3",
    "U+D8+2A1",
    "U+2D4+2A1",
    "U+E8+A2",
    "U+E6+A2",
    "U+D4+A3",
    "U+D5+A2",
    "U+D4+A2",
    "U+A4+A2",
    "U+A3+A2",
}


HIGH_RANK_ROWS = (
    (11, 12, "U+E8+A1"),
    (11, 14, "U+D8+A1"),
    (11, 22, "U+2D4+A1"),
    (11, 90, "U+D4+5A1"),
    (12, 13, "U+E8+A2"),
    (12, 14, "U+E8+2A1"),
    (12, 19, "U+D8+2A1"),
    (12, 59, "U+2D4+2A1"),
    (13, 14, "U+E8+A3"),
    (13, 17, "U+E8+3A1"),
    (13, 39, "U+D8+3A1"),
    (14, 15, "U+E8+D4"),
    (14, 20, "U+D8+D4"),
    (14, 27, "U+E8+4A1"),
    (15, 21, "U+E8+D4+A1"),
    (16, 19, "U+E8+D6"),
    (17, 19, "U+E8+E7"),
    (18, 19, "U+2E8"),
    (19, 24, "U+2E8+A1"),
)


def analyze() -> dict[str, object]:
    total = sum(RANK_COUNTS.values())
    low_rank = sum(RANK_COUNTS[r] for r in (3, 4, 5))
    correction_range = sum(v for r, v in RANK_COUNTS.items() if r >= 6)

    high = [
        {
            "rank": rank,
            "hyperbolic_dimension": rank - 1,
            "facets": facets,
            "lattice": lattice,
            "automorphic_correction": lattice not in NO_CORRECTION,
        }
        for rank, facets, lattice in HIGH_RANK_ROWS
    ]
    high_corrected = [row for row in high if row["automorphic_correction"]]
    high_uncorrected = [row for row in high if not row["automorphic_correction"]]

    checks = {
        "rank_counts_sum_to_118": total == 118,
        "rank_3_to_5_sum_to_49": low_rank == 49,
        "wang_range_sum_to_69": correction_range == 69,
        "wang_split_is_51_plus_18": 51 + len(NO_CORRECTION) == correction_range,
        "no_rank_20_or_above": max(RANK_COUNTS) == 19,
        "compact_plus_noncompact_is_118": 8 + 110 == total,
        "nineteen_high_rank_rows": len(high) == 19,
        "high_rank_split_is_7_plus_12": (
            len(high_corrected) == 7 and len(high_uncorrected) == 12
        ),
        "corrected_high_dimensions": sorted(
            {row["hyperbolic_dimension"] for row in high_corrected}
        )
        == [10, 11, 17, 18],
        "uncorrected_strip_H12_through_H16": all(
            not row["automorphic_correction"]
            for row in high
            if 12 <= row["hyperbolic_dimension"] <= 16
        ),
    }

    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ArithmeticError(f"global K3-root transfer checks failed: {failed}")

    return {
        "object": "global complete integral norm-two chamber return",
        "rank_counts": RANK_COUNTS,
        "marked_finite_volume_types": total,
        "compact_types": 8,
        "finite_volume_noncompact_types": 110,
        "maximum_hyperbolic_dimension": 18,
        "automorphic_correction_range": {
            "minimum_hyperbolic_dimension": 5,
            "maximum_hyperbolic_dimension": 18,
            "types": correction_range,
            "corrected": 51,
            "uncorrected": len(NO_CORRECTION),
        },
        "high_rank_rows": high,
        "high_rank_corrected": len(high_corrected),
        "high_rank_uncorrected": len(high_uncorrected),
        "checks": checks,
        "scope": (
            "Published classifications are inputs; the script checks counts and "
            "the 7/12 intersection, not the classification proofs or wall Grams."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))

