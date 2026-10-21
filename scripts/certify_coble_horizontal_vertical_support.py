#!/usr/bin/env python3
"""Freeze the characteristic-support consequence of the Coble classification.

This is intentionally a tiny, low-memory certificate.  It does not prove
Kondo's classification; it checks the exact bookkeeping consequence used by
the horizontal/vertical incidence-stack theorem.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROWS = (
    ("I", "all_odd", 1),
    ("I", "all_odd", 2),
    ("II", "all_odd", 1),
    ("V", "3", 2),
    ("VI", "5", 1),
    ("VI", "3", 5),
    ("VII", "5", 1),
    ("MI", "3", 2),
    ("MII", "3", 8),
)


def _rows_in_characteristic(characteristic: str) -> list[dict[str, object]]:
    return [
        {"type": kind, "boundary_count": boundary_count}
        for kind, support, boundary_count in ROWS
        if support == "all_odd" or support == characteristic
    ]


def build_certificate() -> dict[str, object]:
    rows = [
        {
            "type": kind,
            "characteristic_support": support,
            "boundary_count": boundary_count,
            "picard_rank": 10 + boundary_count,
            "hyperbolic_dimension": 9 + boundary_count,
        }
        for kind, support, boundary_count in ROWS
    ]

    generic_rows = [row for row in rows if row["characteristic_support"] == "all_odd"]
    all_boundary_counts = sorted({int(row["boundary_count"]) for row in rows})
    generic_boundary_counts = sorted(
        {int(row["boundary_count"]) for row in generic_rows}
    )

    support_by_boundary_count: dict[str, object] = {}
    for boundary_count in all_boundary_counts:
        supports = {
            str(row["characteristic_support"])
            for row in rows
            if row["boundary_count"] == boundary_count
        }
        support_by_boundary_count[str(boundary_count)] = (
            "all_odd" if "all_odd" in supports else sorted(supports)
        )

    certificate = {
        "schema": "coxeter.coble_horizontal_vertical_support.v1",
        "scope": {
            "base": "Spec Z[1/2]",
            "surface_species": "finite-automorphism strong Coble full-nef carriers",
            "classification_input": "Kondo nine-class theorem",
            "classification_is_external": True,
        },
        "classified_rows": rows,
        "allowed_rows_by_characteristic": {
            "3": _rows_in_characteristic("3"),
            "5": _rows_in_characteristic("5"),
            "generic_characteristic_zero": _rows_in_characteristic("generic"),
            "other_odd": _rows_in_characteristic("other"),
        },
        "support_by_boundary_count": support_by_boundary_count,
        "derived": {
            "all_boundary_counts": all_boundary_counts,
            "all_hyperbolic_dimensions": sorted(
                {int(row["hyperbolic_dimension"]) for row in rows}
            ),
            "horizontal_boundary_counts": generic_boundary_counts,
            "horizontal_hyperbolic_dimensions": sorted(
                {int(row["hyperbolic_dimension"]) for row in generic_rows}
            ),
            "horizontal_dimension_ceiling": max(
                int(row["hyperbolic_dimension"]) for row in generic_rows
            ),
            "vertical_only_high_boundary_support": {"5": [3], "8": [3]},
            "unclassified_boundary_counts_have_empty_carrier_stack": True,
        },
        "theorem": {
            "horizontal_component": (
                "If a Coble-Vinberg incidence-stack component dominates "
                "Spec Z[1/2], its generic characteristic-zero point and exact "
                "same-P return force b in {1,2} and dim(P) in {10,11}."
            ),
            "vertical_component": (
                "For b=5 or b=8 every geometric carrier point is supported "
                "over characteristic 3; for b outside {1,2,5,8} the carrier "
                "stack is empty."
            ),
            "design_consequence": (
                "A characteristic-uniform free/evaluation-surjective branch "
                "cannot produce the H14 or H17 Coble endpoints."
            ),
        },
        "limitations": {
            "does_not_prove_kondo_classification": True,
            "does_not_prove_carrier_nonemptiness_for_raw_D_CobVin": True,
            "does_not_make_carrier_success_a_surface_free_input_condition": True,
            "strict_flagship_status": "support theorem and target-class correction, not L4",
        },
    }

    assert len(rows) == 9
    assert all_boundary_counts == [1, 2, 5, 8]
    assert generic_boundary_counts == [1, 2]
    assert certificate["derived"]["all_hyperbolic_dimensions"] == [10, 11, 14, 17]
    assert certificate["derived"]["horizontal_hyperbolic_dimensions"] == [10, 11]
    assert certificate["derived"]["horizontal_dimension_ceiling"] == 11
    assert support_by_boundary_count == {
        "1": "all_odd",
        "2": "all_odd",
        "5": ["3"],
        "8": ["3"],
    }
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-output", type=Path)
    args = parser.parse_args()

    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_output is not None:
        archived = args.check_output.read_text()
        if archived != rendered:
            raise SystemExit(f"archive mismatch: {args.check_output}")
    print(rendered, end="")


if __name__ == "__main__":
    main()
