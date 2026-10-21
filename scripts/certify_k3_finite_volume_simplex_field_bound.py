#!/usr/bin/env python3
"""Exact constant replay for the finite-volume K3 simplex field theorem."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_finite_volume_simplex_field_bound.md"
ARCHIVE = ROOT / "data/enumerations/k3_finite_volume_simplex_field_bound.json"


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rank_rows() -> list[dict[str, int | bool]]:
    rows: list[dict[str, int | bool]] = []
    for module_rank in range(3, 21):
        degree_ceiling = 20 // module_rank
        if degree_ceiling < 2:
            rows.append(
                {
                    "module_rank": module_rank,
                    "hyperbolic_dimension": module_rank - 1,
                    "field_degree_ceiling": 1,
                    "rational_only": True,
                    "binary_discriminant_exponent": 0,
                }
            )
            continue
        exponent = (module_rank + 5) * degree_ceiling * (degree_ceiling - 1)
        rows.append(
            {
                "module_rank": module_rank,
                "hyperbolic_dimension": module_rank - 1,
                "field_degree_ceiling": degree_ceiling,
                "rational_only": False,
                "binary_discriminant_exponent": exponent,
            }
        )
    return rows


def build_certificate() -> dict[str, Any]:
    rows = rank_rows()
    maximum_exponent = max(row["binary_discriminant_exponent"] for row in rows)
    if maximum_exponent != 240:
        raise AssertionError("the uniform simplex discriminant exponent changed")
    return {
        "schema": "coxeter.k3-finite-volume-simplex-field-bound.v1",
        "report_sha256": sha256_file(REPORT),
        "complex_K3_picard_ceiling": 20,
        "rank_rows": rows,
        "uniform_simple_cycle_bound": "2^module_rank",
        "selected_generator_count_upper_bound": 2,
        "primitive_grid_maximum_coefficient": 15,
        "uniform_binary_discriminant_exponent": maximum_exponent,
        "uniform_field_discriminant_bound": str(2**maximum_exponent),
        "uniform_field_discriminant_bound_expression": "256^30",
        "conclusions": {
            "finite_field_set": True,
            "compact_simplices_covered": True,
            "ideal_finite_volume_simplices_covered": True,
            "imposed_mixed_height_needed": False,
            "angle_classification_used": False,
            "reflection_maximality_used_for_bound": False,
            "carrier_existence_proved": False,
            "essential_surjectivity_proved": False,
            "general_higher_facet_polytope_closed": False,
        },
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-output", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()
    encoded = encoded_certificate()
    if args.write_output:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_output:
        if not ARCHIVE.exists() or ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("archived simplex certificate differs from replay")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
