#!/usr/bin/env python3
"""Exact constant replay for the nondotted K3 field-generating-core theorem."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_nondotted_field_generating_core_bound.md"
ARCHIVE = (
    ROOT / "data/enumerations/k3_nondotted_field_generating_core_bound.json"
)


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
        rational_only = degree_ceiling < 2
        exponent = 0
        if not rational_only:
            exponent = (
                (module_rank + 7)
                * degree_ceiling
                * (degree_ceiling - 1)
            )
        rows.append(
            {
                "module_rank": module_rank,
                "hyperbolic_dimension": module_rank - 1,
                "field_degree_ceiling": max(1, degree_ceiling),
                "rational_only": rational_only,
                "chordless_cycle_length_bound": module_rank + 2,
                "binary_discriminant_exponent": exponent,
            }
        )
    return rows


def build_certificate() -> dict[str, Any]:
    rows = rank_rows()
    maximum = max(row["binary_discriminant_exponent"] for row in rows)
    maximizers = [
        {
            "module_rank": row["module_rank"],
            "field_degree": row["field_degree_ceiling"],
        }
        for row in rows
        if row["binary_discriminant_exponent"] == maximum
    ]
    if maximum != 300 or maximizers != [
        {"module_rank": 3, "field_degree": 6}
    ]:
        raise AssertionError("the uniform nondotted-core exponent changed")
    return {
        "schema": "coxeter.k3-nondotted-field-generating-core-bound.v1",
        "report_sha256": sha256_file(REPORT),
        "complex_K3_picard_ceiling": 20,
        "minimum_module_rank": 3,
        "rank_rows": rows,
        "nonrational_degree_ceiling": 6,
        "chordless_cycle_length_bound": "module_rank + 2",
        "generating_cycle_conjugate_bound": "2^(module_rank + 2)",
        "selected_generator_count_upper_bound": 2,
        "primitive_grid": {"coefficient_minimum": 0, "coefficient_maximum": 15},
        "primitive_conjugate_bound": "2^(module_rank + 6)",
        "primitive_difference_bound": "2^(module_rank + 7)",
        "rankwise_binary_discriminant_formula": (
            "(module_rank + 7) * D * (D - 1), "
            "D=floor(20/module_rank)"
        ),
        "uniform_binary_discriminant_exponent": maximum,
        "uniform_field_discriminant_bound_expression": "2^300",
        "maximizers": maximizers,
        "scope": {
            "primitive_exact_full_trace_complex_K3_required": True,
            "minimal_Vinberg_field_generated_by_returned_cycles_required": True,
            "generating_pairs_and_cycle_perimeter_edges_at_most_four": True,
            "dotted_edges_outside_generating_core_allowed": True,
            "facet_count_bound_required": False,
            "finite_volume_required": False,
            "reflection_maximality_required_for_bound": False,
        },
        "conclusions": {
            "finite_field_set": True,
            "carrier_existence_proved": False,
            "essential_surjectivity_proved": False,
            "canonical_mixed_height_bounded": False,
            "all_higher_facet_polytopes_covered": False,
            "chamber_finiteness_proved": False,
            "optimal_numerical_discriminant_constant_claimed": False,
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
            raise SystemExit("archived nondotted-core certificate differs from replay")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
