#!/usr/bin/env python3
"""Exact constant and scope replay for compact full-trace K3 field finiteness."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_compact_narrow_facet_field_finiteness.md"
ARCHIVE = (
    ROOT
    / "data/enumerations/k3_compact_narrow_facet_field_finiteness.json"
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
        trace_degree_ceiling = 20 // module_rank
        rational_only = trace_degree_ceiling < 2
        selected_generators = 0
        exponent = 0
        cycle_height = 0
        if not rational_only:
            selected_generators = trace_degree_ceiling.bit_length() - 1
            exponent = trace_degree_ceiling * (trace_degree_ceiling - 1)
            cycle_height = 14**module_rank
        rows.append(
            {
                "module_rank": module_rank,
                "hyperbolic_dimension": module_rank - 1,
                "full_trace_degree_ceiling": trace_degree_ceiling,
                "field_degree_bound": max(1, trace_degree_ceiling),
                "rational_only": rational_only,
                "selected_generator_count_bound": selected_generators,
                "cycle_conjugate_height_bound": cycle_height,
                "discriminant_exponent": exponent,
            }
        )
    return rows


def build_certificate() -> dict[str, Any]:
    rows = rank_rows()
    nonrational = [row for row in rows if not row["rational_only"]]
    if max(row["field_degree_bound"] for row in nonrational) != 6:
        raise AssertionError("K3 nonrational degree ceiling changed")
    if max(row["module_rank"] for row in nonrational) != 10:
        raise AssertionError("K3 nonrational rank ceiling changed")
    if max(row["selected_generator_count_bound"] for row in nonrational) != 2:
        raise AssertionError("two-generator reduction changed")
    if max(row["discriminant_exponent"] for row in nonrational) != 30:
        raise AssertionError("discriminant exponent changed")

    uniform_cycle_height = 14**10
    discriminant_base = 32 * uniform_cycle_height
    return {
        "schema": "coxeter.k3-compact-narrow-facet-field-finiteness.v1",
        "report_sha256": sha256_file(REPORT),
        "complex_K3_Picard_ceiling": 20,
        "minimum_module_rank": 3,
        "minimum_hyperbolic_dimension": 2,
        "rank_rows": rows,
        "narrow_facet": {
            "sine_ratio_strict_upper_bound": 2,
            "normalized_star_pairing_strict_upper_bound": 14,
            "constant_replay": "4*2*2-2=14",
        },
        "field_descent": {
            "star_spans_full_Lorentz_space": True,
            "actual_star_normals_contain_a_k_basis": True,
            "one_connected_Lorentz_Cartan_block": True,
            "cycle_field_equals_coefficient_field": True,
            "separate_minimal_field_assumption_required": False,
            "nonidentity_forms_positive_definite": True,
        },
        "integrality": {
            "finite_projective_order_module_allowed": True,
            "primitive_normals_required": False,
            "individual_Cartan_entries_integral_required": False,
            "reflection_word_trace_makes_cycles_integral": True,
        },
        "primitive_element": {
            "degree_ceiling": 6,
            "embedding_pair_count_ceiling": 15,
            "coefficient_minimum": 0,
            "coefficient_maximum": 15,
            "grid_size": 16,
            "selected_generator_count_ceiling": 2,
        },
        "uniform_cycle_conjugate_height": uniform_cycle_height,
        "uniform_discriminant_base": discriminant_base,
        "uniform_discriminant_exponent": 30,
        "uniform_field_discriminant_bound_expression": (
            "(32*14^10)^30"
        ),
        "scope": {
            "compact_Coxeter_chamber_required": True,
            "exact_admissible_k_form_required": True,
            "k_linear_integral_facet_reflections_required": True,
            "primitive_full_trace_complex_K3_required": True,
            "hyperbolic_dimension_one_excluded": True,
            "reflection_maximality_required": False,
            "saturation_required": False,
            "facet_count_bound_required": False,
            "imposed_correspondence_height_required": False,
        },
        "conclusions": {
            "finite_field_set_for_compact_exact_K3_carriers": True,
            "with_noncompact_cusp_theorem_all_finite_volume_closed": True,
            "carrier_existence_proved": False,
            "property4_essential_surjectivity_proved": False,
            "finite_host_set_proved": False,
            "finite_chamber_set_proved": False,
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
            raise SystemExit("archived compact-field certificate differs from replay")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
