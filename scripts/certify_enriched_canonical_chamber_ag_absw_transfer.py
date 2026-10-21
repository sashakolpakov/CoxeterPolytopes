#!/usr/bin/env python3
"""Exact dependency and constant replay for the enriched AG/ABSW transfer."""

from __future__ import annotations

import argparse
import hashlib
import json
from math import floor, log2
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/enriched_canonical_chamber_ag_absw_transfer.md"
ARCHIVE = (
    ROOT
    / "data/enumerations/enriched_canonical_chamber_ag_absw_transfer.json"
)

DEPENDENCIES = {
    "results/canonical_chamber_rigidity_absw_transfer_theorem.md": (
        "850d2b5832c3f7bd9fa85c39a6e84449c8f1fb5ede6c98de931f966218da3767"
    ),
    "results/rational_direct_full_ns_n11_dimension_theorem.md": (
        "a5dbd6e2fbbd7cd55025f5af324de6b778a196cf4f72cf06903a1dce2e4ae27f"
    ),
    "results/pfirst_weak_log_complete_carrier_transfer.md": (
        "bff338b5e2f4dde2c918742a853656661c2ac074684938031ebbc95df8fd14b2"
    ),
    "results/equivariant_minimalization_trace_descent_theorem.md": (
        "ff8a89ae37a24ad2a10515f5acd931b15a280ed8f7fafb77e1db0c3d6b61ae5d"
    ),
    "results/rational_pfirst_ordinary_dimension_bound_and_exact_carrier.md": (
        "8ea1ede6cc1717ede57dc902e459ceb831682f9f6cd6328c70fc7657452920af"
    ),
    "results/k3_chordless_cycle_height_field_bound.md": (
        "3372cddd560471830fc0a9cbaac02560ce4a682781db1cdc877b8cd108ea471c"
    ),
    "results/k3_chordless_cycle_height_field_bound_independent_hostile_audit.md": (
        "eb9e8b40ce646e7c7806eb157fbcd6a80393fd9d52a6bb2fb06061fb3fcd3914"
    ),
    "results/rank4_k3_exact_maximal_eigen_interval_infinite_fields.md": (
        "9351875964fbc73da81b55ec3ff68db54af3dbe1ec90272555cff1938ab13c95"
    ),
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def field_record(
    dimension: int,
    correspondence_height: int = 256,
    square_lower_bound: int = 2,
) -> dict[str, Any]:
    if dimension < 2:
        raise ValueError("the enriched field theorem is scoped to n>=2")
    if correspondence_height <= 0 or square_lower_bound <= 0:
        raise ValueError("height and square bounds must be positive")
    trace_degree_ceiling = 20 // (dimension + 1)
    if trace_degree_ceiling < 2:
        return {
            "dimension": dimension,
            "full_trace_K3_degree_ceiling": trace_degree_ceiling,
            "field_degree_bound": 1,
            "ground_field": "Q",
            "nonrational_field_possible": False,
        }
    degree = trace_degree_ceiling
    length = dimension + 3
    selected = floor(log2(degree))
    primitive_grid = degree * (degree - 1) // 2 + 1
    small_bound = 2 ** (length + 1)
    height_ratio_ceiling = (
        correspondence_height + square_lower_bound - 1
    ) // square_lower_bound
    coordinate_bound = max(small_bound, height_ratio_ceiling)
    primitive_bound = selected * primitive_grid * coordinate_bound
    exponent = degree * (degree - 1)
    return {
        "dimension": dimension,
        "full_trace_K3_degree_ceiling": trace_degree_ceiling,
        "field_degree_bound": degree,
        "chordless_cycle_length_bound": length,
        "selected_generator_count_bound": selected,
        "primitive_combination_coefficient_bound": primitive_grid,
        "small_affine_cycle_conjugate_bound": small_bound,
        "height_ratio_ceiling": height_ratio_ceiling,
        "primitive_generator_conjugate_bound": primitive_bound,
        "field_discriminant_exponent": exponent,
        "field_discriminant_bound": str(primitive_bound**exponent),
        "nonrational_field_possible": True,
    }


def build_certificate() -> dict[str, Any]:
    live_dependencies = {
        relative: sha256_file(ROOT / relative)
        for relative in sorted(DEPENDENCIES)
    }
    rows = [field_record(n) for n in range(2, 26)]
    coarse_base = 32 * max(8192, 256 // 2)
    return {
        "schema": "coxeter.enriched-canonical-chamber-ag-absw-transfer.v1",
        "report_sha256": sha256_file(REPORT),
        "dependency_sha256": live_dependencies,
        "dependency_pins_match": live_dependencies == DEPENDENCIES,
        "dimension_branches": {
            "direct_rational_full_NS": {
                "hyperbolic_dimension_bound": 11,
                "minimal_Vinberg_field": "Q",
            },
            "weak_log_rational": {
                "hyperbolic_dimension_bound": 25,
                "minimal_Vinberg_field": "Q",
            },
            "complex_full_trace": {
                "rank_inequality": "[k:Q]*(n+1)<=20",
                "nonrational_dimension_bound": 9,
            },
            "ordinary_property4": {
                "hyperbolic_dimension_bound": 9,
                "exact_carriers_for_dimensions": list(range(3, 10)),
            },
            "global_dimension_bound": 25,
        },
        "field_rows_C256_s2": rows,
        "n2_exact_formula": "(32*max(64,C/s0))^30",
        "n2_C256_s2_discriminant_bound": field_record(2)[
            "field_discriminant_bound"
        ],
        "dimension_uniform_nonrational_formula": (
            "(32*max(8192,C/s0))^30"
        ),
        "dimension_uniform_C256_s2_discriminant_bound": str(
            coarse_base**30
        ),
        "five_gates": [
            "full_trace_fidelity",
            "full_nef_or_exact_eigen_nef_fidelity",
            "canonical_log_MMP_fidelity",
            "raw_essential_surjectivity",
            "uniform_polarized_spectral_height_and_finite_hosts",
        ],
        "conclusions": {
            "conditional_global_dimension_at_most_25": True,
            "nonrational_complex_dimension_at_most_9": True,
            "bounded_large_cycle_bidegrees_give_finite_fields": True,
            "bounded_complete_Hilbert_Chow_grammar_supplies_gate5": True,
            "finite_fields_follow_from_rank_alone": False,
            "finite_fields_imply_finite_integral_hosts": False,
            "finite_hosts_plus_saturation_give_finite_chambers": True,
            "all_raw_maximal_arithmetic_chambers_have_carriers": False,
            "maximality_forces_spectral_height": False,
        },
        "remaining_gates": {
            "rational": "twice-primitive (-4) adapted-host property 4",
            "nonrational": (
                "canonical trace polarization and bounded large-cycle "
                "correspondence degrees"
            ),
        },
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    envelope = {
        "payload": payload,
        "payload_sha256": digest(payload),
    }
    return json.dumps(envelope, indent=2, sort_keys=True) + "\n"


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
            raise SystemExit("archived certificate differs from live replay")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
