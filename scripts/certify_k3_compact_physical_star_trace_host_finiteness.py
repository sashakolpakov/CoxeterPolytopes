#!/usr/bin/env python3
"""Replay constants and finite gates for the physical-star trace-host theorem."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_compact_physical_star_trace_host_finiteness.md"
ARCHIVE = (
    ROOT
    / "data/enumerations/k3_compact_physical_star_trace_host_finiteness.json"
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
        rows.append(
            {
                "module_rank": module_rank,
                "hyperbolic_dimension": module_rank - 1,
                "field_degree_ceiling": degree_ceiling,
                "nonrational_possible": degree_ceiling >= 2,
                "cycle_height_bound": 14**module_rank,
            }
        )
    return rows


def golden_codifferent_box() -> dict[str, int | str]:
    """Enumerate a toy codifferent box for O_K of Q(sqrt(5)).

    O_K^vee=(1/sqrt(5)) O_K.  Write z=(u+v*phi)/sqrt(5).  The
    inequalities |z_0|<=14 and |z_1|<=2 imply |v|<=16 and |u|<58,
    so the displayed finite loop is exhaustive.
    """

    sqrt5 = math.sqrt(5)
    phi0 = (1 + sqrt5) / 2
    phi1 = (1 - sqrt5) / 2
    points = 0
    for u in range(-58, 59):
        for v in range(-16, 17):
            z0 = (u + v * phi0) / sqrt5
            z1 = -(u + v * phi1) / sqrt5
            if abs(z0) <= 14 and abs(z1) <= 2:
                points += 1
    return {
        "field": "Q(sqrt(5))",
        "order": "Z[(1+sqrt(5))/2]",
        "codifferent": "(1/sqrt(5))*order",
        "u_min": -58,
        "u_max": 58,
        "v_min": -16,
        "v_max": 16,
        "bounded_box_point_count": points,
    }


def toy_overorder_interval() -> dict[str, int | list[str]]:
    """The index-two interval Z[sqrt(5)] subset O_K."""

    return {
        "lower_order_discriminant": 20,
        "maximal_order_discriminant": 5,
        "maximal_over_lower_index": 2,
        "additive_quotient_order": 2,
        "intermediate_orders": ["Z[sqrt(5)]", "Z[(1+sqrt(5))/2]"],
        "intermediate_order_count": 2,
    }


def build_certificate() -> dict[str, Any]:
    rows = rank_rows()
    nonrational = [row for row in rows if row["nonrational_possible"]]
    if max(row["field_degree_ceiling"] for row in nonrational) != 6:
        raise AssertionError("nonrational K3 degree ceiling changed")
    if max(row["module_rank"] for row in nonrational) != 10:
        raise AssertionError("nonrational K3 module-rank ceiling changed")

    codifferent_toy = golden_codifferent_box()
    if codifferent_toy["bounded_box_point_count"] <= 0:
        raise AssertionError("toy codifferent box unexpectedly empty")

    return {
        "schema": "coxeter.k3-compact-physical-star-trace-host-finiteness.v1",
        "report_sha256": sha256_file(REPORT),
        "rank": {
            "complex_K3_Picard_ceiling": 20,
            "minimum_module_rank": 3,
            "minimum_hyperbolic_dimension": 2,
            "rows": rows,
        },
        "physical_energy": {
            "physical_root_square": -2,
            "sum_of_positive_place_energies": 2,
            "individual_place_energy_upper_bound": 2,
            "normalized_star_pairing_strict_upper_bound": 14,
            "distinguished_field_Gram_entry_strict_upper_bound": 14,
            "definite_place_field_Gram_entry_upper_bound": 2,
            "integral_trace_pairing_strict_upper_bound": 16,
        },
        "order_gate": {
            "Cartan_cycles_are_order_integral": True,
            "Cartan_cycles_generate_full_field": True,
            "cycle_order_is_full_rank": True,
            "coefficient_order_is_an_overorder_of_cycle_order": True,
            "overorders_correspond_to_subgroups_of_finite_quotient": True,
            "toy_interval": toy_overorder_interval(),
        },
        "codifferent_gate": {
            "field_Gram_entries_lie_in_order_codifferent": True,
            "bounded_Minkowski_box_intersection_is_finite": True,
            "golden_toy": codifferent_toy,
        },
        "index_gate": {
            "root_frame": "F=sum_i O*delta_i",
            "lower_containment": "F subset M",
            "upper_containment": "M subset F#",
            "discriminant_group_order": "abs(det(trace Gram of F))",
            "index_square_divides_root_frame_discriminant": True,
            "even_overlattices_are_isotropic_discriminant_subgroups": True,
            "order_and_reflection_stability_are_finite_filters": True,
        },
        "trace_scalar": {
            "K3_sign_relation": "g=-alpha*h with alpha totally positive",
            "native_reflections_determined_by_g": True,
            "alpha_and_h_separately_enumerated": False,
        },
        "scope": {
            "compact_exact_full_trace_K3_required": True,
            "physical_minus_two_same_vector_for_minimal_star_required": True,
            "projected_proportional_wall_only_is_insufficient": True,
            "primitive_full_trace_required": True,
            "ambient_NS_extensions_bounded": False,
            "property4_proved": False,
        },
        "maximality": {
            "needed_for_finite_trace_hosts": False,
            "complete_oriented_host": "(k,tau0,O,M,g,C+)",
            "distinguished_embedding_choices_per_field_are_finite": True,
            "time_cone_component_choices_at_most": 2,
            "same_marked_k_reflection_host_needed_for_chamber_corollary": True,
            "ordinary_rank_one_K3_reflection_maximality_substitutes": False,
            "k_reflection_Q_minus_eigenspace_dimension": "[k:Q]",
            "complete_k_mirror_arrangement_one_chamber_orbit": True,
        },
        "conclusions": {
            "finite_coefficient_orders": True,
            "finite_field_Gram_matrices": True,
            "finite_integral_trace_modules_and_saturations": True,
            "finite_marked_trace_reflection_hosts": True,
            "finite_maximal_chambers_with_same_host_k_maximality": True,
            "finite_ambient_Neron_Severi_extensions": False,
            "carrier_existence": False,
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
            raise SystemExit("archived physical-star host certificate differs")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
