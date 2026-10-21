#!/usr/bin/env python3
"""Freeze exact silver A1^24/Golay packet and symmetric-prefix audits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.search.silver_a1_golay import (
    SUPPORTED_SILVER_A1_GOLAY_SEEDS,
    first_packet_cube_geometry,
    unit_dense_prefix_through_h20,
    unit_symmetric_slope_upper_bound,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/enumerations/silver_a1_golay_h24_h100.json"


def _silver(value: object) -> object:
    if hasattr(value, "to_json"):
        return value.to_json()  # type: ignore[no-any-return,union-attr]
    if isinstance(value, tuple):
        return [_silver(entry) for entry in value]
    if isinstance(value, list):
        return [_silver(entry) for entry in value]
    return value


def build_payload() -> dict[str, object]:
    first_packets = []
    for seed in SUPPORTED_SILVER_A1_GOLAY_SEEDS:
        geometry = first_packet_cube_geometry(seed)
        first_packets.append(
            {
                "seed_id": seed.seed_id,
                "timelike_coefficient": _silver(seed.timelike_coefficient),
                "canonical_mirror_norms": _silver(seed.canonical_mirror_norms),
                "complete_norm_class_proof": seed.complete_norm_class_proof,
                "first_packet": {
                    "height": _silver(geometry.height),
                    "norm": _silver(geometry.norm),
                    "control": _silver(geometry.control),
                    "cap_coordinate": _silver(geometry.cap_coordinate),
                    "accepted_wall_count": geometry.accepted_wall_count,
                    "total_wall_count": geometry.total_wall_count,
                    "extreme_ray_count": geometry.extreme_ray_count,
                    "finite_vertex_count": geometry.finite_vertex_count,
                    "cusp_count": geometry.cusp_count,
                    "hyperideal_ray_count": geometry.hyperideal_ray_count,
                    "uniform_ray_norm": _silver(geometry.uniform_ray_norm),
                    "finite_volume": geometry.finite_volume,
                    "compact": geometry.compact,
                },
            }
        )
    dense = unit_dense_prefix_through_h20()
    slope = unit_symmetric_slope_upper_bound(100)
    return {
        "schema_version": "0.1.0",
        "field": "Q(sqrt(2))",
        "positive_lattice": "Niemeier N(A1^24), extended Golay Construction A",
        "first_packets": first_packets,
        "unit_seed_complete_dense_prefix": {
            "height_bound": dense.complete_height_bound,
            "first_height": _silver(dense.first_height),
            "first_family_count": dense.first_family_count,
            "second_height": _silver(dense.second_height),
            "second_control": _silver(dense.second_control),
            "second_pattern": [
                {"value": _silver(value), "multiplicity": multiplicity}
                for value, multiplicity in dense.second_pattern
            ],
            "second_family_count": dense.second_family_count,
            "second_family_rank": dense.second_family_rank,
            "same_height_pairing_half_margins": _silver(
                dense.second_same_height_margins
            ),
            "symmetric_coordinate": _silver(dense.symmetric_coordinate),
            "symmetric_norm": _silver(dense.symmetric_norm),
            "symmetric_extreme": dense.symmetric_extreme,
            "symmetric_hyperideal": dense.symmetric_hyperideal,
        },
        "unit_seed_relaxed_h100_slope_certificate": {
            "height_bound": slope.height_bound,
            "relaxed_pattern_count": slope.relaxed_pattern_count,
            "maximum_spatial_sum_over_control": _silver(slope.maximum_slope),
            "maximizing_height": _silver(slope.maximizing_height),
            "maximizing_control": _silver(slope.maximizing_control),
            "critical_slope_square": _silver(slope.critical_slope_square),
            "slope_square_below_critical": slope.slope_square_below_critical,
            "feasible_uniform_coordinate": _silver(
                slope.feasible_uniform_coordinate
            ),
            "feasible_uniform_norm": _silver(slope.feasible_uniform_norm),
            "hyperideal_prefix_witness": slope.hyperideal_prefix_witness,
            "scope": (
                "Later-wall acuteness is deliberately omitted. This enlarges "
                "the root set, so the maximum is a certified upper bound for "
                "the actual Vinberg prefix and the feasible uniform vector is "
                "a rigorous hyperideal witness through H=100."
            ),
        },
        "conclusion": (
            "None of the three silver A1^24 first packets has finite volume. "
            "The unit-coefficient prefix remains provably non-finite-volume "
            "through exact height 100. No global nonreflectivity claim is made."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
