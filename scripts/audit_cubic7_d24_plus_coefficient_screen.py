#!/usr/bin/env python3
"""Freeze the exact H=100 small-coefficient screen for cubic-seven D24+."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_d24_plus_coefficients import (
    P7_D24_SEED,
    SUPPORTED_D24_COEFFICIENT_SEEDS,
    TWO_P7_D24_SEED,
    TWO_THETA_D24_SEED,
    UNIT_ISOMETRY_MULTIPLIER,
    analyze_coefficient_prefix,
    enumerate_coefficient_control_levels,
    run_coefficient_prefix,
    unit_coefficient_isometry_audit,
)


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _vector(vector) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def _search_record(seed, search, geometry) -> dict[str, object]:
    levels = enumerate_coefficient_control_levels(seed, search.height_bound)
    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    return {
        "seed_id": seed.seed_id,
        "timelike_coefficient": _value(seed.timelike_coefficient),
        "ideal_support": seed.support,
        "canonical_mirror_norms": [_value(norm) for norm in seed.canonical_mirror_norms],
        "complete_mirror_norm_proof": seed.complete_norm_class_proof,
        "height_bound": search.height_bound,
        "control_level_count": search.control_level_count,
        "control_levels": [
            {
                "height": _value(level.height),
                "norm": _value(level.norm),
                "control": _value(level.control),
                "spatial_square_target": _value(level.target),
            }
            for level in levels
        ],
        "nonempty_packet_count": len(nonempty),
        "accepted_root_count": len(search.accepted_roots),
        "packets": [
            {
                "height": _value(packet.height),
                "candidate_count": packet.candidate_count,
                "accepted_count": len(packet.accepted_roots),
                "accepted_roots": [
                    {
                        "norm": _value(root.norm),
                        "control": _value(root.control),
                        "target": _value(root.target),
                        "coset": root.coset,
                        "vector_orthonormal": _vector(root.vector),
                    }
                    for root in packet.accepted_roots
                ],
            }
            for packet in search.packets
            if packet.candidate_count or packet.accepted_roots
        ],
        "geometry_scope": {
            "packet_limit": geometry.packet_limit,
            "nonempty_packets_inserted": len(geometry.snapshots),
            "engine": "exact double description; adjacency by common-active rank 23",
        },
        "packet_snapshots": [
            {
                "height": _value(snapshot.height),
                "wall_count": snapshot.wall_count,
                "extreme_ray_count": snapshot.extreme_ray_count,
                "finite_vertex_count": snapshot.finite_vertex_count,
                "cusp_count": snapshot.cusp_count,
                "hyperideal_ray_count": snapshot.hyperideal_ray_count,
            }
            for snapshot in geometry.snapshots
        ],
        "insertion_profiles": [
            {
                "facet_index": insertion.facet_index,
                "negative_ray_count": insertion.negative_ray_count,
                "zero_ray_count": insertion.zero_ray_count,
                "positive_ray_count": insertion.positive_ray_count,
                "output_ray_count": insertion.output_ray_count,
            }
            for insertion in geometry.insertions
        ],
        "geometry_final": {
            "extreme_ray_count": len(geometry.rays),
            "finite_vertex_count": geometry.finite_vertex_count,
            "cusp_count": geometry.cusp_count,
            "hyperideal_ray_count": geometry.hyperideal_ray_count,
            "finite_volume": geometry.finite_volume,
            "compact": geometry.compact,
            "all_active_sets_rank_24": True,
            "ray_stream_sha256": geometry.ray_stream_sha256,
        },
        "spinor_accepted_count": sum(
            root.coset != "D" for root in search.accepted_roots
        ),
    }


def build_audit() -> dict[str, object]:
    if not unit_coefficient_isometry_audit():
        raise AssertionError("the theta/(theta-1) unit isometry changed")
    searches = {
        seed: run_coefficient_prefix(seed, 100)
        for seed in SUPPORTED_D24_COEFFICIENT_SEEDS
    }
    # The 2-theta and 2*(2-theta) squareclasses are the new lanes.  The
    # already-known 2-theta control is carried only through its second packet
    # geometrically, while its complete H=100 wall census is retained.
    geometries = {
        TWO_THETA_D24_SEED: analyze_coefficient_prefix(
            searches[TWO_THETA_D24_SEED], 2
        ),
        P7_D24_SEED: analyze_coefficient_prefix(searches[P7_D24_SEED]),
        TWO_P7_D24_SEED: analyze_coefficient_prefix(searches[TWO_P7_D24_SEED]),
    }
    records = [
        _search_record(seed, searches[seed], geometries[seed])
        for seed in SUPPORTED_D24_COEFFICIENT_SEEDS
    ]
    ranking = sorted(
        (
            {
                "seed_id": record["seed_id"],
                "nonempty_packets_compared": len(record["packet_snapshots"]),
                "finite_vertex_count": record["geometry_final"]["finite_vertex_count"],
                "cusp_count": record["geometry_final"]["cusp_count"],
                "hyperideal_ray_count": record["geometry_final"]["hyperideal_ray_count"],
                "finite_volume": record["geometry_final"]["finite_volume"],
            }
            for record in records
        ),
        key=lambda item: (
            item["hyperideal_ray_count"],
            item["cusp_count"],
            -item["finite_vertex_count"],
        ),
    )
    return {
        "schema_version": "0.1.0",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "spatial_lattice": "D24+ (even unimodular Niemeier lattice)",
        "hyperbolic_dimension": 24,
        "unit_squareclass_deduplication": {
            "theta_and_theta_minus_one_are_isometric": True,
            "time_coordinate_multiplier": _value(UNIT_ISOMETRY_MULTIPLIER),
            "identity": "(theta^2-2)^2*theta=theta-1",
            "consequence": "theta was not redundantly rescreened",
        },
        "seeds": records,
        "ranking_by_final_compared_hyperideal_count": ranking,
        "best_exact_survivor": ranking[0]["seed_id"],
        "discovery": None,
        "scope": (
            "Every primitive mirror norm, control, D24 root coset, and all "
            "seven nonzero spinor residues are exhausted through H=100. Cone "
            "geometry is exact through all nonempty packets for p7 and 2*p7, "
            "and through the first two packets for the already-known 2 lane. "
            "Open H=100 prefixes do not prove global nonreflectivity."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/cubic7_d24_plus_coefficient_screen_h100.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["ranking_by_final_compared_hyperideal_count"], indent=2))


if __name__ == "__main__":
    main()
