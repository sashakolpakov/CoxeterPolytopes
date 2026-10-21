#!/usr/bin/env python3
"""Freeze the exact H=300 prefix for ``<-theta(2-theta)> + D24+``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_d24_plus_coefficients import (
    DELTA7,
    P7_D24_SEED,
    analyze_coefficient_prefix,
    d24_coefficient_pair,
    enumerate_coefficient_control_levels,
    run_coefficient_prefix,
)


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _vector(vector: tuple[Cubic7, ...]) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def build_audit() -> dict[str, object]:
    seed = P7_D24_SEED
    search = run_coefficient_prefix(seed, 300)
    geometry = analyze_coefficient_prefix(search)
    levels = enumerate_coefficient_control_levels(seed, 300)
    accepted_packets = tuple(
        packet for packet in search.packets if packet.accepted_roots
    )
    candidate_packets = tuple(
        packet for packet in search.packets if packet.candidate_count
    )
    if len(accepted_packets) != 1 or len(search.accepted_roots) != 1:
        raise AssertionError("the frozen p7 H=300 wall prefix changed")
    first_root = search.accepted_roots[0]
    if d24_coefficient_pair(seed, first_root.vector, first_root.vector) != first_root.norm:
        raise AssertionError("the frozen affine wall has the wrong norm")
    if geometry.hyperideal_ray_count != 23 or len(geometry.rays) != 25:
        raise AssertionError("the frozen p7 H=300 cone changed")

    level_by_key = {
        (level.height, level.norm, level.control, level.target): level
        for level in levels
    }
    if len(level_by_key) != len(levels):
        raise AssertionError("duplicate control level in exact enumeration")

    return {
        "schema_version": "0.1.0",
        "seed_id": seed.seed_id,
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "lattice": "<-theta(2-theta)> + D24+",
        "hyperbolic_dimension": 24,
        "spatial_lattice": "D24+ (even unimodular Niemeier lattice)",
        "height_bound": 300,
        "timelike_coefficient": _value(seed.timelike_coefficient),
        "ideal_support": seed.support,
        "canonical_mirror_norms": [
            _value(norm) for norm in seed.canonical_mirror_norms
        ],
        "complete_mirror_norm_proof": seed.complete_norm_class_proof,
        "delta_divisible_level_reduction": {
            "delta": _value(DELTA7),
            "statement": (
                "For q=delta or 2*delta, reflection integrality and the "
                "self-duality of D24+ force y in delta*D24+. Hence y=delta*z "
                "is enumerated exactly at z^2=target/delta^2; a nonintegral "
                "reduced target proves that level empty without a spatial search."
            ),
        },
        "control_level_count": search.control_level_count,
        "control_levels": [
            {
                "height": _value(level.height),
                "norm": _value(level.norm),
                "control": _value(level.control),
                "spatial_square_target": _value(level.target),
                "delta_reduced_target": (
                    _value(level.target / (DELTA7 * DELTA7))
                    if level.norm in (DELTA7, 2 * DELTA7)
                    else None
                ),
            }
            for level in levels
        ],
        "height_packets": [
            {
                "height": _value(packet.height),
                "candidate_count": packet.candidate_count,
                "accepted_count": len(packet.accepted_roots),
            }
            for packet in search.packets
        ],
        "candidate_packet_count": len(candidate_packets),
        "nonempty_packet_count": len(accepted_packets),
        "accepted_root_count": len(search.accepted_roots),
        "spinor_accepted_count": sum(
            root.coset != "D" for root in search.accepted_roots
        ),
        "accepted_roots": [
            {
                "height": _value(root.height),
                "norm": _value(root.norm),
                "control": _value(root.control),
                "spatial_square_target": _value(root.target),
                "coset": root.coset,
                "vector_orthonormal": _vector(root.vector),
            }
            for root in search.accepted_roots
        ],
        "same_height_audit": {
            "selection": "all survivors at one exact height are inserted simultaneously",
            "pairwise_positive_pairings": 0,
            "passed": True,
        },
        "geometry": {
            "engine": "exact double description; adjacency by common-active rank 23",
            "wall_count": 24 + len(search.accepted_roots),
            "extreme_ray_count": len(geometry.rays),
            "finite_vertex_count": geometry.finite_vertex_count,
            "cusp_count": geometry.cusp_count,
            "hyperideal_ray_count": geometry.hyperideal_ray_count,
            "finite_volume": geometry.finite_volume,
            "compact": geometry.compact,
            "all_active_sets_rank_24": True,
            "ray_stream_sha256": geometry.ray_stream_sha256,
        },
        "frontier": {
            "last_enumerated_height": _value(levels[-1].height),
            "last_enumerated_norm": _value(levels[-1].norm),
            "last_enumerated_control": _value(levels[-1].control),
            "last_enumerated_target": _value(levels[-1].target),
            "all_controls_at_height_at_most_300_exhausted": True,
        },
        "discovery": None,
        "scope": (
            "Every primitive mirror-norm class, compatible control, D24 root "
            "coset, and all seven nonzero D24+ spinor glue residues are exhausted "
            "through exact Vinberg height H=300. The unchanged open prefix does "
            "not imply global nonreflectivity."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/cubic7_p7_d24_plus_h300.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "controls": payload["control_level_count"],
                "accepted_roots": payload["accepted_root_count"],
                "geometry": payload["geometry"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
