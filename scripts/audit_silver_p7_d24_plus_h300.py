#!/usr/bin/env python3
"""Freeze the exact H=300 prime-seven silver D24+ wall census."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_p7_d24_plus import (
    SILVER_P7_D24_SEED,
    analyze_silver_p7_d24_prefix,
    audit_first_packet_hyperideal_survival,
    audit_frozen_h100_witness,
    certify_h300_replacement_hyperideal_ray,
    enumerate_p7_d24_control_levels,
    run_silver_p7_d24_prefix,
    silver_p7_d24_pair,
)


V1_H100_SHA256 = {
    "backend": "e9d9359209e980eb7e1d76371c3de344e2ac8d1a9276ee5bb84bcb73f367cd81",
    "tests": "719e37b48b2edd2a55a99ff21e3f28fd5efa0af2f581c0cc12caf0a8c8650533",
    "audit": "67997a7ba00214b573b5574ebc74b3cd6ff3a40d7c262e8707d84c1fe36da2c8",
    "artifact": "ec01a747418b4fb0ecbfa25040cbae6d3f6369746f533c5b5447ff9ac11c5fbd",
    "report": "5d88f4454ff5a6d23c51d89ddb04034b6c8be760b0bea621df06d6131740afea",
}


def _value(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in (0, 1)],
    }


def _vector(vector) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def build_audit() -> dict[str, object]:
    search = run_silver_p7_d24_prefix(300)
    levels = enumerate_p7_d24_control_levels(300)
    enumeration_by_level = {
        enumeration.level: enumeration for enumeration in search.level_enumerations
    }
    first_geometry = analyze_silver_p7_d24_prefix(search, 1)
    witness_audit = audit_frozen_h100_witness(search)
    survival = audit_first_packet_hyperideal_survival(search)
    replacement = certify_h300_replacement_hyperideal_ray(search)
    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    first_post_h100 = next(
        packet for packet in nonempty if (packet.height - Silver(100)).sign(0) > 0
    )

    witness_pairings = []
    for packet in nonempty[1:]:
        pairings = tuple(
            silver_p7_d24_pair(witness_audit.vector, root.vector)
            for root in packet.accepted_roots
        )
        witness_pairings.append(
            {
                "height": _value(packet.height),
                "pairings": [_value(value) for value in pairings],
            }
        )

    return {
        "schema_version": "0.3.0",
        "field": "Q(sqrt(2)), ring of integers Z[sqrt(2)]",
        "lattice": "<-(2*sqrt(2)-1)> + D24+",
        "hyperbolic_dimension": 24,
        "height_bound": 300,
        "seed_id": SILVER_P7_D24_SEED.seed_id,
        "v1_h100_lineage": {
            "description": (
                "SHA-256 pins of the independently reverified H=100 freeze "
                "before the H=300 backend/test extension."
            ),
            "sha256": V1_H100_SHA256,
        },
        "canonical_mirror_norms": [
            _value(norm) for norm in SILVER_P7_D24_SEED.canonical_mirror_norms
        ],
        "complete_mirror_norm_proof": (
            SILVER_P7_D24_SEED.complete_norm_class_proof
        ),
        "enumeration": {
            "control_level_count": len(levels),
            "height_packet_count": len(search.packets),
            "accepted_root_count": len(search.accepted_roots),
            "final_wall_count": len(search.walls),
            "glue_classes": [
                "D24",
                "spinor-1",
                "spinor-sqrt(2)",
                "spinor-(1+sqrt(2))",
            ],
            "control_levels": [
                {
                    "norm": _value(level.norm),
                    "height": _value(level.height),
                    "control": _value(level.control),
                    "spatial_square_target": _value(level.target),
                    "reduced_target_after_p7_divisibility": (
                        None
                        if enumeration_by_level[level].reduced_target is None
                        else _value(enumeration_by_level[level].reduced_target)
                    ),
                    "glue_pattern_count": enumeration_by_level[
                        level
                    ].glue_pattern_count,
                    "reflection_integral_count_before_primitivity": (
                        enumeration_by_level[level].reflection_integral_count
                    ),
                    "imprimitive_sqrt2_count": enumeration_by_level[
                        level
                    ].imprimitive_sqrt2_count,
                    "primitive_candidate_count": len(
                        enumeration_by_level[level].candidates
                    ),
                }
                for level in levels
            ],
            "packets": [
                {
                    "height": _value(packet.height),
                    "raw_reflective_count": packet.raw_reflective_count,
                    "imprimitive_sqrt2_count": packet.imprimitive_sqrt2_count,
                    "primitive_candidate_count": packet.candidate_count,
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
            ],
            "nonempty_packet_summary": [
                {
                    "height": _value(packet.height),
                    "raw_reflective_count": packet.raw_reflective_count,
                    "imprimitive_sqrt2_count": packet.imprimitive_sqrt2_count,
                    "primitive_candidate_count": packet.candidate_count,
                    "accepted_count": len(packet.accepted_roots),
                    "accepted_spinor_count": sum(
                        root.coset != "D24" for root in packet.accepted_roots
                    ),
                }
                for packet in nonempty
            ],
            "first_post_h100_nonempty_packet": {
                "height": _value(first_post_h100.height),
                "accepted_count": len(first_post_h100.accepted_roots),
            },
        },
        "exact_rank_certification": {
            "method": (
                "For full-row-rank adjacency sets, specialize sqrt(2) to "
                "both roots modulo the split primes 17, 23, and 31. A "
                "nonzero modular minor proves the same Q(sqrt(2)) minor is "
                "nonzero. Inconclusive specializations fall back to exact "
                "Silver Gaussian elimination."
            ),
            "specializations": [
                {"prime": 17, "roots_of_two": [6, 11]},
                {"prime": 23, "roots_of_two": [5, 18]},
                {"prime": 31, "roots_of_two": [8, 23]},
            ],
            "first_packet_ray_stream_sha256": first_geometry.ray_stream_sha256,
        },
        "h100_witness_extension_audit": {
            "witness_vector": _vector(witness_audit.vector),
            "witness_norm": _value(witness_audit.norm),
            "packet_summaries": [
                {
                    "height": _value(packet.height),
                    "wall_count": packet.wall_count,
                    "negative_pairing_count": packet.negative_pairing_count,
                    "zero_pairing_count": packet.zero_pairing_count,
                    "positive_pairing_count": packet.positive_pairing_count,
                    "minimum_pairing": _value(packet.minimum_pairing),
                    "maximum_pairing": _value(packet.maximum_pairing),
                }
                for packet in witness_audit.packets
            ],
            "exact_pairings_by_packet": witness_pairings,
            "first_cut_height": (
                None
                if witness_audit.first_cut_height is None
                else _value(witness_audit.first_cut_height)
            ),
            "first_cut_packet_positive_indices": list(
                witness_audit.first_cut_packet_positive_indices
            ),
        },
        "first_packet_hyperideal_ray_survival": [
            {
                "height": _value(item.height),
                "cumulative_later_wall_count": item.cumulative_later_wall_count,
                "surviving_hyperideal_ray_count": (
                    item.surviving_hyperideal_ray_count
                ),
            }
            for item in survival
        ],
        "h300_replacement_hyperideal_ray": {
            "vector": _vector(replacement.vector),
            "norm": _value(replacement.norm),
            "wall_count": replacement.wall_count,
            "all_wall_pairings_nonpositive": all(
                value.sign(0) <= 0 for value in replacement.pairings
            ),
            "zero_pairing_count": sum(not value for value in replacement.pairings),
            "active_facets": list(replacement.active_facets),
            "active_rank": replacement.active_rank,
            "minimum_pairing": _value(min(replacement.pairings)),
            "maximum_pairing": _value(max(replacement.pairings)),
            "finite_volume": False,
            "compact": False,
        },
        "geometry_status": {
            "h100_nonfinite_volume_certificate": (
                "retained in the pinned v1 artifact and report"
            ),
            "old_witness_survives_through_height": _value(Silver(56, 40)),
            "all_109_first_packet_hyperideal_rays_cut_by_height": _value(
                Silver(62, 44)
            ),
            "replacement_extreme_ray_search": (
                "a direct exact rank-24 replacement ray was certified "
                "against all 348 walls without completing double description"
            ),
            "h300_finite_volume": False,
            "h300_compact": False,
            "eventual_chamber_status": "unresolved",
        },
        "discovery": None,
        "scope": (
            "Every primitive norm class, control, D24 root coset, and all "
            "three spinor residues are exhausted through algebraic height "
            "300, yielding an exact 348-wall census. The witness-cut audit "
            "records how the old H=100 rays are cut. A separate direct "
            "rank-24 hyperideal ray proves that the complete H=300 prefix "
            "still has infinite volume. Later walls may cut this ray, so the "
            "eventual chamber remains unresolved."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/silver_p7_d24_plus_h300.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "control_level_count": payload["enumeration"][
                    "control_level_count"
                ],
                "accepted_root_count": payload["enumeration"][
                    "accepted_root_count"
                ],
                "final_wall_count": payload["enumeration"]["final_wall_count"],
                "first_cut_height": payload["h100_witness_extension_audit"][
                    "first_cut_height"
                ],
                "eventual_chamber_status": payload["geometry_status"][
                    "eventual_chamber_status"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
