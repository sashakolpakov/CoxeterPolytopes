#!/usr/bin/env python3
"""Freeze the exact H=100 prime-seven silver D24+ search."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_p7_d24_plus import (
    P7_DELTA,
    SILVER_P7_D24_SEED,
    analyze_silver_p7_d24_prefix,
    certify_full_prefix_hyperideal_ray,
    enumerate_p7_d24_control_levels,
    p7_factorization_audit,
    run_silver_p7_d24_prefix,
)


def _value(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in (0, 1)],
    }


def _vector(vector) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def build_audit() -> dict[str, object]:
    if not p7_factorization_audit():
        raise AssertionError("the prime-seven coefficient factorization changed")
    search = run_silver_p7_d24_prefix(100)
    geometry = analyze_silver_p7_d24_prefix(search, 1)
    witness = certify_full_prefix_hyperideal_ray(search)
    levels = enumerate_p7_d24_control_levels(100)
    enumeration_by_level = {
        enumeration.level: enumeration for enumeration in search.level_enumerations
    }
    return {
        "schema_version": "0.1.0",
        "field": "Q(sqrt(2)), ring of integers Z[sqrt(2)]",
        "lattice": "<-(2*sqrt(2)-1)> + D24+",
        "spatial_lattice": "D24+ (even unimodular Niemeier lattice)",
        "hyperbolic_dimension": 24,
        "height_bound": 100,
        "seed_id": SILVER_P7_D24_SEED.seed_id,
        "timelike_coefficient": _value(SILVER_P7_D24_SEED.timelike_coefficient),
        "prime_seven_generator": _value(P7_DELTA),
        "coefficient_factorization": "(1+sqrt(2))*(2*sqrt(2)-1)=3+sqrt(2)",
        "canonical_mirror_norms": [
            _value(norm) for norm in SILVER_P7_D24_SEED.canonical_mirror_norms
        ],
        "complete_mirror_norm_proof": (
            SILVER_P7_D24_SEED.complete_norm_class_proof
        ),
        "enumeration": {
            "control_level_count": len(levels),
            "height_packet_count": len(search.packets),
            "glue_classes": [
                "D24",
                "spinor-1",
                "spinor-sqrt(2)",
                "spinor-(1+sqrt(2))",
            ],
            "dominance": "y1>=...>=y23>=abs(y24)",
            "prime_seven_divisibility": (
                "For q divisible by delta, reflection integrality forces "
                "y in delta*D24+; enumerate y/delta at target y^2/delta^2."
            ),
            "primitivity": (
                "Test divisibility by sqrt(2) in the integral D24+ glue basis; "
                "discard norm-2 vectors whose quotient is a norm-1 mirror."
            ),
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
                    "imprimitive_vectors": [
                        {
                            "raw_norm": _value(level.norm),
                            "raw_vector": _vector(vector),
                            "primitive_quotient": _vector(
                                tuple(value / Silver(0, 1) for value in vector)
                            ),
                        }
                        for vector in enumeration_by_level[level].imprimitive_vectors
                    ],
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
            "accepted_root_count": len(search.accepted_roots),
            "accepted_spinor_count": sum(
                root.coset != "D24" for root in search.accepted_roots
            ),
            "final_wall_count": len(search.walls),
        },
        "first_complete_packet_geometry": {
            "scope": "all three simultaneous walls at height 5+4*sqrt(2)",
            "engine": "exact double description; common-active adjacency rank 23",
            "snapshots": [
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
            "ray_stream_sha256": geometry.ray_stream_sha256,
            "rays": [
                {
                    "vector": _vector(ray.vector),
                    "active_facets": list(ray.active_facets),
                    "norm": _value(ray.norm),
                    "kind": ray.kind,
                }
                for ray in geometry.rays
            ],
        },
        "full_h100_nonfinite_volume_certificate": {
            "method": (
                "A rank-24 hyperideal extreme ray of the first complete "
                "packet satisfies every later accepted half-space. It "
                "therefore remains an extreme hyperideal ray of the full "
                "55-wall H=100 prefix."
            ),
            "surviving_first_packet_hyperideal_ray_count": (
                witness.surviving_hyperideal_ray_count
            ),
            "chosen_ray": {
                "vector": _vector(witness.ray.vector),
                "norm": _value(witness.ray.norm),
                "old_active_facets": list(witness.ray.active_facets),
                "new_zero_facets": list(witness.new_zero_facets),
                "final_active_rank": witness.final_active_rank,
            },
            "later_wall_count": witness.later_wall_count,
            "later_wall_pairings": [_value(value) for value in witness.later_pairings],
            "all_later_pairings_nonpositive": all(
                value.sign(0) <= 0 for value in witness.later_pairings
            ),
            "finite_volume": False,
            "compact": False,
        },
        "discovery": None,
        "scope": (
            "Every canonical primitive mirror norm, control coordinate, D24 "
            "root coset, and all three nonzero spinor residues are exhausted "
            "through algebraic height 100. The exact surviving hyperideal ray "
            "proves the full H=100 prefix has infinite volume. This bounded "
            "prefix is not a global nonreflectivity theorem."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/silver_p7_d24_plus_h100.json"),
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
                "first_packet_snapshots": payload[
                    "first_complete_packet_geometry"
                ]["snapshots"],
                "full_h100_finite_volume": payload[
                    "full_h100_nonfinite_volume_certificate"
                ]["finite_volume"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
