#!/usr/bin/env python3
"""Archive the glue-complete cubic-seven unit D24+ H=100 prefix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_a1_golay import (
    COMPLETE_MIRROR_NORM_PROOF,
    UNIT_A1_GOLAY_SEED,
    enumerate_control_levels,
)
from coxeter_polytopes.search.cubic7_d24_plus_unit import (
    analyze_unit_d24_plus_prefix,
    run_unit_d24_plus_prefix,
)


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _vector(vector) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def build_audit() -> dict[str, object]:
    search = run_unit_d24_plus_prefix(100)
    geometry = analyze_unit_d24_plus_prefix(search)
    levels = enumerate_control_levels(UNIT_A1_GOLAY_SEED, 100)
    packet_by_height = {packet.height: packet for packet in search.packets}
    return {
        "schema_version": "0.1.0",
        "lattice": "<-(theta-1)> + D24+",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "height_bound": 100,
        "canonical_mirror_norms": ["1", "2"],
        "complete_mirror_norm_proof": COMPLETE_MIRROR_NORM_PROOF,
        "enumeration": {
            "control_level_count": search.control_level_count,
            "glue_classes": ["D24", *[f"spinor-parity-{code}" for code in range(1, 8)]],
            "dominance": "y1>=...>=y23>=abs(y24)",
            "control_levels": [
                {
                    "norm": _value(level.norm),
                    "height": _value(level.height),
                    "control": _value(level.control),
                    "spatial_square_target": _value(2 * level.coordinate_square_target),
                    "candidate_count": packet_by_height[level.height].candidate_count,
                    "accepted_count": len(packet_by_height[level.height].accepted_roots),
                }
                for level in levels
            ],
            "accepted_root_count": len(search.accepted_roots),
            "accepted_roots": [
                {
                    "norm": _value(root.norm),
                    "height": _value(root.height),
                    "control": _value(root.control),
                    "coset": root.coset,
                    "vector_orthonormal": _vector(root.vector),
                }
                for root in search.accepted_roots
            ],
        },
        "geometry": {
            "engine": "exact double description; adjacency by common-active rank 23",
            "packet_snapshots": [
                {
                    "height": _value(packet.height),
                    "wall_count": packet.wall_count,
                    "extreme_ray_count": packet.extreme_ray_count,
                    "finite_vertex_count": packet.finite_vertex_count,
                    "cusp_count": packet.cusp_count,
                    "hyperideal_ray_count": packet.hyperideal_ray_count,
                }
                for packet in geometry.packet_geometries
            ],
            "insertion_profiles": [
                {
                    "facet_index": profile.facet_index,
                    "negative_ray_count": profile.negative_ray_count,
                    "zero_ray_count": profile.zero_ray_count,
                    "positive_ray_count": profile.positive_ray_count,
                    "output_ray_count": profile.output_ray_count,
                }
                for profile in geometry.insertion_profiles
            ],
            "final": {
                "extreme_ray_count": len(geometry.extreme_rays),
                "finite_vertex_count": geometry.finite_vertex_count,
                "cusp_count": geometry.cusp_count,
                "hyperideal_ray_count": geometry.hyperideal_ray_count,
                "finite_volume": geometry.finite_volume,
                "compact": geometry.compact,
                "rays": [
                    {
                        "vector": _vector(ray.vector),
                        "active_facets": list(ray.active_facets),
                        "norm": _value(ray.norm),
                        "kind": ray.kind,
                    }
                    for ray in geometry.extreme_rays
                ],
            },
        },
        "scope": (
            "This exhausts every primitive mirror norm, control, D24 root coset, "
            "and all seven spinor residue classes through H=100. It proves "
            "nonclosure only for that prefix; it is not a global nonreflectivity claim."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/cubic7_unit_d24_plus_h100.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["geometry"]["packet_snapshots"], indent=2))


if __name__ == "__main__":
    main()
