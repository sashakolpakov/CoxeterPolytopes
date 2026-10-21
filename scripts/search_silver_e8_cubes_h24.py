#!/usr/bin/env python3
"""Freeze the exact H<=10 silver ``<-a> + E8^3`` discovery screen."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_e8_cubes import (
    INADMISSIBLE_TWO_PLUS_SQRT2,
    SUPPORTED_SILVER_E8_CUBE_SEEDS,
    coefficient_is_arithmetically_admissible,
    first_packet_product_geometry,
    height_zero_dominant_roots,
    run_silver_h24_e8_cube_vinberg,
)


def _silver(value: Silver) -> dict[str, str]:
    return value.to_json()


def _vector(vector: tuple[Silver, ...]) -> list[dict[str, str]]:
    return [_silver(value) for value in vector]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=10)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/silver_e8_cubes_h24_h10.json"),
    )
    args = parser.parse_args()
    records = []
    for seed in SUPPORTED_SILVER_E8_CUBE_SEEDS:
        search = run_silver_h24_e8_cube_vinberg(seed, args.height)
        geometry = first_packet_product_geometry(search)
        records.append(
            {
                "seed_id": seed.seed_id,
                "timelike_coefficient": _silver(seed.timelike_coefficient),
                "coefficient_signs": [
                    seed.timelike_coefficient.sign(0),
                    seed.timelike_coefficient.sign(1),
                ],
                "arithmetically_admissible": seed.admissible,
                "canonical_mirror_norms": [
                    _silver(norm) for norm in seed.canonical_mirror_norms
                ],
                "complete_norm_class_proof": seed.complete_norm_class_proof,
                "height_bound": args.height,
                "height_zero": {
                    "type": "E8^3",
                    "all_root_count": 720,
                    "dominant_orbit_count": len(height_zero_dominant_roots(seed)),
                    "dominant_orbits": [
                        {"norm": _silver(root.norm), "vector": _vector(root.vector)}
                        for root in height_zero_dominant_roots(seed)
                    ],
                },
                "control_level_count": len(search.control_levels),
                "packets": [
                    {
                        "height": _silver(packet.height),
                        "controls": [
                            {
                                "norm": _silver(level.norm),
                                "control": _silver(level.control),
                                "target": _silver(level.target),
                            }
                            for level in packet.controls
                        ],
                        "candidate_count": packet.candidate_count,
                        "accepted_root_count": len(packet.roots),
                        "accepted_roots": [
                            {
                                "norm": _silver(root.norm),
                                "vector": _vector(root.vector),
                            }
                            for root in packet.roots
                        ],
                    }
                    for packet in search.packets
                ],
                "candidate_count": sum(packet.candidate_count for packet in search.packets),
                "accepted_positive_root_count": len(search.accepted_roots),
                "rejected_positive_root_count": len(search.rejected_roots),
                "wall_count": len(search.walls),
                "walls": [_vector(wall) for wall in search.walls],
                "first_packet_geometry": {
                    "structure": "E8-simplex^3",
                    "cap": _silver(geometry.cap),
                    "highest_root": _vector(geometry.highest_root),
                    "extreme_ray_count": geometry.extreme_ray_count,
                    "finite_vertex_count": geometry.finite_vertex_count,
                    "cusp_count": geometry.cusp_count,
                    "hyperideal_ray_count": geometry.hyperideal_ray_count,
                    "finite_volume": geometry.finite_volume,
                    "compact": geometry.compact,
                    "exhaustive_argument": geometry.exhaustive_argument,
                },
                "scope": (
                    "Complete primitive mirror-norm classes and exact Vinberg "
                    "walls through H; hyperideal first-packet rays are not a "
                    "nonreflectivity proof because later walls may cut them."
                ),
            }
        )
    payload = {
        "schema": "silver-h24-e8-cubes-screen-0.1.0",
        "field": "Q(sqrt(2))",
        "ring_of_integers": "Z[sqrt(2)]",
        "coefficient_gate": [
            {
                "coefficient": _silver(value),
                "signs": [value.sign(0), value.sign(1)],
                "admissible": coefficient_is_arithmetically_admissible(value),
                "reason": (
                    "-a is negative at the distinguished embedding and positive "
                    "at the conjugate embedding"
                    if coefficient_is_arithmetically_admissible(value)
                    else "-a is negative at both embeddings, so the conjugate form is not positive definite"
                ),
            }
            for value in (SQRT2, 2 * SQRT2, INADMISSIBLE_TWO_PLUS_SQRT2)
        ],
        "records": records,
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
