#!/usr/bin/env python3
"""Generate exact cubic-seven ``B_(n-1) x A1`` Vinberg-prefix records."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_bna1 import (
    L2_BNA1_SEED,
    L7_BNA1_SEED,
    Cubic7BnA1Seed,
    enumerate_height_zero_patterns,
    run_bna1_vinberg_search,
)


def _parse_dimensions(text: str) -> tuple[int, ...]:
    dimensions: set[int] = set()
    for field in text.split(","):
        field = field.strip()
        if not field:
            continue
        if "-" in field:
            lower, upper = (int(value) for value in field.split("-", 1))
            dimensions.update(range(lower, upper + 1))
        else:
            dimensions.add(int(field))
    if not dimensions or min(dimensions) < 3:
        raise argparse.ArgumentTypeError("dimensions must be integers at least three")
    return tuple(sorted(dimensions))


def _parse_heights(text: str) -> tuple[int, ...]:
    try:
        heights = tuple(sorted({int(field.strip()) for field in text.split(",") if field.strip()}))
    except ValueError as error:
        raise argparse.ArgumentTypeError("heights must be comma-separated integers") from error
    if not heights or heights[0] < 0:
        raise argparse.ArgumentTypeError("heights must be nonnegative integers")
    return heights


def _element(value: Cubic7) -> dict[str, str]:
    return value.to_json()


def _vector(vector: tuple[Cubic7, ...]) -> list[dict[str, str]]:
    return [_element(value) for value in vector]


def _norm_counts(values) -> list[dict[str, object]]:
    counts = Counter(value.norm for value in values)
    return [
        {"norm": _element(norm), "count": count}
        for norm, count in sorted(counts.items())
    ]


def _record(
    seed: Cubic7BnA1Seed,
    dimension: int,
    height: int,
    analyze_geometry: bool,
    geometry_engine: str,
    geometry_ray_limit: int | None,
    geometry_time_limit_seconds: float | None,
) -> dict[str, object]:
    search = run_bna1_vinberg_search(
        seed,
        dimension,
        height,
        analyze_geometry=analyze_geometry,
        geometry_engine=geometry_engine,
        geometry_ray_limit=geometry_ray_limit,
        geometry_time_limit_seconds=geometry_time_limit_seconds,
    )
    patterns = enumerate_height_zero_patterns(seed, dimension)
    geometry = search.geometry
    return {
        "seed": seed.seed_id,
        "hyperbolic_dimension": dimension,
        "height_bound": height,
        "diagonal_form": [_element(value) for value in seed.coefficients(dimension)],
        "complete_canonical_mirror_norms": [
            _element(value) for value in seed.canonical_mirror_norms
        ],
        "local_norm_class_proof": seed.complete_norm_class_proof,
        "height_zero": {
            "weyl_type": f"B{dimension - 1}xA1",
            "orbit_pattern_count": len(patterns),
            "root_count": sum(pattern.orbit_size for pattern in patterns),
            "patterns": [
                {
                    "norm": _element(pattern.norm),
                    "last_coordinate": _element(pattern.last_coordinate),
                    "ordinary_coordinates": _vector(pattern.ordinary_coordinates),
                    "orbit_size": pattern.orbit_size,
                }
                for pattern in patterns
            ],
        },
        "positive_height": {
            "control_counts_by_norm": [
                {"norm": _element(norm), "count": count}
                for norm, count in search.pool.control_counts_by_norm
            ],
            "root_pool_count": len(search.pool.roots),
            "root_pool_counts_by_norm": _norm_counts(search.pool.roots),
            "accepted_root_count": len(search.accepted_roots),
            "accepted_counts_by_norm": _norm_counts(search.accepted_roots),
            "rejected_root_count": len(search.rejected_roots),
            "accepted_roots": [
                {
                    "norm": _element(root.norm),
                    "height": _element(root.height),
                    "vector": _vector(root.vector),
                }
                for root in search.accepted_roots
            ],
        },
        "wall_count": len(search.walls),
        "walls": [_vector(wall) for wall in search.walls],
        "geometry": None
        if geometry is None
        else {
            "finite_vertex_count": geometry.finite_vertex_count,
            "cusp_count": geometry.cusp_count,
            "hyperideal_ray_count": geometry.hyperideal_ray_count,
            "irredundant": geometry.irredundant,
            "consistent_time_sheet": geometry.consistent_time_sheet,
            "finite_volume": geometry.finite_volume,
            "compact": geometry.compact,
            "extreme_rays": [
                {
                    "kind": ray.kind,
                    "norm": _element(ray.norm),
                    "active_facets": list(ray.active_facets),
                    "vector": _vector(ray.vector),
                }
                for ray in geometry.extreme_rays
            ],
        },
        "geometry_status": (
            "exhaustive"
            if geometry is not None
            else "not_run_by_requested_height_policy"
        ),
        "scope": {
            "root_prefix_exhaustive_through_height": height,
            "all_primitive_mirror_norm_classes_exhausted": True,
            "finite_volume_stopping_reached": bool(geometry and geometry.finite_volume),
            "nonreflectivity_proved": False,
            "warning": (
                "A hyperideal ray in a bounded Vinberg prefix proves only that this prefix "
                "has not closed; later exact roots can cut it off."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dimensions", type=_parse_dimensions, default=(5,))
    parser.add_argument("--height", type=int)
    parser.add_argument(
        "--heights",
        type=_parse_heights,
        help="comma-separated adaptive height checkpoints (overrides --height)",
    )
    parser.add_argument("--seeds", choices=("L2", "L7", "both"), default="both")
    parser.add_argument("--no-geometry", action="store_true")
    parser.add_argument(
        "--geometry-through-height",
        type=int,
        help="run geometry only for checkpoints at or below this height",
    )
    parser.add_argument(
        "--geometry-engine", choices=("support", "incremental"), default="support"
    )
    parser.add_argument("--geometry-ray-limit", type=int)
    parser.add_argument("--geometry-time-limit-seconds", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    heights = args.heights or ((args.height if args.height is not None else 5),)
    seeds = (
        (L2_BNA1_SEED, L7_BNA1_SEED)
        if args.seeds == "both"
        else (L2_BNA1_SEED if args.seeds == "L2" else L7_BNA1_SEED,)
    )
    payload = {
        "schema": "cubic7-bna1-vinberg-prefixes-0.2.0",
        "field": {
            "generator": "theta=2*cos(2*pi/7)",
            "minimal_polynomial": "theta^3+theta^2-2theta-1",
            "ring_of_integers": "Z[theta]",
        },
        "enumeration": {
            "ordinary_coordinates": (
                "dominant ordered square partitions; zero padding records stabilization"
            ),
            "arithmetic": "exact Cubic7 arithmetic and certified real-embedding signs",
            "height_zero_orbit_quotient": (
                "full B_(d-1)xA1 Weyl quotient: signed permutations of ordinary "
                "coordinates and the independent sign of the final coordinate"
            ),
            "cone_rays": (
                "exact exhaustive support enumeration"
                if args.geometry_engine == "support"
                else "exact incremental double description with rank-certified adjacency"
            ),
        },
        "records": [
            _record(
                seed,
                dimension,
                height,
                not args.no_geometry
                and (
                    args.geometry_through_height is None
                    or height <= args.geometry_through_height
                ),
                args.geometry_engine,
                args.geometry_ray_limit,
                args.geometry_time_limit_seconds,
            )
            for dimension in args.dimensions
            for seed in seeds
            for height in heights
        ],
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
