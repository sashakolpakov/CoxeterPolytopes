#!/usr/bin/env python3
"""Emit a replayable exact checkpoint for the rank-five L2/L7 search."""

from __future__ import annotations

import argparse
import json
from hashlib import sha256
from pathlib import Path

from coxeter_polytopes.search.cubic7_nonunimodular import (
    L2_SEED,
    L7_SEED,
    P13_R10_SEED,
    P13_R7_SEED,
    P13_R8_SEED,
    P13_SEEDS,
    enumerate_cone_rays_incremental,
    enumerate_height_zero_mirrors,
    nonunimodular_diagonal_pair,
    run_nonunimodular_vinberg_search,
)


def _element(value):
    return value.to_json()


def _vector(value):
    return [_element(entry) for entry in value]


def _root(value):
    return {
        "norm": _element(value.norm),
        "height": _element(value.height),
        "vector": _vector(value.vector),
    }


def checkpoint(seed, height_bound):
    search = run_nonunimodular_vinberg_search(seed, height_bound)
    ray_census = enumerate_cone_rays_incremental(search)
    rays = ray_census.rays
    by_norm = []
    for norm in seed.canonical_mirror_norms:
        by_norm.append(
            {
                "norm": _element(norm),
                "control_count": dict(search.pool.control_counts_by_norm)[norm],
                "candidate_count": sum(root.norm == norm for root in search.pool.roots),
                "accepted_count": sum(root.norm == norm for root in search.accepted_roots),
            }
        )
    return {
        "seed_id": seed.seed_id,
        "diagonal_form": [_element(entry) for entry in seed.coefficients],
        "determinant_field_norm_abs": abs(int(seed.last_coefficient.norm())),
        "canonical_mirror_norms": [_element(norm) for norm in seed.canonical_mirror_norms],
        "local_norm_exhaustiveness_proof": seed.local_norm_proof,
        "height_zero": {
            "mirror_count": len(enumerate_height_zero_mirrors(seed)),
            "weyl_type": "B3xA1",
            "simple_roots": [_vector(root) for root in search.pool.cone_roots],
        },
        "height_bound": height_bound,
        "counts_by_norm": by_norm,
        "positive_height_candidate_count": len(search.pool.roots),
        "accepted_positive_height_count": len(search.accepted_roots),
        "rejected_positive_height_count": len(search.rejected_roots),
        "wall_count": len(search.walls),
        "accepted_roots": [_root(root) for root in search.accepted_roots],
        "wall_norms": [
            _element(nonunimodular_diagonal_pair(seed, root, root)) for root in search.walls
        ],
        "geometry": {
            "finite_ray_count": sum(ray.kind == "finite" for ray in rays),
            "cusp_count": sum(ray.kind == "cusp" for ray in rays),
            "hyperideal_ray_count": sum(ray.kind == "hyperideal" for ray in rays),
            # One exact obstruction is enough to certify that this prefix has
            # not closed; all rays are cheaply recomputed from the archived
            # simple roots by the independent verifier.
            "hyperideal_witnesses": [
                {
                    "active_facets": list(ray.active_facets),
                    "norm": _element(ray.norm),
                    "vector": _vector(ray.vector),
                }
                for ray in rays
                if ray.kind == "hyperideal"
            ][:1],
        },
        "status": (
            "vinberg_stopping_fired"
            if not any(ray.kind == "hyperideal" for ray in rays)
            else "exhaustive_prefix_not_closed"
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=100)
    parser.add_argument(
        "--heights",
        help="comma-separated adaptive checkpoints; overrides --height",
    )
    parser.add_argument(
        "--summary-input",
        type=Path,
        help="freeze externally completed exact summary records without rerunning",
    )
    parser.add_argument(
        "--seed",
        choices=("L2", "L7", "P13-r7", "P13-r10", "P13-r8", "p13", "both", "all"),
        default="both",
    )
    args = parser.parse_args()
    seed_lookup = {
        "L2": (L2_SEED,),
        "L7": (L7_SEED,),
        "P13-r7": (P13_R7_SEED,),
        "P13-r10": (P13_R10_SEED,),
        "P13-r8": (P13_R8_SEED,),
        "p13": P13_SEEDS,
        "both": (L2_SEED, L7_SEED),
        "all": (L2_SEED, L7_SEED, *P13_SEEDS),
    }
    if args.summary_input is not None:
        summaries = json.loads(args.summary_input.read_text(encoding="utf-8"))
        if not isinstance(summaries, list):
            raise ValueError("--summary-input must contain a JSON list")
        payload = {
            "schema": "cubic7-nonunimodular-vinberg-summary-0.1.0",
            "field": {
                "generator": "theta=2*cos(2*pi/7)",
                "minimal_polynomial": "theta^3+theta^2-2*theta-1",
                "ring_of_integers": "Z[theta]",
            },
            "scope": (
                "exact run summaries retained after the full in-memory root pools were "
                "released; replay the listed command to regenerate literal roots"
            ),
            "runs": summaries,
        }
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        payload["payload_sha256"] = sha256(canonical.encode()).hexdigest()
        print(json.dumps(payload, indent=2, sort_keys=True))
        return

    seeds = seed_lookup[args.seed]
    heights = (
        tuple(int(value) for value in args.heights.split(","))
        if args.heights
        else (args.height,)
    )
    payload = {
        "schema": "cubic7-nonunimodular-vinberg-checkpoint-0.1.0",
        "field": {
            "generator": "theta=2*cos(2*pi/7)",
            "minimal_polynomial": "theta^3+theta^2-2*theta-1",
            "ring_of_integers": "Z[theta]",
        },
        "normalization": (
            "primitive reflecting hyperplanes modulo multiplication by units; "
            "totally positive units are squares"
        ),
        "runs": [checkpoint(seed, height) for seed in seeds for height in heights],
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["payload_sha256"] = sha256(canonical.encode()).hexdigest()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
