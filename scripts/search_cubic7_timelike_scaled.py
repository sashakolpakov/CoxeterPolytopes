#!/usr/bin/env python3
"""Freeze replayable exact prefixes for the T2/T7 timelike-scaled search."""

from __future__ import annotations

import argparse
import json
from hashlib import sha256
from pathlib import Path

from coxeter_polytopes.search.cubic7_timelike import (
    MINKOWSKI_BOUND_DENOMINATOR,
    MINKOWSKI_BOUND_NUMERATOR,
    T2_TIMELIKE_SEED,
    T7_TIMELIKE_SEED,
    UNIT_SIGNATURE_BITS,
    UNIT_SIGNATURE_GENERATORS,
    enumerate_height_zero_mirrors,
    enumerate_timelike_cone_rays_incremental,
    run_timelike_vinberg_search,
    timelike_diagonal_pair,
    unit_signature_rank,
)
from coxeter_polytopes.verify.cubic7_alvin import parse_alvin_roots


WORKSPACE = Path(__file__).resolve().parents[1]
LIVE_TRANSCRIPTS = {
    "T2": WORKSPACE
    / "data/enumerations/k7_alvin/pinned_rc7_d92006a_fractionfix_20260812/"
    "timelike-minus-2theta-rank5/stdout.txt",
    "T7": WORKSPACE
    / "data/enumerations/k7_alvin/pinned_rc7_d92006a_fractionfix_20260812/"
    "timelike-minus-theta-times-2-minus-theta-rank5/stdout.txt",
}


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


def _canonical_hash(value) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def checkpoint(seed, height_bound):
    search = run_timelike_vinberg_search(seed, height_bound)
    geometry = enumerate_timelike_cone_rays_incremental(search)
    counts_by_norm = []
    control_counts = dict(search.pool.control_counts_by_norm)
    spatial_counts = dict(search.pool.spatial_alphabet_counts_by_norm)
    for norm in seed.canonical_mirror_norms:
        counts_by_norm.append(
            {
                "norm": _element(norm),
                "control_count": control_counts[norm],
                "spatial_alphabet_count": spatial_counts[norm],
                "candidate_count": sum(root.norm == norm for root in search.pool.roots),
                "accepted_count": sum(root.norm == norm for root in search.accepted_roots),
            }
        )
    record = {
        "seed_id": seed.seed_id,
        "diagonal_form": [_element(entry) for entry in seed.coefficients],
        "determinant_ideal_norm": abs(int(seed.timelike_scale.norm())),
        "canonical_mirror_norms": [_element(norm) for norm in seed.canonical_mirror_norms],
        "complete_norm_class_proof": seed.complete_norm_class_proof,
        "height_zero": {
            "mirror_count": len(enumerate_height_zero_mirrors(seed)),
            "weyl_type": "B4",
            "simple_roots": [_vector(root) for root in search.pool.cone_roots],
        },
        "height_bound": height_bound,
        "counts_by_norm": counts_by_norm,
        "positive_height_candidate_count": len(search.pool.roots),
        "accepted_positive_height_count": len(search.accepted_roots),
        "rejected_positive_height_count": len(search.rejected_roots),
        "wall_count": len(search.walls),
        "accepted_roots": [_root(root) for root in search.accepted_roots],
        "wall_norms": [
            _element(timelike_diagonal_pair(seed, wall, wall)) for wall in search.walls
        ],
        "geometry": {
            "extreme_ray_count": len(geometry.rays),
            "finite_ray_count": geometry.finite_count,
            "cusp_count": geometry.cusp_count,
            "hyperideal_ray_count": geometry.hyperideal_count,
            "irredundant": geometry.irredundant,
            "consistent_time_sheet": geometry.consistent_time_sheet,
            "finite_volume": geometry.finite_volume,
            "compact": geometry.compact,
            "hyperideal_witnesses": [
                {
                    "active_facets": list(ray.active_facets),
                    "norm": _element(ray.norm),
                    "vector": _vector(ray.vector),
                }
                for ray in geometry.rays
                if ray.kind == "hyperideal"
            ][:1],
        },
        "status": (
            "vinberg_stopping_fired"
            if geometry.finite_volume
            else "exhaustive_prefix_not_closed"
        ),
    }
    record["checkpoint_sha256"] = _canonical_hash(record)
    return record, search


def live_crosscheck(seed, search):
    transcript = LIVE_TRANSCRIPTS[seed.seed_id]
    if not transcript.exists():
        return {
            "available": False,
            "path": str(transcript.relative_to(WORKSPACE)),
        }
    live_roots = parse_alvin_roots(transcript)
    live_packet = tuple(live_roots[4:9])
    exact_packet = tuple(root.vector for root in search.accepted_roots[:5])
    return {
        "available": True,
        "path": str(transcript.relative_to(WORKSPACE)),
        "stdout_sha256": sha256(transcript.read_bytes()).hexdigest(),
        "live_root_indices": [5, 6, 7, 8, 9],
        "comparison": "unordered exact hyperplane packet after RC7-to-power-basis conversion",
        "exact_match": set(live_packet) == set(exact_packet),
        "live_packet": [_vector(root) for root in live_packet],
        "note": (
            "The live job was externally interrupted. This certifies only the "
            "early-wall agreement, not AlVin termination or completeness."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--heights", default="20,50,100,200,300")
    parser.add_argument("--seed", choices=("T2", "T7", "both"), default="both")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    heights = tuple(int(value) for value in args.heights.split(","))
    if not heights or any(height < 0 for height in heights):
        raise SystemExit("--heights must be a nonempty comma-separated list of nonnegative integers")
    seeds = {
        "T2": (T2_TIMELIKE_SEED,),
        "T7": (T7_TIMELIKE_SEED,),
        "both": (T2_TIMELIKE_SEED, T7_TIMELIKE_SEED),
    }[args.seed]

    runs = []
    searches = {}
    for height in heights:
        for seed in seeds:
            record, search = checkpoint(seed, height)
            runs.append(record)
            searches[(seed.seed_id, height)] = search

    comparison_height = min((height for height in heights if height >= 50), default=max(heights))
    payload = {
        "schema": "coxeter.cubic7.timelike_scaled_exact_prefixes.v1",
        "field": {
            "generator": "theta=2*cos(2*pi/7)",
            "minimal_polynomial": "theta^3+theta^2-2*theta-1",
            "ring_of_integers": "Z[theta]",
            "discriminant": 49,
        },
        "normalization_certificate": {
            "minkowski_bound": {
                "numerator": MINKOWSKI_BOUND_NUMERATOR,
                "denominator": MINKOWSKI_BOUND_DENOMINATOR,
                "less_than_two": MINKOWSKI_BOUND_NUMERATOR
                < 2 * MINKOWSKI_BOUND_DENOMINATOR,
                "conclusion": "ordinary class number one",
            },
            "unit_signature_generators": [_element(unit) for unit in UNIT_SIGNATURE_GENERATORS],
            "unit_signature_bits": [list(bits) for bits in UNIT_SIGNATURE_BITS],
            "unit_signature_rank_over_f2": unit_signature_rank(),
            "conclusion": (
                "narrow class number one; every totally positive unit is a square"
            ),
        },
        "height_semantics": (
            "Every primitive reflecting hyperplane, normalized modulo unit squares "
            "to a canonical mirror norm and having Vinberg height <= H, is enumerated."
        ),
        "same_height_rule": (
            "Filter a full exact height packet against lower-height walls, prove all "
            "survivors pairwise acute, then adjoin the packet simultaneously."
        ),
        "runs": runs,
        "live_alvin_early_wall_crosschecks": [
            live_crosscheck(seed, searches[(seed.seed_id, comparison_height)]) for seed in seeds
        ],
        "source_sha256": {
            "backend": sha256(
                (WORKSPACE / "src/coxeter_polytopes/search/cubic7_timelike.py").read_bytes()
            ).hexdigest(),
            "script": sha256(Path(__file__).read_bytes()).hexdigest(),
            "tests": sha256(
                (WORKSPACE / "tests/test_cubic7_timelike_search.py").read_bytes()
            ).hexdigest(),
        },
    }
    payload["payload_sha256"] = _canonical_hash(payload)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
