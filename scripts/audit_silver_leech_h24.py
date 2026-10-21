#!/usr/bin/env python3
"""Freeze exact early-packet audits for three silver Leech H24 seeds."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from fractions import Fraction
import argparse
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from coxeter_polytopes.exact.silver import Silver  # noqa: E402
from coxeter_polytopes.search.silver_leech import (  # noqa: E402
    SQRT2_LEECH_SEED,
    TWO_SQRT2_LEECH_SEED,
    UNIT_LEECH_SEED,
    P7_LEECH_SEED,
    nonunit_leech_second_packet_audit,
    p7_leech_first_packet_audit,
    p7_leech_post_first_prefix_h82,
    p7_leech_minimal_direction_prefix_h82,
    p7_leech_second_packet_audit,
    sqrt2_leech_first_packet_audit,
    two_sqrt2_leech_first_packet_audit,
    unit_leech_first_packet_audit,
    unit_leech_minimal_direction_prefix_h82,
    unit_leech_second_packet_audit,
)


def encode(value: object) -> object:
    if isinstance(value, Silver):
        return value.to_json()
    if isinstance(value, Fraction):
        return str(value)
    if is_dataclass(value):
        return encode(asdict(value))
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def seed_record(seed: object) -> dict[str, object]:
    return {
        "seed_id": seed.seed_id,
        "timelike_coefficient": seed.timelike_coefficient,
        "canonical_mirror_norms": seed.canonical_mirror_norms,
        "complete_norm_class_proof": seed.complete_norm_class_proof,
        "admissible": seed.admissible,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/silver_leech_h24_h82.json"),
    )
    args = parser.parse_args()
    payload = {
        "schema": "coxeter.silver_leech_h24_early_packets.v1",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 24,
        "unit_seed": {
            "seed": seed_record(UNIT_LEECH_SEED),
            "first_packet": unit_leech_first_packet_audit(),
            "second_packet": unit_leech_second_packet_audit(),
            "minimal_direction_bounds": unit_leech_minimal_direction_prefix_h82(),
            "complete_height_bound": "82+58*sqrt(2)",
            "conclusion": (
                "An explicit hyperideal minimal-vector direction survives every "
                "arithmetically possible root through the complete H82 prefix."
            ),
        },
        "sqrt2_seed": {
            "seed": seed_record(SQRT2_LEECH_SEED),
            "first_packet": sqrt2_leech_first_packet_audit(),
            "second_packet": nonunit_leech_second_packet_audit(SQRT2_LEECH_SEED),
        },
        "two_sqrt2_seed": {
            "seed": seed_record(TWO_SQRT2_LEECH_SEED),
            "first_packet": two_sqrt2_leech_first_packet_audit(),
            "second_packet": nonunit_leech_second_packet_audit(TWO_SQRT2_LEECH_SEED),
            "local_gate": (
                "For norm 4+2sqrt(2), reflection integrality requires "
                "the spatial component to lie in sqrt(2)*Lambda."
            ),
        },
        "p7_seed": {
            "seed": seed_record(P7_LEECH_SEED),
            "first_packet": p7_leech_first_packet_audit(),
            "post_first_prefix": p7_leech_post_first_prefix_h82(),
            "second_packet": p7_leech_second_packet_audit(),
            "norm_six_direction_bounds": p7_leech_minimal_direction_prefix_h82(),
            "complete_height_bound": "82+58*sqrt(2)",
            "conclusion": (
                "A strict hyperideal norm-six direction survives every "
                "arithmetically possible wall through the complete H82 prefix."
            ),
            "contact_polar_source": {
                "title": "The contact polytope of the Leech lattice (complete version)",
                "authors": "Mathieu Dutour Sikiric, Achill Schuermann, Frank Vallentin",
                "arxiv": "0906.1427v2",
                "doi": "10.1007/s00454-010-9266-z",
                "used_result": (
                    "232 Co0 vertex orbits in the scaled contact polar; unique "
                    "maximum squared norm 8/3, exceptional Co3 orbit with 552 incidences"
                ),
            },
            "conclusion": (
                "An explicit H=31 root cuts the exceptional contact-polar "
                "vertex orbit, but z=(37/67)d for every archived norm-six "
                "direction is a strict hyperideal witness through every "
                "necessary control level up to H=82+58sqrt(2)."
            ),
        },
        "scope": (
            "Exact mirror-norm proofs and early shell/persistent-ray certificates; "
            "no reflectivity, termination, or novelty claim."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(encode(payload), indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
