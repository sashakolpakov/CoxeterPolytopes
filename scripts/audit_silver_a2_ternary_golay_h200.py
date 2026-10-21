#!/usr/bin/env python3
"""Freeze the silver ``N(A2^12)`` necessary frontier through real H=200."""

from __future__ import annotations

import argparse
import json
from dataclasses import fields
from pathlib import Path
from typing import Any

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a2_ternary_golay import (
    SUPPORTED_A2_GOLAY_SEEDS,
    first_h24_cutter,
    frontier_censuses_through_h200,
    h200_uniform_witness,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/enumerations/silver_a2_ternary_golay_h200.json"


def _json(value: Any) -> Any:
    if isinstance(value, Silver):
        return value.to_json()
    if isinstance(value, tuple):
        return [_json(entry) for entry in value]
    if isinstance(value, list):
        return [_json(entry) for entry in value]
    if isinstance(value, dict):
        return {str(key): _json(entry) for key, entry in value.items()}
    if hasattr(value, "__dataclass_fields__"):
        return {
            field.name: _json(getattr(value, field.name))
            for field in fields(value)
        }
    return value


def build_payload() -> dict[str, Any]:
    records = []
    for seed in SUPPORTED_A2_GOLAY_SEEDS:
        censuses = frontier_censuses_through_h200(seed)
        records.append(
            {
                "seed_id": seed.seed_id,
                "timelike_coefficient": _json(seed.timelike_coefficient),
                "complete_mirror_norms": _json(seed.canonical_mirror_norms),
                "control_count_including_first_packet": len(censuses) + 1,
                "post_affine_censuses": _json(censuses),
                "nonempty_packet_count": sum(
                    census.candidate_count > 0 for census in censuses
                ),
                "explicit_h24_accepted_cutter": _json(first_h24_cutter(seed)),
                "strict_uniform_witness": _json(h200_uniform_witness(seed)),
            }
        )
    return {
        "schema": "silver-a2-ternary-golay-h24-h200-v1",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 24,
        "spatial_lattice": "Niemeier N(A2^12)",
        "real_height_bound": 200,
        "records": records,
        "outcome": {
            "compact_survivor_count": 0,
            "finite_volume_survivor_count": 0,
            "statement": (
                "Exactly eight post-affine necessary levels are nonempty in "
                "each host through real H=200. A strict uniform hyperideal "
                "point survives every vector in every complete candidate "
                "superset, hence also all chronologically accepted walls."
            ),
            "scope": (
                "This is a complete bounded necessary-candidate frontier, not "
                "an enumeration of every accepted packet and not a global "
                "nonreflectivity theorem."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = build_payload()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
