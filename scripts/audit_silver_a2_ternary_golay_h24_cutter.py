#!/usr/bin/env python3
"""Freeze the silver ``N(A2^12)`` frontier through H=24+17*sqrt(2)."""

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
    frontier_censuses_through_h24,
    h24_candidate_census,
    h24_uniform_witness,
    ternary_golay_pair_contingencies,
    ternary_golay_two_pointed_pair_contingencies,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/silver_a2_ternary_golay_h24_cutter.json"
)


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
    ordinary_contingencies = ternary_golay_pair_contingencies()
    pointed_contingencies = ternary_golay_two_pointed_pair_contingencies()
    records = []
    for seed in SUPPORTED_A2_GOLAY_SEEDS:
        records.append(
            {
                "seed_id": seed.seed_id,
                "timelike_coefficient": _json(seed.timelike_coefficient),
                "complete_mirror_norms": _json(seed.canonical_mirror_norms),
                "post_affine_control_censuses": _json(
                    frontier_censuses_through_h24(seed)
                ),
                "h24_candidate_superset": _json(h24_candidate_census(seed)),
                "explicit_accepted_cutter": _json(first_h24_cutter(seed)),
                "strict_uniform_witness": _json(h24_uniform_witness(seed)),
            }
        )
    return {
        "schema": "silver-a2-ternary-golay-h24-cutter-v1",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 24,
        "spatial_lattice": "Niemeier N(A2^12)",
        "glue_dynamic_program": {
            "ordinary_contingency_type_count": len(ordinary_contingencies),
            "ordinary_total_multiplicity": sum(
                multiplicity for _, multiplicity in ordinary_contingencies
            ),
            "two_pointed_contingency_type_count": len(pointed_contingencies),
            "two_pointed_total_multiplicity": sum(
                multiplicity for _, _, multiplicity in pointed_contingencies
            ),
            "h8_same_height_pairing_census": {
                "unordered_pair_count": 7_839_180,
                "negative": 7_839_114,
                "zero": 66,
                "positive": 0,
                "maximum": Silver(0).to_json(),
                "method": (
                    "independent exact 3960-by-3960 integral coefficient replay"
                ),
            },
        },
        "records": records,
        "outcome": {
            "compact_survivor_count": 0,
            "finite_volume_survivor_count": 0,
            "statement": (
                "Both hosts have an accepted H=24+17sqrt(2) cutter, but an "
                "exact strict hyperideal point survives every member of the "
                "complete 4,733,256-vector first-cap candidate superset."
            ),
            "scope": (
                "This certifies the complete bounded prefix only. It does not "
                "enumerate the final accepted subset at the endpoint and makes "
                "no global nonreflectivity claim."
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
