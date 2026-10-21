#!/usr/bin/env python3
"""Freeze the exact silver ``N(A2^12)`` H24 screen through real height 17."""

from __future__ import annotations

import argparse
import json
from dataclasses import fields
from pathlib import Path
from typing import Any

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a2_ternary_golay import (
    TERNARY_GOLAY_GENERATOR,
    all_seed_h17_searches,
    niemeier_construction_is_even_unimodular,
    ternary_golay_is_self_dual,
    ternary_golay_pair_contingencies,
    ternary_golay_weight_enumerator,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/silver_a2_ternary_golay_h24_h17.json"
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
    searches = all_seed_h17_searches()
    contingencies = ternary_golay_pair_contingencies()
    return {
        "schema": "silver-a2-ternary-golay-h24-h17-v1",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 24,
        "spatial_lattice": "Niemeier N(A2^12)",
        "construction": {
            "root_lattice": "A2^12",
            "root_lattice_determinant": 3**12,
            "glue_code": "extended ternary Golay [12,6,6]",
            "generator_matrix": [list(row) for row in TERNARY_GOLAY_GENERATOR],
            "code_size": 3**6,
            "weight_enumerator": {
                str(weight): count
                for weight, count in ternary_golay_weight_enumerator().items()
            },
            "self_dual": ternary_golay_is_self_dual(),
            "even_unimodular_lattice": niemeier_construction_is_even_unimodular(),
            "height_zero_root_system": "A2^12",
            "height_zero_root_count": 72,
            "height_zero_simple_walls": 24,
        },
        "glue_pair_contingencies": {
            "type_count": len(contingencies),
            "total_multiplicity": sum(
                multiplicity for _, multiplicity in contingencies
            ),
            "types": [
                {"histogram": list(histogram), "multiplicity": multiplicity}
                for histogram, multiplicity in contingencies
            ],
        },
        "records": [_json(search) for search in searches],
        "outcome": {
            "compact_survivor_count": 0,
            "finite_volume_survivor_count": 0,
            "statement": (
                "Each exact H<=17 prefix has 24 height-zero walls, 12 affine "
                "caps, and a complete 3960-wall norm-two packet at "
                "H=8+6sqrt(2). A strict uniform hyperideal point survives."
            ),
            "scope": (
                "This is exhaustive only for the two displayed hosts through "
                "real height 17. Later mirrors may cut the witness; no global "
                "nonreflectivity claim is made."
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
