#!/usr/bin/env python3
"""Freeze the exact silver ``N(A2^12)`` H24 screen through height ten."""

from __future__ import annotations

import argparse
import json
from dataclasses import fields
from pathlib import Path
from typing import Any

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a2_ternary_golay import (
    TERNARY_GOLAY_GENERATOR,
    all_seed_searches,
    niemeier_construction_is_even_unimodular,
    ternary_golay_is_self_dual,
    ternary_golay_weight_enumerator,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/enumerations/silver_a2_ternary_golay_h24_h10.json"


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
        return {field.name: _json(getattr(value, field.name)) for field in fields(value)}
    return value


def build_payload() -> dict[str, Any]:
    searches = all_seed_searches()
    return {
        "schema": "silver-a2-ternary-golay-h24-h10-v1",
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
        "records": [_json(search) for search in searches],
        "outcome": {
            "compact_survivor_count": sum(search.geometry.compact for search in searches),
            "finite_volume_survivor_count": sum(
                search.geometry.finite_volume for search in searches
            ),
            "statement": (
                "Both complete H<=10 prefixes have 36 walls and exact product-of-"
                "triangles geometry with 531416 hyperideal rays. No compact or "
                "finite-volume survivor occurs at this bounded cutoff."
            ),
            "scope": (
                "This is exhaustive only for the two displayed hosts through H=10. "
                "Later mirrors may cut every current hyperideal ray; no global "
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
