#!/usr/bin/env python3
"""Freeze the exact H=10000 necessary-wall audit for a=3+8sqrt(2)."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from coxeter_polytopes.exact.silver import Silver  # noqa: E402
from coxeter_polytopes.search.silver_leech_norm119 import (  # noqa: E402
    DEFAULT_HEIGHT_BOUND,
    norm119_leech_prefix_audit,
    strongest_norm119_prefix_witness,
)


SOURCE_PATHS = (
    "src/coxeter_polytopes/exact/silver.py",
    "src/coxeter_polytopes/search/silver_e8_cubes.py",
    "src/coxeter_polytopes/search/silver_leech.py",
    "src/coxeter_polytopes/search/silver_leech_manifest.py",
    "src/coxeter_polytopes/search/silver_leech_norm119.py",
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


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--height-bound", type=int, default=DEFAULT_HEIGHT_BOUND
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/silver_leech_norm119_h10000.json"),
    )
    args = parser.parse_args()
    audit = norm119_leech_prefix_audit(args.height_bound)
    strongest = strongest_norm119_prefix_witness(args.height_bound)
    payload = {
        "schema": "coxeter.silver_leech_norm119_necessary_prefix.v1",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 24,
        "height_convention": "H=a*x^2/q at the distinguished embedding",
        "source_sha256": {
            relative: file_sha256(REPOSITORY / relative)
            for relative in SOURCE_PATHS
        },
        "audit": audit,
        "strongest_fixed_direction_witness": strongest,
        "bounded_conclusion": (
            "A fixed norm-six Leech direction, scaled by the archived exact "
            "witness coordinate, is strictly hyperideal and satisfies every "
            "necessary mirror inequality in the complete H<=10000 arithmetic "
            "control prefix."
        ),
        "scope": (
            "The control packets are necessary supersets: Gram triples need "
            "not be realized Leech pairs, primitive roots, or Vinberg-accepted "
            "walls. Surviving a necessary superset certifies only this bounded "
            "prefix obstruction; it is not a global nonreflectivity theorem."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(encode(payload), indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
