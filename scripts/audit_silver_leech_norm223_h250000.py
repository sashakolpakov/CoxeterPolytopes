#!/usr/bin/env python3
"""Freeze the exact H=250000 necessary prefix for a=-13+14sqrt(2)."""

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
from coxeter_polytopes.search.silver_leech_norm223 import (  # noqa: E402
    DEFAULT_HEIGHT_BOUND,
    norm223_leech_prefix_audit,
)


SOURCE_PATHS = (
    "src/coxeter_polytopes/exact/silver.py",
    "src/coxeter_polytopes/search/silver_leech.py",
    "src/coxeter_polytopes/search/silver_leech_manifest.py",
    "src/coxeter_polytopes/search/silver_leech_norm223.py",
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
    parser.add_argument("--height-bound", type=int, default=DEFAULT_HEIGHT_BOUND)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/silver_leech_norm223_h250000.json"),
    )
    args = parser.parse_args()
    audit = norm223_leech_prefix_audit(args.height_bound)
    payload = {
        "schema": "coxeter.silver_leech_norm223_necessary_prefix.v1",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 24,
        "height_convention": "H=a*x^2/q at the distinguished embedding",
        "source_sha256": {
            relative: file_sha256(REPOSITORY / relative)
            for relative in SOURCE_PATHS
        },
        "audit": audit,
        "bounded_conclusion": (
            "A fixed norm-four Leech direction, scaled by the archived exact "
            "witness coordinate, is strictly hyperideal and satisfies every "
            "necessary mirror inequality in the complete H<=250000 arithmetic "
            "control prefix."
        ),
        "scope": (
            "The controls form a necessary super-set. Binary Gram triples are "
            "not asserted to be realized Leech pairs, primitive roots, or "
            "Vinberg-accepted walls. The limiting H=114998+81316*sqrt(2) "
            "component is therefore a necessary-superset limiter/threat, not "
            "a realized cutter. This is a bounded-prefix obstruction, not a "
            "compactness, reflectivity, or global nonreflectivity theorem."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(encode(payload), indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
