#!/usr/bin/env python3
"""Write the exact norm-128 coefficient-direction first-prefix manifest."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from coxeter_polytopes.exact.silver import Silver  # noqa: E402
from coxeter_polytopes.search.silver_leech_manifest import (  # noqa: E402
    full_manifest_first_packet_obstructions,
    norm89_priority_obstructions,
    scalar_shell_priority_obstructions,
    silver_ideal_atoms,
    silver_ideal_squareclasses,
    silver_leech_ideal_squareclass_manifest,
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/silver_leech_low_norm_manifest.json"),
    )
    args = parser.parse_args()
    atoms = silver_ideal_atoms(128)
    squareclasses = silver_ideal_squareclasses(128)
    manifest = silver_leech_ideal_squareclass_manifest(128)
    payload = {
        "schema": "coxeter.silver_leech_coefficient_direction_manifest.v3",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 24,
        "coefficient_unit_normalization": "multiplication by (1+sqrt(2))^(2k)",
        "coefficient_direction_equivalence": (
            "multiplication by a field square, as a diagonal quadratic-space "
            "search direction only; nonunit squares need not preserve the "
            "integral lattice"
        ),
        "coefficient_norm_bound": 128,
        "coefficient_scope": (
            "all squarefree coefficient-ideal directions of norm <=128, "
            "including ramified, inert, and both split-prime choices; not all "
            "integral lattice isometry classes with nonunit square content"
        ),
        "prime_ideal_atoms": atoms,
        "ideal_squareclasses": squareclasses,
        "coefficient_direction_count": len(manifest),
        "rows": manifest,
        "norm89_priority_obstructions": norm89_priority_obstructions(),
        "new_scalar_shell_obstructions": scalar_shell_priority_obstructions(),
        "all_first_packet_obstructions": full_manifest_first_packet_obstructions(128),
        "scope": (
            "Exact squarefree-ideal coefficient-direction and mirror-class census "
            "plus first necessary Leech Gram packets. Necessary triples need not "
            "be realized or accepted Vinberg roots. Every archived hyperideal "
            "witness concerns this finite first-packet prefix only; later mirrors "
            "may cut it, so no final finite-volume verdict is made."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(encode(payload), indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
