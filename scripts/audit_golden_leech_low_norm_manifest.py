#!/usr/bin/env python3
"""Write the exact norm-128 golden Leech ideal-squareclass manifest."""

from __future__ import annotations

import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from coxeter_polytopes.exact.golden import Golden  # noqa: E402
from coxeter_polytopes.search.golden_leech_manifest import (  # noqa: E402
    CLASS_AND_UNIT_PROOF,
    DEFAULT_HEIGHT_BOUND,
    HEIGHT_ZERO_MINIMUM_TRACE,
    HEIGHT_ZERO_PROOF,
    MIRROR_CLASS_PROOF,
    full_manifest_first_packet_obstructions,
    golden_ideal_atoms,
    golden_ideal_squareclasses,
    golden_leech_ideal_squareclass_manifest,
)


def encode(value: object) -> object:
    if isinstance(value, Golden):
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
        default=Path("data/enumerations/golden_leech_low_norm_manifest.json"),
    )
    args = parser.parse_args()
    atoms = golden_ideal_atoms(128)
    squareclasses = golden_ideal_squareclasses(128)
    manifest = golden_leech_ideal_squareclass_manifest(128)
    obstructions = full_manifest_first_packet_obstructions(128)
    tightest = min(obstructions, key=lambda row: row.limiting_hyperideal_excess)
    payload = {
        "schema": "coxeter.golden_leech_ideal_squareclass_manifest.v1",
        "field": "Q(sqrt(5)) = Q(phi)",
        "ring_of_integers": "Z[phi]",
        "hyperbolic_dimension": 24,
        "host_family": "<-a> direct-sum Lambda_24 over Z[phi]",
        "coefficient_equivalence": "multiplication by phi^(2k)",
        "direction_interpretation": (
            "the squarefree ideal representatives parameterize K-square "
            "coefficient directions; removing a nonunit square changes the "
            "integral lattice in general, so this is not an integral-host "
            "classification"
        ),
        "coefficient_norm_bound": 128,
        "control_height_cap": DEFAULT_HEIGHT_BOUND,
        "coefficient_scope": (
            "all squarefree prime-ideal representatives of K-square directions "
            "of norm <=128, including the inert dyadic ideal, the ramified "
            "prime over 5, inert odd primes, and both choices over every split "
            "rational prime"
        ),
        "class_and_unit_proof": CLASS_AND_UNIT_PROOF,
        "height_zero_minimum_trace": HEIGHT_ZERO_MINIMUM_TRACE,
        "height_zero_proof": HEIGHT_ZERO_PROOF,
        "complete_mirror_class_proof": MIRROR_CLASS_PROOF,
        "prime_ideal_atoms": atoms,
        "ideal_squareclasses": squareclasses,
        "coefficient_squareclass_count": len(manifest),
        "rows": manifest,
        "all_first_packet_obstructions": obstructions,
        "tightest_limiting_obstruction": tightest,
        "scope": (
            "Exact bounded census of squarefree K-square coefficient directions, "
            "complete primitive local mirror-norm classes, and first necessary "
            "Leech Gram packets. Necessary triples need not be realised Leech "
            "pairs or accepted Vinberg roots. Nonunit square multiples can give "
            "different integral hosts, and later mirrors can cut every displayed "
            "patch. This is not a census of integral hosts and does not prove "
            "global nonreflectivity."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(encode(payload), indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
