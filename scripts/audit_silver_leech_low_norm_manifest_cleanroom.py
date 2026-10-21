#!/usr/bin/env python3
"""Write the independent silver-Leech norm-128 first-prefix audit."""

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
from coxeter_polytopes.search.silver_leech_manifest_cleanroom import (  # noqa: E402
    EXPLICIT_DIRECTION_VECTORS_SCALED,
    HEIGHT_BOUND,
    clean_atoms,
    clean_directions,
    clean_has_short_combination,
    clean_manifest,
    clean_prefix_witnesses,
    mog_scaled_square,
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
        default=Path(
            "data/enumerations/silver_leech_low_norm_manifest_cleanroom.json"
        ),
    )
    args = parser.parse_args()

    atoms = clean_atoms()
    directions = clean_directions()
    rows = clean_manifest()
    witnesses = clean_prefix_witnesses()

    if (len(atoms), len(directions), len(rows), len(witnesses)) != (32, 57, 57, 57):
        raise AssertionError("the independent norm-128 census changed size")
    if any(row.first_necessary_height > HEIGHT_BOUND for row in rows):
        raise AssertionError("a first packet escaped the exact control cutoff")
    if any(
        clean_has_short_combination(triple)
        for row in rows
        for component in row.first_packet_components
        for triple in component.surviving_triples
    ):
        raise AssertionError("a forbidden short binary combination survived")
    if any(witness.minimum_wall_margin.sign(0) <= 0 for witness in witnesses):
        raise AssertionError("a first-prefix wall margin is not strict")
    if any(witness.witness_hyperideal_excess.sign(0) <= 0 for witness in witnesses):
        raise AssertionError("a first-prefix witness is not hyperideal")

    payload = {
        "schema": "coxeter.silver_leech_manifest_cleanroom.v1",
        "field": "Q(sqrt(2))",
        "hyperbolic_dimension": 24,
        "owner_backend_imported": False,
        "coefficient_direction_scope": (
            "squarefree coefficient-ideal directions of ideal norm <=128, "
            "modulo field squares as diagonal quadratic spaces; this does not "
            "identify integral lattices differing by nonunit square content"
        ),
        "local_scope": (
            "complete necessary mirror-ideal classes and odd spatial gates for "
            "the squarefree directions, including the ramified dyadic atom, "
            "inert rational-content atoms, and both split-prime atoms"
        ),
        "packet_scope": (
            "complete first necessary binary Leech Gram packet through exact "
            f"height bound {HEIGHT_BOUND}; triples need not be Leech-realized "
            "or accepted as Vinberg walls"
        ),
        "volume_scope": (
            "each finite prefix through its first necessary packet has a strict "
            "open hyperideal patch; later mirrors may cut the patch, so this is "
            "not a final finite-volume or reflectivity verdict"
        ),
        "prime_ideal_atoms": atoms,
        "directions": directions,
        "rows": rows,
        "first_prefix_witnesses": witnesses,
        "explicit_mog_direction_vectors_scaled": EXPLICIT_DIRECTION_VECTORS_SCALED,
        "explicit_mog_direction_squares": tuple(
            mog_scaled_square(vector) for vector in EXPLICIT_DIRECTION_VECTORS_SCALED
        ),
        "counts": {
            "prime_ideal_atoms": len(atoms),
            "coefficient_directions": len(directions),
            "first_packets": len(rows),
            "strict_first_prefix_witnesses": len(witnesses),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(encode(payload), indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
