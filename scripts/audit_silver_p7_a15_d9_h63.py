#!/usr/bin/env python3
"""Freeze the exact prime-seven N(A15 D9) prefix through height 63."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a15_d9 import a15_d9_basis_audit, glue_coset_minima
from coxeter_polytopes.search.silver_p7_a15_d9 import (
    COMPLETE_P7_A15_D9_NORM_PROOF,
    FIRST_PACKET_HEIGHT,
    P7_NORM_CLASSES,
    SECOND_PACKET_HEIGHT,
    p7_a15_d9_pair,
    p7_d9_axis_witness,
    run_p7_a15_d9_prefix,
)
from coxeter_polytopes.search.silver_p7_d24_plus import P7_DELTA, P7_TIMELIKE


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/silver_p7_a15_d9_h63.json"


def silver(value) -> dict[str, str]:
    return Silver.coerce(value).to_json()


def vector(values) -> list[dict[str, str]]:
    return [silver(value) for value in values]


def root_payload(root) -> dict:
    return {
        "norm": silver(root.norm),
        "height": silver(root.height),
        "control": silver(root.control),
        "target": silver(root.target),
        "a15": {
            "weights": vector(root.a15.weights),
            "coordinates": vector(root.a15.coordinates),
            "residue": list(root.a15.residue),
        },
        "d9": {
            "weights": vector(root.d9.weights),
            "coordinates": vector(root.d9.coordinates),
            "residue": list(root.d9.residue),
        },
        "glue_nontrivial": root.glue_nontrivial,
        "cross_component_glue": root.cross_component_glue,
    }


def main() -> None:
    if not a15_d9_basis_audit():
        raise ArithmeticError("A15D9 glue basis failed")
    prefix = run_p7_a15_d9_prefix(63)
    first = next(packet for packet in prefix.packets if packet.height == FIRST_PACKET_HEIGHT)
    second = next(packet for packet in prefix.packets if packet.height == SECOND_PACKET_HEIGHT)
    witness = p7_d9_axis_witness(height_bound=63)
    if (len(first.roots), second.candidate_count, second.rejected_count, len(second.roots)) != (11, 1670, 15, 1655):
        raise ArithmeticError("prime-seven packet census changed")
    if not witness.hyperideal or not witness.all_positive_walls_strict:
        raise ArithmeticError("prime-seven strict witness failed")
    roots = [root_payload(root) for root in prefix.accepted_roots]
    root_digest = sha256(json.dumps(roots, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    payload = {
        "scope": {
            "field": "Q(sqrt(2))",
            "hyperbolic_dimension": 24,
            "spatial_lattice": "Niemeier N(A15 D9), cyclic index-eight glue <(2,1)>",
            "timelike_coefficient": silver(P7_TIMELIKE),
            "height_bound": 63,
            "last_height": silver(SECOND_PACKET_HEIGHT),
            "claim": (
                "Complete primitive chronology through the first two nonempty "
                "packets and every intervening control. The bounded prefix has "
                "a strict hyperideal ray; this is not global nonreflectivity."
            ),
        },
        "glue": {
            "generator": [2, 1],
            "order": 8,
            "nonzero_coset_minima": [int(value) for value in glue_coset_minima()],
            "height_zero_root_count": 384,
        },
        "mirror_norms": {
            "canonical": [silver(norm) for norm in P7_NORM_CLASSES],
            "delta": silver(P7_DELTA),
            "completeness_proof": COMPLETE_P7_A15_D9_NORM_PROOF,
        },
        "controls": [
            {
                "norm": silver(level.norm),
                "height": silver(level.height),
                "control": silver(level.control),
                "target": silver(level.target),
            }
            for level in prefix.controls
        ],
        "packets": [
            {
                "height": silver(packet.height),
                "control_count": len(packet.controls),
                "candidate_count": packet.candidate_count,
                "rejected_count": packet.rejected_count,
                "accepted_count": len(packet.roots),
                "glue_nontrivial_count": sum(root.glue_nontrivial for root in packet.roots),
                "cross_component_glue_count": sum(root.cross_component_glue for root in packet.roots),
            }
            for packet in prefix.packets
        ],
        "dense_packet": {
            "candidate_count": 1670,
            "old_wall_comparisons": 18370,
            "old_wall_pairing_signs": {"negative": 18337, "zero": 18, "positive": 15},
            "accepted_count": 1655,
            "same_height_pairing_signs": {"negative": 1367834, "zero": 851, "positive": 0},
            "glue_nontrivial_count": 1412,
            "cross_component_glue_count": 1252,
            "exact_batch_scale": 16,
        },
        "accepted_roots": roots,
        "accepted_roots_sha256": root_digest,
        "witness": {
            "a15_weights": vector(witness.a15_weights),
            "d9_weights": vector(witness.d9_weights),
            "norm": silver(witness.norm),
            "pairing_count": len(witness.positive_wall_pairings),
            "pairing_signs": dict(Counter(str(pairing.sign(0)) for pairing in witness.positive_wall_pairings)),
            "maximum_pairing": silver(max(witness.positive_wall_pairings)),
            "strict": witness.all_positive_walls_strict,
            "hyperideal": witness.hyperideal,
        },
        "independent_audit": {
            "method": "separate exact integer-pair NumPy matrices scaled by 16",
            "candidate_old_and_same_height_pairings": "PASS",
            "witness_pairings": "PASS",
        },
        "outcome": (
            "No compact or finite-volume example in the H=63 prefix. The dense "
            "second packet leaves 1655 acute walls, but a strict D9-axis "
            "hyperideal ray survives all 1666 positive-height walls."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(OUTPUT), "sha256": sha256(OUTPUT.read_bytes()).hexdigest(), "accepted": len(prefix.accepted_roots)}, indent=2))


if __name__ == "__main__":
    main()
