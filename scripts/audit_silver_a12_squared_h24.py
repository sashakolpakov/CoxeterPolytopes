#!/usr/bin/env python3
"""Freeze the exact silver-unit N(A12^2) prefix through H=150."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a12_squared import (
    COMPLETE_UNIT_NORM_PROOF,
    EPSILON,
    a12_squared_basis_audit,
    first_packet_geometry,
    height_zero_dominant_roots,
    pair,
    run_prefix,
    strict_hyperideal_witness,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/silver_a12_squared_h24_h150.json"


def silver(value):
    return Silver.coerce(value).to_json()


def vector(values):
    return [silver(value) for value in values]


def root_payload(root):
    return {
        "norm": silver(root.norm),
        "height": silver(root.height),
        "control": silver(root.control),
        "target": silver(root.target),
        "left": {"square": silver(root.left.square), "weights": vector(root.left.weights), "coordinates": vector(root.left.coordinates), "residue": list(root.left.residue)},
        "right": {"square": silver(root.right.square), "weights": vector(root.right.weights), "coordinates": vector(root.right.coordinates), "residue": list(root.right.residue)},
    }


def main():
    if not a12_squared_basis_audit():
        raise ArithmeticError("A12^2 basis audit failed")
    prefix = run_prefix(150, None)
    nonempty = tuple(packet for packet in prefix.packets if packet.roots)
    counts = tuple(len(packet.roots) for packet in nonempty)
    if counts != (2, 4, 20, 30, 164, 28, 32):
        raise ArithmeticError(f"A12^2 packet counts changed: {counts}")
    if any(root.glue_nontrivial for root in prefix.accepted_roots):
        raise ArithmeticError("an unexpected prime-13 glue root appeared")
    if any(pair(left, right).sign(0) > 0 for i,left in enumerate(prefix.accepted_roots) for right in prefix.accepted_roots[i+1:]):
        raise ArithmeticError("accepted A12^2 roots are not mutually acute")
    witness = strict_hyperideal_witness()
    if not witness.hyperideal or not witness.all_positive_walls_strict:
        raise ArithmeticError("A12^2 strict witness failed")
    roots = [root_payload(root) for root in prefix.accepted_roots]
    geometry = first_packet_geometry()
    payload = {
        "scope": {
            "field": "Q(sqrt(2))",
            "hyperbolic_dimension": 24,
            "spatial_lattice": "Niemeier N(A12^2), prime-13 graph glue {(k,8k)}",
            "timelike_coefficient": silver(EPSILON),
            "height_bound": 150,
            "claim": "Complete primitive chronological prefix through real H=150; strict hyperideal witness; not global nonreflectivity.",
        },
        "glue": {
            "root_lattice_discriminant": 169,
            "glue_index": 13,
            "generator": ["omega_5", "omega_1"],
            "residue_graph_multiplier": 8,
            "generator_square": "4",
            "explicit_basis_even_unimodular": True,
            "height_zero_dominant_orbits": len(height_zero_dominant_roots()),
            "height_zero_root_count": 312,
        },
        "mirror_norms": {"canonical": [silver(1), silver(2)], "completeness_proof": COMPLETE_UNIT_NORM_PROOF},
        "packets": [
            {"height": silver(packet.height), "control_count": len(packet.controls), "candidate_count": packet.candidate_count, "rejected_count": packet.rejected_count, "accepted_count": len(packet.roots), "prime13_glue_count": sum(root.glue_nontrivial for root in packet.roots)}
            for packet in prefix.packets
        ],
        "accepted_roots": roots,
        "accepted_roots_sha256": sha256(json.dumps(roots,sort_keys=True,separators=(",",":")).encode()).hexdigest(),
        "first_affine_product_geometry": {"rays": geometry.extreme_ray_count, "finite": geometry.finite_count, "cusps": geometry.cusp_count, "hyperideal": geometry.hyperideal_count},
        "strict_hyperideal_witness": {"dynkin_index_zero_based": witness.dynkin_index, "weight": str(witness.weight), "norm": silver(witness.norm), "all_positive_walls_strict": witness.all_positive_walls_strict, "hyperideal": witness.hyperideal, "closest_pairing": silver(max(witness.pairings)), "pairing_histogram": dict(sorted(Counter(str(value) for value in witness.pairings).items()))},
        "outcome": "No compact/FV example: 280 accepted walls through H=150, all in root-lattice residue; a strict hyperideal ray survives.",
    }
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"output":str(OUTPUT),"sha256":sha256(OUTPUT.read_bytes()).hexdigest(),"accepted":len(prefix.accepted_roots)},indent=2))


if __name__ == "__main__":
    main()
