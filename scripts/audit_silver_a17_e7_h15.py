#!/usr/bin/env python3
"""Freeze the exact silver-unit N(A17 E7) two-packet prefix."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a17_e7 import (
    COMPLETE_UNIT_NORM_PROOF,
    SECOND_PACKET_HEIGHT,
    a17_e7_basis_audit,
    a17_e7_pair,
    enumerate_dominant_a17_weights,
    enumerate_dominant_e7_weights,
    first_packet_geometry,
    height_zero_dominant_roots,
    run_a17_e7_prefix,
    second_packet_hyperideal_witness,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/silver_a17_e7_h15.json"


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
        "a17": {
            "square": silver(root.a17.square),
            "weights": vector(root.a17.weights),
            "coordinates": vector(root.a17.coordinates),
            "residue": list(root.a17.residue),
        },
        "e7": {
            "square": silver(root.e7.square),
            "weights": vector(root.e7.weights),
            "coordinates": vector(root.e7.coordinates),
            "residue": list(root.e7.residue),
        },
        "glue_nontrivial": root.glue_nontrivial,
        "cross_component_glue": root.cross_component_glue,
    }


def main() -> None:
    if not a17_e7_basis_audit():
        raise ArithmeticError("the A17E7 glue basis audit failed")
    prefix = run_a17_e7_prefix(15)
    geometry = first_packet_geometry()
    witness = second_packet_hyperideal_witness()
    nonempty = tuple(packet for packet in prefix.packets if packet.roots)
    if [len(packet.roots) for packet in nonempty] != [2, 12]:
        raise ArithmeticError("the frozen A17E7 packets changed")
    if not witness.hyperideal or not witness.all_positive_walls_strict:
        raise ArithmeticError("the strict A17E7 hyperideal witness failed")
    if any(a17_e7_pair(root, root) != root.norm for root in prefix.accepted_roots):
        raise ArithmeticError("an accepted root has the wrong norm")
    if any(
        a17_e7_pair(left, right).sign(0) > 0
        for index, left in enumerate(prefix.accepted_roots)
        for right in prefix.accepted_roots[index + 1 :]
    ):
        raise ArithmeticError("the accepted A17E7 walls are not mutually acute")

    second = next(packet for packet in prefix.packets if packet.height == SECOND_PACKET_HEIGHT)
    q2_level = next(level for level in second.controls if level.norm == 2)
    e7 = enumerate_dominant_e7_weights(q2_level.target, q2_level.norm, q2_level.control)
    e7_vectors = tuple(
        weight
        for weight in e7.vectors
        if all((q2_level.target - weight.square).sign(place) >= 0 for place in (0, 1))
    )
    targets = tuple(sorted({q2_level.target - weight.square for weight in e7_vectors}))
    a17 = enumerate_dominant_a17_weights(targets, q2_level.norm, q2_level.control)
    roots = [root_payload(root) for root in prefix.accepted_roots]
    root_digest = sha256(
        json.dumps(roots, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    residue_histogram = Counter(
        f"A{root.a17.residue};E{root.e7.residue}" for root in second.roots
    )
    pairing_histogram = Counter(str(value) for value in witness.positive_wall_pairings)
    payload = {
        "scope": {
            "field": "Q(sqrt(2))",
            "hyperbolic_dimension": 24,
            "spatial_lattice": "Niemeier N(A17 E7), cyclic index-six glue <(3,1)>",
            "timelike_coefficient": silver(Silver(1, 1)),
            "height_bound": 15,
            "claim": (
                "Complete primitive chronological enumeration through the first two "
                "nonempty packets. The resulting prefix has a strict hyperideal ray; "
                "this is not a global nonreflectivity theorem."
            ),
        },
        "glue": {
            "root_lattice_discriminants": {"A17": 18, "E7": 2},
            "generator": {"A17": "omega_3", "E7": "omega_6"},
            "generator_square": "4",
            "residue_generator": [3, 1],
            "glue_order": 6,
            "explicit_basis_even_unimodular": True,
            "height_zero_dominant_orbits": len(height_zero_dominant_roots()),
            "height_zero_root_count": 432,
        },
        "mirror_norms": {
            "canonical": [silver(1), silver(2)],
            "completeness_proof": COMPLETE_UNIT_NORM_PROOF,
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
        "second_packet_residue_histogram": dict(sorted(residue_histogram.items())),
        "rank17_dynamic_program": {
            "e7_weight_count": len(e7.vectors),
            "a17_exact_energy_weight_count": len(a17.vectors),
            "a17_leaf_count": a17.leaf_count,
            "a17_visited_count": sum(a17.visited_by_depth),
            "a17_cap_pruned_count": sum(a17.cap_pruned_by_depth),
            "a17_energy_pruned_count": sum(a17.energy_pruned_by_depth),
        },
        "accepted_roots": roots,
        "accepted_roots_sha256": root_digest,
        "first_packet_geometry": {
            "extreme_ray_count": geometry.extreme_ray_count,
            "finite_count": geometry.finite_count,
            "cusp_count": geometry.cusp_count,
            "hyperideal_count": geometry.hyperideal_count,
        },
        "second_packet_witness": {
            "a17_weights": vector(witness.a17_weights),
            "e7_weights": vector(witness.e7_weights),
            "norm": silver(witness.norm),
            "height_zero_nonpositive": witness.height_zero_nonpositive,
            "all_positive_walls_strict": witness.all_positive_walls_strict,
            "hyperideal": witness.hyperideal,
            "positive_wall_pairing_histogram": dict(sorted(pairing_histogram.items())),
        },
        "outcome": (
            "No compact or finite-volume example in the certified prefix: the "
            "index-six glue contributes six non-root-lattice walls at the second "
            "packet, four genuinely coupling A17 and E7, but a strict hyperideal "
            "middle-A17-axis ray survives all 14 positive-height walls."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "sha256": sha256(OUTPUT.read_bytes()).hexdigest(),
                "accepted": len(prefix.accepted_roots),
                "cross_component_glue": sum(root.cross_component_glue for root in second.roots),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

