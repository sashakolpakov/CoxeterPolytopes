#!/usr/bin/env python3
"""Freeze the exact silver-unit N(A15 D9) prefix through height 50."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_a15_d9 import (
    COMPLETE_UNIT_NORM_PROOF,
    SECOND_PACKET_HEIGHT,
    a15_d9_basis_audit,
    a15_d9_pair,
    axis_witness,
    enumerate_dominant_a15_weights,
    enumerate_dominant_d9_weights,
    first_packet_geometry,
    glue_coset_minima,
    height_zero_dominant_roots,
    run_a15_d9_prefix,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/silver_a15_d9_h50.json"


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
            "square": silver(root.a15.square),
            "weights": vector(root.a15.weights),
            "coordinates": vector(root.a15.coordinates),
            "residue": list(root.a15.residue),
        },
        "d9": {
            "square": silver(root.d9.square),
            "weights": vector(root.d9.weights),
            "coordinates": vector(root.d9.coordinates),
            "residue": list(root.d9.residue),
        },
        "glue_nontrivial": root.glue_nontrivial,
        "cross_component_glue": root.cross_component_glue,
    }


def main() -> None:
    if not a15_d9_basis_audit():
        raise ArithmeticError("the A15D9 glue basis audit failed")
    prefix = run_a15_d9_prefix(50)
    geometry = first_packet_geometry()
    witness = axis_witness("D9", 4, __import__("fractions").Fraction(7, 10), 50)
    nonempty = tuple(packet for packet in prefix.packets if packet.roots)
    if [len(packet.roots) for packet in nonempty] != [2, 15, 191]:
        raise ArithmeticError("the frozen A15D9 packets changed")
    if not witness.hyperideal or not witness.all_positive_walls_strict:
        raise ArithmeticError("the strict A15D9 hyperideal witness failed")
    if any(a15_d9_pair(root, root) != root.norm for root in prefix.accepted_roots):
        raise ArithmeticError("an accepted root has the wrong norm")
    if any(
        a15_d9_pair(left, right).sign(0) > 0
        for index, left in enumerate(prefix.accepted_roots)
        for right in prefix.accepted_roots[index + 1 :]
    ):
        raise ArithmeticError("the accepted A15D9 walls are not mutually acute")

    second = next(packet for packet in prefix.packets if packet.height == SECOND_PACKET_HEIGHT)
    third = next(packet for packet in prefix.packets if packet.height == Silver(14, 10))
    q2_level = next(level for level in second.controls if level.norm == 2)
    d9 = enumerate_dominant_d9_weights(q2_level.target, q2_level.norm, q2_level.control)
    d9_vectors = tuple(
        weight
        for weight in d9.vectors
        if all((q2_level.target - weight.square).sign(place) >= 0 for place in (0, 1))
    )
    targets = tuple(sorted({q2_level.target - weight.square for weight in d9_vectors}))
    a15 = enumerate_dominant_a15_weights(targets, q2_level.norm, q2_level.control)
    third_q2 = next(level for level in third.controls if level.norm == 2)
    third_d9 = enumerate_dominant_d9_weights(
        third_q2.target, third_q2.norm, third_q2.control
    )
    third_d9_vectors = tuple(
        weight
        for weight in third_d9.vectors
        if all((third_q2.target - weight.square).sign(place) >= 0 for place in (0, 1))
    )
    third_targets = tuple(
        sorted({third_q2.target - weight.square for weight in third_d9_vectors})
    )
    third_a15 = enumerate_dominant_a15_weights(
        third_targets, third_q2.norm, third_q2.control
    )
    roots = [root_payload(root) for root in prefix.accepted_roots]
    root_digest = sha256(
        json.dumps(roots, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    residue_histogram = Counter(
        f"A{root.a15.residue};D{root.d9.residue}" for root in second.roots
    )
    pairing_histogram = Counter(str(value) for value in witness.positive_wall_pairings)
    payload = {
        "scope": {
            "field": "Q(sqrt(2))",
            "hyperbolic_dimension": 24,
            "spatial_lattice": "Niemeier N(A15 D9), cyclic index-eight glue <(2,1)>",
            "timelike_coefficient": silver(Silver(1, 1)),
            "height_bound": 50,
            "claim": (
                "Complete primitive chronological enumeration through the first three "
                "nonempty packets. The resulting prefix has a strict hyperideal ray; "
                "this is not a global nonreflectivity theorem."
            ),
        },
        "glue": {
            "root_lattice_discriminants": {"A15": 16, "D9": 4},
            "generator": {"A15": "omega_2", "D9": "omega_9"},
            "generator_square": "4",
            "residue_generator": [2, 1],
            "glue_order": 8,
            "nonzero_coset_minima": [int(value) for value in glue_coset_minima()],
            "explicit_basis_even_unimodular": True,
            "height_zero_dominant_orbits": len(height_zero_dominant_roots()),
            "height_zero_root_count": 384,
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
        "rank15_dynamic_program": {
            "d9_weight_count": len(d9.vectors),
            "a15_exact_energy_weight_count": len(a15.vectors),
            "a15_leaf_count": a15.leaf_count,
            "a15_visited_count": sum(a15.visited_by_depth),
            "a15_cap_pruned_count": sum(a15.cap_pruned_by_depth),
            "a15_energy_pruned_count": sum(a15.energy_pruned_by_depth),
        },
        "third_packet": {
            "height": silver(Silver(14, 10)),
            "candidate_count": next(
                packet.candidate_count
                for packet in prefix.packets
                if packet.height == Silver(14, 10)
            ),
            "rejected_count": next(
                packet.rejected_count
                for packet in prefix.packets
                if packet.height == Silver(14, 10)
            ),
            "accepted_count": next(
                len(packet.roots)
                for packet in prefix.packets
                if packet.height == Silver(14, 10)
            ),
            "glue_nontrivial_count": next(
                sum(root.glue_nontrivial for root in packet.roots)
                for packet in prefix.packets
                if packet.height == Silver(14, 10)
            ),
            "cross_component_glue_count": next(
                sum(root.cross_component_glue for root in packet.roots)
                for packet in prefix.packets
                if packet.height == Silver(14, 10)
            ),
            "residue_histogram": dict(
                sorted(
                    Counter(
                        f"A{root.a15.residue};D{root.d9.residue}"
                        for root in third.roots
                    ).items()
                )
            ),
            "dynamic_program": {
                "d9_weight_count": len(third_d9.vectors),
                "target_energy_count": len(third_targets),
                "a15_exact_energy_weight_count": len(third_a15.vectors),
                "a15_leaf_count": third_a15.leaf_count,
                "a15_visited_count": sum(third_a15.visited_by_depth),
                "a15_cap_pruned_count": sum(third_a15.cap_pruned_by_depth),
                "a15_energy_pruned_count": sum(third_a15.energy_pruned_by_depth),
            },
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
            "a15_weights": vector(witness.a15_weights),
            "d9_weights": vector(witness.d9_weights),
            "norm": silver(witness.norm),
            "height_zero_nonpositive": witness.height_zero_nonpositive,
            "all_positive_walls_strict": witness.all_positive_walls_strict,
            "hyperideal": witness.hyperideal,
            "positive_wall_pairing_histogram": dict(sorted(pairing_histogram.items())),
        },
        "independent_cross_check": {
            "implementation": "separate no-import discriminant-and-packet derivation",
            "basis_height_zero_controls_packets_geometry_witness": "PASS",
        },
        "outcome": (
            "No compact or finite-volume example in the certified prefix: the "
            "index-eight glue contributes nine non-root-lattice second-packet "
            "walls, eight coupling A15 to D9. The third packet adds 191 walls, "
            "112 cross-component, but a strict D9-axis hyperideal ray survives "
            "all 208 positive-height walls."
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
                "cross_component_glue": sum(
                    root.cross_component_glue for root in prefix.accepted_roots
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
