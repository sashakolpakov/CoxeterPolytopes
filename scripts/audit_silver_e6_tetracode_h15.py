#!/usr/bin/env python3
"""Freeze the E6^4/tetracode H15 prefix and Niemeier triage table."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.niemeier_host_triage import niemeier_triage
from coxeter_polytopes.search.silver_e6_tetracode import (
    COMPLETE_UNIT_NORM_PROOF,
    first_packet_geometry,
    height_zero_dominant_roots,
    run_e6_tetracode_prefix,
    second_packet_hyperideal_witness,
    silver_e6_pair,
    ternary_tetracode,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/silver_e6_tetracode_h15.json"


def silver(value):
    return Silver.coerce(value).to_json()


def vector_payload(vector):
    return [silver(value) for value in vector]


def vector_digest(vectors) -> str:
    payload = json.dumps(
        [vector_payload(vector) for vector in vectors],
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return sha256(payload).hexdigest()


def main() -> None:
    prefix = run_e6_tetracode_prefix(15)
    geometry = first_packet_geometry()
    witness = second_packet_hyperideal_witness()
    nonempty = tuple(packet for packet in prefix.packets if packet.roots)
    if [len(packet.roots) for packet in nonempty] != [4, 44]:
        raise ArithmeticError("the two frozen E6^4 packets changed")
    if not witness.hyperideal or not witness.all_positive_walls_strict:
        raise ArithmeticError("the E6^4 hyperideal witness failed")
    if any(
        silver_e6_pair(root.vector, root.vector) != root.norm
        for root in prefix.accepted_roots
    ):
        raise ArithmeticError("an accepted E6^4 root has the wrong norm")

    pairing_histogram = Counter(
        str(pairing) for pairing in witness.positive_wall_pairings
    )
    rows = niemeier_triage()
    payload = {
        "scope": {
            "field": "Q(sqrt(2))",
            "hyperbolic_dimension": 24,
            "spatial_lattice": "Niemeier N(E6^4), exact ternary-tetracode glue",
            "timelike_coefficient": silver(prefix.accepted_roots[0].height / 2),
            "height_bound": 15,
            "claim": (
                "Complete primitive Vinberg enumeration through H=15 for this host; "
                "the prefix has an explicit hyperideal ray. This is not a global "
                "nonreflectivity theorem."
            ),
        },
        "mirror_norms": {
            "canonical": [silver(1), silver(2)],
            "completeness_proof": COMPLETE_UNIT_NORM_PROOF,
        },
        "tetracode": {
            "words": [list(word) for word in sorted(ternary_tetracode())],
            "weight_histogram": dict(
                sorted(
                    Counter(
                        sum(value != 0 for value in word)
                        for word in ternary_tetracode()
                    ).items()
                )
            ),
            "height_zero_dominant_orbits": len(height_zero_dominant_roots()),
            "height_zero_roots": 288,
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
                "controls": len(packet.controls),
                "candidate_count": packet.candidate_count,
                "rejected_count": packet.rejected_count,
                "accepted_count": len(packet.roots),
                "glue_coupled_count": sum(root.glue_coupled for root in packet.roots),
                "root_lattice_count": sum(not root.glue_coupled for root in packet.roots),
            }
            for packet in prefix.packets
        ],
        "accepted_roots": {
            "count": len(prefix.accepted_roots),
            "vector_sha256": vector_digest(root.vector for root in prefix.accepted_roots),
            "all_norms_correct": True,
            "all_same_packet_pairings_nonpositive": True,
        },
        "first_packet_geometry": {
            "cap": silver(geometry.cap),
            "extreme_ray_count": geometry.extreme_ray_count,
            "finite_vertex_count": geometry.finite_vertex_count,
            "cusp_count": geometry.cusp_count,
            "hyperideal_ray_count": geometry.hyperideal_ray_count,
            "compact": geometry.compact,
        },
        "second_packet_witness": {
            "dynkin_weights": [vector_payload(block) for block in witness.dynkin_weights],
            "norm": silver(witness.norm),
            "height_zero_nonpositive": witness.height_zero_nonpositive,
            "all_positive_walls_strict": witness.all_positive_walls_strict,
            "hyperideal": witness.hyperideal,
            "positive_wall_pairing_histogram": dict(sorted(pairing_histogram.items())),
        },
        "unsearched_niemeier_triage": [
            {
                **asdict(row),
                "components": [list(component) for component in row.components],
            }
            for row in rows
        ],
        "outcome": (
            "No compact or finite-volume example: 24 genuine tetracode-coupled "
            "walls appear in the second packet, but the exact hyperideal witness "
            "survives all 48 positive-height walls through H=15."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(OUTPUT), "sha256": sha256(OUTPUT.read_bytes()).hexdigest(), "accepted": len(prefix.accepted_roots)}, indent=2))


if __name__ == "__main__":
    main()
