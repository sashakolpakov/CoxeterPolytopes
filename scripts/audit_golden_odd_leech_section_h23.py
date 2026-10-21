#!/usr/bin/env python3
"""Freeze the two exact Q(sqrt(5)) first-prefix screens in H23."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.search.golden_odd_leech_section import (
    COMPLETE_MIRROR_NORM_PROOF,
    HEIGHT_ZERO_PROOF,
    NORM_FIVE_SHELL_PROOF,
    NORM_FOUR_SHELL_PROOF,
    RADICAL_NORM_EXCLUSION,
    GoldenSectionControlLevel,
    section_norm_five_shell_audit,
    section_norm_four_shell_audit,
    twice_unit_section_first_packet_audit,
    unit_section_first_packet_audit,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/golden_odd_leech_section_h23_first_packets.json"


def golden(value: Golden | int) -> dict[str, str]:
    return Golden.coerce(value).to_json()


def triple_payload(triple) -> dict[str, int]:
    return {
        "u_square": triple.u_square,
        "v_square": triple.v_square,
        "u_dot_v": triple.u_dot_v,
    }


def control_payload(level: GoldenSectionControlLevel) -> dict:
    return {
        "height": golden(level.height),
        "norm": golden(level.norm),
        "control": golden(level.control),
        "target": golden(level.target),
        "pairing_modulus": level.pairing_modulus,
        "coefficient_minimum": level.coefficient_minimum,
        "necessary_gram_triples": [
            triple_payload(triple) for triple in level.necessary_gram_triples
        ],
    }


def shell_payload(shell, proof: str) -> dict:
    return {
        "norm": shell.norm,
        "even_coset_count": shell.even_coset_count,
        "odd_coset_count": shell.odd_coset_count,
        "component_counts": [
            {"component": label, "count": count}
            for label, count in shell.component_counts
        ],
        "total_count": shell.total_count,
        "proof": proof,
    }


def prefix_payload(audit) -> dict:
    return {
        "coefficient": golden(audit.coefficient),
        "complete_controls": [control_payload(level) for level in audit.controls],
        "first_packet": {
            "height": golden(audit.first_height),
            "norm": golden(audit.first_norm),
            "control": golden(audit.first_control),
            "target": golden(audit.first_target),
            "spatial_scale": golden(audit.first_scale),
            "spatial_norm": audit.first_spatial_norm,
            "gram_type": triple_payload(audit.first_gram_type),
            "root_count": audit.first_packet_root_count,
            "pairing_upper_bound": golden(audit.packet_pairing_upper_bound),
        },
        "witness": {
            "vector": f"(1, {audit.witness_spatial_scale}*e)",
            "spatial_vector": audit.witness_label,
            "square": golden(audit.witness_square),
            "closest_wall_pairing": golden(audit.closest_wall_pairing),
            "strict": audit.strict_hyperideal,
            "hyperideal": audit.witness_square.sign(0) > 0,
        },
        "finite_volume": audit.finite_volume,
        "compact": audit.compact,
    }


def build_payload() -> dict:
    unit = unit_section_first_packet_audit()
    twice = twice_unit_section_first_packet_audit()
    shell4 = section_norm_four_shell_audit()
    shell5 = section_norm_five_shell_audit()
    model = unit.model
    if (
        model.determinant != 3
        or model.smith_invariants != (1,) * 22 + (3,)
        or model.minimum != 3
        or model.minimal_shell_count != 2576
        or shell4.total_count != 53_682
        or shell5.total_count != 553_104
        or unit.first_packet_root_count != 553_104
        or twice.first_packet_root_count != 53_682
        or not unit.strict_hyperideal
        or not twice.strict_hyperideal
    ):
        raise ArithmeticError("the frozen golden S23 certificate changed")

    return {
        "scope": {
            "field": "Q(sqrt(5))",
            "dimension": 23,
            "spatial_lattice": "S23=(1^24)^perp in the catalogue odd Leech lattice",
            "coefficients": [golden(unit.coefficient), golden(twice.coefficient)],
            "claim": (
                "Complete local norm alphabet, height-zero census, and exact "
                "chronological prefix through the first nonempty packet for "
                "a=phi and a=2phi; both prefixes have strict hyperideal points."
            ),
            "bounded_claim": (
                "These two first prefixes are neither finite-volume nor compact. "
                "Later mirrors and other timelike coefficient ideals were not "
                "enumerated and can change the final chamber."
            ),
        },
        "section_model": {
            "rank": 23,
            "determinant": model.determinant,
            "smith_invariants": list(model.smith_invariants),
            "minimum": model.minimum,
            "minimal_shell_count": model.minimal_shell_count,
            "minimal_shell_sha256": model.minimal_shell_digest,
            "minimal_shell_span_rank": model.minimal_shell_span_rank,
            "gram_rank_mod_three": model.radical_mod_three_rank,
            "pairing_radical_dimension_mod_three": 23 - model.radical_mod_three_rank,
            "radical_minimum_lower_bound": model.radical_minimum_lower_bound,
        },
        "norm_alphabet": {
            "canonical": [golden(norm) for norm in unit.canonical_mirror_norms],
            "proof": COMPLETE_MIRROR_NORM_PROOF,
            "three_gate_radical_exclusion": RADICAL_NORM_EXCLUSION,
        },
        "height_zero": {
            "wall_count": 0,
            "ambient_necessary_type_counts": [
                {"norm": golden(norm), "count": len(types)}
                for norm, types in unit.height_zero_ambient_types
            ],
            "gated_type_counts": [
                {"norm": golden(norm), "count": len(types)}
                for norm, types in unit.height_zero_gated_types
            ],
            "proof": HEIGHT_ZERO_PROOF,
        },
        "section_shells": {
            "norm_4": shell_payload(shell4, NORM_FOUR_SHELL_PROOF),
            "norm_5": shell_payload(shell5, NORM_FIVE_SHELL_PROOF),
        },
        "prefixes": {
            "a=phi": prefix_payload(unit),
            "a=2phi": prefix_payload(twice),
        },
        "outcome": (
            "No compact or finite-volume first-prefix survivor for either "
            "golden coefficient: both exact finite prefix cones contain a "
            "strict spacelike interior point."
        ),
    }


def encoded_payload() -> bytes:
    return (json.dumps(build_payload(), indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    encoded = encoded_payload()
    if args.verify:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != encoded:
            raise SystemExit("FAIL: frozen golden S23 JSON differs from recomputation")
        status = "PASS"
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(encoded)
        status = "WROTE"
    print(
        json.dumps(
            {
                "status": status,
                "output": str(OUTPUT),
                "sha256": sha256(encoded).hexdigest(),
                "a=phi_walls": build_payload()["prefixes"]["a=phi"]["first_packet"]["root_count"],
                "a=2phi_walls": build_payload()["prefixes"]["a=2phi"]["first_packet"]["root_count"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
