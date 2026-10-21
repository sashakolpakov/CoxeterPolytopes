#!/usr/bin/env python3
"""Freeze the two exact Q(sqrt(2)) first-prefix screens in H23."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_odd_leech_section import (
    COMPLETE_MIRROR_NORM_PROOF,
    HEIGHT_ZERO_PROOF,
    NORM_SIX_SHELL_PROOF,
    RADICAL_NORM_EXCLUSION,
    SilverSectionControlLevel,
    norm_six_section_shell_audit,
    sqrt2_section_first_packet_audit,
    unit_section_first_packet_audit,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/silver_odd_leech_section_h23_first_packets.json"


def silver(value: Silver | int) -> dict[str, str]:
    return Silver.coerce(value).to_json()


def triple_payload(triple) -> dict[str, int]:
    return {
        "u_square": triple.u_square,
        "v_square": triple.v_square,
        "u_dot_v": triple.u_dot_v,
    }


def control_payload(level: SilverSectionControlLevel) -> dict:
    return {
        "height": silver(level.height),
        "norm": silver(level.norm),
        "control": silver(level.control),
        "target": silver(level.target),
        "coefficient_minimum": level.coefficient_minimum,
        "necessary_gram_triples": [
            triple_payload(triple) for triple in level.necessary_gram_triples
        ],
    }


def prefix_payload(audit) -> dict:
    return {
        "coefficient": silver(audit.coefficient),
        "complete_controls": [control_payload(level) for level in audit.controls],
        "first_height": silver(audit.first_height),
        "first_packet_norm_counts": [
            {"norm": silver(norm), "primitive_walls": count}
            for norm, count in audit.first_packet_norm_counts
        ],
        "imprimitive_candidates_rejected": audit.imprimitive_candidate_count,
        "accepted_walls": audit.accepted_wall_count,
        "packet_pairing_upper_bounds": [
            {"pair": label, "upper_bound": silver(value)}
            for label, value in audit.packet_pairing_upper_bounds
        ],
        "witness": {
            "vector": "(1, t*(4,-4,0^22)/sqrt(8))",
            "t": str(audit.witness_spatial_scale),
            "square": silver(audit.witness_square),
            "closest_wall_pairing": silver(audit.closest_wall_pairing),
            "strict": audit.strict_hyperideal,
            "hyperideal": audit.witness_square.sign(0) > 0,
        },
        "finite_volume": audit.finite_volume,
        "compact": audit.compact,
    }


def main() -> None:
    unit = unit_section_first_packet_audit()
    dyadic = sqrt2_section_first_packet_audit()
    shell6 = norm_six_section_shell_audit()
    model = unit.model
    if (
        model.determinant != 3
        or model.smith_invariants != (1,) * 22 + (3,)
        or model.minimum != 3
        or model.minimal_shell_count != 2576
        or shell6.total_count != 3_740_352
        or unit.accepted_wall_count != 2_576
        or dyadic.accepted_wall_count != 3_742_928
        or not unit.strict_hyperideal
        or not dyadic.strict_hyperideal
    ):
        raise ArithmeticError("the frozen silver S23 certificate changed")

    payload = {
        "scope": {
            "field": "Q(sqrt(2))",
            "dimension": 23,
            "spatial_lattice": "S23=(1^24)^perp in the catalogue odd Leech lattice",
            "coefficients": [silver(unit.coefficient), silver(dyadic.coefficient)],
            "claim": (
                "Complete local norm alphabet, height-zero census, and exact "
                "chronological prefix through the first nonempty packet for "
                "both coefficients; strict hyperideal witnesses."
            ),
            "bounded_claim": (
                "The two first prefixes are neither finite-volume nor compact. "
                "Later mirrors were not enumerated and can cut the witnesses."
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
            "radical_mod_three_rank": model.radical_mod_three_rank,
            "radical_minimum_lower_bound": model.radical_minimum_lower_bound,
            "norm_four_witness_pairing_histogram": [
                [pairing, count]
                for pairing, count in model.witness_shell_pairing_histogram
            ],
        },
        "norm_alphabet": {
            "canonical": [silver(norm) for norm in unit.canonical_mirror_norms],
            "proof": COMPLETE_MIRROR_NORM_PROOF,
            "three_gate_radical_exclusion": RADICAL_NORM_EXCLUSION,
        },
        "height_zero": {
            "wall_count": 0,
            "ambient_necessary_type_counts": [
                {"norm": silver(norm), "count": len(types)}
                for norm, types in unit.height_zero_ambient_types
            ],
            "gated_type_counts": [
                {"norm": silver(norm), "count": len(types)}
                for norm, types in unit.height_zero_gated_types
            ],
            "proof": HEIGHT_ZERO_PROOF,
        },
        "norm_six_shell": {
            "proof": NORM_SIX_SHELL_PROOF,
            "odd_coset_count": shell6.odd_coset_count,
            "even_weight_zero_count": shell6.even_weight_zero_count,
            "even_weight_eight_count": shell6.even_weight_eight_count,
            "even_weight_twelve_count": shell6.even_weight_twelve_count,
            "total_count": shell6.total_count,
            "norm_four_witness_pairing_histogram": [
                [pairing, count]
                for pairing, count in shell6.witness_pairing_histogram
            ],
        },
        "prefixes": {
            "a=1+sqrt(2)": prefix_payload(unit),
            "a=sqrt(2)": prefix_payload(dyadic),
        },
        "outcome": (
            "No compact or finite-volume first-prefix survivor: both exact "
            "prefix cones contain a strict spacelike interior point."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "sha256": sha256(OUTPUT.read_bytes()).hexdigest(),
                "unit_walls": unit.accepted_wall_count,
                "sqrt2_walls": dyadic.accepted_wall_count,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
