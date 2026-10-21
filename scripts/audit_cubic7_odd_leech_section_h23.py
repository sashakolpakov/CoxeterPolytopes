#!/usr/bin/env python3
"""Archive the exact first packet for the cubic-seven odd-Leech section."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_odd_leech_section import (
    COMPLETE_MIRROR_NORM_PROOF,
    RADICAL_NORM_EXCLUSION,
    odd_leech_section_first_packet_audit,
)


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _gram(value) -> dict[str, object]:
    return {
        "tuple": [
            value.u_square,
            value.v_square,
            value.w_square,
            value.u_dot_v,
            value.u_dot_w,
            value.v_dot_w,
        ],
        "determinant": value.determinant,
        "rank": value.rank,
    }


def build_audit() -> dict[str, object]:
    audit = odd_leech_section_first_packet_audit()
    model = audit.model
    return {
        "schema_version": "0.1.0",
        "status": "FIRST_PACKET_INFINITE_VOLUME; EVENTUAL_CHAMBER_UNRESOLVED",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "hyperbolic_dimension": 23,
        "host": "<-theta> + ((1^24)^perp in O24)",
        "spatial_lattice": {
            "construction": (
                "primitive orthogonal complement of d=(1^24)/sqrt(8) "
                "inside the Nebe-Sloane odd-Leech catalogue lattice"
            ),
            "catalogue_source": (
                "https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/oddLeech.html"
            ),
            "rank": 23,
            "determinant": model.determinant,
            "smith_invariants": list(model.smith_invariants),
            "leading_principal_minors": list(model.leading_principal_minors),
            "minimum": model.minimum,
            "minimal_shell_count": model.minimal_shell_count,
            "minimal_shell_description": (
                "the 2576 weight-12 extended-Golay sign vectors"
            ),
            "minimal_shell_digest": model.minimal_shell_digest,
            "minimal_shell_span_rank": model.minimal_shell_span_rank,
            "pairing_radical_mod_3": {
                "gram_rank_mod_3": model.radical_mod_three_rank,
                "generator_in_section_basis": list(
                    model.radical_generator_coordinates
                ),
                "generator_raw_residue_mod_3": list(
                    model.radical_generator_raw_residue
                ),
                "minimum_strict_lower_bound": model.radical_minimum_lower_bound,
                "norm_3_and_6_exclusion_proof": RADICAL_NORM_EXCLUSION,
            },
        },
        "mirror_norm_theorem": {
            "canonical_norms": [
                _value(norm) for norm in audit.canonical_mirror_norms
            ],
            "complete_local_proof": COMPLETE_MIRROR_NORM_PROOF,
            "norm_3_and_6_spatial_gate": (
                "each of u,v,w lies in the mod-three pairing radical"
            ),
        },
        "height_zero": {
            "ambient_gram_types_by_norm": [
                {
                    "norm": _value(norm),
                    "types": [_gram(gram) for gram in types],
                }
                for norm, types in audit.height_zero_ambient_types
            ],
            "gated_gram_types_by_norm": [
                {
                    "norm": _value(norm),
                    "types": [_gram(gram) for gram in types],
                }
                for norm, types in audit.height_zero_gated_types
            ],
            "reflecting_wall_count": audit.height_zero_wall_count,
        },
        "complete_control_census": {
            "safe_integer_height_cap": 7,
            "exact_endpoint": _value(audit.first_height),
            "control_count": len(audit.controls),
            "levels": [
                {
                    "height": _value(level.height),
                    "norm": _value(level.norm),
                    "control": _value(level.control),
                    "target": _value(level.target),
                    "pairing_modulus": level.pairing_modulus,
                    "coefficient_minimum": level.coefficient_minimum,
                    "necessary_gram_types": [
                        _gram(gram) for gram in level.necessary_gram_types
                    ],
                }
                for level in audit.controls
            ],
        },
        "first_packet": {
            "height": _value(audit.first_height),
            "norm": _value(audit.first_norm),
            "control": _value(audit.first_control),
            "spatial_scale": _value(audit.first_scale),
            "target": _value(audit.first_target),
            "unique_gram_type": _gram(audit.first_gram_type),
            "root_formula": (
                "r_b=(2+theta,(1+theta)b), b in S23, b^2=3"
            ),
            "primitive": True,
            "reflection_integral": True,
            "root_count": audit.first_packet_root_count,
            "distinct_spatial_pairings": list(audit.distinct_spatial_pairings),
            "distinct_root_pairings": [
                _value(pairing) for pairing in audit.distinct_root_pairings
            ],
            "simultaneously_vinberg_acute": True,
            "accepted_root_count": audit.first_packet_root_count,
        },
        "strict_hyperideal_certificate": {
            "affine_point": (
                "z=(1,e), e=(4,-4,0^22)/sqrt(8) in S23"
            ),
            "spatial_square": audit.witness_square,
            "shell_pairing_histogram": [
                [pairing, count]
                for pairing, count in model.witness_shell_pairing_histogram
            ],
            "maximum_spatial_pairing": audit.witness_maximum_pairing,
            "minimum_wall_margin": _value(audit.witness_minimum_wall_margin),
            "lorentzian_excess": _value(audit.witness_hyperideal_excess),
            "strict": True,
            "finite_volume": audit.finite_volume,
            "compact": audit.compact,
        },
        "scope": (
            "Every primitive mirror norm, height-zero Gram type, control, and "
            "necessary coefficient Gram type is exhausted through the first "
            "nonempty packet. The strict spacelike point proves that exact "
            "prefix has infinite volume. It does not prove that the eventual "
            "Vinberg chamber has infinite volume, because later mirrors may "
            "cut the patch. No compact or finite-volume example is claimed."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "data/enumerations/cubic7_odd_leech_section_h23_first_packet.json"
        ),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": payload["status"],
                "spatial_lattice": payload["spatial_lattice"],
                "first_packet": payload["first_packet"],
                "strict_hyperideal_certificate": payload[
                    "strict_hyperideal_certificate"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
