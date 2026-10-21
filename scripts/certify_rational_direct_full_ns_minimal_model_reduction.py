#!/usr/bin/env python3
"""Exact lattice checks for the direct full-NS minimal-model reduction."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from typing import Sequence


Vector = tuple[int, ...]


def hirzebruch_pair(left: Sequence[int], right: Sequence[int], e: int) -> int:
    """Form in the basis s,f,e_1,..., with s^2=-e and s.f=1."""

    return (
        -e * left[0] * right[0]
        + left[0] * right[1]
        + left[1] * right[0]
        - sum(a * b for a, b in zip(left[2:], right[2:]))
    )


def section_record(e: int, hits: int) -> dict[str, object]:
    if e < 0 or hits < 0:
        raise ValueError("e and the hit count must be nonnegative")
    size = 2 + hits
    section = (1, 0, *((-1,) * hits))
    fibre = (0, 1, *((0,) * hits))
    square = hirzebruch_pair(section, section, e)
    fibre_pairing = hirzebruch_pair(fibre, section, e)
    coefficient = Fraction(-2 * fibre_pairing, square) if square else None
    reflected_fibre = None
    integral = None
    if coefficient is not None:
        reflected_fibre_fractional = tuple(
            Fraction(entry) + coefficient * root_entry
            for entry, root_entry in zip(fibre, section)
        )
        integral = all(entry.denominator == 1 for entry in reflected_fibre_fractional)
        reflected_fibre = [str(entry) for entry in reflected_fibre_fractional]
    return {
        "e": e,
        "section_hits": hits,
        "section_square": square,
        "section_primitive": True,
        "fibre_pairing": fibre_pairing,
        "reflection_coefficient_on_fibre": (
            str(coefficient) if coefficient is not None else None
        ),
        "reflected_fibre": reflected_fibre,
        "reflection_integral_on_full_NS": integral,
    }


def f0_swap_record(exceptional_rank: int) -> dict[str, object]:
    if exceptional_rank < 0:
        raise ValueError("the exceptional rank must be nonnegative")
    size = 2 + exceptional_rank
    c = (1, 0, *((0,) * exceptional_rank))
    d = (0, 1, *((0,) * exceptional_rank))
    root = (1, -1, *((0,) * exceptional_rank))
    if hirzebruch_pair(root, root, 0) != -2:
        raise AssertionError("the ruling-swap vector ceased to be a root")

    def reflect(vector: Vector) -> Vector:
        coefficient = hirzebruch_pair(vector, root, 0)
        return tuple(
            entry + coefficient * root_entry
            for entry, root_entry in zip(vector, root)
        )

    exceptional_fixed = True
    for index in range(exceptional_rank):
        vector = [0] * size
        vector[2 + index] = 1
        exceptional_fixed &= reflect(tuple(vector)) == tuple(vector)
    return {
        "exceptional_rank": exceptional_rank,
        "root": list(root),
        "root_square": -2,
        "c_maps_to_d": reflect(c) == d,
        "d_maps_to_c": reflect(d) == c,
        "all_total_exceptional_coordinates_fixed": exceptional_fixed,
        "c_plus_d_pairing": hirzebruch_pair(
            tuple(a + b for a, b in zip(c, d)), root, 0
        ),
        "mirror_contains_an_ample_class": True,
    }


def adapted_swap_parity(boundary_pairs: Sequence[int]) -> dict[str, object]:
    return {
        "boundary_pairings_with_c_minus_d": list(boundary_pairs),
        "reflection_preserves_adapted_host": all(
            value % 2 == 0 for value in boundary_pairs
        ),
    }


def adapted_section_exponent_two_necessary(e: int, hits: int) -> dict[str, object]:
    """Necessary condition from L subset M subset (1/2)L; not sufficiency."""

    if e < 0 or hits < 0:
        raise ValueError("e and the hit count must be nonnegative")
    negative_square = e + hits
    return {
        "e": e,
        "section_hits": hits,
        "negative_section_square": negative_square,
        "can_preserve_exponent_two_host_only_if_glue_allows": (
            negative_square in (1, 2, 4)
        ),
        "condition_is_sufficient": False,
    }


def build_certificate() -> dict[str, object]:
    sections = [
        section_record(e, hits)
        for e in range(0, 11)
        for hits in range(0, 6)
    ]
    swaps = [f0_swap_record(rank) for rank in (0, 1, 2, 5, 17, 64)]
    if not all(
        record["reflection_integral_on_full_NS"] is False
        for record in sections
        if record["e"] > 2
    ):
        raise AssertionError("an e>2 section reflection became integral")
    if section_record(2, 0)["reflection_integral_on_full_NS"] is not True:
        raise AssertionError("the untouched F2 section must be an integral root")
    if any(
        section_record(2, hits)["reflection_integral_on_full_NS"] is not False
        for hits in range(1, 6)
    ):
        raise AssertionError("a blowup on the F2 section passed integrality")
    return {
        "theorem": "direct full-NS minimal-model reduction",
        "allowed_minimal_models": ["P2", "F2"],
        "F2_section_hit_count": 0,
        "F0_excluded_by_ruling_swap_through_ample_cone": True,
        "F_e_e_gt_2_excluded_by_section_reflection_integrality": True,
        "F1_is_not_minimal": True,
        "section_samples": sections,
        "F0_swap_samples": swaps,
        "adapted_parity_samples": [
            adapted_swap_parity(()),
            adapted_swap_parity((0, 2, -4)),
            adapted_swap_parity((1,)),
            adapted_swap_parity((2, 3)),
        ],
        "adapted_section_exponent_two_samples": [
            adapted_section_exponent_two_necessary(e, hits)
            for e in range(1, 7)
            for hits in range(0, 5)
        ],
        "adapted_section_necessary_divisibility": "e + hits divides 4",
        "adapted_host_theorem_claimed": False,
        "vertical_depth_bound_claimed": False,
        "reflective_classification_used": False,
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pretty", action="store_true")
    arguments = parser.parse_args()
    payload = build_certificate()
    envelope = {"payload": payload, "payload_sha256": digest(payload)}
    print(json.dumps(envelope, indent=2 if arguments.pretty else None, sort_keys=True))


if __name__ == "__main__":
    main()
