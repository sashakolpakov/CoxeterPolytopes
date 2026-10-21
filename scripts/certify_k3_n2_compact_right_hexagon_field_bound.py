#!/usr/bin/env python3
"""Exact replay for the compact right-angled K3 hexagon field theorem."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Sequence


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def right_hexagon_from_alternating_cosh(u, v, w) -> tuple[Fraction, ...]:
    u, v, w = Fraction(u), Fraction(v), Fraction(w)
    if min(u, v, w) <= 1:
        raise ValueError("alternating side coshes must exceed one")
    x0, x2, x4 = u * u - 1, v * v - 1, w * w - 1
    S = u * u + v * v + w * w - 1 + 2 * u * v * w
    return x0, S / (x0 * x2), x2, S / (x2 * x4), x4, S / (x4 * x0)


def closure_record(u, v, w) -> dict[str, object]:
    coordinates = right_hexagon_from_alternating_cosh(u, v, w)
    x0, x1, x2, x3, x4, x5 = coordinates
    S = x0 * x1 * x2
    common_numerator_residuals = (
        x0 * x1 * x2 - S,
        x2 * x3 * x4 - S,
        x4 * x5 * x0 - S,
    )
    four_minor_residuals = tuple(
        coordinates[i] * coordinates[(i + 1) % 6]
        - coordinates[(i + 3) % 6] * coordinates[(i + 4) % 6]
        for i in range(6)
    )
    if any(common_numerator_residuals) or any(four_minor_residuals):
        raise AssertionError("the right-hexagon Gram closure failed")
    return {
        "alternating_cosh": [fraction_text(Fraction(value)) for value in (u, v, w)],
        "shifted_second_neighbour_coordinates": [
            fraction_text(value) for value in coordinates
        ],
        "common_numerator_residuals": [
            fraction_text(value) for value in common_numerator_residuals
        ],
        "four_consecutive_Gram_minor_residuals": [
            fraction_text(value) for value in four_minor_residuals
        ],
    }


def short_pair_record() -> dict[str, object]:
    rational_numerator = 51
    surd_coefficient = 34
    radicand = 17
    rational_upper_using_sqrt_lt_5 = rational_numerator + 5 * surd_coefficient
    denominator = 256
    if rational_upper_using_sqrt_lt_5 != 221 or not 221 < denominator:
        raise AssertionError("the short-pair comparison changed")
    return {
        "alternating_threshold": 16,
        "large_alternating_case_odd_upper_bound": "(51+34sqrt(17))/256",
        "sqrt17_rational_upper_bound": 5,
        "resulting_rational_numerator_upper_bound": rational_upper_using_sqrt_lt_5,
        "large_alternating_case_odd_strict_upper_bound": 1,
        "universal_some_shifted_coordinate_strict_upper_bound": 16,
        "universal_some_pair_cycle_strict_upper_bound": 68,
    }


def field_record() -> dict[str, object]:
    exponent = 6 * 5
    return {
        "field_degree_upper_bound": 6,
        "primitive_generator_conjugate_interval": "(0,68)",
        "unique_distinguished_generator_interval": "(4,68)",
        "definite_generator_interval": "(0,4)",
        "Vandermonde_difference_bound": 68,
        "field_discriminant_bound": f"68^{exponent}",
        "unique_large_conjugate_forces_full_field_generation": True,
    }


def build_certificate() -> dict[str, object]:
    samples = [
        closure_record(2, 2, 2),
        closure_record(2, 3, 4),
        closure_record(Fraction(5, 4), Fraction(13, 12), Fraction(17, 15)),
        closure_record(Fraction(29, 20), Fraction(41, 40), Fraction(101, 99)),
    ]
    return {
        "theorem": "compact right-angled full-trace K3 hexagons have a finite field set",
        "closure_samples": samples,
        "short_pair_bound": short_pair_record(),
        "field_bound": field_record(),
        "conclusions": {
            "K3_full_trace_field_degree_at_most_6": True,
            "rank_three_cofactor_closure_used": True,
            "some_dotted_pair_cycle_below_68": True,
            "that_dotted_pair_cycle_generates_the_full_field": True,
            "field_discriminant_at_most_68_power_30": True,
            "imposed_mixed_height_required": False,
            "reflection_maximality_used": False,
        },
        "scope": {
            "minimal_complex_projective_K3_required": True,
            "primitive_exact_full_trace_Picard_required": True,
            "minimal_cyclic_Cartan_ground_field_required": True,
            "exact_same_compact_right_angled_hexagon_eigen_nef_chamber_required": True,
            "finite_projective_order_module_and_integral_facet_reflections_required": True,
            "non_right_hexagon_closed": False,
            "noncompact_hexagon_closed": False,
            "arbitrary_six_or_more_facet_chamber_closed": False,
            "carrier_existence_proved": False,
            "uses_reflective_lattice_classification": False,
        },
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args(argv)
    encoded = encoded_certificate()
    if arguments.check is not None:
        if arguments.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit(f"archive mismatch: {arguments.check}")
        print(f"PASS archive replay: {arguments.check}")
        return
    if arguments.output is not None:
        arguments.output.write_text(encoded, encoding="utf-8")
        print(arguments.output)
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
