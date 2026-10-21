#!/usr/bin/env python3
"""Exact replay for the general rank-six K3 right-pentagon field bound."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
from typing import Sequence


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def shifted_pentagon(
    first: int | Fraction,
    third: int | Fraction,
) -> tuple[Fraction, ...]:
    x0 = Fraction(first)
    x2 = Fraction(third)
    x1 = x0 * x2 - 1
    if x0 <= 0 or x2 <= 0 or x1 <= 0:
        raise ValueError("the distinguished free coordinates need x0,x2,x0*x2-1>0")
    x3 = (x2 + 1) / x1
    x4 = (x0 + 1) / x1
    return x0, x1, x2, x3, x4


def recurrence_residuals(coordinates: Sequence[Fraction]) -> tuple[Fraction, ...]:
    if len(coordinates) != 5:
        raise ValueError("an A2 pentagon has five coordinates")
    return tuple(
        coordinates[(index - 1) % 5] * coordinates[(index + 1) % 5]
        - coordinates[index]
        - 1
        for index in range(5)
    )


def closure_record(
    first: int | Fraction,
    third: int | Fraction,
) -> dict[str, object]:
    coordinates = shifted_pentagon(first, third)
    residuals = recurrence_residuals(coordinates)
    if any(residuals):
        raise AssertionError("the five A2 recurrence identities failed")
    A, B, C, D, E = (coordinate + 1 for coordinate in coordinates)
    minors = (
        C * E - C - D - E + 1,
        -A + B * E - B - E + 1,
        B * D - B - C - D + 1,
        A * D - A - D - E + 1,
        A * C - A - B - C + 1,
    )
    if any(minors):
        raise AssertionError("a principal four-by-four Gram minor survived")
    return {
        "shifted_coordinates": [fraction_text(value) for value in coordinates],
        "squared_pair_coordinates": [
            fraction_text(value) for value in (A, B, C, D, E)
        ],
        "A2_recurrence_residuals": [fraction_text(value) for value in residuals],
        "principal_Gram_minor_residuals": [fraction_text(value) for value in minors],
        "moving_denominator": fraction_text(coordinates[1]),
    }


def quarter_integral_bound_record(scale: int = 4) -> dict[str, object]:
    if scale <= 0:
        raise ValueError("the integral scale must be positive")
    lower = Fraction(1, scale * scale)
    upper = Fraction(1, lower) + (1 + lower) / (lower * lower)
    # Equivalently 1/eps + (1+eps)/(eps*z), maximized as z descends to eps.
    if scale == 4 and upper != 288:
        raise AssertionError("the quarter-integral A2 constant changed")
    return {
        "integral_scale": scale,
        "conjugate_scaled_interval": f"(-{scale},0)",
        "strict_shifted_coordinate_lower_bound": fraction_text(lower),
        "strict_shifted_coordinate_upper_bound": fraction_text(upper),
        "moving_denominator_upper_bound": fraction_text(upper),
        "proof_identity": "x[i+3]=(z+1)/(x[i]*z-1)",
    }


def general_field_bound_record() -> dict[str, object]:
    shifted_upper = 288
    pair_upper = 4 * (shifted_upper + 1)
    conjugate_pair_upper = 4
    trace_upper = pair_upper + conjugate_pair_upper - 1
    norm_upper = pair_upper * conjugate_pair_upper - 1
    discriminant_upper = trace_upper * trace_upper - 4
    if (
        pair_upper,
        trace_upper,
        norm_upper,
        discriminant_upper,
    ) != (1156, 1159, 4623, 1_343_277):
        raise AssertionError("the general field constants changed")
    return {
        "shifted_coordinate_upper_bound": shifted_upper,
        "pair_cycle_distinguished_strict_upper_bound": pair_upper,
        "pair_cycle_conjugate_strict_upper_bound": conjugate_pair_upper,
        "integral_trace_upper_bound": trace_upper,
        "integral_norm_upper_bound": norm_upper,
        "polynomial_and_field_discriminant_upper_bound": discriminant_upper,
        "each_pair_cycle_generates_quadratic_field": True,
    }


def squarefree_part(value: int) -> int:
    if value <= 0:
        raise ValueError("value must be positive")
    result = 1
    remaining = value
    prime = 2
    while prime * prime <= remaining:
        parity = 0
        while remaining % prime == 0:
            remaining //= prime
            parity ^= 1
        if parity:
            result *= prime
        prime += 1
    if remaining > 1:
        result *= remaining
    return result


def integral_coordinate_candidate_packet() -> dict[str, object]:
    rows: list[dict[str, int]] = []
    radicands: set[int] = set()
    for trace in range(3, 35):
        for norm in range(1, 34):
            discriminant = trace * trace - 4 * norm
            if discriminant <= 0 or isqrt(discriminant) ** 2 == discriminant:
                continue
            # The smaller root is in (0,1): sqrt(D)>T-2.
            if discriminant <= (trace - 2) ** 2:
                continue
            radicand = squarefree_part(discriminant)
            rows.append(
                {
                    "trace": trace,
                    "norm": norm,
                    "polynomial_discriminant": discriminant,
                    "field_squarefree_radicand": radicand,
                }
            )
            radicands.add(radicand)
    if len(rows) != 528 or len(radicands) != 336:
        raise AssertionError("the integral-coordinate candidate count changed")
    if max(row["polynomial_discriminant"] for row in rows) != 1152:
        raise AssertionError("the sharp candidate discriminant ceiling changed")
    return {
        "trace_range": [3, 34],
        "norm_range": [1, 33],
        "candidate_trace_norm_row_count": len(rows),
        "candidate_squarefree_radicand_count": len(radicands),
        "candidate_squarefree_radicands": sorted(radicands),
        "candidate_rows": rows,
        "candidate_not_realization_count": True,
    }


# Arithmetic in Q(phi), phi^2=phi+1, represented by a+b*phi.
PhiElement = tuple[Fraction, Fraction]


def phi_add(left: PhiElement, right: PhiElement) -> PhiElement:
    return left[0] + right[0], left[1] + right[1]


def phi_subtract(left: PhiElement, right: PhiElement) -> PhiElement:
    return left[0] - right[0], left[1] - right[1]


def phi_multiply(left: PhiElement, right: PhiElement) -> PhiElement:
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def phi_scale(scalar: int | Fraction, value: PhiElement) -> PhiElement:
    scalar = Fraction(scalar)
    return scalar * value[0], scalar * value[1]


def phi_conjugate(value: PhiElement) -> PhiElement:
    # phi'=1-phi.
    return value[0] + value[1], -value[1]


def phi_norm(value: PhiElement) -> Fraction:
    product = phi_multiply(value, phi_conjugate(value))
    if product[1]:
        raise AssertionError("quadratic norm was not rational")
    return product[0]


def phi_divide(numerator: PhiElement, denominator: PhiElement) -> PhiElement:
    norm = phi_norm(denominator)
    if norm == 0:
        raise ZeroDivisionError
    return phi_scale(1 / norm, phi_multiply(numerator, phi_conjugate(denominator)))


def phi_text(value: PhiElement) -> list[str]:
    return [fraction_text(value[0]), fraction_text(value[1])]


def poison_packet_record() -> dict[str, object]:
    one = (Fraction(1), Fraction(0))
    coordinates: tuple[PhiElement, ...] = (
        (Fraction(-1, 2), Fraction(1, 2)),
        (Fraction(-1, 4), Fraction(1, 4)),
        (Fraction(1, 2), Fraction(2)),
        (Fraction(8), Fraction(14)),
        (Fraction(2), Fraction(4)),
    )
    residuals = tuple(
        phi_subtract(
            phi_subtract(
                phi_multiply(coordinates[(index - 1) % 5], coordinates[(index + 1) % 5]),
                coordinates[index],
            ),
            one,
        )
        for index in range(5)
    )
    if any(value != (0, 0) for value in residuals):
        raise AssertionError("the poison A2 packet failed its recurrence")
    scaled_integral = tuple(
        all(coefficient.denominator == 1 for coefficient in phi_scale(4, value))
        for value in coordinates
    )
    if not all(scaled_integral):
        raise AssertionError("the poison packet lost quarter-integrality")
    integral = tuple(
        all(coefficient.denominator == 1 for coefficient in value)
        for value in coordinates
    )
    x0, x1, x2, _, _ = coordinates
    numerator = phi_scale(
        -32,
        phi_multiply(
            phi_multiply(phi_multiply(x0, x2), phi_add(x0, one)),
            phi_add(x2, one),
        ),
    )
    cycle = phi_divide(numerator, x1)
    if cycle != (Fraction(-408), Fraction(-664)):
        raise AssertionError("the poison length-five cycle changed")
    return {
        "field": "Q(sqrt(5))",
        "basis_relation": "phi^2=phi+1",
        "shifted_coordinates_in_1_phi_basis": [phi_text(value) for value in coordinates],
        "recurrence_residuals": [phi_text(value) for value in residuals],
        "four_times_each_coordinate_integral": list(scaled_integral),
        "each_coordinate_integral": list(integral),
        "directed_acute_five_cycle_in_1_phi_basis": phi_text(cycle),
        "first_three_coordinates_not_integral": integral[:3] == (False, False, False),
        "poisons_division_by_four_only": True,
        "asserts_exact_K3_or_maximal_carrier": False,
    }


def build_certificate() -> dict[str, object]:
    samples = []
    for first, third in (
        (Fraction(1, 8), 9),
        (Fraction(1, 2), 3),
        (1, 2),
        (2, 2),
        (Fraction(7, 3), Fraction(11, 5)),
        (17, Fraction(2, 17)),
    ):
        samples.append(closure_record(first, third))
    return {
        "theorem": "quarter-integral A2 closure bounds every compact rank-six K3 right-pentagon field",
        "closure_samples": samples,
        "quarter_integral_bound": quarter_integral_bound_record(),
        "general_field_bound": general_field_bound_record(),
        "integral_coordinate_corollary": integral_coordinate_candidate_packet(),
        "division_by_four_poison_packet": poison_packet_record(),
        "conclusions": {
            "all_five_Gram_identities_form_cyclic_A2_recurrence": True,
            "automatic_integrality_is_four_times_shifted_coordinate": True,
            "every_distinguished_shifted_coordinate_between_1_over_16_and_288": True,
            "moving_denominator_AC_minus_A_minus_C_below_288": True,
            "every_pair_cycle_generates_quadratic_field": True,
            "general_field_discriminant_at_most_1343277": True,
            "power_coupling_required": False,
            "bounded_mixed_height_premise_required": False,
            "integral_A_C_corollary_field_discriminant_at_most_1152": True,
            "integral_A_C_corollary_candidate_squareclasses": 336,
        },
        "scope": {
            "minimal_complex_projective_K3_required": True,
            "primitive_exact_full_trace_Picard_rank_six_required": True,
            "exact_same_compact_right_pentagon_eigen_nef_chamber_required": True,
            "finite_projective_order_module_and_integral_facet_reflections_required": True,
            "complete_saturated_reflection_packet_retained": True,
            "reflection_maximality_in_declared_order_host_retained": True,
            "reflection_maximality_used_in_bound": False,
            "A_and_C_algebraic_integers_required_for_general_theorem": False,
            "A_and_C_algebraic_integers_required_for_1152_corollary": True,
            "poison_packet_is_exact_K3_carrier": False,
            "arbitrary_n2_chamber_closed": False,
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
