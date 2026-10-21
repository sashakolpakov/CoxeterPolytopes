#!/usr/bin/env python3
"""Exact replay for the all-degree compact K3 right-pentagon field bound."""

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


def shifted_pentagon(first, third) -> tuple[Fraction, ...]:
    x0, x2 = Fraction(first), Fraction(third)
    x1 = x0 * x2 - 1
    if x0 <= 0 or x2 <= 0 or x1 <= 0:
        raise ValueError("positive pentagon coordinates require x0,x2,x0*x2-1>0")
    return x0, x1, x2, (x2 + 1) / x1, (x0 + 1) / x1


def recurrence_residuals(coordinates: Sequence[Fraction]) -> tuple[Fraction, ...]:
    if len(coordinates) != 5:
        raise ValueError("an A2 pentagon has five coordinates")
    return tuple(
        coordinates[(i - 1) % 5] * coordinates[(i + 1) % 5]
        - coordinates[i]
        - 1
        for i in range(5)
    )


def closure_record(first, third) -> dict[str, object]:
    coordinates = shifted_pentagon(first, third)
    residuals = recurrence_residuals(coordinates)
    if any(residuals):
        raise AssertionError("the cyclic A2 recurrence failed")
    A, B, C, D, E = (coordinate + 1 for coordinate in coordinates)
    minors = (
        C * E - C - D - E + 1,
        -A + B * E - B - E + 1,
        B * D - B - C - D + 1,
        A * D - A - D - E + 1,
        A * C - A - B - C + 1,
    )
    if any(minors):
        raise AssertionError("a principal four-normal Gram identity failed")
    return {
        "shifted_coordinates": [fraction_text(value) for value in coordinates],
        "recurrence_residuals": [fraction_text(value) for value in residuals],
        "principal_Gram_minor_residuals": [fraction_text(value) for value in minors],
    }


def degree_bound_record(degree: int) -> dict[str, object]:
    if degree < 2 or degree > 6:
        raise ValueError("the nonrational full-trace K3 degree lies in [2,6]")
    reciprocal_floor = 4**degree
    shifted_upper = reciprocal_floor**2 + 2 * reciprocal_floor
    squared_pair_upper = shifted_upper + 1
    if squared_pair_upper != (reciprocal_floor + 1) ** 2:
        raise AssertionError("the A2 upper-bound square identity failed")
    five_cycle_bound = 32 * (reciprocal_floor + 1) ** 5
    return {
        "field_degree": degree,
        "strict_shifted_coordinate_lower_bound": f"1/{reciprocal_floor}",
        "strict_shifted_coordinate_upper_bound": shifted_upper,
        "strict_squared_pair_coordinate_upper_bound": squared_pair_upper,
        "strict_pair_cycle_upper_bound": 4 * squared_pair_upper,
        "strict_length_five_cycle_absolute_upper_bound": five_cycle_bound,
        "proof_identity": "x[i+3]=(z+1)/(x[i]*z-1)",
    }


def all_degree_field_record() -> dict[str, object]:
    record = degree_bound_record(6)
    H = record["strict_length_five_cycle_absolute_upper_bound"]
    primitive_bound = 16 * H
    difference_bound = 2 * primitive_bound
    expected = (
        16_785_408,
        16_785_409,
        36_938_546_139_294_728_224,
        1_182_033_476_457_431_303_168,
    )
    actual = (
        record["strict_shifted_coordinate_upper_bound"],
        record["strict_squared_pair_coordinate_upper_bound"],
        H,
        difference_bound,
    )
    if actual != expected:
        raise AssertionError("the all-degree pentagon constants changed")
    return {
        "field_degree_upper_bound": 6,
        "uniform_shifted_coordinate_strict_upper_bound": actual[0],
        "uniform_squared_pair_coordinate_strict_upper_bound": actual[1],
        "simple_cyclic_generator_conjugate_bound": H,
        "selected_generator_count_upper_bound": 2,
        "bad_primitive_combination_coefficient_count_upper_bound": 15,
        "primitive_element_conjugate_bound": primitive_bound,
        "Vandermonde_difference_bound": difference_bound,
        "field_discriminant_bound": f"{difference_bound}^30",
    }


def quadratic_record() -> dict[str, object]:
    record = degree_bound_record(2)
    pair_upper = record["strict_pair_cycle_upper_bound"]
    trace_upper = pair_upper + 4 - 1
    norm_upper = pair_upper * 4 - 1
    discriminant_upper = trace_upper**2 - 4
    if (
        record["strict_shifted_coordinate_upper_bound"],
        pair_upper,
        trace_upper,
        norm_upper,
        discriminant_upper,
    ) != (288, 1156, 1159, 4623, 1_343_277):
        raise AssertionError("the quadratic compatibility constants changed")
    return {
        "strict_shifted_coordinate_lower_bound": "1/16",
        "strict_shifted_coordinate_upper_bound": 288,
        "pair_cycle_generates_quadratic_field": True,
        "quadratic_field_discriminant_upper_bound": discriminant_upper,
    }


def build_certificate() -> dict[str, object]:
    samples = [
        closure_record(Fraction(1, 8), 9),
        closure_record(Fraction(1, 2), 3),
        closure_record(1, 2),
        closure_record(2, 2),
        closure_record(Fraction(7, 3), Fraction(11, 5)),
    ]
    return {
        "theorem": "all full-trace K3 compact right-pentagon cyclic ground fields form a finite set",
        "closure_samples": samples,
        "degree_bounds": [degree_bound_record(degree) for degree in range(2, 7)],
        "all_degree_field_bound": all_degree_field_record(),
        "quadratic_compatibility": quadratic_record(),
        "conclusions": {
            "K3_full_trace_field_degree_at_most_6": True,
            "automatic_shift_floor_is_4_to_minus_d": True,
            "all_degree_every_shift_below_16785408": True,
            "all_degree_length_five_cycle_below_36938546139294728224": True,
            "bounded_cyclic_packet_gives_finite_field_set": True,
            "suggested_69648_bound_is_not_claimed": True,
            "bounded_mixed_height_premise_required": False,
            "power_coupling_required": False,
        },
        "scope": {
            "minimal_complex_projective_K3_required": True,
            "primitive_exact_full_trace_Picard_required": True,
            "returned_field_is_cyclic_Cartan_ground_field_required": True,
            "exact_same_compact_right_pentagon_eigen_nef_chamber_required": True,
            "finite_projective_order_module_and_integral_facet_reflections_required": True,
            "complete_saturated_reflection_packet_retained": True,
            "reflection_maximality_in_declared_order_host_retained": True,
            "reflection_maximality_used_in_bound": False,
            "non_right_or_noncompact_pentagon_closed": False,
            "arbitrary_n2_polygon_closed": False,
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
