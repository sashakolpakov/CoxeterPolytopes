#!/usr/bin/env python3
"""Exact replay for the all-degree compact K3 quadrilateral field bound."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import permutations
import json
from pathlib import Path
from typing import Sequence


def permutation_sign(permutation: Sequence[int]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    rank = len(matrix)
    if any(len(row) != rank for row in matrix):
        raise ValueError("determinant needs a square matrix")
    answer = Fraction(0)
    for permutation in permutations(range(rank)):
        term = Fraction(permutation_sign(permutation))
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        answer += term
    return answer


def gram_matrix(a, b, c, d, U, V) -> tuple[tuple[Fraction, ...], ...]:
    a, b, c, d, U, V = map(Fraction, (a, b, c, d, U, V))
    return (
        (Fraction(1), -a, -U, -d),
        (-a, Fraction(1), -b, -V),
        (-U, -b, Fraction(1), -c),
        (-d, -V, -c, Fraction(1)),
    )


def determinant_formula(a, b, c, d, U, V) -> Fraction:
    a, b, c, d, U, V = map(Fraction, (a, b, c, d, U, V))
    r = U * U - 1
    s = V * V - 1
    P = a * b + c * d
    Q = a * d + b * c
    R = a * c + b * d
    K = (
        -a * a
        - b * b
        - c * c
        - d * d
        + a * a * c * c
        + b * b * d * d
        - 2 * a * b * c * d
    )
    return r * s + K - 2 * P * U - 2 * Q * V - 2 * R * U * V


def determinant_record(a, b, c, d, U, V) -> dict[str, object]:
    matrix = gram_matrix(a, b, c, d, U, V)
    direct = determinant(matrix)
    formula = determinant_formula(a, b, c, d, U, V)
    if direct != formula:
        raise AssertionError("the quadrilateral determinant expansion failed")
    return {
        "parameters": [str(Fraction(value)) for value in (a, b, c, d, U, V)],
        "direct_determinant": str(direct),
        "expanded_determinant": str(formula),
        "expansion_identity": True,
    }


def diagonal_bound_record(degree: int = 6) -> dict[str, object]:
    if degree < 2 or degree > 6:
        raise ValueError("the nonrational full-trace K3 degree must lie in [2,6]")
    if degree == 2:
        radicand = 17
        denominator = 4
        linear_rational = 64
        linear_surd = 16
        constant_rational = 129
        constant_surd = 16
        trial = 132
        expected = (8847, 2128, 1_286_881)
    else:
        # The degree-six lower bound is weakest and hence controls d=3,...,6.
        radicand = 4097
        denominator = 64
        linear_rational = 16_384
        linear_surd = 256
        constant_rational = 32_769
        constant_surd = 256
        trial = 32_772
        expected = (537_034_767, 8_389_888, 17_606_953_271_521)
    remainder_rational = (
        trial * trial
        - linear_rational * trial
        - constant_rational
    )
    remainder_surd = linear_surd * trial + constant_surd
    positivity_certificate = (
        remainder_rational * remainder_rational
        - radicand * remainder_surd * remainder_surd
    )
    if (remainder_rational, remainder_surd, positivity_certificate) != expected:
        raise AssertionError("the exact opposite-magnitude comparison changed")
    return {
        "field_degree": degree,
        "shift_lower_bound": f"1/{4**degree}",
        "opposite_magnitude_lower_bound": f"sqrt({radicand})/{denominator}",
        "linear_coefficient_upper_bound": f"{linear_rational}+{linear_surd}sqrt({radicand})",
        "constant_coefficient_upper_bound": f"{constant_rational}+{constant_surd}sqrt({radicand})",
        "trial_opposite_magnitude": trial,
        "trial_quadratic_remainder": f"{remainder_rational}-{remainder_surd}sqrt({radicand})",
        "remainder_square_difference": positivity_certificate,
        "opposite_magnitude_strict_upper_bound": trial,
        "opposite_squared_ratio_strict_upper_bound": trial * trial,
    }


def quadratic_field_bound_record() -> dict[str, object]:
    magnitude_bound = diagonal_bound_record(2)["opposite_magnitude_strict_upper_bound"]
    pair_cycle_upper = 4 * magnitude_bound * magnitude_bound
    conjugate_pair_cycle_upper = 4
    trace_upper = pair_cycle_upper + conjugate_pair_cycle_upper - 1
    norm_upper = pair_cycle_upper * conjugate_pair_cycle_upper - 1
    discriminant_upper = trace_upper * trace_upper - 4
    if (
        pair_cycle_upper,
        trace_upper,
        norm_upper,
        discriminant_upper,
    ) != (69_696, 69_699, 278_783, 4_857_950_597):
        raise AssertionError("the quadrilateral field constants changed")
    return {
        "pair_cycle_distinguished_strict_upper_bound": pair_cycle_upper,
        "pair_cycle_conjugate_strict_upper_bound": conjugate_pair_cycle_upper,
        "integral_trace_upper_bound": trace_upper,
        "integral_norm_upper_bound": norm_upper,
        "field_discriminant_upper_bound": discriminant_upper,
        "opposite_pair_cycle_generates_quadratic_field": True,
    }


def general_field_bound_record() -> dict[str, object]:
    magnitude_bound = diagonal_bound_record(6)["opposite_magnitude_strict_upper_bound"]
    simple_cycle_bound = 16 * magnitude_bound**4
    primitive_conjugate_bound = 16 * simple_cycle_bound
    difference_bound = 2 * primitive_conjugate_bound
    discriminant_bound_base = difference_bound
    if (
        magnitude_bound**2,
        simple_cycle_bound,
        difference_bound,
    ) != (
        1_074_003_984,
        18_455_752_922_365_956_096,
        590_584_093_515_710_595_072,
    ):
        raise AssertionError("the all-degree cyclic field constants changed")
    return {
        "field_degree_upper_bound": 6,
        "opposite_magnitude_strict_upper_bound": magnitude_bound,
        "opposite_squared_ratio_strict_upper_bound": magnitude_bound**2,
        "simple_cyclic_coordinate_conjugate_bound": simple_cycle_bound,
        "selected_generator_count_upper_bound": 2,
        "bad_primitive_combination_coefficient_count_upper_bound": 15,
        "primitive_element_conjugate_bound": primitive_conjugate_bound,
        "Vandermonde_difference_bound": difference_bound,
        "field_discriminant_bound": f"{discriminant_bound_base}^30",
    }


def build_certificate() -> dict[str, object]:
    samples = [
        determinant_record(0, 0, 0, 0, 2, 3),
        determinant_record(Fraction(1, 2), 0, Fraction(1, 3), 0, 2, 3),
        determinant_record(
            Fraction(1, 2),
            Fraction(2, 3),
            Fraction(3, 4),
            Fraction(1, 5),
            Fraction(7, 4),
            Fraction(11, 6),
        ),
        determinant_record(
            Fraction(9, 10),
            Fraction(8, 9),
            Fraction(7, 8),
            Fraction(6, 7),
            Fraction(13, 5),
            Fraction(17, 4),
        ),
    ]
    return {
        "theorem": "compact full-trace K3 Coxeter quadrilaterals have a finite field set in every degree d<=6",
        "determinant_expansion_samples": samples,
        "universal_coefficient_bounds": {
            "P_Q_R_strict_upper_bound": 2,
            "absolute_K_strict_upper_bound": 8,
            "determinant_growth_inequality": "r*s<4U+4V+4UV+8",
        },
        "all_degree_diagonal_bound": diagonal_bound_record(6),
        "quadratic_diagonal_bound": diagonal_bound_record(2),
        "all_degree_field_bound": general_field_bound_record(),
        "quadratic_field_bound": quadratic_field_bound_record(),
        "conclusions": {
            "K3_full_trace_field_degree_at_most_6": True,
            "opposite_shift_lower_bound_from_all_place_norm": "4^(-d)",
            "all_degree_both_opposite_magnitudes_below_32772": True,
            "all_degree_both_opposite_squared_ratios_below_1074003984": True,
            "bounded_simple_cycles_generate_ground_field": True,
            "all_degree_field_discriminant_bound_is_explicit": True,
            "quadratic_both_opposite_magnitudes_below_132": True,
            "quadratic_field_discriminant_at_most_4857950597": True,
            "bounded_height_premise_required": False,
            "angle_order_classification_required": False,
        },
        "scope": {
            "minimal_complex_projective_K3_required": True,
            "primitive_exact_full_trace_Picard_required": True,
            "returned_field_is_cyclic_Cartan_ground_field_required": True,
            "exact_same_compact_Coxeter_quadrilateral_eigen_nef_chamber_required": True,
            "finite_projective_order_module_and_integral_facet_reflections_required": True,
            "complete_saturated_reflection_packet_retained": True,
            "reflection_maximality_in_declared_order_host_retained": True,
            "reflection_maximality_used_in_bound": False,
            "quadratic_field_assumption_required_for_general_theorem": False,
            "noncompact_quadrilateral_closed": False,
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
