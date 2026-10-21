#!/usr/bin/env python3
"""Exact replay for the dominant-correspondence field-height theorem.

The general proof is symbolic.  This executable records its finite algebra,
the primitive-positive-combination bound, and the rank-four escape family
which shows why the mixed correspondence degree cannot be omitted.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, isqrt
from pathlib import Path
from typing import Sequence


def as_fraction(value: int | Fraction) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def spectral_identity_record(
    positive_square: int | Fraction,
    negative_place_weights: Sequence[int | Fraction],
    eigenvalues: Sequence[int | Fraction],
) -> dict[str, object]:
    """Check b-lambda_0*s=sum_i(lambda_i-lambda_0)w_i exactly.

    The first eigenvalue is distinguished.  The non-distinguished quadratic
    weights must be nonpositive and their eigenvalues positive and no larger
    than the distinguished one.
    """

    square = as_fraction(positive_square)
    weights = tuple(as_fraction(value) for value in negative_place_weights)
    lambdas = tuple(as_fraction(value) for value in eigenvalues)
    if square <= 0:
        raise ValueError("positive_square must be positive")
    if len(lambdas) != len(weights) + 1:
        raise ValueError("one eigenvalue is required for every real place")
    if any(weight > 0 for weight in weights):
        raise ValueError("non-distinguished weights must be nonpositive")
    if any(value <= 0 for value in lambdas):
        raise ValueError("the generator must be totally positive")
    if any(value > lambdas[0] for value in lambdas[1:]):
        raise ValueError("the distinguished eigenvalue must be dominant")

    distinguished_weight = square - sum(weights, Fraction(0))
    mixed_height = lambdas[0] * distinguished_weight + sum(
        value * weight for value, weight in zip(lambdas[1:], weights)
    )
    gap = mixed_height - lambdas[0] * square
    expanded_gap = sum(
        (value - lambdas[0]) * weight
        for value, weight in zip(lambdas[1:], weights)
    )
    if gap != expanded_gap or gap < 0:
        raise AssertionError("signed spectral identity failed")
    return {
        "positive_square": str(square),
        "distinguished_weight": str(distinguished_weight),
        "negative_place_weights": [str(value) for value in weights],
        "eigenvalues": [str(value) for value in lambdas],
        "mixed_height": str(mixed_height),
        "spectral_gap": str(gap),
        "expanded_gap": str(expanded_gap),
        "dominant_spectral_inequality": True,
    }


def primitive_combination_bound(degree: int, generator_count: int) -> dict[str, int]:
    """Grid bound for a positive primitive combination of joint generators."""

    if degree < 2:
        raise ValueError("degree must be at least two")
    if generator_count < 1:
        raise ValueError("generator_count must be positive")
    embedding_pairs = degree * (degree - 1) // 2
    coefficient_bound = embedding_pairs + 1
    grid_size = coefficient_bound**generator_count
    union_bound = embedding_pairs * coefficient_bound ** (generator_count - 1)
    if union_bound >= grid_size:
        raise AssertionError("the hyperplane-avoidance grid did not close")
    return {
        "degree": degree,
        "generator_count": generator_count,
        "embedding_pairs": embedding_pairs,
        "positive_coefficient_bound": coefficient_bound,
        "grid_size": grid_size,
        "bad_grid_union_bound": union_bound,
        "good_grid_point_exists": True,
    }


def coefficient_bounds(degree: int, root_bound: int | Fraction) -> list[int]:
    if degree < 1:
        raise ValueError("degree must be positive")
    bound = as_fraction(root_bound)
    if bound < 1:
        raise ValueError("root_bound must be at least one")
    return [
        floor_fraction(Fraction(comb(degree, index)) * bound**index)
        for index in range(1, degree + 1)
    ]


def polynomial_search_count(degree_bound: int, root_bound: int | Fraction) -> int:
    """Crude finite count of monic integer coefficient packets."""

    if degree_bound < 1:
        raise ValueError("degree_bound must be positive")
    total = 0
    for degree in range(1, degree_bound + 1):
        count = 1
        for bound in coefficient_bounds(degree, root_bound):
            count *= 2 * bound + 1
        total += count
    return total


def squarefree_part(value: int) -> tuple[int, int]:
    if value <= 0:
        raise ValueError("value must be positive")
    # The largest square divisor leaves a squarefree quotient.
    for candidate in range(isqrt(value), 0, -1):
        if value % (candidate * candidate) == 0:
            return value // (candidate * candidate), candidate
    raise AssertionError("squarefree decomposition failed")


def fundamental_quadratic_discriminant(squarefree: int) -> int:
    if squarefree <= 1:
        raise ValueError("a nontrivial squarefree parameter is required")
    return squarefree if squarefree % 4 == 1 else 4 * squarefree


def quadratic_escape_record(a: int) -> dict[str, object]:
    """Replay the fixed-square rank-four K3 family against the new gate."""

    if a < 2:
        raise ValueError("a must be at least two")
    radicand = a * a - 1
    squarefree, pell_y = squarefree_part(radicand)
    if squarefree * pell_y * pell_y != radicand:
        raise AssertionError("Pell decomposition failed")
    square = 14
    height = 32 * a - 10
    ratio = Fraction(height, square)

    # ratio >= a+sqrt(a^2-1): both sides after subtracting a are positive,
    # and this is the exact numerator after squaring with denominator 49.
    dominance_margin_numerator = 32 * a * a - 90 * a + 74
    if ratio <= a or dominance_margin_numerator <= 0:
        raise AssertionError("rank-four spectral inequality failed")

    polynomial_discriminant = 4 * radicand
    discriminant_margin_numerator = height * height - square * square * polynomial_discriminant
    if discriminant_margin_numerator <= 0:
        raise AssertionError("quadratic discriminant envelope failed")

    return {
        "a": a,
        "field_squarefree_parameter": squarefree,
        "Pell_y": pell_y,
        "polarization_square": square,
        "mixed_correspondence_height": height,
        "height_over_square": str(ratio),
        "dominance_square_margin_numerator_over_49": dominance_margin_numerator,
        "minimal_polynomial": [1, -2 * a, 1],
        "polynomial_discriminant": polynomial_discriminant,
        "field_discriminant": fundamental_quadratic_discriminant(squarefree),
        "height_ratio_squared_discriminant_margin_numerator": discriminant_margin_numerator,
        "fixed_square_but_unbounded_mixed_height": True,
    }


def k3_rank_record(
    hyperbolic_dimension: int,
    generator_count: int,
    correspondence_height_bound: int,
    square_lower_bound: int = 2,
) -> dict[str, object]:
    if hyperbolic_dimension < 2:
        raise ValueError("this record is for hyperbolic dimension at least two")
    if correspondence_height_bound < square_lower_bound:
        raise ValueError("height bound cannot be smaller than the positive square gate")
    degree_bound = 20 // (hyperbolic_dimension + 1)
    if degree_bound < 2:
        raise ValueError("no non-rational full-trace K3 field fits this dimension")
    if generator_count < 1:
        raise ValueError("generator_count must be positive")
    pair_count = degree_bound * (degree_bound - 1) // 2
    coefficient_bound = pair_count + 1
    coordinate_root_bound = max(
        Fraction(4), Fraction(correspondence_height_bound, square_lower_bound)
    )
    primitive_root_bound = (
        generator_count * coefficient_bound * coordinate_root_bound
    )
    discriminant_bound = primitive_root_bound ** (
        degree_bound * (degree_bound - 1)
    )
    if discriminant_bound.denominator != 1:
        displayed_discriminant_bound = str(discriminant_bound)
    else:
        displayed_discriminant_bound = str(discriminant_bound.numerator)
    return {
        "hyperbolic_dimension": hyperbolic_dimension,
        "K3_Picard_rank_ceiling": 20,
        "field_degree_bound": degree_bound,
        "generator_count": generator_count,
        "positive_combination_coefficient_bound": coefficient_bound,
        "positive_square_lower_bound": square_lower_bound,
        "mixed_correspondence_height_bound": correspondence_height_bound,
        "individual_conjugate_bound": str(coordinate_root_bound),
        "primitive_generator_conjugate_bound": str(primitive_root_bound),
        "field_discriminant_bound": displayed_discriminant_bound,
        "finite_monic_polynomial_search_count": polynomial_search_count(
            degree_bound, primitive_root_bound
        ),
    }


def build_certificate(max_a: int = 128) -> dict[str, object]:
    if max_a < 2:
        raise ValueError("max_a must be at least two")
    samples = tuple(
        dict.fromkeys(value for value in (2, 3, 4, 7, 19, 64, max_a) if value <= max_a)
    )
    spectral_samples = (
        spectral_identity_record(2, (-1,), (7, 3)),
        spectral_identity_record(14, (-3, -5), (31, 7, 2)),
        spectral_identity_record(
            Fraction(9, 2),
            (Fraction(-2, 3), Fraction(-7, 5), Fraction(0)),
            (Fraction(23, 3), Fraction(5, 2), Fraction(1, 7), Fraction(23, 3)),
        ),
    )
    return {
        "theorem": "dominant correspondence height bounds the full-trace field",
        "symbolic_identity": "H.tau(H)-lambda0*H^2=sum_{sigma!=0}(lambda_sigma-lambda0)q(H_sigma)",
        "spectral_identity_records": list(spectral_samples),
        "primitive_combination_records": [
            primitive_combination_bound(2, 1),
            primitive_combination_bound(4, 2),
            primitive_combination_bound(6, 3),
        ],
        "K3_n_at_least_2_records": [
            k3_rank_record(2, 1, 32),
            k3_rank_record(2, 3, 32),
            k3_rank_record(3, 2, 32),
            k3_rank_record(5, 1, 32),
        ],
        "rank_four_escape_records": [quadratic_escape_record(a) for a in samples],
        "conclusions": {
            "one_positive_place_spectral_inequality": True,
            "bounded_positive_mixed_degree_bounds_dominant_eigenvalue": True,
            "bounded_rank_and_generator_packet_give_finite_field_set": True,
            "primitive_generator_gives_explicit_field_discriminant_bound": True,
            "pair_Cartan_packet_uses_conjugate_Cauchy_bound_four": True,
            "K3_n_at_least_two_field_degree_at_most_six": True,
            "finite_field_set_not_only_bounded_degree": True,
            "bounded_polarization_square_alone_suffices": False,
        },
        "scope": {
            "requires_exact_full_Picard_or_positive_H_in_trace_block": True,
            "arbitrary_field_generators_assumed_dominant": False,
            "pair_coordinate_integrality_from_diagram_alone": False,
            "proves_maximality_forces_the_height_gate": False,
            "constructs_nonorthogonal_rank_six_counterfamily": False,
            "rules_out_counterfamilies_with_unbounded_correspondence_degree": False,
            "uses_reflective_lattice_classification": False,
            "height_hypothesis_is_project_native": True,
        },
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def encoded_certificate(max_a: int = 128) -> str:
    payload = build_certificate(max_a=max_a)
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-a", type=int, default=128)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args()
    encoded = encoded_certificate(arguments.max_a)
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
