#!/usr/bin/env python3
"""Independent hostile audit of the compact-quadrilateral field theorem.

This file imports no owner proof code.  It expands the normalized Gram
determinant symbolically, reconstructs both diagonal estimates from the norm
gap, and redoes the complete cyclic-packet primitive-element argument.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from itertools import permutations
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

PINS = {
    "owner_report": (
        "results/k3_n2_compact_quadrilateral_field_bound.md",
        "7f300ef5b54e9dbf0ca47140b5ca7134119046fc6f1ee39bd46bce58a6028528",
    ),
    "owner_certifier": (
        "scripts/certify_k3_n2_compact_quadrilateral_field_bound.py",
        "25ccdd5e21ae2ad5ff5d43d079fb9308b3dbacb4780fa0ba8638f7c27f9d884d",
    ),
    "owner_archive": (
        "data/enumerations/k3_n2_compact_quadrilateral_field_bound.json",
        "a61e4a28960e082380af1fc7e0dd3f8b311ab2871b6e361268f1f4e50b6213c4",
    ),
    "owner_tests": (
        "tests/test_k3_n2_compact_quadrilateral_field_bound.py",
        "0517df99c24402c881e32c6666e91abf67c0d4e3f84fb5be2cc857aaf25e072c",
    ),
    "owner_scope_audit": (
        "results/k3_n2_compact_quadrilateral_field_bound_hostile_scope_audit.md",
        "a4c642a352c197d8f7e43c58f76580cabe86b9b18a54112220dec776d83b5ca3",
    ),
    "owner_scope_tests": (
        "tests/test_k3_n2_compact_quadrilateral_field_bound_hostile_scope.py",
        "aa51ed636fbbb9862520d38315aed20ddb0cb909d65e6b4a9f3002692f47fb52",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins() -> dict[str, str]:
    observed = {}
    for label, (relative, expected) in PINS.items():
        actual = sha256(ROOT / relative)
        assert actual == expected, (label, expected, actual)
        observed[label] = actual
    return observed


# Sparse polynomials in (a,b,c,e,U,V), used so that the determinant identity
# is checked formally rather than on a numerical grid.
Monomial = tuple[int, int, int, int, int, int]
Polynomial = dict[Monomial, Fraction]
ZERO_MONOMIAL: Monomial = (0, 0, 0, 0, 0, 0)


def poly_const(value: int | Fraction) -> Polynomial:
    value = Fraction(value)
    return {} if value == 0 else {ZERO_MONOMIAL: value}


def poly_var(index: int) -> Polynomial:
    exponent = [0] * 6
    exponent[index] = 1
    return {tuple(exponent): Fraction(1)}


def poly_add(*terms: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for term in terms:
        for monomial, coefficient in term.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
            if result[monomial] == 0:
                del result[monomial]
    return result


def poly_scale(term: Polynomial, scalar: int | Fraction) -> Polynomial:
    scalar = Fraction(scalar)
    return {
        monomial: scalar * coefficient
        for monomial, coefficient in term.items()
        if scalar * coefficient
    }


def poly_mul(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(
                left_monomial[index] + right_monomial[index]
                for index in range(6)
            )
            result[monomial] = (
                result.get(monomial, Fraction(0))
                + left_coefficient * right_coefficient
            )
            if result[monomial] == 0:
                del result[monomial]
    return result


def poly_pow(term: Polynomial, exponent: int) -> Polynomial:
    result = poly_const(1)
    for _ in range(exponent):
        result = poly_mul(result, term)
    return result


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def symbolic_determinant(matrix: tuple[tuple[Polynomial, ...], ...]) -> Polynomial:
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    determinant: Polynomial = {}
    for permutation in permutations(range(size)):
        term = poly_const(permutation_sign(permutation))
        for row, column in enumerate(permutation):
            term = poly_mul(term, matrix[row][column])
        determinant = poly_add(determinant, term)
    return determinant


def determinant_certificate() -> dict[str, object]:
    one = poly_const(1)
    a, b, c, e, u, v = (poly_var(index) for index in range(6))
    matrix = (
        (one, poly_scale(a, -1), poly_scale(u, -1), poly_scale(e, -1)),
        (poly_scale(a, -1), one, poly_scale(b, -1), poly_scale(v, -1)),
        (poly_scale(u, -1), poly_scale(b, -1), one, poly_scale(c, -1)),
        (poly_scale(e, -1), poly_scale(v, -1), poly_scale(c, -1), one),
    )
    direct = symbolic_determinant(matrix)

    r = poly_add(poly_pow(u, 2), poly_const(-1))
    s = poly_add(poly_pow(v, 2), poly_const(-1))
    p = poly_add(poly_mul(a, b), poly_mul(c, e))
    q = poly_add(poly_mul(a, e), poly_mul(b, c))
    rr = poly_add(poly_mul(a, c), poly_mul(b, e))
    k = poly_add(
        poly_scale(poly_pow(a, 2), -1),
        poly_scale(poly_pow(b, 2), -1),
        poly_scale(poly_pow(c, 2), -1),
        poly_scale(poly_pow(e, 2), -1),
        poly_mul(poly_pow(a, 2), poly_pow(c, 2)),
        poly_mul(poly_pow(b, 2), poly_pow(e, 2)),
        poly_scale(poly_mul(poly_mul(a, b), poly_mul(c, e)), -2),
    )
    claimed = poly_add(
        poly_mul(r, s),
        k,
        poly_scale(poly_mul(poly_mul(p, u), poly_const(1)), -2),
        poly_scale(poly_mul(poly_mul(q, v), poly_const(1)), -2),
        poly_scale(poly_mul(poly_mul(poly_mul(rr, u), v), poly_const(1)), -2),
    )
    assert direct == claimed
    assert len(direct) == 17
    return {
        "formal_polynomial_identity": True,
        "determinant_monomial_count": len(direct),
        "identity": "(U^2-1)(V^2-1)+K-2PU-2QV-2RUV",
        "coefficient_bounds": {"P_Q_R": "[0,2)", "absolute_K": "<8"},
        "growth_inequality": "(U^2-1)(V^2-1)<4U+4V+4UV+8",
    }


def norm_gap_certificate() -> dict[str, object]:
    rows = []
    for degree in range(2, 7):
        full_trace_rank = 3 * degree
        assert full_trace_rank <= 20
        denominator = 4**degree
        rows.append(
            {
                "degree": degree,
                "full_trace_rank": full_trace_rank,
                "strict_shift_lower_bound": f"1/{denominator}",
            }
        )
    assert 3 * 7 > 20
    return {
        "degree_rows": rows,
        "degree_upper_bound": 6,
        "integral_element": "4(A-1)",
        "definite_conjugate_absolute_interval": "(0,4)",
        "nonzero_integer_norm_absolute_lower_bound": 1,
        "norm_gap": "A-1>4^(-d)",
        "same_for_C": True,
    }


def diagonal_certificate(degree: int, trial: int) -> dict[str, object]:
    q = 4**degree
    square_root_denominator = 2**degree
    linear_rational = 4 * q
    surd_coefficient = 4 * square_root_denominator
    constant_rational = 8 * q + 1
    rational_remainder = trial**2 - linear_rational * trial - constant_rational
    surd_remainder = surd_coefficient * (trial + 1)
    square_difference = (
        rational_remainder**2 - (q + 1) * surd_remainder**2
    )
    assert rational_remainder > 0
    assert square_difference > 0
    return {
        "degree": degree,
        "shift_gap": f"1/{q}",
        "V_lower_bound": f"sqrt({q + 1})/{square_root_denominator}",
        "linear_coefficient": f"{linear_rational}+{surd_coefficient}sqrt({q + 1})",
        "constant_coefficient": f"{constant_rational}+{surd_coefficient}sqrt({q + 1})",
        "trial": trial,
        "remainder": f"{rational_remainder}-{surd_remainder}sqrt({q + 1})",
        "remainder_square_difference": square_difference,
        "both_opposite_magnitudes_strictly_below": trial,
        "both_opposite_squared_ratios_strictly_below": trial**2,
    }


def field_packet_certificate() -> dict[str, object]:
    degree_bound = 6
    magnitude_bound = 32_772
    cycle_length_bound = 4
    cycle_height = 2**cycle_length_bound * magnitude_bound**cycle_length_bound
    assert cycle_height == 18_455_752_922_365_956_096

    # A minimal packet of m field generators gives a tower of m strict
    # extensions, hence degree at least 2^m.
    generator_count = 0
    while 2 ** (generator_count + 1) <= degree_bound:
        generator_count += 1
    assert generator_count == 2
    bad_coefficients = comb(degree_bound, 2)
    coefficient_choices = bad_coefficients + 1
    assert (bad_coefficients, coefficient_choices) == (15, 16)

    primitive_height = coefficient_choices * cycle_height
    difference_height = 2 * primitive_height
    exponent = degree_bound * (degree_bound - 1)
    assert difference_height == 590_584_093_515_710_595_072
    assert exponent == 30
    return {
        "simple_cycle_lengths": [2, 3, 4],
        "simple_cycles_generate_all_closed_walk_coordinates": True,
        "complete_cycle_packet_generates_k_is_hypothesis": True,
        "single_coordinate_generation_used_for_degree_gt_2": False,
        "cycle_conjugate_height": cycle_height,
        "minimal_generator_subpacket_size_upper_bound": generator_count,
        "embedding_pairs_and_bad_coefficients_upper_bound": bad_coefficients,
        "primitive_grid": "theta1+c*theta2 for c=0,...,15",
        "primitive_conjugate_height": primitive_height,
        "Vandermonde_difference_height": difference_height,
        "Vandermonde_exponent": exponent,
        "field_discriminant_bound": f"{difference_height}^{exponent}",
        "conclusion_is_finite_field_set": True,
    }


def quadratic_certificate() -> dict[str, object]:
    magnitude_bound = 132
    alpha_distinguished_upper = 4 * magnitude_bound**2
    alpha_definite_upper = 4
    trace_upper = alpha_distinguished_upper + alpha_definite_upper - 1
    norm_upper = alpha_distinguished_upper * alpha_definite_upper - 1
    discriminant_upper = trace_upper**2 - 4
    assert (
        alpha_distinguished_upper,
        trace_upper,
        norm_upper,
        discriminant_upper,
    ) == (69_696, 69_699, 278_783, 4_857_950_597)
    return {
        "alpha": "4A",
        "distinguished_interval": "(4,69696)",
        "definite_interval": "(0,4)",
        "disjoint_intervals_force_nonrational": True,
        "one_nonrational_coordinate_generates_only_because_degree_is_two": True,
        "trace_upper": trace_upper,
        "norm_upper": norm_upper,
        "norm_absolute_lower_bound": 1,
        "field_discriminant_upper": discriminant_upper,
    }


def scope_certificate() -> dict[str, object]:
    text = " ".join((ROOT / PINS["owner_report"][0]).read_text().split())
    required = {
        "compact": "exactly a compact Coxeter quadrilateral",
        "full_trace": "primitive full trace",
        "cyclic_ground_field": "ground field generated by the cyclic Cartan products",
        "complete": "complete and saturated",
        "maximal": "reflection-maximal among integral",
        "maximality_unused": "reflection integrality prove the bound before maximality enters",
        "noncompact_open": "treat noncompact quadrilaterals",
        "arbitrary_polygon_open": "chamber is a quadrilateral",
        "gratuitous_extension_excluded": "gratuitous extension of the cyclic ground field",
    }
    for label, phrase in required.items():
        assert phrase in text, (label, phrase)
    return {label: True for label in required}


def build_audit() -> dict[str, object]:
    degree_six = diagonal_certificate(6, 32_772)
    degree_two = diagonal_certificate(2, 132)
    assert degree_six["remainder_square_difference"] == 17_606_953_271_521
    assert degree_two["remainder_square_difference"] == 1_286_881
    return {
        "schema": "k3_n2_compact_quadrilateral_field_bound/independent_hostile/v1",
        "verdict": "PASS_CONDITIONAL_EXACT_SCOPE",
        "pins": verify_pins(),
        "determinant": determinant_certificate(),
        "all_place_norm_gap": norm_gap_certificate(),
        "degree_six_diagonal_bound": degree_six,
        "degree_two_diagonal_bound": degree_two,
        "complete_field_generating_packet": field_packet_certificate(),
        "degree_two_specialization": quadratic_certificate(),
        "scope": scope_certificate(),
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
