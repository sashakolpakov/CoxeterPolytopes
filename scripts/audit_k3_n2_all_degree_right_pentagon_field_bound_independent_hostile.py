#!/usr/bin/env python3
"""Independent hostile audit of the all-degree right-pentagon field bound."""

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
        "results/k3_n2_all_degree_right_pentagon_field_bound.md",
        "4177f65e13e4fdc19333d5b5ab6097799cd293d961a76aca685630d62cbf2bfc",
    ),
    "owner_certifier": (
        "scripts/certify_k3_n2_all_degree_right_pentagon_field_bound.py",
        "499a1ff95ea7f40ad789e83b8533c2fffcc064c0927a08d7fd6464698190907c",
    ),
    "owner_archive": (
        "data/enumerations/k3_n2_all_degree_right_pentagon_field_bound.json",
        "5c8c05c10a5f1c5c79dac0abc4df2be18447f743986a4bd7783f347c4b1d669f",
    ),
    "owner_tests": (
        "tests/test_k3_n2_all_degree_right_pentagon_field_bound.py",
        "fd16417cbfbce39bd467e1ea488cbdded18f435bba23569650b2e347ee12c88e",
    ),
    "local_scope_report": (
        "results/k3_n2_all_degree_right_pentagon_field_bound_hostile_scope_audit.md",
        "56b0781c6a8626fc2562ea93281084c42e5083c8b11a01344d539fc0a7bc678b",
    ),
    "local_scope_tests": (
        "tests/test_k3_n2_all_degree_right_pentagon_field_bound_hostile_scope.py",
        "bfe7e0b89c909b2124ce8564028a8edab0929fbe89943072e57669dca3a1fece",
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


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    size = len(matrix)
    assert all(len(row) == size for row in matrix)
    result = Fraction(0)
    for permutation in permutations(range(size)):
        term = Fraction(permutation_sign(permutation))
        for row, column in enumerate(permutation):
            term *= matrix[row][column]
        result += term
    return result


def tree_minor(first, middle, last) -> Fraction:
    first, middle, last = map(Fraction, (first, middle, last))
    matrix = (
        (Fraction(1), -first, Fraction(0), Fraction(0)),
        (-first, Fraction(1), -middle, Fraction(0)),
        (Fraction(0), -middle, Fraction(1), -last),
        (Fraction(0), Fraction(0), -last, Fraction(1)),
    )
    direct = determinant(matrix)
    claimed = 1 - first**2 - middle**2 - last**2 + first**2 * last**2
    assert direct == claimed
    return direct


def gram_recurrence_certificate() -> dict[str, object]:
    tested = 0
    for first in (Fraction(1, 3), Fraction(1, 2), Fraction(3, 2), Fraction(7, 3)):
        for middle in (Fraction(1, 4), Fraction(2, 3), Fraction(5, 4)):
            for last in (Fraction(1, 5), Fraction(4, 5), Fraction(9, 5)):
                tree_minor(first, middle, last)
                tested += 1
    assert tested == 36

    # For a vanishing minor put A_-=first^2, A=middle^2, A_+=last^2.
    # The determinant reads (A_--1)(A_+-1)=A.
    samples = (
        (Fraction(1, 8), Fraction(9)),
        (Fraction(1, 2), Fraction(3)),
        (Fraction(1), Fraction(2)),
        (Fraction(7, 3), Fraction(11, 5)),
    )
    closure_rows = []
    for x0, x2 in samples:
        x1 = x0 * x2 - 1
        assert x1 > 0
        coordinates = (
            x0,
            x1,
            x2,
            (x2 + 1) / x1,
            (x0 + 1) / x1,
        )
        residuals = tuple(
            coordinates[(index - 1) % 5] * coordinates[(index + 1) % 5]
            - coordinates[index]
            - 1
            for index in range(5)
        )
        assert residuals == (0,) * 5
        closure_rows.append([str(value) for value in coordinates])
    return {
        "independent_tree_determinants_checked": tested,
        "tree_determinant": "1-A_left-A_mid-A_right+A_left*A_right",
        "cyclic_recurrence": "x[i-1]*x[i+1]=x[i]+1",
        "all_five_relations_checked": True,
        "positive_closure_samples": closure_rows,
    }


def norm_and_upper_bound_certificate() -> dict[str, object]:
    rows = []
    for degree in range(2, 7):
        q = 4**degree
        epsilon = Fraction(1, q)
        upper = q**2 + 2 * q
        # In x < eps^-1 + (1+eps)/(eps*z), the strict z>eps
        # supplies the second occurrence of the same all-degree floor.
        substituted = Fraction(1, epsilon) + (1 + epsilon) / (epsilon**2)
        assert substituted == upper
        rows.append(
            {
                "degree": degree,
                "norm_floor": f"1/{q}",
                "floor_used_for_x_i_plus_3": f"1/{q}",
                "floor_used_for_z": f"1/{q}",
                "strict_shift_upper": upper,
                "strict_A_upper": (q + 1) ** 2,
            }
        )
    terminal = rows[-1]
    assert terminal["strict_shift_upper"] == 16_785_408
    assert terminal["strict_A_upper"] == 16_785_409
    return {
        "degree_bound_from_3d_le_20": 6,
        "integral_norm_element": "4*x_i",
        "definite_conjugate_interval": "(-4,0)",
        "rows": rows,
        "degree_six_shift_upper": terminal["strict_shift_upper"],
        "degree_six_A_upper": terminal["strict_A_upper"],
        "invalid_mixed_floor_value_not_used": 69_648,
        "same_epsilon_used_twice": True,
    }


def cyclic_field_packet_certificate() -> dict[str, object]:
    q = 4**6
    A_upper = (q + 1) ** 2
    five_cycle_height = 32 * (q + 1) ** 5
    assert A_upper == 16_785_409
    assert five_cycle_height == 36_938_546_139_294_728_224

    # Symmetry identifies the two directed five-cycles and gives
    # c5^2=product of the five pair cycles=4^5 product(A_i).
    minimal_generators = 2
    bad_values = comb(6, 2)
    choices = bad_values + 1
    primitive_height = choices * five_cycle_height
    difference_height = 2 * primitive_height
    exponent = 6 * 5
    assert (bad_values, choices, difference_height, exponent) == (
        15,
        16,
        1_182_033_476_457_431_303_168,
        30,
    )
    return {
        "nonzero_Cartan_graph": "C5",
        "simple_cycle_packet": "five pair cycles plus one directed length-five cycle",
        "cyclic_integrality_source": "traces of products of T_i=1-s_i in End_O(M)",
        "individual_directed_entries_assumed_integral": False,
        "closed_walks_factor_into_simple_cycles": True,
        "five_cycle_square_identity": "c5^2=4^5*product(A_i)",
        "uniform_packet_height": five_cycle_height,
        "complete_packet_generates_k_is_hypothesis": True,
        "single_pair_generates_for_degree_gt_2": False,
        "minimal_generating_subpacket_size_upper_bound": minimal_generators,
        "primitive_grid": "theta1+c*theta2, c=0,...,15",
        "Vandermonde_difference_height": difference_height,
        "Vandermonde_exponent": exponent,
        "field_discriminant_bound": f"{difference_height}^{exponent}",
        "finite_field_set": True,
    }


def quadratic_certificate() -> dict[str, object]:
    q = 16
    shift_upper = q**2 + 2 * q
    pair_upper = 4 * ((q + 1) ** 2)
    trace_upper = pair_upper + 4 - 1
    discriminant_upper = trace_upper**2 - 4
    assert (shift_upper, pair_upper, trace_upper, discriminant_upper) == (
        288,
        1156,
        1159,
        1_343_277,
    )
    return {
        "shift_interval": "(1/16,288)",
        "pair_distinguished_interval": "(4,1156)",
        "pair_definite_interval": "(0,4)",
        "disjoint_intervals_force_quadratic_generation": True,
        "field_discriminant_upper_bound": discriminant_upper,
    }


def scope_certificate() -> dict[str, object]:
    text = " ".join((ROOT / PINS["owner_report"][0]).read_text().split())
    required = {
        "right_compact": "exactly a compact right-angled Coxeter pentagon",
        "full_trace": "primitive full trace",
        "cyclic_ground_field": "ground field generated by the cyclic Cartan products",
        "complete": "complete and saturated",
        "maximal": "reflection-maximal among integral",
        "mixed_floor_warning": r"16+17\cdot4^6=69,648",
        "nonright_open": "treat non-right or noncompact pentagons",
        "arbitrary_open": "chamber is a right pentagon",
        "gratuitous_extension_open": "gratuitous declared field extension",
    }
    for label, phrase in required.items():
        assert phrase in text, (label, phrase)
    return {label: True for label in required}


def build_audit() -> dict[str, object]:
    return {
        "schema": "k3_n2_all_degree_right_pentagon_field_bound/independent_hostile/v1",
        "verdict": "PASS_CONDITIONAL_EXACT_SCOPE",
        "pins": verify_pins(),
        "five_Gram_closure": gram_recurrence_certificate(),
        "all_place_floor_and_recurrence_bound": norm_and_upper_bound_certificate(),
        "complete_cyclic_ground_field_packet": cyclic_field_packet_certificate(),
        "quadratic_specialization": quadratic_certificate(),
        "scope": scope_certificate(),
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
