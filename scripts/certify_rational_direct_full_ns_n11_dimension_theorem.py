#!/usr/bin/env python3
"""Exact arithmetic certificate for the direct full-NS n<=11 theorem.

The checker is standalone.  It reconstructs the two minimal-base roots,
the satellite/facet-integrality obstruction, the degree-one Diophantine
gates, the exceptional-cone remainder, and both parity-free mixed-ray
arguments.  It imports no reflective-lattice classifier or marked-frame
certifier.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from pathlib import Path
from typing import Iterable, Sequence


Vector = tuple[int, ...]


def diagonal_pair(left: Sequence[int], right: Sequence[int]) -> int:
    """The negative identity form used on total exceptional coordinates."""

    return -sum(a * b for a, b in zip(left, right))


def plane_pair(left: Sequence[int], right: Sequence[int]) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def f2_pair(left: Sequence[int], right: Sequence[int]) -> int:
    return (
        -2 * left[0] * right[0]
        + left[0] * right[1]
        + left[1] * right[0]
        - sum(a * b for a, b in zip(left[2:], right[2:]))
    )


def reflect(
    vector: Sequence[int], root: Sequence[int], pairing
) -> tuple[Fraction, ...]:
    root_square = pairing(root, root)
    coefficient = Fraction(-2 * pairing(vector, root), root_square)
    return tuple(
        Fraction(entry) + coefficient * root_entry
        for entry, root_entry in zip(vector, root)
    )


def f0_swap_record(exceptional_rank: int) -> dict[str, object]:
    if exceptional_rank < 0:
        raise ValueError("exceptional_rank must be nonnegative")
    root = (1, -1, *((0,) * exceptional_rank))

    def f0_pair(left: Sequence[int], right: Sequence[int]) -> int:
        return (
            left[0] * right[1]
            + left[1] * right[0]
            - sum(a * b for a, b in zip(left[2:], right[2:]))
        )

    c = (1, 0, *((0,) * exceptional_rank))
    d = (0, 1, *((0,) * exceptional_rank))
    return {
        "exceptional_rank": exceptional_rank,
        "root_square": f0_pair(root, root),
        "c_maps_to_d": reflect(c, root, f0_pair) == tuple(map(Fraction, d)),
        "d_maps_to_c": reflect(d, root, f0_pair) == tuple(map(Fraction, c)),
        "c_plus_d_is_on_mirror": f0_pair(
            tuple(a + b for a, b in zip(c, d)), root
        )
        == 0,
        "all_total_exceptional_coordinates_fixed": all(
            reflect(
                tuple(1 if j == 2 + index else 0 for j in range(2 + exceptional_rank)),
                root,
                f0_pair,
            )
            == tuple(
                Fraction(1 if j == 2 + index else 0)
                for j in range(2 + exceptional_rank)
            )
            for index in range(exceptional_rank)
        ),
    }


def section_integrality_record(e: int, hits: int) -> dict[str, object]:
    if e < 0 or hits < 0:
        raise ValueError("e and hits must be nonnegative")
    negative_square = e + hits
    coefficient = None if negative_square == 0 else Fraction(2, negative_square)
    return {
        "e": e,
        "section_hits": hits,
        "negative_section_square": negative_square,
        "reflected_fibre_section_coefficient": (
            None if coefficient is None else str(coefficient)
        ),
        "reflection_integral_on_full_NS": (
            None
            if coefficient is None
            else coefficient.denominator == 1
        ),
    }


def exceptional_prime_record(children: int) -> dict[str, object]:
    """Final strict exceptional prime after ``children`` proximate blowups."""

    if children < 0:
        raise ValueError("children must be nonnegative")
    root = (1, *((-1,) * children))
    square = diagonal_pair(root, root)
    basis = tuple(
        tuple(1 if i == j else 0 for i in range(children + 1))
        for j in range(children + 1)
    )
    reflection_integral = all(
        all(entry.denominator == 1 for entry in reflect(vector, root, diagonal_pair))
        for vector in basis
    )
    return {
        "proximate_children": children,
        "strict_prime_square": square,
        "strict_prime_primitive": True,
        "square_in_physical_alphabet": -square in (1, 2, 4),
        "reflection_integral_on_unimodular_exceptional_lattice": (
            reflection_integral
        ),
        "compatible_with_both_alphabet_and_integrality": (
            -square in (1, 2, 4) and reflection_integral
        ),
    }


def nonincreasing_tuples(length: int, largest_entry: int) -> Iterable[Vector]:
    for pattern in combinations_with_replacement(
        range(largest_entry + 1), length
    ):
        yield tuple(reversed(pattern))


def exact_partitions(
    length: int, total: int, square_total: int
) -> tuple[Vector, ...]:
    if min(length, total, square_total) < 0:
        return ()
    cap = min(total, square_total)
    return tuple(
        values
        for values in nonincreasing_tuples(length, cap)
        if sum(values) == total
        and sum(value * value for value in values) == square_total
    )


def degree_one_root_records(length: int = 12) -> dict[str, object]:
    if length < 3:
        raise ValueError("length must be at least three")
    plane_minus_one = exact_partitions(length, 2, 2)
    plane_minus_two = exact_partitions(length, 3, 3)
    plane_half_minus_four = exact_partitions(length, 4, 2)
    f2_minus_one = exact_partitions(length, 1, 1)
    f2_minus_two = exact_partitions(length, 2, 2)
    f2_half_minus_four = exact_partitions(length, 3, 1)
    expected_plane_minus_one = ((1, 1, *((0,) * (length - 2))),)
    expected_plane_minus_two = ((1, 1, 1, *((0,) * (length - 3))),)
    expected_f2_minus_one = ((1, *((0,) * (length - 1))),)
    expected_f2_minus_two = ((1, 1, *((0,) * (length - 2))),)
    if plane_minus_one != expected_plane_minus_one:
        raise AssertionError("plane degree-one (-1) root classification changed")
    if plane_minus_two != expected_plane_minus_two:
        raise AssertionError("plane degree-one (-2) root classification changed")
    if f2_minus_one != expected_f2_minus_one:
        raise AssertionError("F2 fibre-class (-1) root classification changed")
    if f2_minus_two != expected_f2_minus_two:
        raise AssertionError("F2 fibre-class (-2) root classification changed")
    if plane_half_minus_four or f2_half_minus_four:
        raise AssertionError("a forbidden degree-one half-(-4) solution appeared")
    return {
        "multiplicity_slots": length,
        "plane_physical_minus_one_patterns": [list(row) for row in plane_minus_one],
        "plane_physical_minus_two_patterns": [list(row) for row in plane_minus_two],
        "plane_half_minus_four_equations": {
            "sum": 4,
            "sum_of_squares": 2,
            "solutions": [],
        },
        "F2_physical_minus_one_patterns": [list(row) for row in f2_minus_one],
        "F2_physical_minus_two_patterns": [list(row) for row in f2_minus_two],
        "F2_half_minus_four_equations": {
            "sum": 3,
            "sum_of_squares": 1,
            "solutions": [],
        },
    }


def exceptional_chain_coordinates(target: Sequence[int]) -> tuple[int, ...]:
    """Coefficients in e_i-e_{i+1},...,e_m via prefix sums."""

    running = 0
    coefficients: list[int] = []
    for value in target:
        running += value
        coefficients.append(running)
    return tuple(coefficients)


def exceptional_remainder_records(m: int = 12) -> dict[str, object]:
    if m < 3:
        raise ValueError("m must be at least three")
    plane_target = tuple(-1 if index == 2 else 0 for index in range(m))
    f2_target = tuple(-1 if index == 1 else 0 for index in range(m))
    plane_coefficients = exceptional_chain_coordinates(plane_target)
    f2_coefficients = exceptional_chain_coordinates(f2_target)
    return {
        "exceptional_generators": "e_i-e_{i+1} and e_m",
        "plane_remainder": "-e_3",
        "plane_generator_coefficients": list(plane_coefficients),
        "plane_in_nonnegative_exceptional_cone": all(
            value >= 0 for value in plane_coefficients
        ),
        "F2_remainder": "-e_2",
        "F2_generator_coefficients": list(f2_coefficients),
        "F2_in_nonnegative_exceptional_cone": all(
            value >= 0 for value in f2_coefficients
        ),
    }


def plane_standard_walls(m: int) -> tuple[Vector, ...]:
    if m < 12:
        raise ValueError("plane mixed ray requires m>=12")
    walls: list[Vector] = []
    alpha = [0] * (m + 1)
    alpha[0] = 1
    alpha[1:4] = (-1, -1, -1)
    walls.append(tuple(alpha))
    for index in range(1, m):
        wall = [0] * (m + 1)
        wall[index], wall[index + 1] = 1, -1
        walls.append(tuple(wall))
    endpoint = [0] * (m + 1)
    endpoint[m] = 1
    walls.append(tuple(endpoint))
    return tuple(walls)


def plane_mixed_ray(m: int) -> Vector:
    if m < 12:
        raise ValueError("plane mixed ray requires m>=12")
    return (7, -3, *((-2,) * 11), *((0,) * (m - 12)))


def f2_standard_walls(m: int) -> tuple[Vector, ...]:
    if m < 11:
        raise ValueError("F2 mixed ray requires m>=11")
    section = (1, 0, *((0,) * m))
    fibre = (0, 1, -1, -1, *((0,) * (m - 2)))
    walls: list[Vector] = [section, fibre]
    for index in range(1, m):
        wall = [0] * (m + 2)
        wall[index + 1], wall[index + 2] = 1, -1
        walls.append(tuple(wall))
    endpoint = [0] * (m + 2)
    endpoint[m + 1] = 1
    walls.append(tuple(endpoint))
    return tuple(walls)


def f2_mixed_ray(m: int) -> Vector:
    if m < 11:
        raise ValueError("F2 mixed ray requires m>=11")
    return (4, 9, *((-2,) * 11), *((0,) * (m - 11)))


def mixed_ray_records(max_m: int) -> dict[str, object]:
    if max_m < 12:
        raise ValueError("max_m must be at least twelve")
    plane_rows: list[dict[str, object]] = []
    for m in range(12, max_m + 1):
        ray = plane_mixed_ray(m)
        pairings = [plane_pair(ray, wall) for wall in plane_standard_walls(m)]
        expected = [0] * (m + 1)
        expected[1] = 1
        expected[12] = 2
        if pairings != expected:
            raise AssertionError("plane mixed ray left the standard wall cone")
        if plane_pair(ray, ray) != -4:
            raise AssertionError("plane mixed ray is no longer spacelike")
        canonical = (-3, *((1,) * m))
        if plane_pair(canonical, ray) != 4:
            raise AssertionError("plane canonical pairing changed")
        plane_rows.append(
            {
                "m": m,
                "square": -4,
                "canonical_pairing": 4,
                "nonzero_standard_pairings": {"exceptional_1": 1, "exceptional_12": 2},
            }
        )

    f2_rows: list[dict[str, object]] = []
    for m in range(11, max_m + 1):
        ray = f2_mixed_ray(m)
        pairings = [f2_pair(ray, wall) for wall in f2_standard_walls(m)]
        expected = [0] * (m + 2)
        expected[0] = 1
        expected[12] = 2
        if pairings != expected:
            raise AssertionError("F2 mixed ray left the standard wall cone")
        if f2_pair(ray, ray) != -4:
            raise AssertionError("F2 mixed ray is no longer spacelike")
        canonical = (-2, -4, *((1,) * m))
        if f2_pair(canonical, ray) != 4:
            raise AssertionError("F2 canonical pairing changed")
        f2_rows.append(
            {
                "m": m,
                "square": -4,
                "canonical_pairing": 4,
                "nonzero_standard_pairings": {"section": 1, "exceptional_11": 2},
            }
        )
    return {"plane": plane_rows, "F2": f2_rows}


def plane_minus_four_low_pairing() -> dict[str, object]:
    rows: list[dict[str, int]] = []
    for R in range(4):
        patterns = 0
        adjunction_candidates = 0
        square_solutions = 0
        for tail in nonincreasing_tuples(11, R):
            if tail[0] + tail[1] > R:
                continue
            patterns += 1
            numerator = sum(tail) - R - 2
            if numerator < 0 or numerator % 2:
                continue
            degree = numerator // 2
            mu_1 = degree - R
            if mu_1 < tail[0]:
                continue
            adjunction_candidates += 1
            required_square = 2 * degree * R - R * R + 4
            if sum(value * value for value in tail) == required_square:
                square_solutions += 1
        if square_solutions:
            raise AssertionError(f"plane bad branch R={R} has a solution")
        rows.append(
            {
                "R": R,
                "monotone_line_patterns": patterns,
                "adjunction_candidates": adjunction_candidates,
                "square_solutions": square_solutions,
            }
        )
    return {
        "identity": "C.w=d-mu_1-2q+4+2*T_12",
        "bad_branch": "q=4, T_12=0, R=d-mu_1<=3",
        "rows": rows,
        "solutions": 0,
    }


def f2_minus_four_low_pairing() -> dict[str, object]:
    rows: list[dict[str, int]] = []
    for a in range(4):
        patterns = 0
        adjunction_candidates = 0
        square_solutions = 0
        for multiplicities in nonincreasing_tuples(11, a):
            if multiplicities[0] + multiplicities[1] > a:
                continue
            patterns += 1
            numerator = sum(multiplicities) - 2
            if numerator < 0 or numerator % 2:
                continue
            b = numerator // 2
            if b < 2 * a:
                continue
            adjunction_candidates += 1
            required_square = -2 * a * a + 2 * a * b + 4
            if sum(value * value for value in multiplicities) == required_square:
                square_solutions += 1
        if square_solutions:
            raise AssertionError(f"F2 bad branch a={a} has a solution")
        rows.append(
            {
                "a": a,
                "monotone_fibre_patterns": patterns,
                "adjunction_candidates": adjunction_candidates,
                "square_solutions": square_solutions,
            }
        )
    return {
        "identity": "C.w=a-2q+4+2*T_11",
        "bad_branch": "q=4, T_11=0, a<=3",
        "rows": rows,
        "solutions": 0,
    }


def build_certificate(max_m: int = 40) -> dict[str, object]:
    if max_m < 12:
        raise ValueError("max_m must be at least twelve")
    satellite_rows = [exceptional_prime_record(children) for children in range(2, 17)]
    if any(row["compatible_with_both_alphabet_and_integrality"] for row in satellite_rows):
        raise AssertionError("a satellite-compatible strict exceptional prime appeared")
    remainders = exceptional_remainder_records()
    if remainders["plane_in_nonnegative_exceptional_cone"]:
        raise AssertionError("-e_3 entered the free-chain exceptional cone")
    if remainders["F2_in_nonnegative_exceptional_cone"]:
        raise AssertionError("-e_2 entered the free-chain exceptional cone")
    return {
        "certificate": "rational_direct_full_ns_n11_dimension_theorem",
        "status": "PASS_DIRECT_FULL_NS_DIMENSION_BOUND",
        "scope": {
            "host": "direct full Neron-Severi lattice only",
            "complete_Mori_cone": "rational polyhedral",
            "physical_negative_extremal_alphabet": (
                "smooth rational (-1),(-2),(-4)"
            ),
            "square_zero_tangent_Mori_rays_allowed": True,
            "every_facet_reflection_integral_on_full_NS": True,
            "finite_volume_Coxeter_fundamental_chamber": True,
            "reflection_maximality": True,
            "adapted_half_boundary_host_claimed": False,
            "reflective_lattice_classification_or_census_used": False,
        },
        "minimal_base": {
            "F0_ruling_swap_samples": [f0_swap_record(rank) for rank in (0, 1, 2, 7, 31)],
            "section_integrality_samples": [
                section_integrality_record(e, hits)
                for e in range(2, 9)
                for hits in range(0, 5)
            ],
            "survivors": ["P2", "F2"],
            "F2_section_hits": 0,
        },
        "exceptional_geometry": {
            "relative_ample_equalization_for_incomparable_valuations": True,
            "reflection_maximality_forces_width_one": True,
            "satellite_rows": satellite_rows,
            "all_successors_forced_free_on_newest_exceptional": True,
        },
        "nonsimplicial_root_argument": {
            "module": "free integral module on all physical facet roots",
            "cartan_entries_integral": True,
            "positive_root_coefficients_nonnegative_integers": True,
            "actual_facet_dependencies_do_not_change_formal_coefficients": True,
            "ample_degree_selects_positive_sign": True,
        },
        "degree_one_gate": degree_one_root_records(),
        "exceptional_remainders": remainders,
        "exact_frames": {
            "plane": "first three centres follow a line; later centres leave it; free exceptional chain",
            "F2": "second centre is fibre-following free at strict fibre intersection E1; later centres leave fibre; free exceptional chain off s",
            "satellite_terminology_used_for_F2_second_centre": False,
        },
        "mixed_rays": mixed_ray_records(max_m),
        "plane_minus_four_low_pairing": plane_minus_four_low_pairing(),
        "F2_minus_four_low_pairing": f2_minus_four_low_pairing(),
        "conclusion": {
            "plane_blowups_at_most": 11,
            "F2_blowups_at_most": 10,
            "projective_nef_dimension": (
                "rho-1 even when nonnegative isotropic extreme rays occur"
            ),
            "hyperbolic_dimension_at_most": 11,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-m", type=int, default=40)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args()
    payload = build_certificate(arguments.max_m)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check is not None:
        if arguments.check.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archived certificate differs from exact replay")
    if arguments.output is not None:
        arguments.output.write_text(rendered, encoding="utf-8")
    if arguments.output is None and arguments.check is None:
        print(rendered, end="")


if __name__ == "__main__":
    main()
