#!/usr/bin/env python3
"""Symbolic certificate for the first-mixed-ray standard-frame bound.

This checker uses only the two minimal-base intersection forms, adjunction,
and the closed-form mixed rays.  It imports no root enumerator or reflective
lattice classification.
"""

from __future__ import annotations

import argparse
from itertools import combinations_with_replacement
import json
from pathlib import Path
from typing import Sequence


Vector = tuple[int, ...]


def nonincreasing_tuples(length: int, largest_entry: int) -> tuple[Vector, ...]:
    """Return the finite monotone patterns used in the low-pairing proof."""
    return tuple(
        tuple(reversed(pattern))
        for pattern in combinations_with_replacement(
            range(largest_entry + 1), length
        )
    )


def plane_pair(left: Sequence[int], right: Sequence[int]) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def plane_standard_walls(m: int) -> tuple[Vector, ...]:
    if m < 12:
        raise ValueError("the first plane mixed ray requires m>=12")
    walls: list[Vector] = []
    alpha_0 = [0] * (m + 1)
    alpha_0[0] = 1
    alpha_0[1:4] = (-1, -1, -1)
    walls.append(tuple(alpha_0))
    for index in range(1, m):
        wall = [0] * (m + 1)
        wall[index] = 1
        wall[index + 1] = -1
        walls.append(tuple(wall))
    endpoint = [0] * (m + 1)
    endpoint[m] = 1
    walls.append(tuple(endpoint))
    return tuple(walls)


def plane_first_mixed_ray(m: int) -> Vector:
    if m < 12:
        raise ValueError("the first plane mixed ray requires m>=12")
    ray = [0] * (m + 1)
    ray[0] = 7
    ray[1] = -3
    for index in range(2, 13):
        ray[index] = -2
    return tuple(ray)


def plane_canonical(m: int) -> Vector:
    return (-3, *((1,) * m))


def f2_pair(left: Sequence[int], right: Sequence[int]) -> int:
    return (
        -2 * left[0] * right[0]
        + left[0] * right[1]
        + left[1] * right[0]
        - sum(a * b for a, b in zip(left[2:], right[2:]))
    )


def f2_standard_walls(m: int) -> tuple[Vector, ...]:
    if m < 11:
        raise ValueError("the first F2 mixed ray requires m>=11")
    section = [0] * (m + 2)
    section[0] = 1
    walls: list[Vector] = [tuple(section)]
    delta_0 = [0] * (m + 2)
    delta_0[1] = 1
    delta_0[2:4] = (-1, -1)
    walls.append(tuple(delta_0))
    for index in range(1, m):
        wall = [0] * (m + 2)
        wall[index + 1] = 1
        wall[index + 2] = -1
        walls.append(tuple(wall))
    endpoint = [0] * (m + 2)
    endpoint[m + 1] = 1
    walls.append(tuple(endpoint))
    return tuple(walls)


def f2_first_mixed_ray(m: int) -> Vector:
    if m < 11:
        raise ValueError("the first F2 mixed ray requires m>=11")
    ray = [0] * (m + 2)
    ray[0:2] = (4, 9)
    for index in range(2, 13):
        ray[index] = -2
    return tuple(ray)


def f2_canonical(m: int) -> Vector:
    return (-2, -4, *((1,) * m))


def standard_ray_records(max_m: int) -> dict[str, object]:
    plane_rows: list[dict[str, object]] = []
    for m in range(12, max_m + 1):
        ray = plane_first_mixed_ray(m)
        pairings = [plane_pair(ray, wall) for wall in plane_standard_walls(m)]
        expected = [0] * (m + 1)
        expected[1] = 1
        expected[12] = 2
        if pairings != expected:
            raise AssertionError("plane mixed ray stopped lying in the wall cone")
        if plane_pair(ray, ray) != -4 or plane_pair(plane_canonical(m), ray) != 4:
            raise AssertionError("plane mixed-ray norm or canonical pairing changed")
        plane_rows.append(
            {
                "m": m,
                "square": plane_pair(ray, ray),
                "canonical_pairing": plane_pair(plane_canonical(m), ray),
                "nonzero_standard_wall_pairings": {"alpha_1": 1, "alpha_12": 2},
            }
        )

    f2_rows: list[dict[str, object]] = []
    for m in range(11, max_m + 1):
        ray = f2_first_mixed_ray(m)
        pairings = [f2_pair(ray, wall) for wall in f2_standard_walls(m)]
        expected = [0] * (m + 2)
        expected[0] = 1
        expected[12] = 2
        if pairings != expected:
            raise AssertionError("F2 mixed ray stopped lying in the wall cone")
        if f2_pair(ray, ray) != -4 or f2_pair(f2_canonical(m), ray) != 4:
            raise AssertionError("F2 mixed-ray norm or canonical pairing changed")
        f2_rows.append(
            {
                "m": m,
                "square": f2_pair(ray, ray),
                "canonical_pairing": f2_pair(f2_canonical(m), ray),
                "nonzero_standard_wall_pairings": {"section_s": 1, "delta_11": 2},
            }
        )
    return {"plane": plane_rows, "F2": f2_rows}


def plane_symbolic_proof() -> dict[str, object]:
    # q=-C^2.  Adjunction gives total multiplicity 3d+q-2.
    species = {
        "minus_one": {
            "q": 1,
            "pairing_lower_bound_from_d_ge_mu1": 2,
        },
        "minus_two": {
            "q": 2,
            "pairing_lower_bound_from_d_ge_mu1": 0,
        },
        "minus_four": {
            "q": 4,
            "pairing_reparameterization": [
                "R=d-mu_1>=mu_2+mu_3",
                "T=sum_{i>12}mu_i",
                "C.w_P2=R-4+2T",
            ],
            "T_at_least_2": "pairing nonnegative",
            "T_equals_1": (
                "proximity gives mu_2,mu_3>=1, hence R>=2 and pairing>=0"
            ),
            "T_equals_0_R_at_most_3_cases": {
                "R=0_or_1": "total multiplicity sum is immediately too large",
                "R=2": (
                    "the only large-sum pattern forces d=3 and ten remaining "
                    "ones; remaining square sum 10 contradicts required 12"
                ),
                "R=3": (
                    "eleven remaining positions have total at most 12, while "
                    "the adjunction total is 2d+5 and the small d cases also fail"
                ),
            },
            "conclusion": "R>=4, so C.w_P2>=0",
        },
    }
    return {
        "adjunction_total_multiplicity": "sum mu_i=3d+q-2",
        "mixed_pairing_identity": (
            "C.w_P2=d-mu_1-2q+4+2*sum_{i>12}mu_i"
        ),
        "proximity": "mu_1>=mu_2>=...>=0",
        "line_nonintersection": "d>=mu_1+mu_2+mu_3",
        "species": species,
        "all_physical_pairings_nonnegative": True,
    }


def certify_plane_minus_four_low_R() -> dict[str, object]:
    """Exclude the T=0, R=0,1,2,3 branches by exact integer arithmetic."""
    rows: list[dict[str, int]] = []
    for R in range(4):
        patterns = 0
        sum_and_proximity_candidates = 0
        square_solutions = 0
        for tail in nonincreasing_tuples(11, R):
            if tail[0] + tail[1] > R:
                continue
            patterns += 1
            numerator = sum(tail) - R - 2
            if numerator < 0 or numerator % 2:
                continue
            d = numerator // 2
            mu_1 = d - R
            if mu_1 < tail[0]:
                continue
            sum_and_proximity_candidates += 1
            required_square_sum = 2 * d * R - R * R + 4
            if sum(value * value for value in tail) == required_square_sum:
                square_solutions += 1
        if square_solutions:
            raise AssertionError(f"plane low-R branch R={R} has a solution")
        rows.append(
            {
                "R": R,
                "monotone_line_patterns": patterns,
                "sum_and_proximity_candidates": sum_and_proximity_candidates,
                "square_solutions": square_solutions,
            }
        )
    return {
        "branch": "q=4, T_12=0, C.w_P2<0",
        "tested_R_values": [0, 1, 2, 3],
        "rows": rows,
        "solutions": 0,
    }


def f2_symbolic_proof() -> dict[str, object]:
    # q=-C^2.  Adjunction gives total multiplicity 2b+q-2.
    species = {
        "minus_one": {
            "q": 1,
            "pairing_lower_bound_from_a_ge_0": 2,
        },
        "minus_two": {
            "q": 2,
            "pairing_lower_bound_from_a_ge_0": 0,
        },
        "minus_four": {
            "q": 4,
            "pairing_reparameterization": [
                "T=sum_{i>11}mu_i",
                "C.w_F2=a-4+2T",
            ],
            "T_at_least_2": "pairing nonnegative",
            "T_equals_1": (
                "proximity and fibre nonintersection give a>=mu_1+mu_2>=2"
            ),
            "T_equals_0_a_at_most_3_cases": {
                "a=0_or_1": "total multiplicity sum contradicts b>=2a",
                "a=2": (
                    "only b=4 and ten ones survive the sum bound; square sum "
                    "10 contradicts required 12"
                ),
                "a=3": "total multiplicity is at most 12 but at least 14",
            },
            "conclusion": "a>=4, so C.w_F2>=0",
        },
    }
    return {
        "adjunction_total_multiplicity": "sum mu_i=2b+q-2",
        "mixed_pairing_identity": "C.w_F2=a-2q+4+2*sum_{i>11}mu_i",
        "proximity": (
            "mu_1>=mu_2>=...>=0 along the fibre-following free chain"
        ),
        "fibre_nonintersection": "a>=mu_1+mu_2",
        "section_nonintersection": "b>=2a",
        "species": species,
        "all_physical_pairings_nonnegative": True,
    }


def certify_f2_minus_four_low_a() -> dict[str, object]:
    """Exclude the T=0, a=0,1,2,3 branches by exact integer arithmetic."""
    rows: list[dict[str, int]] = []
    for a in range(4):
        patterns = 0
        sum_section_candidates = 0
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
            sum_section_candidates += 1
            required_square_sum = -2 * a * a + 2 * a * b + 4
            if (
                sum(value * value for value in multiplicities)
                == required_square_sum
            ):
                square_solutions += 1
        if square_solutions:
            raise AssertionError(f"F2 low-a branch a={a} has a solution")
        rows.append(
            {
                "a": a,
                "monotone_fibre_patterns": patterns,
                "sum_and_section_candidates": sum_section_candidates,
                "square_solutions": square_solutions,
            }
        )
    return {
        "branch": "q=4, T_11=0, C.w_F2<0",
        "tested_a_values": [0, 1, 2, 3],
        "rows": rows,
        "solutions": 0,
    }


def build_certificate(max_m: int = 32) -> dict[str, object]:
    if max_m < 12:
        raise ValueError("max_m must be at least 12")
    return {
        "certificate": "rational_single_chain_first_mixed_ray_bound",
        "status": "PASS_EXACT_CONDITIONAL_DEPTH_BOUND",
        "scope": {
            "host": "direct or adapted; no parity/doubling is used",
            "complete_physical_alphabet": "smooth rational (-1),(-2),(-4) facets",
            "plane_frame": (
                "alpha_0 and the free exceptional chain are actual unblown "
                "irreducible walls"
            ),
            "F2_frame": (
                "s and delta_0 are actual unblown irreducible walls; the "
                "second centre is the fibre-following free infinitely-near "
                "point (only E1 exceptional), and later centres are free off s"
            ),
            "adapted_half_boundary_host_covered_on_the_exact_frames": True,
            "arbitrary_single_chain_covered": False,
            "reflection_classification_or_census_used": False,
        },
        "standard_mixed_rays": standard_ray_records(max_m),
        "plane_symbolic_proof": plane_symbolic_proof(),
        "plane_low_pairing_certificate": certify_plane_minus_four_low_R(),
        "F2_symbolic_proof": f2_symbolic_proof(),
        "F2_low_pairing_certificate": certify_f2_minus_four_low_a(),
        "finite_volume_conclusion": {
            "plane_chain_length_at_most": 11,
            "F2_chain_length_at_most": 10,
            "reason": (
                "the spacelike first mixed ray lies in every physical wall "
                "halfspace; the homogeneous cone of a finite-volume "
                "finite-sided hyperbolic polyhedron cannot contain it"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-m", type=int, default=32)
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
