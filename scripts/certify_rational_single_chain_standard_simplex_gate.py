#!/usr/bin/env python3
"""Exact certificate for the marked P2 single-chain standard-simplex gate.

The surface convention is ``H^2=1`` and ``e_i^2=-1``.  Vectors are stored
as ``(a,c_1,...,c_m)`` for ``a H + sum c_i e_i``.  No reflective-lattice
classification or root enumeration is imported.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence


Vector = tuple[int, ...]


def pair(left: Sequence[int], right: Sequence[int]) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def basis_vector(m: int, index: int) -> Vector:
    values = [0] * (m + 1)
    values[index] = 1
    return tuple(values)


def standard_walls(m: int) -> tuple[Vector, ...]:
    if m < 3:
        raise ValueError("the plane standard frame requires m>=3")
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
    walls.append(basis_vector(m, m))
    return tuple(walls)


def dual_rays(m: int) -> tuple[Vector, ...]:
    rays: list[Vector] = [basis_vector(m, 0)]
    for index in range(1, m + 1):
        ray = [0] * (m + 1)
        ray[0] = min(index, 3)
        for exceptional in range(1, index + 1):
            ray[exceptional] = -1
        rays.append(tuple(ray))
    return tuple(rays)


def canonical_class(m: int) -> Vector:
    return (-3, *((1,) * m))


def arithmetic_genus(m: int, divisor: Vector) -> int:
    numerator = pair(divisor, divisor) + pair(canonical_class(m), divisor)
    if numerator % 2:
        raise ArithmeticError("adjunction numerator is odd")
    return 1 + numerator // 2


def visible_record(m: int) -> dict[str, object]:
    walls = standard_walls(m)
    rays = dual_rays(m)
    kronecker = [
        [pair(ray, wall) for wall in walls]
        for ray in rays
    ]
    expected = [
        [int(row == column) for column in range(m + 1)]
        for row in range(m + 1)
    ]
    if kronecker != expected:
        raise AssertionError("closed-form vectors are not the dual basis")
    squares = [pair(ray, ray) for ray in rays]
    formula = [1] + [min(index, 3) ** 2 - index for index in range(1, m + 1)]
    if squares != formula:
        raise AssertionError("inverse-diagonal formula changed")
    return {
        "m": m,
        "dual_squares": squares,
        "all_dual_squares_nonnegative": all(value >= 0 for value in squares),
        "first_negative_index": next(
            (index for index, value in enumerate(squares) if value < 0), None
        ),
        "kronecker_pairing": True,
    }


def tail_ray(m: int, index: int) -> Vector:
    if not 3 <= index <= m:
        raise ValueError("tail-ray index must satisfy 3<=index<=m")
    return dual_rays(m)[index]


def doubled_stage_ten_wall(m: int) -> Vector:
    if m < 10:
        raise ValueError("the doubled stage-ten wall requires m>=10")
    return tuple(2 * value for value in tail_ray(m, 10))


def augmented_extreme_rays(m: int) -> tuple[Vector, ...]:
    """Extreme rays after adjoining the half-root v_10 as one extra wall.

    In coordinates dual to the standard walls, the extra inequality has
    coefficients p_i=(v_i,v_10).  They are positive for i<9, zero at 9,
    and -1 from 10 onward.  The extreme rays are therefore the surviving
    v_i and the two-coordinate rays v_i+p_i v_j.
    """

    if m < 10:
        raise ValueError("the augmented cone requires m>=10")
    rays = dual_rays(m)
    root = tail_ray(m, 10)
    coefficients = [pair(ray, root) for ray in rays]
    if coefficients[:10] != [3, 2, 4, 6, 5, 4, 3, 2, 1, 0]:
        raise AssertionError("positive stage-ten coefficients changed")
    if any(value != -1 for value in coefficients[10:]):
        raise AssertionError("stage-ten inequality coefficients changed")
    positive = [index for index, value in enumerate(coefficients) if value > 0]
    nonnegative = [
        rays[index] for index, value in enumerate(coefficients) if value >= 0
    ]
    mixed: list[Vector] = []
    for left in positive:
        p_left = coefficients[left]
        for right in range(10, m + 1):
            mixed.append(
                tuple(
                    rays[left][coordinate] + p_left * rays[right][coordinate]
                    for coordinate in range(m + 1)
                )
            )
    return tuple((*nonnegative, *mixed))


def f2_pair(left: Sequence[int], right: Sequence[int]) -> int:
    """Intersection on ``NS(F_2) + <-1>^m`` in the basis ``(s,f,e_i)``."""

    return (
        -2 * left[0] * right[0]
        + left[0] * right[1]
        + left[1] * right[0]
        - sum(a * b for a, b in zip(left[2:], right[2:]))
    )


def f2_basis_vector(m: int, index: int) -> Vector:
    values = [0] * (m + 2)
    values[index] = 1
    return tuple(values)


def f2_standard_walls(m: int) -> tuple[Vector, ...]:
    if m < 2:
        raise ValueError("the F2 affine-fibre frame requires m>=2")
    walls: list[Vector] = [f2_basis_vector(m, 0)]
    delta_0 = [0] * (m + 2)
    delta_0[1] = 1
    delta_0[2:4] = (-1, -1)
    walls.append(tuple(delta_0))
    for index in range(1, m):
        wall = [0] * (m + 2)
        wall[index + 1] = 1
        wall[index + 2] = -1
        walls.append(tuple(wall))
    walls.append(f2_basis_vector(m, m + 1))
    return tuple(walls)


def f2_dual_rays(m: int) -> tuple[Vector, ...]:
    rays: list[Vector] = [f2_basis_vector(m, 1)]
    section_dual = [0] * (m + 2)
    section_dual[0:2] = (1, 2)
    rays.append(tuple(section_dual))
    for index in range(1, m + 1):
        coefficient = min(index, 2)
        ray = [0] * (m + 2)
        ray[0:2] = (coefficient, 2 * coefficient)
        for exceptional in range(1, index + 1):
            ray[exceptional + 1] = -1
        rays.append(tuple(ray))
    return tuple(rays)


def f2_canonical_class(m: int) -> Vector:
    return (-2, -4, *((1,) * m))


def f2_arithmetic_genus(m: int, divisor: Vector) -> int:
    numerator = f2_pair(divisor, divisor) + f2_pair(
        f2_canonical_class(m), divisor
    )
    if numerator % 2:
        raise ArithmeticError("F2 adjunction numerator is odd")
    return 1 + numerator // 2


def f2_visible_record(m: int) -> dict[str, object]:
    walls = f2_standard_walls(m)
    rays = f2_dual_rays(m)
    kronecker = [
        [f2_pair(ray, wall) for wall in walls]
        for ray in rays
    ]
    expected = [
        [int(row == column) for column in range(m + 2)]
        for row in range(m + 2)
    ]
    if kronecker != expected:
        raise AssertionError("closed-form F2 vectors are not the dual basis")
    squares = [f2_pair(ray, ray) for ray in rays]
    formula = [0, 2] + [
        2 * min(index, 2) ** 2 - index for index in range(1, m + 1)
    ]
    if squares != formula:
        raise AssertionError("F2 inverse-diagonal formula changed")
    return {
        "m": m,
        "dual_squares": squares,
        "all_dual_squares_nonnegative": all(value >= 0 for value in squares),
        "first_negative_wall_position": next(
            (index for index, value in enumerate(squares) if value < 0), None
        ),
        "kronecker_pairing": True,
    }


def f2_tail_ray(m: int, index: int) -> Vector:
    if not 2 <= index <= m:
        raise ValueError("F2 tail-ray index must satisfy 2<=index<=m")
    return f2_dual_rays(m)[index + 1]


def f2_doubled_stage_nine_wall(m: int) -> Vector:
    if m < 9:
        raise ValueError("the doubled F2 stage-nine wall requires m>=9")
    return tuple(2 * value for value in f2_tail_ray(m, 9))


def f2_augmented_extreme_rays(m: int) -> tuple[Vector, ...]:
    """Extreme rays after adjoining the half-root u_9 to the F2 frame."""

    if m < 9:
        raise ValueError("the augmented F2 cone requires m>=9")
    rays = f2_dual_rays(m)
    root = f2_tail_ray(m, 9)
    coefficients = [f2_pair(ray, root) for ray in rays]
    if coefficients[:10] != [2, 4, 3, 6, 5, 4, 3, 2, 1, 0]:
        raise AssertionError("positive F2 stage-nine coefficients changed")
    if any(value != -1 for value in coefficients[10:]):
        raise AssertionError("F2 stage-nine inequality coefficients changed")
    positive = [index for index, value in enumerate(coefficients) if value > 0]
    nonnegative = [
        rays[index] for index, value in enumerate(coefficients) if value >= 0
    ]
    mixed: list[Vector] = []
    for left in positive:
        p_left = coefficients[left]
        for right in range(10, m + 2):
            mixed.append(
                tuple(
                    rays[left][coordinate] + p_left * rays[right][coordinate]
                    for coordinate in range(m + 2)
                )
            )
    return tuple((*nonnegative, *mixed))


def build_f2_certificate(max_m: int) -> dict[str, object]:
    if max_m < 11:
        raise ValueError("max_m must expose the first F2 mixed escape at m=11")
    visible = [f2_visible_record(m) for m in range(2, max_m + 1)]
    finite_visible = [
        row["m"] for row in visible if row["all_dual_squares_nonnegative"]
    ]
    if finite_visible != list(range(2, 9)):
        raise AssertionError("F2 visible finite-volume threshold changed")

    later_rays: list[dict[str, int]] = []
    for m in range(9, max_m + 1):
        canonical = f2_canonical_class(m)
        for index in range(9, m + 1):
            ray = f2_tail_ray(m, index)
            later_rays.append(
                {
                    "m": m,
                    "index": index,
                    "square": f2_pair(ray, ray),
                    "canonical_pairing": f2_pair(canonical, ray),
                    "arithmetic_genus": f2_arithmetic_genus(m, ray),
                }
            )
    if any(row["arithmetic_genus"] != 1 for row in later_rays):
        raise AssertionError("an F2 tail ray stopped having genus one")

    doubled_checks: list[dict[str, object]] = []
    for m in range(9, max_m + 1):
        wall = f2_doubled_stage_nine_wall(m)
        half_root = tuple(value // 2 for value in wall)
        pairings = [
            f2_pair(wall, f2_tail_ray(m, index))
            for index in range(9, m + 1)
        ]
        if f2_pair(wall, wall) != -4 or f2_pair(half_root, half_root) != -1:
            raise AssertionError("doubled F2 wall norms changed")
        if (
            f2_pair(f2_canonical_class(m), wall) != 2
            or f2_arithmetic_genus(m, wall) != 0
        ):
            raise AssertionError("doubled F2 wall adjunction changed")
        if pairings != [-2] * (m - 8):
            raise AssertionError("F2 doubled wall stopped cutting old tail rays")
        doubled_checks.append(
            {
                "m": m,
                "wall_square": f2_pair(wall, wall),
                "half_root_square": f2_pair(half_root, half_root),
                "canonical_pairing": f2_pair(f2_canonical_class(m), wall),
                "arithmetic_genus": f2_arithmetic_genus(m, wall),
                "pairing_with_u_9_through_u_m": pairings,
            }
        )

    augmented: list[dict[str, object]] = []
    for m in range(9, max_m + 1):
        rays = f2_augmented_extreme_rays(m)
        squares = [f2_pair(ray, ray) for ray in rays]
        augmented.append(
            {
                "m": m,
                "extreme_ray_count": len(rays),
                "minimum_square": min(squares),
                "all_squares_nonnegative": all(value >= 0 for value in squares),
            }
        )
    finite_augmented = [
        row["m"] for row in augmented if row["all_squares_nonnegative"]
    ]
    if finite_augmented != [9, 10]:
        raise AssertionError("one-wall augmented F2 threshold changed")
    f_ray = f2_dual_rays(11)[0]
    u_11 = f2_tail_ray(11, 11)
    first_mixed = tuple(a + 2 * b for a, b in zip(f_ray, u_11))
    if f2_pair(first_mixed, first_mixed) != -4:
        raise AssertionError("first F2 mixed escape changed")

    return {
        "base": "F2",
        "marking": (
            "negative section plus a fibre-following free point at the "
            "strict-fibre/E1 intersection, followed by a free exceptional "
            "chain"
        ),
        "centres_off_negative_section": True,
        "negative_section_and_affine_roots_physical_is_conditional": True,
        "visible_simplex": {
            "inverse_diagonal_formula": (
                "0,2; 2*min(i,2)^2-i for delta_i duals, 1<=i<=m"
            ),
            "finite_volume_exactly_m_at_most": 8,
            "rows": visible,
        },
        "genus_one_escape_family": later_rays,
        "first_escape": "u_9=-K_(S_9)",
        "doubled_stage_nine_numerical_wall": doubled_checks,
        "one_wall_augmented_cone": {
            "finite_volume_exactly_m": [9, 10],
            "rows": augmented,
            "first_mixed_escape": {
                "m": 11,
                "vector_description": "f+2u_11",
                "square": f2_pair(first_mixed, first_mixed),
            },
        },
        "geometric_effectivity_of_doubled_numerical_wall_claimed": False,
        "uniform_chain_length_bound_claimed": False,
    }


def build_certificate(max_m: int = 32) -> dict[str, object]:
    if max_m < 12:
        raise ValueError("max_m must expose the first mixed escape at m=12")
    visible = [visible_record(m) for m in range(3, max_m + 1)]
    finite_visible = [
        row["m"] for row in visible if row["all_dual_squares_nonnegative"]
    ]
    if finite_visible != list(range(3, 10)):
        raise AssertionError("visible finite-volume threshold changed")

    later_rays: list[dict[str, int]] = []
    for m in range(10, max_m + 1):
        canonical = canonical_class(m)
        for index in range(10, m + 1):
            ray = tail_ray(m, index)
            later_rays.append(
                {
                    "m": m,
                    "index": index,
                    "square": pair(ray, ray),
                    "canonical_pairing": pair(canonical, ray),
                    "arithmetic_genus": arithmetic_genus(m, ray),
                }
            )
    if any(row["arithmetic_genus"] != 1 for row in later_rays):
        raise AssertionError("a standard tail ray stopped having genus one")

    doubled_checks: list[dict[str, object]] = []
    for m in range(10, max_m + 1):
        wall = doubled_stage_ten_wall(m)
        half_root = tuple(value // 2 for value in wall)
        pairings = [pair(wall, tail_ray(m, index)) for index in range(10, m + 1)]
        if pair(wall, wall) != -4 or pair(half_root, half_root) != -1:
            raise AssertionError("doubled wall norms changed")
        if pair(canonical_class(m), wall) != 2 or arithmetic_genus(m, wall) != 0:
            raise AssertionError("doubled wall adjunction changed")
        if pairings != [-2] * (m - 9):
            raise AssertionError("doubled wall no longer cuts all original tail rays")
        doubled_checks.append(
            {
                "m": m,
                "wall_square": pair(wall, wall),
                "half_root_square": pair(half_root, half_root),
                "canonical_pairing": pair(canonical_class(m), wall),
                "arithmetic_genus": arithmetic_genus(m, wall),
                "pairing_with_v_10_through_v_m": pairings,
            }
        )

    augmented: list[dict[str, object]] = []
    for m in range(10, max_m + 1):
        rays = augmented_extreme_rays(m)
        squares = [pair(ray, ray) for ray in rays]
        augmented.append(
            {
                "m": m,
                "extreme_ray_count": len(rays),
                "minimum_square": min(squares),
                "all_squares_nonnegative": all(value >= 0 for value in squares),
            }
        )
    finite_augmented = [
        row["m"] for row in augmented if row["all_squares_nonnegative"]
    ]
    if finite_augmented != [10, 11]:
        raise AssertionError("one-wall augmented threshold changed")
    first_mixed = tuple(
        dual_rays(12)[1][coordinate] + 2 * dual_rays(12)[12][coordinate]
        for coordinate in range(13)
    )
    if pair(first_mixed, first_mixed) != -4:
        raise AssertionError("first mixed escape changed")

    return {
        "certificate": "rational_single_chain_standard_simplex_gate",
        "status": "PASS_EXACT_SCOPED_NO_GO",
        "scope": {
            "base": "P2",
            "marking": "standard free-chain root frame",
            "alpha_0_not_universal": True,
            "geometric_effectivity_of_doubled_numerical_wall_claimed": False,
            "uniform_chain_length_bound_claimed": False,
            "classification_or_census_used": False,
        },
        "visible_simplex": {
            "inverse_diagonal_formula": "1; min(i,3)^2-i for 1<=i<=m",
            "finite_volume_exactly_m_at_most": 9,
            "rows": visible,
        },
        "genus_one_escape_family": later_rays,
        "full_NS_hidden_wall_identity": {
            "formula": "C.v_j=C^2+2+sum_{i>j} mu_i",
            "smooth_rational_plane_curve_multiplicities_nonnegative": True,
            "only_square_minus_four_can_pair_negatively": True,
            "integral_reflection_for_square_minus_four_forces_C_equals_2r": True,
            "tail_multiplicities_even": True,
            "negative_pairing_forces_tail_sum_zero": True,
        },
        "doubled_stage_ten_numerical_wall": doubled_checks,
        "one_wall_augmented_cone": {
            "finite_volume_exactly_m": [10, 11],
            "rows": augmented,
            "first_mixed_escape": {
                "m": 12,
                "vector_description": "v_1+2v_12",
                "square": pair(first_mixed, first_mixed),
            },
        },
        "F2_minimal_base_companion": build_f2_certificate(max_m),
        "conclusion": (
            "on both P2 and F2, the standard inverse diagonal cannot bound "
            "chain depth: one admissible doubled numerical wall cuts every "
            "original tail ray; new mixed extreme rays, not the old diagonal "
            "family, carry the remaining global debt"
        ),
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
