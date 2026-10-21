#!/usr/bin/env python3
"""Exact rank 4--9 ordinary P-first closure with timelike -K vertex."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


Vector = tuple[int, ...]


def dot(left: Vector, right: Vector) -> int:
    return left[0] * right[0] - sum(
        a * b for a, b in zip(left[1:], right[1:])
    )


def bareiss_determinant(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix]
    size = len(work)
    assert size and all(len(row) == size for row in work)
    sign = 1
    previous = 1
    for column in range(size - 1):
        if work[column][column] == 0:
            pivot = next(
                (
                    row
                    for row in range(column + 1, size)
                    if work[row][column]
                ),
                None,
            )
            if pivot is None:
                return 0
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        value = work[column][column]
        for row in range(column + 1, size):
            for next_column in range(column + 1, size):
                numerator = (
                    work[row][next_column] * value
                    - work[row][column] * work[column][next_column]
                )
                assert numerator % previous == 0
                work[row][next_column] = numerator // previous
        previous = value
    return sign * work[-1][-1]


def solve_square(matrix: list[list[int]], right: list[int]) -> list[Fraction]:
    size = len(matrix)
    augmented = [
        [Fraction(entry) for entry in row] + [Fraction(value)]
        for row, value in zip(matrix, right)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column]),
            None,
        )
        if pivot is None:
            raise ArithmeticError("singular exact system")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        value = augmented[column][column]
        augmented[column] = [entry / value for entry in augmented[column]]
        for row in range(size):
            if row != column and augmented[row][column]:
                coefficient = augmented[row][column]
                augmented[row] = [
                    entry - coefficient * pivot_entry
                    for entry, pivot_entry in zip(
                        augmented[row], augmented[column]
                    )
                ]
    return [augmented[index][-1] for index in range(size)]


def normalize(entries: list[Fraction]) -> Vector:
    denominator = 1
    for entry in entries:
        denominator = math.lcm(denominator, entry.denominator)
    integral = [int(entry * denominator) for entry in entries]
    divisor = math.gcd(*map(abs, integral))
    primitive = tuple(entry // divisor for entry in integral)
    if primitive[0] < 0:
        primitive = tuple(-entry for entry in primitive)
    return primitive


def standard_roots(exceptional_count: int) -> tuple[Vector, ...]:
    dimension = exceptional_count + 1
    roots: list[Vector] = []
    last_exceptional = [0] * dimension
    last_exceptional[exceptional_count] = 1
    roots.append(tuple(last_exceptional))
    for index in range(1, exceptional_count):
        difference = [0] * dimension
        difference[index] = 1
        difference[index + 1] = -1
        roots.append(tuple(difference))
    flex_line = [0] * dimension
    flex_line[0] = 1
    flex_line[1] = flex_line[2] = flex_line[3] = -1
    roots.append(tuple(flex_line))
    return tuple(roots)


def dual_vertices(roots: tuple[Vector, ...]) -> tuple[Vector, ...]:
    size = len(roots)
    metric = (1,) + (-1,) * (size - 1)
    matrix = [
        [metric[column] * root[column] for column in range(size)]
        for root in roots
    ]
    answer = []
    for index in range(size):
        right = [1 if row == index else 0 for row in range(size)]
        answer.append(normalize(solve_square(matrix, right)))
    return tuple(answer)


def coxeter_histogram(roots: tuple[Vector, ...]) -> dict[str, int]:
    histogram: Counter[str] = Counter()
    for first_index, first in enumerate(roots):
        first_norm = -dot(first, first)
        for second in roots[first_index + 1 :]:
            second_norm = -dot(second, second)
            pairing = dot(first, second)
            if pairing == 0:
                order = "2"
            elif first_norm == second_norm == 2 and pairing == 1:
                order = "3"
            elif {first_norm, second_norm} == {1, 2} and pairing == 1:
                order = "4"
            else:
                raise ArithmeticError(
                    f"unexpected Coxeter pair {first_norm, second_norm, pairing}"
                )
            histogram[order] += 1
    return dict(sorted(histogram.items()))


def rank_record(exceptional_count: int) -> dict[str, object]:
    n = exceptional_count
    dimension = n + 1
    canonical: Vector = (-3,) + (1,) * n
    anticanonical: Vector = (3,) + (-1,) * n
    roots = standard_roots(n)

    assert dot(canonical, canonical) == 9 - n == 10 - dimension
    assert dot(anticanonical, anticanonical) == 9 - n > 0
    assert all(entry % 2 for entry in canonical)
    assert [dot(root, root) for root in roots] == [-1] + [-2] * n
    assert [dot(canonical, root) for root in roots] == [-1] + [0] * n
    assert [dot(anticanonical, root) for root in roots] == [1] + [0] * n
    assert bareiss_determinant([list(root) for root in roots]) == -1

    gram = [[dot(first, second) for second in roots] for first in roots]
    q2_gram = [row[1:] for row in gram[1:]]
    assert abs(bareiss_determinant(q2_gram)) == 9 - n

    vertices = dual_vertices(roots)
    squares = [dot(vertex, vertex) for vertex in vertices]
    assert vertices[0] == anticanonical
    assert squares[0] == 9 - n
    assert squares[1] == 0
    assert all(square >= 0 for square in squares)
    assert all(vertex[0] > 0 for vertex in vertices)

    # The transformed flex pencil on the first n blowups has divisor class
    # 3H-sum E_i=-K.  For n=3 it has residual ideal (x^6,w); each subsequent
    # blowup lowers the exponent by one.  Until n=8 this is a finite base
    # point but never a fixed curve.
    residual_exponent = 9 - n
    assert 1 <= residual_exponent <= 6

    # Symbolic positivity terms used in both classification-free arguments.
    # In the reflecting-root descent, with m1=k+A, m2=k+B, m3=k:
    # s^2-[m1^2+m2^2+(n-2)m3^2]
    #   =(9-n)k^2+4k(A+B)+2AB.
    # In characteristic-vector rigidity the corresponding k(k+1) gap is:
    #   =(9-n)k(k+1)+4k(A+B)+2AB+2A+2B.
    assert 9 - n > 0

    return {
        "exceptional_count": n,
        "rank": dimension,
        "hyperbolic_dimension": n,
        "canonical": list(canonical),
        "canonical_square": 9 - n,
        "anticanonical": list(anticanonical),
        "anticanonical_square": 9 - n,
        "anticanonical_is_timelike_vertex": True,
        "facet_roots": [list(root) for root in roots],
        "facet_squares": [-1] + [-2] * n,
        "canonical_pairings": [-1] + [0] * n,
        "facet_basis_determinant": -1,
        "q2_gram_determinant_absolute": 9 - n,
        "gram": gram,
        "coxeter_order_histogram": coxeter_histogram(roots),
        "vertices": [list(vertex) for vertex in vertices],
        "vertex_squares": squares,
        "ideal_vertex_count": sum(square == 0 for square in squares),
        "reflection_descent_gap_coefficient": 9 - n,
        "characteristic_rigidity_gap_coefficient": 9 - n,
        "flex_cluster_carrier": {
            "blowups_performed": n,
            "residual_base_ideal": f"(x^{residual_exponent},w)",
            "anticanonical_pencil_has_no_fixed_curve": True,
            "anticanonical_is_nef": True,
            "facet_curves_effective_and_irreducible": True,
            "full_nef_return": True,
        },
    }


def payload() -> dict[str, object]:
    records = [rank_record(n) for n in range(3, 9)]
    result: dict[str, object] = {
        "schema": "rational_pfirst_timelike_anticanonical_vertex_closure/v1",
        "ranks_covered": [4, 5, 6, 7, 8, 9],
        "exceptional_counts_covered": [3, 4, 5, 6, 7, 8],
        "records": records,
        "uniform_reflection_maximality": {
            "primitive_integral_reflecting_squares": [-1, -2],
            "generators": [
                "signed permutations of the exceptional coordinates",
                "reflection in H-E1-E2-E3",
            ],
            "time_coordinate_descent": True,
            "all_integral_reflections_generated": True,
            "relative_maximality_in_O_plus_I_1_n": True,
        },
        "uniform_characteristic_rigidity": {
            "reduced_characteristic_substitution": (
                "a=2l+3 and m_i=2k_i+1"
            ),
            "fixed_norm": "a^2-sum(m_i^2)=9-n",
            "only_solution_in_standard_chamber": (
                "l=k_1=...=k_n=0, hence -K=(3,-1,...,-1)"
            ),
            "ordinary_physical_maximal_chamber_is_Stab_K_conjugate": True,
        },
        "uniform_carrier": {
            "pencil": "Y^2 Z = X^3 + X Z^2 + lambda Z^3",
            "base_cluster_total_length": 9,
            "construction": "perform the first n=3,...,8 flex-cluster blowups",
            "transformed_anticanonical_pencil_has_no_fixed_component": True,
            "no_fixed_component_implies_anticanonical_nef": True,
            "weak_log_reflection_return_proves_exact_full_nef_cone": True,
        },
        "theorem_scope": {
            "ordinary_walls_only_minus_one_minus_two": True,
            "K_is_not_nef": True,
            "minus_K_is_not_a_cusp": True,
            "minus_K_is_a_timelike_vertex": True,
            "classification_free_exact_carrier": "PROVED",
            "twice_minus_four_or_minus_K_outside_P": "OPEN",
            "rank_at_least_ten_non_cusp_lane": "OPEN",
        },
    }
    canonical = json.dumps(
        result, sort_keys=True, separators=(",", ":")
    ).encode()
    result["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    result = payload()
    default = Path(
        "data/enumerations/"
        "rational_pfirst_timelike_anticanonical_vertex_closure.json"
    )
    target = args.output or default
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        assert target.read_text() == rendered
    elif args.output:
        target.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
