#!/usr/bin/env python3
"""Exact hostile replay for the rational H3 star trichotomy.

The replay uses only integer and rational arithmetic.  It certifies the local
intersection, discrepancy, reflection, fibration-core, and wrong-sign
calculations.  It intentionally does not assert semiampleness or an absolute
tail/rank bound.
"""

from __future__ import annotations

import json
from fractions import Fraction
from itertools import product


def mat_vec(matrix: list[list[int]], vector: list[int | Fraction]) -> list[Fraction]:
    return [
        sum(Fraction(entry) * Fraction(value) for entry, value in zip(row, vector))
        for row in matrix
    ]


def dot(left: list[int | Fraction], right: list[int | Fraction]) -> Fraction:
    return sum(Fraction(x) * Fraction(y) for x, y in zip(left, right))


def bilinear(
    left: list[int | Fraction],
    matrix: list[list[int]],
    right: list[int | Fraction],
) -> Fraction:
    return dot(left, mat_vec(matrix, right))


def solve(
    matrix: list[list[int]], right: list[int | Fraction]
) -> list[Fraction]:
    size = len(matrix)
    rows = [
        [Fraction(value) for value in row] + [Fraction(right[i])]
        for i, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [value / scale for value in rows[column]]
        for row in range(size):
            if row == column:
                continue
            scale = rows[row][column]
            if scale:
                rows[row] = [
                    value - scale * pivot_value
                    for value, pivot_value in zip(rows[row], rows[column])
                ]
    return [rows[i][-1] for i in range(size)]


def determinant(matrix: list[list[int]]) -> Fraction:
    size = len(matrix)
    rows = [[Fraction(value) for value in row] for row in matrix]
    answer = Fraction(1)
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if rows[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            answer *= -1
        pivot_value = rows[column][column]
        answer *= pivot_value
        for row in range(column + 1, size):
            scale = rows[row][column] / pivot_value
            for entry in range(column + 1, size):
                rows[row][entry] -= scale * rows[column][entry]
    return answer


def physical_star(profile: tuple[int, ...]) -> list[list[int]]:
    size = 1 + len(profile)
    matrix = [[0] * size for _ in range(size)]
    matrix[0][0] = -1
    for index, incidence in enumerate(profile, start=1):
        matrix[index][index] = -4
        matrix[0][index] = matrix[index][0] = incidence
    return matrix


def host_star(profile: tuple[int, ...]) -> list[list[int]]:
    size = 1 + len(profile)
    matrix = [[0] * size for _ in range(size)]
    for index in range(size):
        matrix[index][index] = -2
    for index, incidence in enumerate(profile, start=1):
        matrix[0][index] = matrix[index][0] = incidence
    return matrix


def reflect(
    vector: list[int | Fraction],
    root: list[int | Fraction],
    matrix: list[list[int]],
) -> list[Fraction]:
    root_square = bilinear(root, matrix, root)
    coefficient = -2 * bilinear(vector, matrix, root) / root_square
    return [
        Fraction(value) + coefficient * Fraction(root_value)
        for value, root_value in zip(vector, root)
    ]


def tail_matrix(length: int) -> list[list[int]]:
    matrix = [[0] * length for _ in range(length)]
    for index in range(length):
        matrix[index][index] = -1 if index == length - 1 else -2
        if index:
            matrix[index - 1][index] = matrix[index][index - 1] = 1
    return matrix


def encode(value: Fraction | int) -> int | str:
    value = Fraction(value)
    if value.denominator == 1:
        return value.numerator
    return f"{value.numerator}/{value.denominator}"


def certificate() -> dict[str, object]:
    low_profiles: set[tuple[int, ...]] = set()
    for length in range(1, 7):
        for profile in product(range(5), repeat=length):
            normalized = tuple(sorted((value for value in profile if value), reverse=True))
            if not normalized or normalized != profile[: len(normalized)]:
                continue
            total = sum(normalized)
            square_total = sum(value * value for value in normalized)
            if total >= 3 and square_total <= 4:
                low_profiles.add(normalized)
    assert low_profiles == {(1, 1, 1), (1, 1, 1, 1)}

    spherical_profile = (1, 1, 1)
    spherical_physical = physical_star(spherical_profile)
    spherical_host = host_star(spherical_profile)
    canonical = [-1, 2, 2, 2]
    fundamental = [3, 1, 1, 1]
    discrepancy = solve(spherical_physical, canonical)
    assert discrepancy == [Fraction(-2), Fraction(-1), Fraction(-1), Fraction(-1)]
    intersections = mat_vec(spherical_physical, fundamental)
    assert intersections == [0, -1, -1, -1]
    z_square = bilinear(fundamental, spherical_physical, fundamental)
    k_dot_z = dot(canonical, fundamental)
    arithmetic_genus = 1 + (z_square + k_dot_z) / 2
    assert (z_square, k_dot_z, arithmetic_genus) == (-3, 3, 1)
    assert [
        encode(determinant([row[:size] for row in spherical_physical[:size]]))
        for size in range(1, 5)
    ] == [-1, 3, -8, 16]

    # Coefficientwise minimality of the fundamental anti-nef cycle.
    anti_nef_full_support = []
    for candidate in product(*(range(1, value + 1) for value in fundamental)):
        if all(value <= 0 for value in mat_vec(spherical_physical, list(candidate))):
            anti_nef_full_support.append(candidate)
    assert anti_nef_full_support == [tuple(fundamental)]

    basis = [[int(i == j) for i in range(4)] for j in range(4)]
    reflected = [Fraction(value) for value in basis[1]]
    for root_index in (0, 2, 3, 0):
        reflected = reflect(reflected, basis[root_index], spherical_host)
    highest_root = [2, 1, 1, 1]
    assert reflected == highest_root
    assert bilinear(highest_root, spherical_host, highest_root) == -2

    # Concrete unimodular rational-surface model:
    # NS(P1xP1 blown up at 13 points) = U + <-1>^13.
    # B_i=f-sum_{j=1}^4 e_ij and E=g-e_0.  Leaves are disjoint
    # square-minus-four classes meeting E once, but B_1-B_2 is odd.
    dimension = 15
    u_gram = [[0] * dimension for _ in range(dimension)]
    u_gram[0][1] = u_gram[1][0] = 1
    for index in range(2, dimension):
        u_gram[index][index] = -1
    e_curve = [0] * dimension
    e_curve[1] = 1
    e_curve[2] = -1
    leaves: list[list[int]] = []
    for leaf_index in range(3):
        leaf = [0] * dimension
        leaf[0] = 1
        start = 3 + 4 * leaf_index
        for exceptional_index in range(start, start + 4):
            leaf[exceptional_index] = -1
        leaves.append(leaf)
    assert bilinear(e_curve, u_gram, e_curve) == -1
    assert all(bilinear(leaf, u_gram, leaf) == -4 for leaf in leaves)
    assert all(bilinear(e_curve, u_gram, leaf) == 1 for leaf in leaves)
    assert all(
        bilinear(leaves[i], u_gram, leaves[j]) == 0
        for i in range(3)
        for j in range(i + 1, 3)
    )
    leaf_difference = [left - right for left, right in zip(leaves[0], leaves[1])]
    parity_witness = [0] * dimension
    parity_witness[3] = 1
    odd_pairing = bilinear(parity_witness, u_gram, leaf_difference)
    assert odd_pairing == 1
    assert any(value % 2 for value in leaf_difference)

    affine_profile = (1, 1, 1, 1)
    affine_physical = physical_star(affine_profile)
    affine_host = host_star(affine_profile)
    affine_f = [4, 1, 1, 1, 1]
    affine_r = [2, 1, 1, 1, 1]
    affine_canonical = [-1, 2, 2, 2, 2]
    assert mat_vec(affine_physical, affine_f) == [0] * 5
    assert mat_vec(affine_host, affine_r) == [0] * 5
    assert dot(affine_canonical, affine_f) == 4
    fibration_cases = []
    for multiplicity in range(1, 9):
        for genus in range(2, 9):
            if multiplicity * (2 * genus - 2) == 4:
                rho_bound = 10 + 8 * (genus - 1)
                fibration_cases.append((multiplicity, genus, rho_bound))
    assert fibration_cases == [(1, 3, 26), (2, 2, 18)]

    hyperbolic_rows = []
    for profile in ((2, 1), (3,), (1, 1, 1, 1, 1), (2, 2)):
        total = sum(profile)
        square_total = sum(value * value for value in profile)
        assert total >= 3 and square_total > 4
        matrix = physical_star(profile)
        f_vector = [4, *profile]
        canonical_vector = [-1, *([2] * len(profile))]
        f_square = bilinear(f_vector, matrix, f_vector)
        k_dot_f = dot(canonical_vector, f_vector)
        assert f_square == 4 * (square_total - 4)
        assert k_dot_f == -4 + 2 * total
        bar_e_square = Fraction(square_total, 4) - 1
        ky_dot_bar_e = Fraction(total, 2) - 1
        assert bar_e_square > 0 and ky_dot_bar_e > 0
        coefficient_left = 1 - Fraction(square_total, 4)
        coefficient_right = Fraction(total, 2) - 1
        assert coefficient_left < 0 < coefficient_right
        hyperbolic_rows.append(
            {
                "profile": list(profile),
                "T": total,
                "Q": square_total,
                "F_square": encode(f_square),
                "K_dot_F": encode(k_dot_f),
                "bar_E_square": encode(bar_e_square),
                "K_Y_dot_bar_E": encode(ky_dot_bar_e),
                "boundary_left_coefficient": encode(coefficient_left),
                "boundary_positive_rhs": encode(coefficient_right),
            }
        )

    tail_determinants = {}
    for length in range(1, 13):
        matrix = tail_matrix(length)
        leading = [
            determinant([row[:size] for row in matrix[:size]])
            for size in range(1, length + 1)
        ]
        expected = [
            Fraction((-1) ** size * (size + 1))
            for size in range(1, length)
        ] + [Fraction((-1) ** length)]
        assert leading == expected
        tail_determinants[str(length)] = encode(leading[-1])

    return {
        "schema": "rational-h3-hostile-incidence-trichotomy-hostile-audit-v1",
        "scope": {
            "claim": "exact local H3 reduction only",
            "semiampleness": "conditional",
            "absolute_picard_bound": False,
            "tail_bound": False,
        },
        "profile_identity": "Q-T=sum_i t_i(t_i-1)>=0",
        "low_profiles": [list(profile) for profile in sorted(low_profiles)],
        "spherical": {
            "physical_gram": spherical_physical,
            "host_gram": spherical_host,
            "leading_principal_determinants": [-1, 3, -8, 16],
            "fundamental_cycle": fundamental,
            "cycle_intersections": [encode(value) for value in intersections],
            "cycle_square": encode(z_square),
            "canonical_degree": encode(k_dot_z),
            "arithmetic_genus": encode(arithmetic_genus),
            "discrepancies": [encode(value) for value in discrepancy],
            "highest_root": highest_root,
            "highest_root_square": -2,
            "leaf_swap_odd_pairing_witness": encode(odd_pairing),
        },
        "affine": {
            "physical_null_cycle": affine_f,
            "host_null_root": affine_r,
            "K_dot_F": 4,
            "conditional_fibration_cases_m_g_rho_core": [
                list(row) for row in fibration_cases
            ],
        },
        "hyperbolic_samples": hyperbolic_rows,
        "tail_chain_determinants": tail_determinants,
        "verdict": "PASS_AFTER_CONTRACTION_REPAIR",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
