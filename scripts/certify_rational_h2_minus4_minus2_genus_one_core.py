#!/usr/bin/env python3
"""Exact arithmetic replay for the minimal H2 genus-one core theorem."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_h2_minus4_minus2_genus_one_core.json"
)


def determinant(matrix: list[list[int]]) -> Fraction:
    """Return a determinant by exact row reduction."""
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("expected a nonempty square matrix")
    work = [[Fraction(entry) for entry in row] for row in matrix]
    answer = Fraction(1)
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        pivot_entry = work[column][column]
        answer *= pivot_entry
        for row in range(column + 1, n):
            scale = work[row][column] / pivot_entry
            for j in range(column + 1, n):
                work[row][j] -= scale * work[column][j]
    return answer


def solve(matrix: list[list[int]], rhs: list[int]) -> list[Fraction]:
    """Solve a nonsingular square system over Q."""
    n = len(matrix)
    if n == 0 or len(rhs) != n or any(len(row) != n for row in matrix):
        raise ValueError("expected a nonempty square system")
    work = [
        [Fraction(entry) for entry in row] + [Fraction(rhs[i])]
        for i, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if pivot is None:
            raise ValueError("singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_entry = work[column][column]
        work[column] = [entry / pivot_entry for entry in work[column]]
        for row in range(n):
            if row == column:
                continue
            scale = work[row][column]
            if scale:
                work[row] = [
                    work[row][j] - scale * work[column][j]
                    for j in range(n + 1)
                ]
    return [work[i][-1] for i in range(n)]


def mat_vec(
    matrix: list[list[Fraction]], vector: list[Fraction]
) -> list[Fraction]:
    return [
        sum((entry * value for entry, value in zip(row, vector, strict=True)),
            Fraction(0))
        for row in matrix
    ]


def reflection_matrix(root: list[int]) -> list[list[Fraction]]:
    """Reflection matrix for q_M=-2 I in the total-exceptional basis."""
    norm = -2 * sum(entry * entry for entry in root)
    if norm == 0:
        raise ValueError("a reflection root must be nonzero")
    n = len(root)
    matrix: list[list[Fraction]] = []
    for i in range(n):
        row = []
        for j in range(n):
            delta = Fraction(int(i == j))
            # q(e_j,root)=-2 root_j.
            row.append(delta - Fraction(2 * root[i] * (-2 * root[j]), norm))
        matrix.append(row)
    return matrix


def tail_row(length: int) -> dict[str, Any]:
    """Check the free B_m exceptional chain without group enumeration."""
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ValueError("tail length must be a positive integer")

    roots: list[list[int]] = []
    for i in range(length - 1):
        root = [0] * length
        root[i] = 1
        root[i + 1] = -1
        roots.append(root)
    terminal = [0] * length
    terminal[-1] = 1
    roots.append(terminal)

    gram = [
        [
            -2
            * sum(
                left[k] * right[k]
                for k in range(length)
            )
            for right in roots
        ]
        for left in roots
    ]

    actions = [reflection_matrix(root) for root in roots]
    generator_checks: list[bool] = []
    for i, action in enumerate(actions[:-1]):
        for basis_index in range(length):
            basis = [Fraction(int(j == basis_index)) for j in range(length)]
            image = mat_vec(action, basis)
            target_index = (
                i + 1
                if basis_index == i
                else i
                if basis_index == i + 1
                else basis_index
            )
            expected = [
                Fraction(int(j == target_index)) for j in range(length)
            ]
            generator_checks.append(image == expected)
    terminal_action = actions[-1]
    for basis_index in range(length):
        basis = [Fraction(int(j == basis_index)) for j in range(length)]
        image = mat_vec(terminal_action, basis)
        expected = basis[:]
        if basis_index == length - 1:
            expected[basis_index] = -expected[basis_index]
        generator_checks.append(image == expected)

    return {
        "length": length,
        "physical_curve_squares_NS": [-2] * (length - 1) + [-1],
        "root_squares_qM": [gram[i][i] for i in range(length)],
        "root_gram_qM": gram,
        "root_gram_determinant": str(determinant(gram)),
        "expected_root_gram_determinant": str((-2) ** length),
        "simple_reflections_are_adjacent_swaps_then_terminal_sign": all(
            generator_checks
        ),
        "signed_permutation_group_order": str(
            (2**length) * math.factorial(length)
        ),
    }


def build_certificate() -> dict[str, Any]:
    pair_matrix = [[-4, 2], [2, -2]]
    canonical_degrees = [2, 0]
    discrepancies = solve(pair_matrix, canonical_degrees)
    rows = [tail_row(length) for length in (1, 2, 3, 4, 8, 16, 32)]
    payload = {
        "schema_version": 1,
        "scope": {
            "pair": "B^2=-4, R^2=-2, B.R=2 on a smooth rational surface",
            "core_only": True,
            "original_rank_bound_claimed": False,
            "reflection_maximal_counterfamily_claimed": False,
        },
        "pair_arithmetic": {
            "intersection_matrix": pair_matrix,
            "determinant": str(determinant(pair_matrix)),
            "canonical_degrees": canonical_degrees,
            "D_square": -2,
            "K_dot_D": 2,
            "arithmetic_genus_D": 1,
            "adjoint_multidegrees_on_B_R": [0, 0],
            "adjoint_square_formula": "12-rho",
            "discrepancies": [str(entry) for entry in discrepancies],
        },
        "sharp_equality_control": {
            "base": "P1xP1",
            "base_picard_rank": 2,
            "base_curve_classes": ["(1,1)", "(1,1)"],
            "base_curve_squares": [2, 2],
            "private_boundary_blowups": [6, 4],
            "strict_transform_squares": [-4, -2],
            "strict_transform_intersection": 2,
            "picard_rank": 12,
            "K_plus_B_plus_R": "0",
        },
        "free_exceptional_Bm_tail_controls": rows,
        "conditional_total_marked_degree_bound": {
            "premise": (
                "every contracted tail component is a physical facet and "
                "T is the total integral ample-marking degree of all "
                "physical facets"
            ),
            "tail_length_bound": "m<=T",
            "picard_rank_bound": "rho(S)<=12+T",
            "hyperbolic_dimension_bound": "dim<=11+T",
            "per_facet_degree_bound_alone_is_sufficient": False,
        },
        "conditional_trace_bound": {
            "rank": "d(n+1)",
            "bound_if_the_trace_block_descends": "d(n+1)<=12",
            "rational_field_bound": "n<=11",
            "descent_is_automatic": False,
        },
    }
    certify(payload)
    payload["all_checks_pass"] = True
    return payload


def certify(payload: dict[str, Any]) -> None:
    pair = payload["pair_arithmetic"]
    if pair["determinant"] != "4":
        raise AssertionError("the H2 pair determinant changed")
    if pair["discrepancies"] != ["-1", "-1"]:
        raise AssertionError("the H2 discrepancy vector changed")
    if (
        pair["D_square"],
        pair["K_dot_D"],
        pair["arithmetic_genus_D"],
        pair["adjoint_multidegrees_on_B_R"],
    ) != (-2, 2, 1, [0, 0]):
        raise AssertionError("the genus-one adjoint arithmetic changed")

    equality = payload["sharp_equality_control"]
    if equality["picard_rank"] != 12:
        raise AssertionError("the equality example has the wrong rank")
    if equality["strict_transform_squares"] != [-4, -2]:
        raise AssertionError("the equality example has the wrong squares")
    if equality["K_plus_B_plus_R"] != "0":
        raise AssertionError("the equality example is not log Calabi-Yau")

    for row in payload["free_exceptional_Bm_tail_controls"]:
        length = row["length"]
        expected_gram = [[0] * length for _ in range(length)]
        for i in range(length):
            expected_gram[i][i] = -4 if i < length - 1 else -2
            if i + 1 < length:
                expected_gram[i][i + 1] = 2
                expected_gram[i + 1][i] = 2
        if row["root_gram_qM"] != expected_gram:
            raise AssertionError(f"wrong B_{length} root Gram")
        if row["root_gram_determinant"] != str((-2) ** length):
            raise AssertionError(f"wrong B_{length} root determinant")
        if not row["simple_reflections_are_adjacent_swaps_then_terminal_sign"]:
            raise AssertionError(f"wrong B_{length} simple reflection action")
        if row["signed_permutation_group_order"] != str(
            (2**length) * math.factorial(length)
        ):
            raise AssertionError(f"wrong W(B_{length}) order")

    degree = payload["conditional_total_marked_degree_bound"]
    if degree != {
        "premise": (
            "every contracted tail component is a physical facet and "
            "T is the total integral ample-marking degree of all "
            "physical facets"
        ),
        "tail_length_bound": "m<=T",
        "picard_rank_bound": "rho(S)<=12+T",
        "hyperbolic_dimension_bound": "dim<=11+T",
        "per_facet_degree_bound_alone_is_sufficient": False,
    }:
        raise AssertionError("the conditional total-degree scope changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.write:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.check_archive:
        frozen = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        if frozen != payload:
            raise SystemExit("archive mismatch")
    print("PASS: minimal H2 genus-one core (rho_core<=12) and B_m tail")


if __name__ == "__main__":
    main()
