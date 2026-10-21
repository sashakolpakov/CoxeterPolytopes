#!/usr/bin/env python3
"""Exact replay for the rational finite-vertex discrepancy/tail gate."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data/enumerations/rational_finite_vertex_discrepancy_tail_gate.json"
)

SPECIES_SQUARE = {"B": -4, "R": -2, "E": -1}
CANONICAL_DEGREE = {"B": 2, "R": 0, "E": -1}


def _solve(matrix: list[list[int]], rhs: list[int]) -> list[Fraction]:
    """Solve a nonsingular rational system by exact row reduction."""
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix) or len(rhs) != n:
        raise ValueError("expected a nonempty square matrix and matching rhs")
    aug = [
        [Fraction(value) for value in row] + [Fraction(rhs[i])]
        for i, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if aug[row][column]), None
        )
        if pivot is None:
            raise ValueError("singular matrix")
        aug[column], aug[pivot] = aug[pivot], aug[column]
        scale = aug[column][column]
        aug[column] = [value / scale for value in aug[column]]
        for row in range(n):
            if row == column:
                continue
            scale = aug[row][column]
            if scale:
                aug[row] = [
                    aug[row][j] - scale * aug[column][j]
                    for j in range(n + 1)
                ]
    return [aug[i][-1] for i in range(n)]


def _determinant(matrix: list[list[int]]) -> Fraction:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("expected a nonempty square matrix")
    work = [[Fraction(value) for value in row] for row in matrix]
    determinant = Fraction(1)
    for column in range(n):
        pivot = next(
            (row for row in range(column, n) if work[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant *= pivot_value
        for row in range(column + 1, n):
            scale = work[row][column] / pivot_value
            for j in range(column + 1, n):
                work[row][j] -= scale * work[column][j]
    return determinant


def discrepancy_vector(
    species: list[str], intersections: dict[tuple[int, int], int]
) -> tuple[list[list[int]], list[Fraction]]:
    """Return the physical intersection matrix and discrepancy vector."""
    if any(item not in SPECIES_SQUARE for item in species):
        raise ValueError("unknown physical species")
    n = len(species)
    matrix = [[0 for _ in range(n)] for _ in range(n)]
    for i, item in enumerate(species):
        matrix[i][i] = SPECIES_SQUARE[item]
    for (i, j), value in intersections.items():
        if not (0 <= i < j < n) or value < 0:
            raise ValueError("invalid intersection entry")
        matrix[i][j] = matrix[j][i] = value
    rhs = [CANONICAL_DEGREE[item] for item in species]
    return matrix, _solve(matrix, rhs)


def weak_barrier_residual(
    species: list[str], intersections: dict[tuple[int, int], int]
) -> list[Fraction]:
    """Compute (-Q)u-k for u_B=1/2 and u_R=u_E=0."""
    matrix, _ = discrepancy_vector(species, intersections)
    u = [Fraction(1, 2) if item == "B" else Fraction(0) for item in species]
    k = [Fraction(CANONICAL_DEGREE[item]) for item in species]
    return [
        sum(Fraction(-matrix[i][j]) * u[j] for j in range(len(species)))
        - k[i]
        for i in range(len(species))
    ]


def _row(
    row_id: str,
    species: list[str],
    intersections: dict[tuple[int, int], int],
) -> dict[str, Any]:
    matrix, discrepancy = discrepancy_vector(species, intersections)
    return {
        "id": row_id,
        "species": species,
        "intersection_matrix": matrix,
        "determinant": str(_determinant(matrix)),
        "discrepancies": [str(value) for value in discrepancy],
        "log_discrepancies": [str(value + 1) for value in discrepancy],
        "weak_barrier_residual": [
            str(value) for value in weak_barrier_residual(species, intersections)
        ],
    }


def hostile_controls() -> list[dict[str, Any]]:
    return [
        _row("H1_BB_intersection_2", ["B", "B"], {(0, 1): 2}),
        _row("H2_BR_intersection_2", ["B", "R"], {(0, 1): 2}),
        _row(
            "H3_E_three_B_leaves",
            ["E", "B", "B", "B"],
            {(0, 1): 1, (0, 2): 1, (0, 3): 1},
        ),
    ]


def safe_controls() -> list[dict[str, Any]]:
    return [
        _row("single_boundary", ["B"], {}),
        _row("B_E_B", ["B", "E", "B"], {(0, 1): 1, (1, 2): 1}),
        _row(
            "B_E_R_R",
            ["B", "E", "R", "R"],
            {(0, 1): 1, (1, 2): 1, (2, 3): 1},
        ),
    ]


def exceptional_tail(length: int) -> dict[str, Any]:
    """The chain R_1-...-R_{m-1}-E_m from m successive blowups."""
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ValueError("tail length must be a positive integer")
    species = ["R"] * (length - 1) + ["E"]
    intersections = {(i, i + 1): 1 for i in range(length - 1)}
    row = _row(f"exceptional_tail_{length}", species, intersections)
    matrix = row["intersection_matrix"]
    row["leading_principal_determinants"] = [
        str(_determinant([values[:size] for values in matrix[:size]]))
        for size in range(1, length + 1)
    ]
    row["expected_discrepancies"] = [str(i) for i in range(1, length + 1)]
    discrepancies = [Fraction(value) for value in row["discrepancies"]]
    correction = -sum(
        discrepancies[i] * Fraction(matrix[i][j]) * discrepancies[j]
        for i in range(length)
        for j in range(length)
    )
    if correction.denominator != 1:
        raise AssertionError("exceptional-tail correction is not integral")
    row["canonical_square_loss_before_contraction"] = length
    row["discrepancy_correction"] = correction.numerator
    row["canonical_square_is_preserved_by_contraction"] = correction == length
    return row


def build_certificate() -> dict[str, Any]:
    safe = safe_controls()
    hostile = hostile_controls()
    tails = [exceptional_tail(length) for length in range(1, 33)]
    payload = {
        "schema_version": 1,
        "theorem": {
            "barrier": (
                "For a negative-definite weak-incidence packet, "
                "a_B>=-1/2 and a_R,a_E>=0."
            ),
            "singularity_scope": (
                "For a normal Q-Gorenstein contraction whose full exceptional "
                "divisor is the SNC packet, all log discrepancies are >=1/2; "
                "under the strict convention this is epsilon-klt for every "
                "0<epsilon<1/2, not necessarily 1/2-klt."
            ),
            "tail_obstruction": (
                "Finite-nef-vertex discrepancy bounds alone cannot bound rank: "
                "successive point blowups of P2 give arbitrary R...R-E tails; "
                "no global full-nef Coxeter property is asserted."
            ),
        },
        "safe_controls": safe,
        "sharp_hostile_controls": hostile,
        "exceptional_tail_controls": tails,
    }
    certify(payload)
    payload["all_checks_pass"] = True
    return payload


def certify(payload: dict[str, Any]) -> None:
    safe = payload["safe_controls"]
    if not all(
        all(value >= Fraction(0) for value in map(Fraction, row["weak_barrier_residual"]))
        and all(value >= Fraction(1, 2) for value in map(Fraction, row["log_discrepancies"]))
        for row in safe
    ):
        raise AssertionError("a safe control violates the discrepancy barrier")
    expected_hostile = {
        "H1_BB_intersection_2": ["-1", "-1"],
        "H2_BR_intersection_2": ["-1", "-1"],
        "H3_E_three_B_leaves": ["-2", "-1", "-1", "-1"],
    }
    for row in payload["sharp_hostile_controls"]:
        if row["discrepancies"] != expected_hostile[row["id"]]:
            raise AssertionError(f"wrong hostile discrepancy row: {row['id']}")
    for length, row in enumerate(payload["exceptional_tail_controls"], start=1):
        expected_leading = [
            str(((-1) ** size) * (size + 1))
            for size in range(1, length)
        ] + [str((-1) ** length)]
        if row["leading_principal_determinants"] != expected_leading:
            raise AssertionError("exceptional tail is not negative definite")
        if row["determinant"] != str((-1) ** length):
            raise AssertionError("exceptional tail is not unimodular")
        if row["discrepancies"] != row["expected_discrepancies"]:
            raise AssertionError("wrong exceptional-tail discrepancy vector")
        if row["discrepancy_correction"] != length:
            raise AssertionError("wrong canonical-square correction")
        if not row["canonical_square_is_preserved_by_contraction"]:
            raise AssertionError("tail contraction failed to restore K^2")


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
    print("PASS: weak finite-vertex discrepancy barrier and arbitrary tail")


if __name__ == "__main__":
    main()
