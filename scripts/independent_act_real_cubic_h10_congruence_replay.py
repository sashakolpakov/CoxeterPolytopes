#!/usr/bin/env python3
"""Independent exact replay of the ACT H10 congruence obstruction.

This file deliberately imports no project search or certificate module.  It
implements Q(omega), omega^2+omega+1=0, with rational pairs, inverts the ACT
fixed frame independently, and checks the full 3^11 residue gate in chunks by
literal multiplication in F3[epsilon]/(epsilon^2).
"""

from __future__ import annotations

from fractions import Fraction
import itertools
import json
from pathlib import Path
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "act_real_cubic_h10_congruence_obstruction.json"
)
N = 11
E = tuple[Fraction, Fraction]
ZERO: E = (Fraction(0), Fraction(0))
ONE: E = (Fraction(1), Fraction(0))

CONSTRAINTS = np.array(
    [
        [0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2],
    ],
    dtype=np.int16,
)


def _e(a: int | Fraction = 0, b: int | Fraction = 0) -> E:
    return (Fraction(a), Fraction(b))


def _add(left: E, right: E) -> E:
    return (left[0] + right[0], left[1] + right[1])


def _neg(value: E) -> E:
    return (-value[0], -value[1])


def _mul(left: E, right: E) -> E:
    a, b = left
    c, d = right
    return (a * c - b * d, a * d + b * c - b * d)


def _scale(value: E, scalar: int | Fraction) -> E:
    scalar = Fraction(scalar)
    return (scalar * value[0], scalar * value[1])


def _conjugate(value: E) -> E:
    # conjugate(omega)=-1-omega.
    a, b = value
    return (a - b, -b)


def _inverse(value: E) -> E:
    a, b = value
    norm = a * a - a * b + b * b
    if norm == 0:
        raise ZeroDivisionError
    return _scale(_conjugate(value), Fraction(1, 1) / norm)


def _matrix(rows: int, columns: int) -> list[list[E]]:
    return [[ZERO for _ in range(columns)] for _ in range(rows)]


def _identity(size: int) -> list[list[E]]:
    answer = _matrix(size, size)
    for index in range(size):
        answer[index][index] = ONE
    return answer


def _matmul(left: list[list[E]], right: list[list[E]]) -> list[list[E]]:
    answer = _matrix(len(left), len(right[0]))
    for row in range(len(left)):
        for column in range(len(right[0])):
            value = ZERO
            for middle in range(len(right)):
                value = _add(value, _mul(left[row][middle], right[middle][column]))
            answer[row][column] = value
    return answer


def _transpose(matrix: list[list[E]]) -> list[list[E]]:
    return [list(row) for row in zip(*matrix)]


def _conjugate_matrix(matrix: list[list[E]]) -> list[list[E]]:
    return [[_conjugate(value) for value in row] for row in matrix]


def _matrix_inverse(matrix: list[list[E]]) -> list[list[E]]:
    augmented = [left[:] + right[:] for left, right in zip(matrix, _identity(len(matrix)))]
    size = len(matrix)
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if augmented[row][column] != ZERO
        )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        multiplier = _inverse(augmented[column][column])
        augmented[column] = [_mul(multiplier, value) for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            multiplier = augmented[row][column]
            if multiplier == ZERO:
                continue
            augmented[row] = [
                _add(value, _neg(_mul(multiplier, pivot_value)))
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [row[size:] for row in augmented]


def _act_frame() -> tuple[list[list[E]], list[list[E]], list[list[E]]]:
    omega = _e(0, 1)
    omega_bar = _e(-1, -1)
    theta = _e(1, 2)

    hermitian = _matrix(N, N)
    hermitian[0][0] = _e(3)
    for offset in (1, 5):
        for index in range(4):
            hermitian[offset + index][offset + index] = _e(3)
        for index in range(3):
            hermitian[offset + index][offset + index + 1] = theta
            hermitian[offset + index + 1][offset + index] = _neg(theta)
    hermitian[9][10] = theta
    hermitian[10][9] = _neg(theta)

    fixed = _matrix(N, N)
    fixed[0][0] = ONE
    column = 1
    for left, right in ((1, 4), (2, 3), (5, 8), (6, 7), (9, 10)):
        fixed[left][column] = fixed[right][column] = ONE
        column += 1
        fixed[left][column] = omega_bar
        fixed[right][column] = omega
        column += 1

    p4 = (
        (-1, -1, -1, 0),
        (-1, -1, 0, 0),
        (0, 1, 0, 0),
        (-1, -1, -1, -1),
    )
    diagonalizer = _matrix(N, N)
    diagonalizer[0][0] = ONE
    for offset in (1, 5):
        for row in range(4):
            for col in range(4):
                diagonalizer[offset + row][offset + col] = _e(p4[row][col])
    diagonalizer[9][9] = ONE
    diagonalizer[10][9] = ONE
    diagonalizer[10][10] = ONE

    fixed_gram = _matmul(
        _transpose(fixed), _matmul(hermitian, _conjugate_matrix(fixed))
    )
    q_fixed = [[_scale(value, Fraction(2, 3)) for value in row] for row in fixed_gram]
    diagonal_gram = _matmul(
        _transpose(diagonalizer), _matmul(q_fixed, diagonalizer)
    )
    expected = _matrix(N, N)
    for index in range(10):
        expected[index][index] = _e(2)
    expected[10][10] = _e(-2)
    if diagonal_gram != expected:
        raise ArithmeticError("independent fixed Gram audit failed")

    frame = _matmul(fixed, diagonalizer)
    inverse = _matrix_inverse(frame)
    if _matmul(frame, inverse) != _identity(N):
        raise ArithmeticError("independent frame inversion failed")
    adjoint = [[_scale(value, 3) for value in row] for row in inverse]
    if any(value[0].denominator != 1 or value[1].denominator != 1 for row in adjoint for value in row):
        raise ArithmeticError("3F^-1 is not Eisenstein integral")
    return frame, adjoint, diagonal_gram


def _dual_coefficients(matrix: list[list[E]]) -> tuple[np.ndarray, np.ndarray]:
    constant = np.empty((N, N), dtype=np.int16)
    epsilon = np.empty((N, N), dtype=np.int16)
    for row in range(N):
        for column in range(N):
            a, b = matrix[row][column]
            if a.denominator != 1 or b.denominator != 1:
                raise ArithmeticError("nonintegral Eisenstein residue")
            constant[row, column] = (int(a) + int(b)) % 3
            epsilon[row, column] = (2 * int(b)) % 3
    return constant, epsilon


def _literal_dual_number_gate(
    residues: np.ndarray,
    frame: list[list[E]],
    adjoint: list[list[E]],
) -> np.ndarray:
    fc, fe = _dual_coefficients(frame)
    jc, je = _dual_coefficients(adjoint)
    accepted = np.empty(len(residues), dtype=bool)
    for start in range(0, len(residues), 512):
        batch = residues[start : start + 512]
        left_c = (batch @ fc.T) % 3
        left_e = (batch @ fe.T) % 3
        q_batch = batch.copy()
        q_batch[:, -1] *= -1
        right_c = (q_batch @ jc) % 3
        right_e = (q_batch @ je) % 3
        constant_products = (
            left_c[:, :, None] * right_c[:, None, :]
        ) % 3
        epsilon_products = (
            left_c[:, :, None] * right_e[:, None, :]
            + left_e[:, :, None] * right_c[:, None, :]
        ) % 3
        accepted[start : start + len(batch)] = np.all(
            (constant_products == 0) & (epsilon_products == 0), axis=(1, 2)
        )
    return accepted


def replay() -> dict[str, object]:
    frame, adjoint, _ = _act_frame()
    residues = np.array(
        list(itertools.product(range(3), repeat=N)), dtype=np.int16
    )
    literal = _literal_dual_number_gate(residues, frame, adjoint)
    linear = np.all((residues @ CONSTRAINTS.T) % 3 == 0, axis=1)
    if not np.array_equal(literal, linear):
        raise ArithmeticError("independent literal and linear residue gates differ")
    code = residues[linear]
    q = (np.sum(code[:, :10] ** 2, axis=1) - code[:, 10] ** 2) % 3
    if not np.array_equal(q, code[:, 0] ** 2 % 3):
        raise ArithmeticError("independent Q=x0^2 identity failed")

    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    if not (
        int(np.sum(linear)) == archive["exact_extension_congruence"]["accepted_residue_count"] == 729
        and archive["primitive_reflection_root_gate"]["norm_two_excluded_by_code_identity"]
        and archive["primitive_reflection_root_gate"]["orthogonal_pair_count_in_the_complete_root_set"] == 0
        and not archive["global_finite_volume_theorem"]["Vinberg_first_finite_volume_gate_passes"]
    ):
        raise ArithmeticError("independent replay disagrees with the frozen theorem")
    return {
        "fixed_Gram": "2*diag(1^10,-1)",
        "literal_residue_count": int(np.sum(literal)),
        "linear_residue_count": int(np.sum(linear)),
        "literal_equals_five_equations": True,
        "Q_equals_x0_squared_on_code": True,
        "norm_two_roots": 0,
        "orthogonal_norm_one_root_pairs": 0,
        "Vinberg_first_finite_volume_gate": False,
        "status": "PASS",
    }


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
