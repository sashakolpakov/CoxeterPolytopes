"""Small exact matrix routines used by certificate verification."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, TypeAlias

from .golden import Golden, RationalLike, ZERO


Scalar: TypeAlias = Golden | RationalLike
Matrix: TypeAlias = Sequence[Sequence[Scalar]]


def _coerce_matrix(matrix: Matrix, *, square: bool = False) -> list[list[Golden]]:
    rows = [[Golden.coerce(entry) for entry in row] for row in matrix]
    if not rows:
        return rows
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("matrix rows must all have the same length")
    if square and len(rows) != width:
        raise ValueError("matrix must be square")
    return rows


def determinant(matrix: Matrix) -> Golden:
    work = _coerce_matrix(matrix, square=True)
    size = len(work)
    result = Golden(1)
    sign = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return ZERO
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, size):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for index in range(column + 1, size):
                work[row][index] -= factor * work[column][index]
    return sign * result


def rank(matrix: Matrix) -> int:
    work = _coerce_matrix(matrix)
    if not work:
        return 0
    row_count = len(work)
    column_count = len(work[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next((row for row in range(pivot_row, row_count) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        for row in range(row_count):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for index in range(column, column_count):
                work[row][index] -= factor * work[pivot_row][index]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def rectangular_pairing(
    left: Matrix,
    gram: Matrix,
    right: Matrix,
) -> list[list[Golden]]:
    """Compute ``left.T * gram * right`` for rectangular column matrices."""

    middle = _coerce_matrix(gram, square=True)
    left_matrix = _coerce_matrix(left)
    right_matrix = _coerce_matrix(right)
    ambient = len(middle)
    if len(left_matrix) != ambient or len(right_matrix) != ambient:
        raise ValueError("column matrices and Gram matrix must have the same row count")
    left_width = len(left_matrix[0]) if left_matrix else 0
    right_width = len(right_matrix[0]) if right_matrix else 0
    return [
        [
            sum(
                left_matrix[k][i] * middle[k][ell] * right_matrix[ell][j]
                for k in range(ambient)
                for ell in range(ambient)
            )
            for j in range(right_width)
        ]
        for i in range(left_width)
    ]


@dataclass(frozen=True, slots=True)
class Inertia:
    positive: int
    negative: int
    zero: int

    @property
    def rank(self) -> int:
        return self.positive + self.negative


def inertia(matrix: Matrix, *, embedding: int = 0) -> Inertia:
    """Return exact inertia of a symmetric matrix over a real embedding.

    Symmetric congruence elimination implements Sylvester's law of inertia.
    A nonzero diagonal pivot contributes its sign.  When all diagonal entries
    vanish, a nonzero off-diagonal 2-by-2 pivot contributes one sign of each.
    """

    work = _coerce_matrix(matrix, square=True)
    if any(work[i][j] != work[j][i] for i in range(len(work)) for j in range(i)):
        raise ValueError("inertia is defined here only for symmetric matrices")
    return _inertia_symmetric(work, embedding)


def _inertia_symmetric(work: list[list[Golden]], embedding: int) -> Inertia:
    size = len(work)
    if size == 0:
        return Inertia(0, 0, 0)

    diagonal_pivot = next((index for index in range(size) if work[index][index]), None)
    if diagonal_pivot is not None:
        if diagonal_pivot != 0:
            order = [diagonal_pivot, *[i for i in range(size) if i != diagonal_pivot]]
            work = [[work[i][j] for j in order] for i in order]
        pivot = work[0][0]
        remainder = [
            [work[i][j] - work[i][0] * work[0][j] / pivot for j in range(1, size)]
            for i in range(1, size)
        ]
        tail = _inertia_symmetric(remainder, embedding)
        pivot_sign = pivot.sign(embedding)
        if pivot_sign > 0:
            return Inertia(tail.positive + 1, tail.negative, tail.zero)
        return Inertia(tail.positive, tail.negative + 1, tail.zero)

    pair = next(
        ((i, j) for i in range(size) for j in range(i + 1, size) if work[i][j]),
        None,
    )
    if pair is None:
        return Inertia(0, 0, size)

    first, second = pair
    order = [first, second, *[i for i in range(size) if i not in pair]]
    work = [[work[i][j] for j in order] for i in order]
    off_diagonal = work[0][1]
    remainder = [
        [
            work[i][j]
            - (work[i][0] * work[1][j] + work[i][1] * work[0][j]) / off_diagonal
            for j in range(2, size)
        ]
        for i in range(2, size)
    ]
    tail = _inertia_symmetric(remainder, embedding)
    return Inertia(tail.positive + 1, tail.negative + 1, tail.zero)
