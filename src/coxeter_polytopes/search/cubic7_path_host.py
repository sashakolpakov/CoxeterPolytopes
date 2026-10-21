"""Exact rank-23 conductor-seven terminal-dotted path host.

Put ``theta=2*cos(2*pi/7)`` and ``c=2+theta``.  The symmetric Gram matrix
has diagonal ``(2*c,2^22)``, first off diagonal ``-c``, and all remaining
path entries ``-1``.  This module freezes host arithmetic and the exact
geometry of the initial 23-wall path simplex; it makes no Vinberg-completeness
claim above height zero.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    LABEL7_PRODUCT,
    Cubic7,
    cubic7_determinant,
    cubic7_inertia,
)
from coxeter_polytopes.search.cubic7_diagonal import bounded_cubic7_integers


Vector: TypeAlias = tuple[Cubic7, ...]
Matrix: TypeAlias = tuple[tuple[Cubic7, ...], ...]

NATIVE_RANK = 23
DIMENSION = 22
SPATIAL_RANK = 22
COEFFICIENT = LABEL7_PRODUCT  # 2+theta = 4*cos(pi/7)^2
DETERMINANT = Cubic7(4, -42, -22)
LORENTZ_THRESHOLD = Cubic7(23) / 11


@lru_cache(maxsize=1)
def lattice_gram() -> Matrix:
    rows = [[CUBIC7_ZERO for _ in range(NATIVE_RANK)] for _ in range(NATIVE_RANK)]
    rows[0][0] = 2 * COEFFICIENT
    for index in range(1, NATIVE_RANK):
        rows[index][index] = Cubic7(2)
    rows[0][1] = rows[1][0] = -COEFFICIENT
    for index in range(1, NATIVE_RANK - 1):
        rows[index][index + 1] = rows[index + 1][index] = -CUBIC7_ONE
    return tuple(tuple(row) for row in rows)


def basis_vector(index: int) -> Vector:
    if not 0 <= index < NATIVE_RANK:
        raise ValueError("basis index is out of range")
    return tuple(CUBIC7_ONE if index == j else CUBIC7_ZERO for j in range(NATIVE_RANK))


def pair(left: Sequence[Cubic7], right: Sequence[Cubic7]) -> Cubic7:
    if len(left) != NATIVE_RANK or len(right) != NATIVE_RANK:
        raise ValueError("cubic-seven path vectors must have rank 23")
    gram = lattice_gram()
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(NATIVE_RANK)
            for j in range(NATIVE_RANK)
            if left[i] and gram[i][j] and right[j]
        ),
        start=CUBIC7_ZERO,
    )


@lru_cache(maxsize=1)
def weyl_vector() -> Vector:
    """Vector orthogonal to the 22 spatial simple roots."""

    return (
        CUBIC7_ONE,
        *(
            COEFFICIENT * Cubic7(NATIVE_RANK - index) / NATIVE_RANK
            for index in range(1, NATIVE_RANK)
        ),
    )


@lru_cache(maxsize=1)
def inverse_gram() -> Matrix:
    gram = lattice_gram()
    size = len(gram)
    work = [
        [*row, *(CUBIC7_ONE if i == j else CUBIC7_ZERO for j in range(size))]
        for i, row in enumerate(gram)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [value / pivot_value for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                work[row][j] - multiplier * work[column][j]
                for j in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


def seed_simplex_vertex_types() -> tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]:
    """Omitted-facet indices of finite, ideal, and hyperideal simplex rays."""

    diagonal = tuple(inverse_gram()[i][i] for i in range(NATIVE_RANK))
    finite = tuple(i for i, value in enumerate(diagonal) if value.sign(0) < 0)
    cusps = tuple(i for i, value in enumerate(diagonal) if value.sign(0) == 0)
    hyperideal = tuple(i for i, value in enumerate(diagonal) if value.sign(0) > 0)
    return finite, cusps, hyperideal


def admissible_terminal_coefficient(value: Cubic7) -> bool:
    """Test the exact rank-23 path signature inequalities.

    The Schur complement of the spatial ``A22`` block is
    ``2*value-(22/23)*value^2``.  It is negative at the selected place
    exactly when ``value>23/11`` and positive at a conjugate place exactly
    when ``0<value<23/11``.
    """

    value = Cubic7.coerce(value)
    return (value - LORENTZ_THRESHOLD).sign(0) > 0 and all(
        value.sign(place) > 0
        and (LORENTZ_THRESHOLD - value).sign(place) > 0
        for place in (1, 2)
    )


@lru_cache(maxsize=None)
def small_admissible_coefficients(coefficient_bound: int = 10) -> tuple[Cubic7, ...]:
    """Exact coefficient-box census used by the frozen discovery screen."""

    if coefficient_bound < 0:
        raise ValueError("coefficient bound must be nonnegative")
    values = (
        Cubic7(a, b, c)
        for a in range(-coefficient_bound, coefficient_bound + 1)
        for b in range(-coefficient_bound, coefficient_bound + 1)
        for c in range(-coefficient_bound, coefficient_bound + 1)
    )
    return tuple(sorted(value for value in values if admissible_terminal_coefficient(value)))


@lru_cache(maxsize=1)
def admissible_coefficients_through_selected_coefficient() -> tuple[Cubic7, ...]:
    """Complete Minkowski census from threshold through ``2+theta``.

    The signature inequalities imply squared embedding caps ``(16,5,5)``;
    the repository's trace-dual enumerator therefore makes this finite box
    complete rather than a coefficient heuristic.
    """

    return tuple(
        value
        for value in bounded_cubic7_integers((16, 5, 5))
        if admissible_terminal_coefficient(value)
        and (COEFFICIENT - value).sign(0) >= 0
    )


def exact_host_checks() -> None:
    gram = lattice_gram()
    if cubic7_determinant(gram) != DETERMINANT:
        raise AssertionError("cubic-seven path determinant changed")
    expected = ((22, 1, 0), (23, 0, 0), (23, 0, 0))
    actual = tuple(
        (item.positive, item.negative, item.zero)
        for item in (cubic7_inertia(gram, embedding=place) for place in range(3))
    )
    if actual != expected:
        raise AssertionError("cubic-seven path inertias changed")
    if seed_simplex_vertex_types() != ((0, 1, 2), (), tuple(range(3, 23))):
        raise AssertionError("cubic-seven path seed geometry changed")


__all__ = [
    "COEFFICIENT",
    "DETERMINANT",
    "DIMENSION",
    "LORENTZ_THRESHOLD",
    "NATIVE_RANK",
    "SPATIAL_RANK",
    "admissible_coefficients_through_selected_coefficient",
    "admissible_terminal_coefficient",
    "basis_vector",
    "exact_host_checks",
    "inverse_gram",
    "lattice_gram",
    "pair",
    "seed_simplex_vertex_types",
    "small_admissible_coefficients",
    "weyl_vector",
]
