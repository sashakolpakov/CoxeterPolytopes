"""Restriction-of-scalars trace lattices for Q(sqrt(5))."""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

from .golden import Golden, ONE, PHI, RationalLike


SQRT5 = 2 * PHI - 1
INV_SQRT5 = SQRT5 / 5
GOLDEN_RING_BASIS = (ONE, PHI)


def is_codifferent(element: Golden | RationalLike) -> bool:
    """Whether ``element`` lies in the codifferent of Z[phi]."""

    return (Golden.coerce(element) * SQRT5).is_integral


def trace_gram(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
    alpha: Golden | RationalLike,
) -> list[list[Fraction]]:
    """Return the Q-trace Gram matrix in the interleaved basis e_i, phi*e_i.

    The resulting pairing is ``Tr(alpha * h(x,y))``.  Coefficients are
    returned as Fractions so nonintegral candidate trace lattices remain
    representable and can be rejected explicitly.
    """

    gram = [[Golden.coerce(entry) for entry in row] for row in field_gram]
    rank = len(gram)
    if any(len(row) != rank for row in gram):
        raise ValueError("field Gram matrix must be square")
    if any(gram[i][j] != gram[j][i] for i in range(rank) for j in range(i)):
        raise ValueError("field Gram matrix must be symmetric")

    alpha = Golden.coerce(alpha)
    result = [[Fraction(0) for _ in range(2 * rank)] for _ in range(2 * rank)]
    for i in range(rank):
        for j in range(rank):
            for basis_i, multiplier_i in enumerate(GOLDEN_RING_BASIS):
                for basis_j, multiplier_j in enumerate(GOLDEN_RING_BASIS):
                    result[2 * i + basis_i][2 * j + basis_j] = (
                        alpha * multiplier_i * gram[i][j] * multiplier_j
                    ).trace()
    return result


def integral_gram(
    gram: Sequence[Sequence[Fraction | int]],
) -> list[list[int]]:
    """Convert an integral rational Gram matrix to integers, or reject it."""

    result: list[list[int]] = []
    for row in gram:
        integer_row: list[int] = []
        for entry in row:
            entry = Fraction(entry)
            if entry.denominator != 1:
                raise ValueError(f"trace pairing is not integral: {entry}")
            integer_row.append(entry.numerator)
        result.append(integer_row)
    return result


def is_even_integral(gram: Sequence[Sequence[Fraction | int]]) -> bool:
    try:
        integer = integral_gram(gram)
    except ValueError:
        return False
    return all(integer[index][index] % 2 == 0 for index in range(len(integer)))


def change_gram_basis(
    gram: Sequence[Sequence[Golden | RationalLike]],
    change: Sequence[Sequence[Golden | RationalLike]],
) -> list[list[Golden]]:
    """Compute ``change.T * gram * change`` exactly."""

    source = [[Golden.coerce(entry) for entry in row] for row in gram]
    transform = [[Golden.coerce(entry) for entry in row] for row in change]
    size = len(source)
    if any(len(row) != size for row in source):
        raise ValueError("Gram matrix must be square")
    if len(transform) != size or any(len(row) != size for row in transform):
        raise ValueError("basis-change matrix must have the same square size")
    return [
        [
            sum(
                transform[k][i] * source[k][ell] * transform[ell][j]
                for k in range(size)
                for ell in range(size)
            )
            for j in range(size)
        ]
        for i in range(size)
    ]

