"""Exact divisor-slice carriers on blow-ups of the projective plane.

The Neron--Severi lattice of the blow-up of ``P^2`` at ``b`` distinct
points is the odd unimodular lattice ``I_(1,b)``.  This module implements a
constructive rational version of the elementary fact that every integral
Lorentzian lattice embeds, after one common square rescaling, in such a
lattice.  It also extends a self-adjoint operator on the embedded subspace to
an integral self-adjoint operator of the ambient divisor lattice.

The construction is deliberately a *weak divisor-slice carrier*: it does not
assert effectivity of the embedded classes, equality with a nef cone, or
realization of the reflections by surface automorphisms.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import lcm
from typing import Sequence

import sympy as sp
from sympy.solvers.diophantine.diophantine import sum_of_four_squares


@dataclass(frozen=True)
class LorentzianEmbedding:
    """An exact embedding with ``E.T * J * E = vector_scale^2 * G``."""

    vector_scale: int
    embedding: sp.ImmutableMatrix
    ambient_gram: sp.ImmutableMatrix
    rational_diagonal: tuple[sp.Rational, ...]

    @property
    def blowup_count(self) -> int:
        return self.ambient_gram.rows - 1


@dataclass(frozen=True)
class IntegralOperatorExtension:
    """An integral self-adjoint extension of a source operator."""

    operator_multiplier: int
    operator: sp.ImmutableMatrix
    rational_projection: sp.ImmutableMatrix
    complement_scalar: sp.Rational


def _matrix(value: Sequence[Sequence[object]] | sp.MatrixBase) -> sp.Matrix:
    matrix = sp.Matrix(value)
    if matrix.rows == 0 or matrix.rows != matrix.cols:
        raise ValueError("the Gram matrix must be nonempty and square")
    return matrix


def _denominator_lcm(matrix: sp.MatrixBase) -> int:
    result = 1
    for entry in matrix:
        result = lcm(result, int(sp.denom(entry)))
    return result


def rational_congruence_diagonalization(
    gram: Sequence[Sequence[object]] | sp.MatrixBase,
) -> tuple[sp.ImmutableMatrix, tuple[sp.Rational, ...]]:
    """Return ``C,d`` with ``C.T*gram*C=diag(d)`` over ``Q``.

    The Gram--Schmidt step also handles a basis whose remaining vectors are
    all isotropic: nondegeneracy then supplies a pair with nonzero pairing,
    and their sum is anisotropic.
    """

    source = _matrix(gram)
    if source != source.T or source.det() == 0:
        raise ValueError("the Gram matrix must be symmetric and nondegenerate")
    if any(not entry.is_Rational for entry in source):
        raise ValueError("the Gram matrix must have rational entries")

    remaining = [sp.eye(source.rows).col(index) for index in range(source.rows)]
    orthogonal_basis: list[sp.Matrix] = []
    while remaining:
        pivot = next(
            (
                index
                for index, vector in enumerate(remaining)
                if (vector.T * source * vector)[0] != 0
            ),
            None,
        )
        if pivot is None:
            paired: tuple[int, int] | None = None
            for left in range(len(remaining)):
                for right in range(left + 1, len(remaining)):
                    if (remaining[left].T * source * remaining[right])[0] != 0:
                        paired = (left, right)
                        break
                if paired is not None:
                    break
            if paired is None:  # pragma: no cover - guarded by nondegeneracy
                raise AssertionError("nondegenerate restricted form has no pivot")
            left, right = paired
            remaining[left] = remaining[left] + remaining[right]
            pivot = left

        vector = remaining.pop(pivot)
        square = (vector.T * source * vector)[0]
        orthogonal_basis.append(vector)
        remaining = [
            other - (vector.T * source * other)[0] / square * vector
            for other in remaining
        ]

    change = sp.Matrix.hstack(*orthogonal_basis)
    diagonalized = change.T * source * change
    if any(
        diagonalized[row, column] != 0
        for row in range(source.rows)
        for column in range(source.cols)
        if row != column
    ):
        raise AssertionError("rational Gram--Schmidt did not diagonalize the form")
    diagonal = tuple(sp.Rational(diagonalized[index, index]) for index in range(source.rows))
    return sp.ImmutableMatrix(change), diagonal


def embed_lorentzian_lattice(
    gram: Sequence[Sequence[object]] | sp.MatrixBase,
) -> LorentzianEmbedding:
    """Embed an integral form of signature ``(1,m)`` in ``I_(1,1+4m)``.

    If ``G`` is the input Gram matrix, the returned integral matrix ``E``
    satisfies ``E.T*J*E = D^2*G``.  The common square factor changes neither
    hyperplanes nor Coxeter angles.

    The proof is constructive.  Rationally diagonalize ``G``.  Represent its
    positive coefficient as a difference of two rational squares and each
    negative coefficient as minus a sum of four rational squares.  Lagrange's
    four-square theorem supplies the latter representations; clearing all
    denominators gives ``E``.
    """

    source = _matrix(gram)
    if any(not entry.is_Integer for entry in source):
        raise ValueError("the source lattice must be integral")
    change_raw, diagonal_raw = rational_congruence_diagonalization(source)
    positive = [index for index, entry in enumerate(diagonal_raw) if entry > 0]
    negative = [index for index, entry in enumerate(diagonal_raw) if entry < 0]
    if len(positive) != 1 or len(positive) + len(negative) != source.rows:
        raise ValueError("the source lattice must have signature (1,m)")

    order = positive + negative
    change = sp.Matrix(change_raw)[:, order]
    diagonal = tuple(diagonal_raw[index] for index in order)
    negative_rank = len(negative)
    diagonal_embedding = sp.zeros(2 + 4 * negative_rank, source.rows)

    positive_square = diagonal[0]
    diagonal_embedding[0, 0] = (positive_square + 1) / 2
    diagonal_embedding[1, 0] = (positive_square - 1) / 2

    for index, entry in enumerate(diagonal[1:]):
        value = -entry
        numerator = int(sp.numer(value))
        denominator = int(sp.denom(value))
        four_squares = sum_of_four_squares(numerator * denominator)
        for offset, square_root in enumerate(four_squares):
            diagonal_embedding[2 + 4 * index + offset, index + 1] = sp.Rational(
                square_root, denominator
            )

    rational_embedding = diagonal_embedding * change.inv()
    vector_scale = _denominator_lcm(rational_embedding)
    embedding = vector_scale * rational_embedding
    ambient = sp.diag(1, *([-1] * (embedding.rows - 1)))
    if any(not entry.is_Integer for entry in embedding):
        raise AssertionError("denominator clearing did not give an integral embedding")
    if embedding.T * ambient * embedding != vector_scale**2 * source:
        raise AssertionError("the Lorentzian embedding identity failed")

    return LorentzianEmbedding(
        vector_scale=vector_scale,
        embedding=sp.ImmutableMatrix(embedding),
        ambient_gram=sp.ImmutableMatrix(ambient),
        rational_diagonal=diagonal,
    )


def extend_self_adjoint_operator(
    embedding: Sequence[Sequence[object]] | sp.MatrixBase,
    ambient_gram: Sequence[Sequence[object]] | sp.MatrixBase,
    source_operator: Sequence[Sequence[object]] | sp.MatrixBase,
    *,
    complement_scalar: object = 0,
) -> IntegralOperatorExtension:
    """Extend a self-adjoint source operator and clear denominators.

    On the rational image of ``E`` the extension is ``T``.  On its orthogonal
    complement it is the chosen rational scalar.  Multiplication by the
    returned positive integer makes the whole ambient matrix integral without
    changing any real eigenspace.
    """

    embed = sp.Matrix(embedding)
    ambient = _matrix(ambient_gram)
    operator = sp.Matrix(source_operator)
    if embed.rows != ambient.rows or operator.rows != embed.cols or operator.cols != embed.cols:
        raise ValueError("incompatible embedding, ambient form, and operator sizes")
    if ambient != ambient.T or ambient.det() == 0:
        raise ValueError("the ambient form must be symmetric and nondegenerate")
    if any(not entry.is_Rational for entry in embed):
        raise ValueError("the embedding must have rational entries")
    if any(not entry.is_Rational for entry in ambient):
        raise ValueError("the ambient form must have rational entries")
    if any(not entry.is_Rational for entry in operator):
        raise ValueError("the source operator must have rational entries")
    source_gram = embed.T * ambient * embed
    if source_gram.det() == 0 or operator.T * source_gram != source_gram * operator:
        raise ValueError("the source operator must be self-adjoint and nondegenerate")

    projection = embed * source_gram.inv() * embed.T * ambient
    scalar = sp.Rational(complement_scalar)
    rational_extension = (
        embed * operator * source_gram.inv() * embed.T * ambient
        + scalar * (sp.eye(ambient.rows) - projection)
    )
    multiplier = _denominator_lcm(rational_extension)
    integral_extension = multiplier * rational_extension
    if any(not entry.is_Integer for entry in integral_extension):
        raise AssertionError("operator denominator clearing failed")
    if integral_extension.T * ambient != ambient * integral_extension:
        raise AssertionError("the ambient operator is not self-adjoint")
    if integral_extension * embed != embed * (multiplier * operator):
        raise AssertionError("the ambient operator has the wrong restriction")

    return IntegralOperatorExtension(
        operator_multiplier=multiplier,
        operator=sp.ImmutableMatrix(integral_extension),
        rational_projection=sp.ImmutableMatrix(projection),
        complement_scalar=scalar,
    )
