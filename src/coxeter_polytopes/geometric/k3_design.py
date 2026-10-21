"""Exact necessary constraints for inverse-designing a K3 Coxeter section.

The compact example track uses an even ``Z[phi]`` lattice whose simple field
roots have norm two.  With the K3 trace scalar ``-phi/sqrt(5)``, their field
Gram matrix has a second interpretation: it is the intersection matrix of
the ambient trace-lattice ``(-2)`` classes.  This module records the resulting
cheap design filters before a lattice embedding or Vinberg enumeration is
attempted.

These tests are necessary, not sufficient.  In particular, they do not prove
that a proposed list is the complete norm-two chamber.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from coxeter_polytopes.exact.coxeter import classify_root_pair
from coxeter_polytopes.exact.golden import Golden, PHI, RationalLike
from coxeter_polytopes.exact.matrix import Inertia, inertia, rank

from .k3 import k3_trace_square


FieldGram = Sequence[Sequence[Golden | RationalLike]]


def _coerce_norm_two_gram(field_gram: FieldGram) -> list[list[Golden]]:
    gram = [[Golden.coerce(entry) for entry in row] for row in field_gram]
    size = len(gram)
    if any(len(row) != size for row in gram):
        raise ValueError("field wall Gram matrix must be square")
    if any(gram[i][j] != gram[j][i] for i in range(size) for j in range(i)):
        raise ValueError("field wall Gram matrix must be symmetric")
    if any(gram[i][i] != 2 for i in range(size)):
        raise ValueError("K3 inverse-design audit requires norm-two field roots")
    return gram


def trace_curve_gram(field_gram: FieldGram) -> tuple[tuple[Fraction, ...], ...]:
    """Convert a norm-two field wall Gram matrix to K3 intersections.

    For ``alpha=-phi/sqrt(5)``, this applies ``Tr(alpha * -)`` entrywise.
    The diagonal is therefore ``-2``.  In the normalized golden case a
    Coxeter 3-edge has field product ``-1`` and a 5-edge has product
    ``-phi``; both become curve intersection ``+1``.
    """

    gram = _coerce_norm_two_gram(field_gram)
    return tuple(tuple(k3_trace_square(entry) for entry in row) for row in gram)


def correspondence_gram(field_gram: FieldGram) -> tuple[tuple[Fraction, ...], ...]:
    """Return ``q(T r_i, r_j)`` for multiplication ``T=phi``.

    Together with :func:`trace_curve_gram`, this integral symmetric matrix
    retains the golden decoration that ordinary curve intersections forget.
    In particular, norm-two 3- and 5-edges give respectively ``(Q,B)=(1,1)``
    and ``(Q,B)=(1,2)``.
    """

    gram = _coerce_norm_two_gram(field_gram)
    return tuple(
        tuple(k3_trace_square(PHI * entry) for entry in row) for row in gram
    )


def paired_trace_curve_gram(
    field_gram: FieldGram,
) -> tuple[tuple[Fraction, ...], ...]:
    """Return the Gram matrix of both K3 roots above every golden wall.

    Put ``u = phi - 1``.  In the even golden trace construction a norm-two
    field root ``r_i`` and ``u*r_i`` are distinct trace-lattice roots of
    square ``-2``.  They are orthogonal and their hyperplanes restrict to the
    same wall in the distinguished eigenspace.  If

    ``Q_ij = q(r_i, r_j)`` and ``B_ij = q(phi*r_i, r_j)``, then the Gram
    matrix in the ordered list ``(r_1,...,r_n,u*r_1,...,u*r_n)`` is

    ``[[Q, B-Q], [B-Q, 2Q-B]]``.

    This is the actual ``2n``-curve lift of a completed eigen-nef facet list;
    it is different from :func:`~coxeter_polytopes.geometric.k3_presentation.doubled_trace_gram`,
    whose second block consists of the generally non-root vectors
    ``phi*r_i``.
    """

    q = trace_curve_gram(field_gram)
    b = correspondence_gram(field_gram)
    size = len(q)
    result = [[Fraction(0) for _ in range(2 * size)] for _ in range(2 * size)]
    for i in range(size):
        for j in range(size):
            result[i][j] = q[i][j]
            result[i][size + j] = b[i][j] - q[i][j]
            result[size + i][j] = b[i][j] - q[i][j]
            result[size + i][size + j] = 2 * q[i][j] - b[i][j]
    return tuple(tuple(row) for row in result)


def recover_field_gram(
    curve_gram: Sequence[Sequence[int | Fraction]],
    rm_gram: Sequence[Sequence[int | Fraction]],
) -> tuple[tuple[Golden, ...], ...]:
    """Recover ``h(r_i,r_j)`` from ``q`` and ``q(T-, -)``.

    If ``h=a+b*phi``, then for the fixed trace scalar

    ``q=-(a+b)`` and ``q(T-, -)=-(a+2b)``.

    Hence ``a=B-2Q`` and ``b=Q-B``.  Integral input is required because the
    intended inverse-design problem lives in the Neron--Severi lattice.
    """

    q = [[Fraction(entry) for entry in row] for row in curve_gram]
    b = [[Fraction(entry) for entry in row] for row in rm_gram]
    size = len(q)
    if any(len(row) != size for row in q) or len(b) != size or any(
        len(row) != size for row in b
    ):
        raise ValueError("the two divisor Gram matrices must be square and equal-sized")
    if any(q[i][j] != q[j][i] or b[i][j] != b[j][i] for i in range(size) for j in range(i)):
        raise ValueError("the two divisor Gram matrices must be symmetric")
    if any(entry.denominator != 1 for matrix in (q, b) for row in matrix for entry in row):
        raise ValueError("divisor and correspondence pairings must be integral")
    return tuple(
        tuple(Golden(b[i][j] - 2 * q[i][j], q[i][j] - b[i][j]) for j in range(size))
        for i in range(size)
    )


def finite_labels(field_gram: FieldGram) -> tuple[int, ...]:
    """Return the distinct finite Coxeter labels in a simple-root Gram matrix."""

    gram = _coerce_norm_two_gram(field_gram)
    labels: set[int] = set()
    for i in range(len(gram)):
        for j in range(i):
            relation = classify_root_pair(2, gram[i][j], 2)
            if relation.kind == "finite_unrecognized":
                raise ValueError(f"unrecognized finite angle at wall pair ({j}, {i})")
            if relation.kind == "finite":
                assert relation.label is not None
                labels.add(relation.label)
    return tuple(sorted(labels))


@dataclass(frozen=True, slots=True)
class K3DesignAudit:
    """Result of the exact necessary K3 lift tests."""

    wall_count: int
    field_rank: int
    finite_labels: tuple[int, ...]
    trace_off_diagonal_values: tuple[Fraction, ...]
    effective_curve_pairings: bool
    paired_trace_off_diagonal_values: tuple[Fraction, ...]
    paired_effective_curve_pairings: bool
    identity_inertia: Inertia
    conjugate_inertia: Inertia
    vertex_trace_negative_definite: bool

    @property
    def passes(self) -> bool:
        return (
            set(self.finite_labels) <= {2, 3, 5}
            and self.effective_curve_pairings
            and self.paired_effective_curve_pairings
            and self.vertex_trace_negative_definite
        )


def audit_k3_coxeter_design(
    field_gram: FieldGram,
    vertices: Sequence[Sequence[int]],
) -> K3DesignAudit:
    """Audit necessary curve and vertex constraints for a K3 eigen-nef lift.

    ``vertices`` lists the active facets at each finite vertex.  Their trace
    intersection matrices must be negative definite: geometrically, the
    corresponding rational curves are simultaneously orthogonal to a positive
    nef class.  When their cross-intersections are only zero or one, this is
    exactly the ADE condition on the undecorated curve graph.
    """

    gram = _coerce_norm_two_gram(field_gram)
    size = len(gram)
    curve_gram = trace_curve_gram(gram)
    paired_curve_gram = paired_trace_curve_gram(gram)
    off_diagonal = tuple(
        sorted({curve_gram[i][j] for i in range(size) for j in range(i)})
    )
    effective = all(
        value.denominator == 1 and value >= 0
        for i in range(size)
        for value in curve_gram[i][:i]
    )
    paired_off_diagonal = tuple(
        sorted(
            {
                paired_curve_gram[i][j]
                for i in range(2 * size)
                for j in range(i)
            }
        )
    )
    paired_effective = all(
        value.denominator == 1 and value >= 0
        for i in range(2 * size)
        for value in paired_curve_gram[i][:i]
    )

    vertex_negative = True
    for vertex in vertices:
        indices = tuple(vertex)
        if len(indices) != len(set(indices)):
            raise ValueError("a vertex facet list contains a repeated index")
        if any(index < 0 or index >= size for index in indices):
            raise IndexError("a vertex facet index is out of range")
        subgram = [[curve_gram[i][j] for j in indices] for i in indices]
        if inertia(subgram) != Inertia(0, len(indices), 0):
            vertex_negative = False

    return K3DesignAudit(
        wall_count=size,
        field_rank=rank(gram),
        finite_labels=finite_labels(gram),
        trace_off_diagonal_values=off_diagonal,
        effective_curve_pairings=effective,
        paired_trace_off_diagonal_values=paired_off_diagonal,
        paired_effective_curve_pairings=paired_effective,
        identity_inertia=inertia(gram, embedding=0),
        conjugate_inertia=inertia(gram, embedding=1),
        vertex_trace_negative_definite=vertex_negative,
    )
