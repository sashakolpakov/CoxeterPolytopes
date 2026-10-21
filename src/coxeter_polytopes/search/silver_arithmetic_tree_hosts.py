"""Exact rank-23 silver tree-host census around finite Dynkin skeletons.

The field-generating edge has Cartan product ``2 + sqrt(2)``.  It is the
finite Coxeter edge of label eight (not a geometric dotted edge).  A second,
separate calculation records the hitting-set obstruction created by two or
more genuine dotted edges (Cartan product greater than four).
"""

from __future__ import annotations

from dataclasses import dataclass

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import (
    SILVER_ONE,
    SILVER_ZERO,
    SQRT2,
    Silver,
    silver_determinant,
    silver_inertia,
)


ARITHMETIC_EDGE_PRODUCT = 2 + SQRT2
GENTLE_DOTTED_PRODUCT = 3 + 2 * SQRT2


SilverGram = tuple[tuple[Silver, ...], ...]


def _zero_matrix(rank: int) -> list[list[Silver]]:
    return [[SILVER_ZERO for _ in range(rank)] for _ in range(rank)]


def finite_skeleton(kind: str, rank: int) -> tuple[SilverGram, tuple[Silver, ...]]:
    """Return an integral symmetric Gram and its Cartan symmetrizer.

    ``A`` is the ordinary path.  ``D`` is numbered with branching node zero,
    short leaves one and two, and long arm ``0-3-...-(rank-1)``.  ``B`` is
    the path ``0-1-...-(rank-1)``, with product two on ``0-1`` and
    symmetrizer ``(1,2,...,2)``.
    """

    if kind not in {"A", "B", "D"}:
        raise ValueError("kind must be A, B, or D")
    if rank < (4 if kind == "D" else 2):
        raise ValueError("the requested Dynkin skeleton is too small")

    gram = _zero_matrix(rank)
    if kind == "B":
        symmetrizer = (SILVER_ONE, *(Silver(2) for _ in range(rank - 1)))
    else:
        symmetrizer = tuple(SILVER_ONE for _ in range(rank))
    for i, scale in enumerate(symmetrizer):
        gram[i][i] = 2 * scale

    if kind in {"A", "B"}:
        edges = [(i, i + 1) for i in range(rank - 1)]
    else:
        edges = [(0, 1), (0, 2), (0, 3), *[(i, i + 1) for i in range(3, rank - 1)]]
    for i, j in edges:
        cross = Silver(-2) if kind == "B" and (i, j) == (0, 1) else -symmetrizer[i]
        # Along a simply-laced B arm both endpoints have symmetrizer two.
        if kind == "B" and i >= 1:
            cross = Silver(-2)
        gram[i][j] = gram[j][i] = cross
    return tuple(tuple(row) for row in gram), symmetrizer


def attach_arithmetic_leaf(kind: str, rank: int, attachment: int) -> SilverGram:
    """Attach one leaf by Cartan pair ``(-c,-1)``, ``c=2+sqrt(2)``."""

    skeleton, symmetrizer = finite_skeleton(kind, rank)
    if not 0 <= attachment < rank:
        raise IndexError("attachment vertex is outside the finite skeleton")
    total = rank + 1
    gram = _zero_matrix(total)
    for i in range(rank):
        for j in range(rank):
            gram[i][j] = skeleton[i][j]
    scale = symmetrizer[attachment] * ARITHMETIC_EDGE_PRODUCT
    gram[rank][rank] = 2 * scale
    gram[rank][attachment] = gram[attachment][rank] = -scale
    return tuple(tuple(row) for row in gram)


def _minor(matrix: SilverGram, omitted: int) -> SilverGram:
    keep = tuple(i for i in range(len(matrix)) if i != omitted)
    return tuple(tuple(matrix[i][j] for j in keep) for i in keep)


def inverse_diagonal(matrix: SilverGram) -> tuple[Silver, ...]:
    """Compute the exact inverse diagonal by principal cofactors."""

    determinant = silver_determinant(matrix)
    if not determinant:
        raise ValueError("inverse diagonal is undefined for a singular Gram")
    return tuple(silver_determinant(_minor(matrix, i)) / determinant for i in range(len(matrix)))


def ray_kind(value: Silver) -> str:
    """Classify a simplex dual ray from its exact selected-place norm."""

    sign = value.sign(0)
    return "finite" if sign < 0 else "ideal" if sign == 0 else "hyperideal"


@dataclass(frozen=True, slots=True)
class TreeHostRow:
    kind: str
    attachment: int
    selected_inertia: Inertia
    conjugate_inertia: Inertia
    inverse_diagonal: tuple[Silver, ...]
    ray_kinds: tuple[str, ...]

    @property
    def finite_count(self) -> int:
        return self.ray_kinds.count("finite")

    @property
    def ideal_count(self) -> int:
        return self.ray_kinds.count("ideal")

    @property
    def hyperideal_count(self) -> int:
        return self.ray_kinds.count("hyperideal")

    @property
    def arithmetic(self) -> bool:
        rank = len(self.inverse_diagonal)
        return (
            self.selected_inertia == Inertia(rank - 1, 1, 0)
            and self.conjugate_inertia == Inertia(rank, 0, 0)
        )


def analyze_tree_host(kind: str, rank: int, attachment: int) -> TreeHostRow:
    gram = attach_arithmetic_leaf(kind, rank, attachment)
    diagonal = inverse_diagonal(gram)
    return TreeHostRow(
        kind=kind,
        attachment=attachment,
        selected_inertia=silver_inertia(gram, embedding=0),
        conjugate_inertia=silver_inertia(gram, embedding=1),
        inverse_diagonal=diagonal,
        ray_kinds=tuple(ray_kind(value) for value in diagonal),
    )


def rank23_tree_host_census() -> tuple[TreeHostRow, ...]:
    """Enumerate every attachment orbit representative before quotienting.

    Returning all labelled attachment vertices keeps the exact inverse data
    transparent.  Symmetric vertices naturally produce identical rows.
    """

    return tuple(
        analyze_tree_host(kind, 22, attachment)
        for kind in ("A", "D", "B")
        for attachment in range(22)
    )


def genuine_dotted_hyperideal_lower_bound(
    rank: int, dotted_edges: tuple[tuple[int, int], ...]
) -> int:
    """Lower-bound simplex hyperideal rays using the dotted-edge cover.

    A vertex ray is classified by deleting one wall.  If the remaining
    principal Gram contains a genuine dotted edge, its two-wall principal
    block is indefinite, so that deletion is hyperideal.  Consequently only
    vertices incident to every dotted edge can possibly be finite or ideal.
    """

    if rank < 3:
        raise ValueError("a Coxeter simplex here has rank at least three")
    if len(dotted_edges) < 2:
        raise ValueError("the obstruction is for at least two dotted edges")
    normalized = tuple(tuple(sorted(edge)) for edge in dotted_edges)
    if any(len(set(edge)) != 2 or min(edge) < 0 or max(edge) >= rank for edge in normalized):
        raise ValueError("a dotted edge has invalid endpoints")
    common = set(normalized[0])
    for edge in normalized[1:]:
        common.intersection_update(edge)
    return rank - len(common)


def adjacent_two_dotted_rank23_host() -> SilverGram:
    """A concrete arithmetic rank-23 tree with two adjacent dotted edges.

    Node zero meets dotted leaves one and two, both of product
    ``p=3+2*sqrt(2)``, and the ordinary label-three arm
    ``0-3-4-...-22``.  Orienting each dotted Cartan pair as ``(-p,-1)``
    gives symmetrizer ``(1,p,p,1,...,1)`` and this integral symmetric Gram.
    """

    rank = 23
    gram = _zero_matrix(rank)
    symmetrizer = (SILVER_ONE, GENTLE_DOTTED_PRODUCT, GENTLE_DOTTED_PRODUCT,
                   *(SILVER_ONE for _ in range(rank - 3)))
    for i, scale in enumerate(symmetrizer):
        gram[i][i] = 2 * scale
    for leaf in (1, 2):
        gram[0][leaf] = gram[leaf][0] = -GENTLE_DOTTED_PRODUCT
    gram[0][3] = gram[3][0] = Silver(-1)
    for i in range(3, rank - 1):
        gram[i][i + 1] = gram[i + 1][i] = Silver(-1)
    return tuple(tuple(row) for row in gram)


__all__ = [
    "ARITHMETIC_EDGE_PRODUCT",
    "GENTLE_DOTTED_PRODUCT",
    "TreeHostRow",
    "analyze_tree_host",
    "adjacent_two_dotted_rank23_host",
    "attach_arithmetic_leaf",
    "finite_skeleton",
    "genuine_dotted_hyperideal_lower_bound",
    "inverse_diagonal",
    "rank23_tree_host_census",
    "ray_kind",
]
