"""Exact recognition primitives for Coxeter subdiagrams.

The functions here use only principal Gram matrices.  Connectivity is read
from nonzero off-diagonal entries, so the routines work with roots of unequal
length and do not require a normalized Coxeter matrix.
"""

from __future__ import annotations

from collections.abc import Sequence

from .golden import Golden
from .matrix import Inertia, Matrix, inertia


def _coxeter_gram(gram: Matrix) -> list[list[Golden]]:
    """Coerce a square symmetric matrix, rejecting malformed input."""

    rows = [[Golden.coerce(entry) for entry in row] for row in gram]
    if any(len(row) != len(rows) for row in rows):
        raise ValueError("a Coxeter Gram matrix must be square")
    if any(rows[i][j] != rows[j][i] for i in range(len(rows)) for j in range(i)):
        raise ValueError("a Coxeter Gram matrix must be symmetric")
    return rows


def _has_coxeter_signs(gram: Sequence[Sequence[Golden]]) -> bool:
    """Whether the matrix has positive diagonal and nonpositive cross terms."""

    return all(gram[i][i].sign() > 0 for i in range(len(gram))) and all(
        gram[i][j].sign() <= 0
        for i in range(len(gram))
        for j in range(i)
    )


def induced_subgram(gram: Matrix, vertices: Sequence[int]) -> list[list[Golden]]:
    """Return the principal Gram matrix induced by ``vertices``.

    The requested order is preserved.  Repeated vertices are rejected because
    an induced subdiagram is a set of vertices, not a multiset.
    """

    matrix = _coxeter_gram(gram)
    indices = list(vertices)
    if not all(isinstance(index, int) for index in indices):
        raise TypeError("subdiagram vertices must be integers")
    if len(indices) != len(set(indices)):
        raise ValueError("subdiagram vertices must be distinct")
    if any(index < 0 or index >= len(matrix) for index in indices):
        raise IndexError("subdiagram vertex is out of range")
    return [[matrix[i][j] for j in indices] for i in indices]


def connected_components(gram: Matrix) -> tuple[tuple[int, ...], ...]:
    """Return connected components of a Coxeter Gram matrix.

    Two vertices are adjacent exactly when their off-diagonal pairing is
    nonzero.  Components and their vertices are returned in increasing order.
    """

    matrix = _coxeter_gram(gram)
    unseen = set(range(len(matrix)))
    components: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            neighbours = {
                other
                for other in unseen
                if matrix[vertex][other]
            }
            unseen.difference_update(neighbours)
            component.update(neighbours)
            frontier.extend(neighbours)
        components.append(tuple(sorted(component)))
    return tuple(components)


def is_elliptic(gram: Matrix) -> bool:
    """Return whether ``gram`` is an elliptic Coxeter diagram.

    The empty diagram is elliptic of rank zero, as is customary for induced
    subdiagram arguments.
    """

    matrix = _coxeter_gram(gram)
    return _has_coxeter_signs(matrix) and inertia(matrix) == Inertia(len(matrix), 0, 0)


def is_affine_parabolic_component(gram: Matrix) -> bool:
    """Return whether ``gram`` is one connected affine Coxeter component.

    Besides exact positive-semidefinite corank one inertia, all one-vertex
    deletions are required to be elliptic.  This is the intrinsic Gram-matrix
    characterization and excludes merely singular or wrongly signed matrices.
    """

    matrix = _coxeter_gram(gram)
    size = len(matrix)
    if size < 2 or not _has_coxeter_signs(matrix):
        return False
    if len(connected_components(matrix)) != 1:
        return False
    if inertia(matrix) != Inertia(size - 1, 0, 1):
        return False
    return all(
        is_elliptic(induced_subgram(matrix, [j for j in range(size) if j != i]))
        for i in range(size)
    )


def is_parabolic(gram: Matrix) -> bool:
    """Return whether every component is affine parabolic.

    Thus a reducible parabolic diagram is a disjoint union of connected affine
    diagrams.  The empty diagram is not parabolic.
    """

    matrix = _coxeter_gram(gram)
    components = connected_components(matrix)
    return bool(components) and all(
        is_affine_parabolic_component(induced_subgram(matrix, component))
        for component in components
    )


def parabolic_rank(gram: Matrix) -> int:
    """Return the Coxeter rank of a componentwise affine diagram.

    A connected affine component on ``n`` vertices has rank ``n-1``.  The
    result therefore equals both the sum of component ranks and the exact
    matrix rank.  Non-parabolic input is rejected instead of assigning it a
    misleading rank.
    """

    matrix = _coxeter_gram(gram)
    components = connected_components(matrix)
    if not components or any(
        not is_affine_parabolic_component(induced_subgram(matrix, component))
        for component in components
    ):
        raise ValueError("parabolic rank is defined only for a parabolic diagram")
    return sum(len(component) - 1 for component in components)
