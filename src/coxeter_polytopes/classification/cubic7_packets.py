"""Small exact discovery slices over ``Q(cos(2*pi/7))``.

The routines here are intentionally bounded and exhaustive in their stated
boxes.  They do two jobs:

* calibrate the cubic backend on equal-root-length simplex packets with
  Coxeter labels ``2, 3, 7``; and
* enumerate the symmetric-base triangular-prism slice in dimension three.

The prism slice is useful as a warning as well as a calibration.  It returns
two arithmetic compact chambers, but both have a visible base-swapping
reflection and are doubles of straight prisms.  They are therefore candidates
only in the literal ``new chamber`` sense, not candidates for a new maximal
reflection group or commensurability class.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations, product
from typing import Iterable, Sequence

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    RHO7,
    Cubic7,
    cubic7_determinant,
    cubic7_inertia,
)
from coxeter_polytopes.exact.matrix import Inertia


Edge = tuple[int, int]
Matrix = tuple[tuple[Cubic7, ...], ...]

EQUAL_NORM_LABELS = (2, 3, 7)


def _edges(rank: int) -> tuple[Edge, ...]:
    return tuple(combinations(range(rank), 2))


def _entry(label: int) -> Cubic7:
    if label == 2:
        return CUBIC7_ZERO
    if label == 3:
        return -CUBIC7_ONE
    if label == 7:
        return -RHO7
    raise ValueError("the equal-norm cubic alphabet contains only labels 2, 3, 7")


def _submatrix(matrix: Matrix, indices: Sequence[int]) -> Matrix:
    return tuple(tuple(matrix[i][j] for j in indices) for i in indices)


def _connected(rank: int, labels: Sequence[int]) -> bool:
    adjacency = [set() for _ in range(rank)]
    for (i, j), label in zip(_edges(rank), labels, strict=True):
        if label != 2:
            adjacency[i].add(j)
            adjacency[j].add(i)
    reached = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for neighbor in adjacency[vertex] - reached:
            reached.add(neighbor)
            frontier.append(neighbor)
    return len(reached) == rank


def equal_norm_gram(rank: int, labels: Sequence[int]) -> Matrix:
    """Return the norm-two symmetric Gram matrix for a complete label tuple."""

    edges = _edges(rank)
    if len(labels) != len(edges):
        raise ValueError(f"rank {rank} needs {len(edges)} edge labels")
    matrix = [[CUBIC7_ZERO for _ in range(rank)] for _ in range(rank)]
    for i in range(rank):
        matrix[i][i] = Cubic7(2)
    for (i, j), label in zip(edges, labels, strict=True):
        matrix[i][j] = matrix[j][i] = _entry(label)
    return tuple(tuple(row) for row in matrix)


def _relabel_complete_graph(rank: int, labels: Sequence[int], old_to_new: Sequence[int]) -> tuple[int, ...]:
    dictionary = {edge: label for edge, label in zip(_edges(rank), labels, strict=True)}
    return tuple(
        dictionary[tuple(sorted((old_to_new[i], old_to_new[j])))]
        for i, j in _edges(rank)
    )


def _canonical_simplex_labels(rank: int, labels: Sequence[int]) -> tuple[int, ...]:
    return min(
        _relabel_complete_graph(rank, labels, permutation)
        for permutation in permutations(range(rank))
    )


@dataclass(frozen=True, slots=True)
class Cubic7SimplexPacket:
    rank: int
    labels: tuple[int, ...]
    gram: Matrix

    @property
    def dimension(self) -> int:
        return self.rank - 1


def _is_compact_arithmetic_simplex(gram: Matrix) -> bool:
    rank = len(gram)
    if cubic7_inertia(gram, embedding=0) != Inertia(rank - 1, 1, 0):
        return False
    for omitted in range(rank):
        indices = tuple(i for i in range(rank) if i != omitted)
        if cubic7_inertia(_submatrix(gram, indices), embedding=0) != Inertia(
            rank - 1, 0, 0
        ):
            return False
    return all(
        cubic7_inertia(gram, embedding=embedding) == Inertia(rank, 0, 0)
        for embedding in (1, 2)
    )


def enumerate_equal_norm_compact_simplices(rank: int) -> tuple[Cubic7SimplexPacket, ...]:
    """Exhaust compact arithmetic simplex packets in the ``{2,3,7}`` box.

    At least one label seven is required, so the returned packet has cubic,
    rather than rational, ground field.  Facet permutations are quotiented.
    """

    if rank < 3 or rank > 5:
        raise ValueError("the bounded simplex calibration supports ranks 3 through 5")
    edge_count = len(_edges(rank))
    representatives: dict[tuple[int, ...], Cubic7SimplexPacket] = {}
    for raw_labels in product(EQUAL_NORM_LABELS, repeat=edge_count):
        if 7 not in raw_labels or not _connected(rank, raw_labels):
            continue
        gram = equal_norm_gram(rank, raw_labels)
        if not _is_compact_arithmetic_simplex(gram):
            continue
        labels = _canonical_simplex_labels(rank, raw_labels)
        representatives.setdefault(
            labels,
            Cubic7SimplexPacket(rank, labels, equal_norm_gram(rank, labels)),
        )
    return tuple(representatives[key] for key in sorted(representatives))


PRISM_EDGES: tuple[Edge, ...] = tuple(
    edge for edge in combinations(range(5), 2) if edge != (0, 1)
)
PRISM_VERTICES: tuple[tuple[int, int, int], ...] = tuple(
    (base, i, j)
    for base in (0, 1)
    for i, j in combinations((2, 3, 4), 2)
)


def symmetric_prism_gram(
    base_labels: Sequence[int],
    lateral_labels: Sequence[int],
    dotted_entry: Cubic7,
) -> Matrix:
    """Build a triangular-prism Gram matrix with interchangeable bases."""

    if len(base_labels) != 3 or len(lateral_labels) != 3:
        raise ValueError("a triangular prism needs three base and three lateral labels")
    matrix = [[CUBIC7_ZERO for _ in range(5)] for _ in range(5)]
    for i in range(5):
        matrix[i][i] = Cubic7(2)
    matrix[0][1] = matrix[1][0] = -dotted_entry
    for offset, lateral in enumerate((2, 3, 4)):
        value = _entry(base_labels[offset])
        matrix[0][lateral] = matrix[lateral][0] = value
        matrix[1][lateral] = matrix[lateral][1] = value
    for (i, j), label in zip(combinations((2, 3, 4), 2), lateral_labels, strict=True):
        matrix[i][j] = matrix[j][i] = _entry(label)
    return tuple(tuple(row) for row in matrix)


def _nontrivial_symmetric_prism_root(
    base_labels: Sequence[int], lateral_labels: Sequence[int]
) -> Cubic7 | None:
    """Solve ``det(G(d))=0`` away from the formal root ``d=-2``.

    Base symmetry forces ``det(G(d))`` to be ``(d+2)`` times a linear
    polynomial.  Three determinant evaluations recover its other root without
    a square-root oracle.
    """

    value_zero = cubic7_determinant(
        symmetric_prism_gram(base_labels, lateral_labels, CUBIC7_ZERO)
    )
    value_one = cubic7_determinant(
        symmetric_prism_gram(base_labels, lateral_labels, CUBIC7_ONE)
    )
    value_minus_one = cubic7_determinant(
        symmetric_prism_gram(base_labels, lateral_labels, -CUBIC7_ONE)
    )
    quadratic = (value_one + value_minus_one) / 2 - value_zero
    if not quadratic:
        return None
    # Product of the two roots is C/A; one root is -2.
    return -value_zero / (2 * quadratic)


def _canonical_prism_labels(
    base_labels: Sequence[int], lateral_labels: Sequence[int]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    lateral_edges = tuple(combinations(range(3), 2))
    lateral_dictionary = {
        edge: label for edge, label in zip(lateral_edges, lateral_labels, strict=True)
    }
    candidates = []
    for permutation in permutations(range(3)):
        new_base = tuple(base_labels[permutation[i]] for i in range(3))
        new_lateral = tuple(
            lateral_dictionary[tuple(sorted((permutation[i], permutation[j])))]
            for i, j in lateral_edges
        )
        candidates.append((new_base, new_lateral))
    return min(candidates)


def _solve_linear(matrix: Matrix, right: Sequence[Cubic7]) -> tuple[Cubic7, ...]:
    size = len(matrix)
    work = [[*matrix[i], Cubic7.coerce(right[i])] for i in range(size)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            raise ValueError("linear system is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                work[row][j] - multiplier * work[column][j]
                for j in range(size + 1)
            ]
    return tuple(work[i][-1] for i in range(size))


def _primitive_radical_vector(gram: Matrix) -> tuple[Cubic7, ...] | None:
    size = len(gram)
    for free in range(size):
        indices = tuple(i for i in range(size) if i != free)
        coefficient_matrix = _submatrix(gram, indices)
        if not cubic7_determinant(coefficient_matrix):
            continue
        solution = _solve_linear(
            coefficient_matrix, tuple(-gram[i][free] for i in indices)
        )
        vector = [CUBIC7_ZERO for _ in range(size)]
        vector[free] = CUBIC7_ONE
        for index, value in zip(indices, solution, strict=True):
            vector[index] = value
        candidate = tuple(vector)
        if not all(value.is_integral for value in candidate):
            continue
        if all(
            sum((gram[i][j] * candidate[j] for j in range(size)), start=CUBIC7_ZERO)
            == 0
            for i in range(size)
        ):
            return candidate
    return None


@dataclass(frozen=True, slots=True)
class Cubic7SymmetricPrismPacket:
    base_labels: tuple[int, int, int]
    lateral_labels: tuple[int, int, int]
    dotted_entry: Cubic7
    gram: Matrix
    radical_vector: tuple[Cubic7, ...]

    @property
    def cosh_base_distance(self) -> Cubic7:
        return self.dotted_entry / 2

    @property
    def base_swap_root(self) -> tuple[Cubic7, ...]:
        """The integral root whose reflection interchanges the two bases."""

        return (CUBIC7_ONE, -CUBIC7_ONE, CUBIC7_ZERO, CUBIC7_ZERO, CUBIC7_ZERO)


def _is_compact_arithmetic_prism(gram: Matrix) -> bool:
    if cubic7_inertia(gram, embedding=0) != Inertia(3, 1, 1):
        return False
    if any(
        cubic7_inertia(_submatrix(gram, vertex), embedding=0) != Inertia(3, 0, 0)
        for vertex in PRISM_VERTICES
    ):
        return False
    return all(
        cubic7_inertia(gram, embedding=embedding) == Inertia(4, 0, 1)
        for embedding in (1, 2)
    )


def enumerate_symmetric_compact_prisms() -> tuple[Cubic7SymmetricPrismPacket, ...]:
    """Exhaust the symmetric-base equal-norm cubic triangular-prism box."""

    representatives: dict[
        tuple[tuple[int, ...], tuple[int, ...]], Cubic7SymmetricPrismPacket
    ] = {}
    for base_labels in product(EQUAL_NORM_LABELS, repeat=3):
        for lateral_labels in product(EQUAL_NORM_LABELS, repeat=3):
            if 7 not in (*base_labels, *lateral_labels):
                continue
            dotted = _nontrivial_symmetric_prism_root(base_labels, lateral_labels)
            if dotted is None or not dotted.is_integral or (dotted - 2).sign(0) <= 0:
                continue
            gram = symmetric_prism_gram(base_labels, lateral_labels, dotted)
            if not _is_compact_arithmetic_prism(gram):
                continue
            radical = _primitive_radical_vector(gram)
            if radical is None:
                continue
            canonical = _canonical_prism_labels(base_labels, lateral_labels)
            canonical_gram = symmetric_prism_gram(*canonical, dotted)
            canonical_radical = _primitive_radical_vector(canonical_gram)
            assert canonical_radical is not None
            representatives.setdefault(
                canonical,
                Cubic7SymmetricPrismPacket(
                    canonical[0],
                    canonical[1],
                    dotted,
                    canonical_gram,
                    canonical_radical,
                ),
            )
    return tuple(representatives[key] for key in sorted(representatives))


__all__ = [
    "Cubic7SimplexPacket",
    "Cubic7SymmetricPrismPacket",
    "EQUAL_NORM_LABELS",
    "PRISM_EDGES",
    "PRISM_VERTICES",
    "enumerate_equal_norm_compact_simplices",
    "enumerate_symmetric_compact_prisms",
    "equal_norm_gram",
    "symmetric_prism_gram",
]
