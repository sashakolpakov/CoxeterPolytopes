"""Exact bounded root search on a cubic-field diagonal Lorentz lattice.

Let ``K=Q(theta)``, ``theta=2*cos(2*pi/7)``, and

``L_n = <-theta> + <1>^n``.

At the distinguished place this has signature ``(n,1)``; at the other two
places it is positive definite.  The lattice is unimodular.  Moreover, ``2``
is inert in ``K`` and the narrow class number is one.  Consequently every
primitive reflection hyperplane has a representative of norm exactly ``1``
or ``2``.  (A primitive root norm divides ``2``; every totally positive unit
is a square.)  This turns the first Vinberg stage into a genuinely exhaustive
finite calculation, rather than a heuristic small-norm scan.

This module enumerates the complete pool through a prescribed Vinberg height.
Selecting irredundant simple walls and proving finite volume are deliberately
separate steps: a later mirror may cut off a wall from an earlier non-Coxeter
dihedral sector.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import isqrt
from typing import Iterable, Sequence

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    THETA7,
    Cubic7,
    cubic7_inertia,
)
from coxeter_polytopes.exact.matrix import Inertia


Vector = tuple[Cubic7, ...]


def diagonal_pair(left: Sequence[Cubic7], right: Sequence[Cubic7]) -> Cubic7:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("diagonal cubic vectors must have the same rank at least two")
    return -THETA7 * left[0] * right[0] + sum(
        (left[i] * right[i] for i in range(1, len(left))),
        start=CUBIC7_ZERO,
    )


def _integer_ceiling_square_root(value: Fraction) -> int:
    if value < 0:
        raise ValueError("a square cap cannot be negative")
    floor = isqrt(value.numerator // value.denominator)
    while Fraction(floor * floor) < value:
        floor += 1
    return floor


def _integer_ceiling(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


# Certified rational upper bounds for the absolute values of the three
# trace-dual basis elements, with columns indexed by the real embeddings.
# Each of the eighteen inequalities ``-bound < sigma(d_i) < bound`` can be
# replayed using Cubic7.sign; the deliberately slack bounds keep enumeration
# boxes tiny without introducing a numerical oracle.
_DUAL_ABSOLUTE_UPPER = (
    (Fraction(1, 6), Fraction(1), Fraction(1, 7)),
    (Fraction(4, 9), Fraction(1, 4), Fraction(1, 5)),
    (Fraction(1, 5), Fraction(4, 9), Fraction(1, 4)),
)


def bounded_cubic7_integers(square_caps: Sequence[Fraction | int]) -> tuple[Cubic7, ...]:
    """Enumerate all ``z in O_K`` with ``sigma_i(z)^2 <= square_caps[i]``.

    The coefficient box is certified by the trace-dual basis.  For
    ``(1,theta,theta^2)`` that basis is

    ``1-theta/7-3theta^2/7``, ``-1/7+2theta/7+theta^2/7``,
    ``-3/7+theta/7+2theta^2/7``.

    Rational upper bounds for the absolute values of these dual elements are
    applied place by place.  Exact sign tests then apply the requested caps;
    no floating-point lattice-point test occurs.
    """

    if len(square_caps) != 3:
        raise ValueError("the cubic field has exactly three real embeddings")
    caps = tuple(Fraction(cap) for cap in square_caps)
    if any(cap < 0 for cap in caps):
        return ()
    embedding_bounds = tuple(_integer_ceiling_square_root(cap) for cap in caps)
    coefficient_bounds = tuple(
        _integer_ceiling(
            sum(
                (
                    Fraction(embedding_bounds[embedding])
                    * _DUAL_ABSOLUTE_UPPER[coefficient][embedding]
                    for embedding in range(3)
                ),
                start=Fraction(0),
            )
        )
        for coefficient in range(3)
    )
    values = []
    for a in range(-coefficient_bounds[0], coefficient_bounds[0] + 1):
        for b in range(-coefficient_bounds[1], coefficient_bounds[1] + 1):
            for c in range(-coefficient_bounds[2], coefficient_bounds[2] + 1):
                value = Cubic7(a, b, c)
                square = value * value
                if all((Cubic7(caps[i]) - square).sign(i) >= 0 for i in range(3)):
                    values.append(value)
    return tuple(sorted(values))


@dataclass(frozen=True, slots=True)
class Cubic7DiagonalRoot:
    norm: int
    height: Cubic7
    vector: Vector

    @property
    def control(self) -> Cubic7:
        return self.vector[0]


@dataclass(frozen=True, slots=True)
class Cubic7DiagonalRootPool:
    dimension: int
    height_bound: int
    cone_roots: tuple[Vector, ...]
    roots: tuple[Cubic7DiagonalRoot, ...]
    control_counts_by_norm: tuple[tuple[int, int], ...]
    spatial_alphabet_counts_by_norm: tuple[tuple[int, int], ...]

    @property
    def all_vectors(self) -> tuple[Vector, ...]:
        return (*self.cone_roots, *(root.vector for root in self.roots))


def _oriented_cone_roots(dimension: int) -> tuple[Vector, ...]:
    # The height-zero roots form B_n, not A1^n.  This chamber is
    # y_1 >= y_2 >= ... >= y_n >= 0.
    roots = []
    for index in range(dimension - 1):
        root = [CUBIC7_ZERO for _ in range(dimension + 1)]
        root[index + 1] = -CUBIC7_ONE
        root[index + 2] = CUBIC7_ONE
        roots.append(tuple(root))
    final_root = [CUBIC7_ZERO for _ in range(dimension + 1)]
    final_root[-1] = -CUBIC7_ONE
    roots.append(tuple(final_root))
    return tuple(roots)


def _sum_of_squares_solutions(
    alphabet: Sequence[Cubic7], dimension: int, target: Cubic7
) -> tuple[tuple[Cubic7, ...], ...]:
    by_square: dict[Cubic7, list[Cubic7]] = {}
    for value in alphabet:
        by_square.setdefault(value * value, []).append(value)

    solutions: list[tuple[Cubic7, ...]] = []

    def extend(prefix: tuple[Cubic7, ...], remainder: Cubic7) -> None:
        remaining = dimension - len(prefix)
        if remaining == 1:
            for value in by_square.get(remainder, ()):
                solutions.append((*prefix, value))
            return
        for value in alphabet:
            next_remainder = remainder - value * value
            if all(next_remainder.sign(embedding) >= 0 for embedding in range(3)):
                extend((*prefix, value), next_remainder)

    extend((), target)
    return tuple(solutions)


def enumerate_diagonal_root_pool(
    dimension: int = 3, height_bound: int = 10
) -> Cubic7DiagonalRootPool:
    """Exhaust norm-one/two walls through ``delta <= height_bound``.

    Roots are oriented by ``x>0`` and by the coordinate chamber
    ``y_i>=0`` at the distinguished embedding.  The Vinberg height is

    ``delta(r) = theta*x^2 / norm(r)``.
    """

    if dimension < 2:
        raise ValueError("hyperbolic dimension must be at least two")
    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")

    roots: set[Cubic7DiagonalRoot] = set()
    control_counts: list[tuple[int, int]] = []
    spatial_counts: list[tuple[int, int]] = []
    for norm in (1, 2):
        # At place zero theta>1, so x^2 <= norm*height_bound.  At the
        # middle place -theta>2/5, hence x^2<3*norm; at the last place
        # -theta>1, hence x^2<=norm.  The exact target test follows.
        controls = bounded_cubic7_integers(
            (norm * height_bound, 3 * norm, norm)
        )
        controls = tuple(
            x
            for x in controls
            if x.sign(0) > 0
            and (Cubic7(norm * height_bound) - THETA7 * x * x).sign(0) >= 0
            and all((Cubic7(norm) + THETA7 * x * x).sign(i) >= 0 for i in (1, 2))
        )
        control_counts.append((norm, len(controls)))

        # Globally valid caps for any spatial coordinate at this height.
        spatial_pool = bounded_cubic7_integers(
            (norm * (height_bound + 1), norm, norm)
        )
        spatial_counts.append((norm, len(spatial_pool)))

        for control in controls:
            control_square = control * control
            target = Cubic7(norm) + THETA7 * control_square
            alphabet = tuple(
                value
                for value in spatial_pool
                if value.sign(0) >= 0
                and all((target - value * value).sign(i) >= 0 for i in range(3))
            )
            for spatial in _sum_of_squares_solutions(alphabet, dimension, target):
                if any(
                    (spatial[i] - spatial[i + 1]).sign(0) < 0
                    for i in range(dimension - 1)
                ):
                    continue
                vector = (control, *spatial)
                if diagonal_pair(vector, vector) != norm:
                    raise AssertionError("sum-of-squares enumerator returned a wrong norm")
                height = THETA7 * control_square / norm
                roots.add(Cubic7DiagonalRoot(norm, height, vector))

    ordered = tuple(
        sorted(
            roots,
            key=lambda root: (
                root.height,
                root.norm,
                root.vector,
            ),
        )
    )
    return Cubic7DiagonalRootPool(
        dimension=dimension,
        height_bound=height_bound,
        cone_roots=_oriented_cone_roots(dimension),
        roots=ordered,
        control_counts_by_norm=tuple(control_counts),
        spatial_alphabet_counts_by_norm=tuple(spatial_counts),
    )


def _wall_gram(walls: Sequence[Vector]) -> tuple[tuple[Cubic7, ...], ...]:
    return tuple(tuple(diagonal_pair(left, right) for right in walls) for left in walls)


def _submatrix(
    matrix: Sequence[Sequence[Cubic7]], indices: Sequence[int]
) -> tuple[tuple[Cubic7, ...], ...]:
    return tuple(tuple(matrix[i][j] for j in indices) for i in indices)


@dataclass(frozen=True, slots=True)
class Cubic7DiagonalVinbergSearch:
    pool: Cubic7DiagonalRootPool
    accepted_roots: tuple[Cubic7DiagonalRoot, ...]
    rejected_roots: tuple[Cubic7DiagonalRoot, ...]
    walls: tuple[Vector, ...]
    wall_gram: tuple[tuple[Cubic7, ...], ...]
    compact_prism_vertices: tuple[tuple[int, int, int], ...]
    compact: bool


def run_diagonal_vinberg_search(
    dimension: int = 3, height_bound: int = 7
) -> Cubic7DiagonalVinbergSearch:
    """Run the exact acute-wall selection on the bounded root pool.

    In dimension three, bound seven reaches a compact chamber.  The returned
    compactness flag is proved by all six elliptic vertex links of its
    triangular-prism combinatorics, so Vinberg stopping makes the prefix the
    complete reflection chamber.
    """

    pool = enumerate_diagonal_root_pool(dimension, height_bound)
    walls = list(pool.cone_roots)
    accepted: list[Cubic7DiagonalRoot] = []
    rejected: list[Cubic7DiagonalRoot] = []
    for candidate in pool.roots:
        if all(diagonal_pair(candidate.vector, wall).sign(0) <= 0 for wall in walls):
            walls.append(candidate.vector)
            accepted.append(candidate)
        else:
            rejected.append(candidate)
    if any(
        diagonal_pair(left, right).sign(0) > 0
        for i, left in enumerate(walls)
        for right in walls[i + 1 :]
    ):
        raise AssertionError("the accepted Vinberg walls are not pairwise acute")

    gram = _wall_gram(walls)
    compact_vertices: tuple[tuple[int, int, int], ...] = ()
    compact = False
    if dimension == 3 and len(walls) == 5:
        # In the resulting order, cone wall 2 and the last positive-height
        # wall are the unique divergent pair.  The other three are lateral
        # facets of a triangular prism.
        bases = (2, 4)
        laterals = (0, 1, 3)
        compact_vertices = tuple(
            (base, *pair)
            for base in bases
            for pair in combinations(laterals, 2)
        )
        compact = (
            cubic7_inertia(gram, embedding=0) == Inertia(3, 1, 1)
            and all(
                cubic7_inertia(gram, embedding=embedding) == Inertia(4, 0, 1)
                for embedding in (1, 2)
            )
            and all(
                cubic7_inertia(_submatrix(gram, vertex), embedding=0)
                == Inertia(3, 0, 0)
                for vertex in compact_vertices
            )
            and (
                4 * gram[bases[0]][bases[1]] ** 2
                - gram[bases[0]][bases[0]] * gram[bases[1]][bases[1]] * 4
            ).sign(0)
            > 0
        )
    return Cubic7DiagonalVinbergSearch(
        pool=pool,
        accepted_roots=tuple(accepted),
        rejected_roots=tuple(rejected),
        walls=tuple(walls),
        wall_gram=gram,
        compact_prism_vertices=compact_vertices,
        compact=compact,
    )


__all__ = [
    "Cubic7DiagonalRoot",
    "Cubic7DiagonalRootPool",
    "Cubic7DiagonalVinbergSearch",
    "bounded_cubic7_integers",
    "diagonal_pair",
    "enumerate_diagonal_root_pool",
    "run_diagonal_vinberg_search",
]
