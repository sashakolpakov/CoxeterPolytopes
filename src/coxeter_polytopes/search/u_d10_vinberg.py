"""Exact Vinberg search on ``U + D10(-1)``.

This even rank-twelve Lorentz lattice has determinant four and is
Scharlau--Kirschmer ``List11/Gram11`` entry 6,
``II_{11,1}([2^2]_2)``.  It is the fixed lattice of the Roulleau
``U+E8(-1)+A3(-1)`` terminal-diamond involution.

For ``r=(a,b;v)``, control vector ``p=e+f``, ``m=-(a+b)>0`` and
``d=a-b``, a root of square ``-q`` satisfies

    d^2 + 2*v_D10^2 = m^2 + 2*q,   d == m (mod 2).

The height-zero reflection system is ``A1+B10``.  The discriminant group has
exponent two, so primitive reflecting roots have norm ``-2`` or ``-4``;
primitivity and reflection integrality are checked directly nonetheless.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt
from typing import Sequence

import sympy as sp

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.verify.chamber import ChamberVerification, diagnose_chamber_cone


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


@lru_cache(maxsize=1)
def d10_cartan() -> Matrix:
    matrix = [[2 if left == right else 0 for right in range(10)] for left in range(10)]
    for left, right in (*((index, index + 1) for index in range(7)), (7, 8), (7, 9)):
        matrix[left][right] = matrix[right][left] = -1
    return tuple(tuple(row) for row in matrix)


@lru_cache(maxsize=1)
def lattice_gram() -> Matrix:
    gram = [[0 for _ in range(12)] for _ in range(12)]
    gram[0][1] = gram[1][0] = 1
    cartan = d10_cartan()
    for left in range(10):
        for right in range(10):
            gram[left + 2][right + 2] = -cartan[left][right]
    return tuple(tuple(row) for row in gram)


def _mat_vec(matrix: Sequence[Sequence[int]], vector: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        sum(row[column] * vector[column] for column in range(len(vector)))
        for row in matrix
    )


def pair(left: Sequence[int], right: Sequence[int]) -> int:
    gram = lattice_gram()
    return sum(
        left[row] * gram[row][column] * right[column]
        for row in range(12)
        for column in range(12)
    )


def divisibility(root: Sequence[int]) -> int:
    result = 0
    for value in _mat_vec(lattice_gram(), root):
        result = gcd(result, abs(value))
    return result


def _primitive(root: Sequence[int]) -> bool:
    result = 0
    for value in root:
        result = gcd(result, abs(value))
    return result == 1


@lru_cache(maxsize=1)
def height_zero_simple_roots() -> tuple[Vector, ...]:
    roots: list[Vector] = [(1, -1, *(0,) * 10)]
    roots.extend(
        (0, 0, *tuple(-1 if index == simple else 0 for index in range(10)))
        for simple in range(9)
    )
    roots.append((0, 0, *(0,) * 8, 1, -1))
    result = tuple(roots)
    if len(result) != 11:
        raise AssertionError("wrong height-zero rank")
    if any(
        pair(left, right) < 0
        for position, left in enumerate(result)
        for right in result[position + 1 :]
    ):
        raise AssertionError("height-zero roots are not pairwise acute")
    if {(pair(root, root), divisibility(root)) for root in result} != {
        (-2, 1),
        (-4, 2),
    }:
        raise AssertionError("unexpected height-zero root types")
    return result


@dataclass(frozen=True, slots=True)
class DominantD10Vector:
    coordinates: tuple[int, ...]
    weights: tuple[int, ...]
    square: int


@lru_cache(maxsize=None)
def dominant_d10_vectors(max_square: int) -> tuple[DominantD10Vector, ...]:
    """Enumerate the ``B10``-dominant vectors of the ``D10`` root lattice.

    In standard orthonormal coordinates the chosen height-zero chamber is
    ``x1 >= ... >= x10 >= 0``.  The ``D10`` integrality condition is that the
    coordinate sum is even.  Enumerating square partitions is dramatically
    smaller than a box of Dynkin weights and remains exact.
    """

    if max_square < 0:
        return ()
    simple_columns = sp.zeros(10, 10)
    for index in range(9):
        simple_columns[index, index] = 1
        simple_columns[index + 1, index] = -1
    simple_columns[8, 9] = simple_columns[9, 9] = 1
    inverse_simple = simple_columns.inv()

    def partitions(total: int, length: int, maximum: int):
        if length == 0:
            if total == 0:
                yield ()
            return
        upper = min(maximum, isqrt(total))
        for value in range(upper, -1, -1):
            remainder = total - value * value
            if remainder < 0 or remainder > (length - 1) * value * value:
                continue
            for suffix in partitions(remainder, length - 1, value):
                yield (value, *suffix)

    vectors: list[DominantD10Vector] = []
    for square in range(max_square + 1):
        for standard in partitions(square, 10, isqrt(square)):
            if sum(standard) % 2:
                continue
            standard_vector = sp.Matrix(standard)
            coordinates = inverse_simple * standard_vector
            if any(value.q != 1 for value in coordinates):
                raise ArithmeticError("even-sum D10 vector left the root lattice")
            weights = simple_columns.T * standard_vector
            vectors.append(
                DominantD10Vector(
                    coordinates=tuple(int(value) for value in coordinates),
                    weights=tuple(int(value) for value in weights),
                    square=square,
                )
            )
    return tuple(sorted(vectors, key=lambda vector: vector.coordinates))


@dataclass(frozen=True, slots=True)
class Candidate:
    root: Vector
    control: int
    norm_abs: int
    root_divisibility: int
    difference: int
    d10_square: int

    @property
    def height(self) -> Fraction:
        return Fraction(self.control * self.control, self.norm_abs)


@lru_cache(maxsize=None)
def enumerate_level_candidates(control: int, norm_abs: int) -> tuple[Candidate, ...]:
    if control <= 0:
        raise ValueError("control must be positive")
    if norm_abs not in (2, 4):
        raise ValueError("supported root norms are -2 and -4")
    radius = control * control + 2 * norm_abs
    dominant: dict[int, list[DominantD10Vector]] = {}
    for vector in dominant_d10_vectors(radius // 2):
        if vector.weights[9] >= vector.weights[8]:
            dominant.setdefault(vector.square, []).append(vector)
    candidates: list[Candidate] = []
    for difference in range(-isqrt(radius), 1):
        if (difference - control) % 2:
            continue
        remainder = radius - difference * difference
        if remainder < 0 or remainder % 2:
            continue
        d10_square = remainder // 2
        a = (-control + difference) // 2
        b = (-control - difference) // 2
        for vector in dominant.get(d10_square, ()):
            root = (a, b, *vector.coordinates)
            if not _primitive(root) or pair(root, root) != -norm_abs:
                continue
            root_divisibility = divisibility(root)
            if (2 * root_divisibility) % norm_abs:
                continue
            if any(pair(root, wall) < 0 for wall in height_zero_simple_roots()):
                continue
            candidates.append(
                Candidate(
                    root=root,
                    control=control,
                    norm_abs=norm_abs,
                    root_divisibility=root_divisibility,
                    difference=difference,
                    d10_square=d10_square,
                )
            )
    return tuple(sorted(candidates, key=lambda candidate: candidate.root))


def ordered_level_keys(height_bound: Fraction | int) -> tuple[tuple[int, int], ...]:
    bound = Fraction(height_bound)
    keys: list[tuple[int, int]] = []
    for norm_abs in (2, 4):
        control = 1
        while Fraction(control * control, norm_abs) <= bound:
            keys.append((control, norm_abs))
            control += 1
    return tuple(sorted(keys, key=lambda key: (Fraction(key[0] ** 2, key[1]), key)))


@dataclass(frozen=True, slots=True)
class HeightLevel:
    height: Fraction
    keys: tuple[tuple[int, int], ...]
    candidate_count: int
    accepted_roots: tuple[Vector, ...]


@dataclass(frozen=True, slots=True)
class StoppingAudit:
    attempted: bool
    fired: bool
    finite_vertex_count: int | None
    cusp_count: int | None
    hyperideal_ray_count: int | None
    opposite_sheet_ray_count: int | None
    geometry: ChamberVerification | None


def _diagnostic_interior(roots: Sequence[Vector]) -> tuple[Golden, ...]:
    height_zero = height_zero_simple_roots()
    finite_gram = sp.Matrix(
        [[pair(left, right) for right in height_zero] for left in height_zero]
    )
    coefficients = finite_gram.inv() * sp.Matrix([-1] * len(height_zero))
    root_matrix = sp.Matrix.hstack(*(sp.Matrix(root) for root in height_zero))
    perturbation = tuple(Fraction(value) for value in root_matrix * coefficients)
    gram = tuple(tuple(Golden(-value) for value in row) for row in lattice_gram())
    normals = tuple(tuple(Golden(-value) for value in root) for root in roots)
    for denominator in (10, 100, 1_000, 10_000):
        point = (
            Golden(1) + Golden(perturbation[0]) / denominator,
            Golden(1) + Golden(perturbation[1]) / denominator,
            *(Golden(value) / denominator for value in perturbation[2:]),
        )
        norm = sum(
            point[row] * gram[row][column] * point[column]
            for row in range(12)
            for column in range(12)
        )
        values = tuple(
            sum(
                point[row] * gram[row][column] * root[column]
                for row in range(12)
                for column in range(12)
            )
            for root in normals
        )
        if norm.sign() < 0 and all(value.sign() < 0 for value in values):
            return point
    raise AssertionError("failed to construct a strict chamber point")


def audit_stopping(roots: Sequence[Vector]) -> StoppingAudit:
    if len(roots) < 12:
        return StoppingAudit(False, False, None, None, None, None, None)
    gram = tuple(tuple(Golden(-value) for value in row) for row in lattice_gram())
    normals = tuple(tuple(Golden(-value) for value in root) for root in roots)
    geometry = diagnose_chamber_cone(
        normals, gram, _diagnostic_interior(roots), dimension=11
    )
    fired = not geometry.hyperideal_rays and not geometry.opposite_sheet_rays
    return StoppingAudit(
        True,
        fired,
        len(geometry.finite_vertices),
        len(geometry.cusps),
        len(geometry.hyperideal_rays),
        len(geometry.opposite_sheet_rays),
        geometry,
    )


@dataclass(frozen=True, slots=True)
class Search:
    height_bound: Fraction
    levels: tuple[HeightLevel, ...]
    roots: tuple[Vector, ...]
    stopping: StoppingAudit


@lru_cache(maxsize=None)
def run_search(height_bound: Fraction | int, *, stop_when_finite: bool = True) -> Search:
    bound = Fraction(height_bound)
    roots = list(height_zero_simple_roots())
    levels: list[HeightLevel] = []
    keys = ordered_level_keys(bound)
    position = 0
    stopping = StoppingAudit(False, False, None, None, None, None, None)
    while position < len(keys):
        height = Fraction(keys[position][0] ** 2, keys[position][1])
        same_height: list[tuple[int, int]] = []
        while position < len(keys) and Fraction(keys[position][0] ** 2, keys[position][1]) == height:
            same_height.append(keys[position])
            position += 1
        candidates = tuple(
            candidate
            for control, norm_abs in same_height
            for candidate in enumerate_level_candidates(control, norm_abs)
        )
        accepted = tuple(
            candidate.root
            for candidate in candidates
            if all(pair(candidate.root, wall) >= 0 for wall in roots)
        )
        if any(
            pair(left, right) < 0
            for offset, left in enumerate(accepted)
            for right in accepted[offset + 1 :]
        ):
            raise AssertionError("same-height accepted roots are not pairwise acute")
        roots.extend(accepted)
        levels.append(HeightLevel(height, tuple(same_height), len(candidates), accepted))
        stopping = audit_stopping(roots)
        if stop_when_finite and stopping.fired:
            break
    return Search(bound, tuple(levels), tuple(roots), stopping)
