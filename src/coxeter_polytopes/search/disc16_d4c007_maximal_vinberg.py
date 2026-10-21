"""Exact Vinberg search on ``U + D6(-1) + A1(-1)^3``.

This is the unique maximal integral overlattice preserved by the 27 mirror
reflections of the Shimada discriminant-16 atlas cast ``D0004_C007``.  Write
a vector as ``r=(a,b;v,z1,z2,z3)`` and use ``p=e+f`` as the control vector.
For ``n=-(p,r)=-(a+b)>0`` and ``d=a-b``, a root of square ``-q`` obeys

    d^2 + 2*v_D6^2 + 4*(z1^2+z2^2+z3^2) = n^2 + 2*q,
    d == n (mod 2).

The full height-zero reflection arrangement is ``B6 + B3 + A1``: norm-four
roots extend the evident ``D6 + A1^3`` factors, while the ``e-f`` root remains
a separate ``A1`` because its divisibility is one.  Since the discriminant
group has exponent two, primitive integral reflection roots can only have
square ``-2`` or ``-4``; integrality is checked directly for every candidate.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import gcd, isqrt
from typing import Sequence

import sympy as sp

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.verify.chamber import ChamberVerification, diagnose_chamber_cone


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


@lru_cache(maxsize=1)
def d6_cartan() -> Matrix:
    matrix = [[2 if left == right else 0 for right in range(6)] for left in range(6)]
    for left, right in ((0, 1), (1, 2), (2, 3), (3, 4), (3, 5)):
        matrix[left][right] = matrix[right][left] = -1
    return tuple(tuple(row) for row in matrix)


@lru_cache(maxsize=1)
def lattice_gram() -> Matrix:
    gram = [[0 for _ in range(11)] for _ in range(11)]
    gram[0][1] = gram[1][0] = 1
    d6 = d6_cartan()
    for left in range(6):
        for right in range(6):
            gram[left + 2][right + 2] = -d6[left][right]
    for index in range(8, 11):
        gram[index][index] = -2
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
        for row in range(11)
        for column in range(11)
    )


def divisibility(root: Sequence[int]) -> int:
    result = 0
    for entry in _mat_vec(lattice_gram(), root):
        result = gcd(result, abs(entry))
    return result


def _primitive(root: Sequence[int]) -> bool:
    content = 0
    for entry in root:
        content = gcd(content, abs(entry))
    return content == 1


@lru_cache(maxsize=1)
def height_zero_simple_roots() -> tuple[Vector, ...]:
    roots: list[Vector] = []
    # B6 in D6 root-basis coordinates: five D6 simple roots followed by the
    # norm-four root alpha_5-alpha_6.
    roots.extend(
        (0, 0, *tuple(-1 if index == simple else 0 for index in range(6)), 0, 0, 0)
        for simple in range(5)
    )
    roots.append((0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0))
    # The divisibility-one e-f root and B3 on z1,z2,z3.
    roots.append((1, -1, *(0,) * 9))
    roots.append((0, 0, *(0,) * 6, 1, -1, 0))
    roots.append((0, 0, *(0,) * 6, 0, 1, -1))
    roots.append((0, 0, *(0,) * 6, 0, 0, 1))
    result = tuple(roots)
    if len(result) != 10:
        raise AssertionError("wrong height-zero rank")
    if any(
        pair(left, right) < 0
        for position, left in enumerate(result)
        for right in result[position + 1 :]
    ):
        raise AssertionError("height-zero roots are not pairwise acute")
    return result


@dataclass(frozen=True, slots=True)
class DominantD6Vector:
    coordinates: tuple[int, ...]
    weights: tuple[int, ...]
    square: int


@lru_cache(maxsize=None)
def dominant_d6_vectors(max_square: int) -> tuple[DominantD6Vector, ...]:
    if max_square < 0:
        return ()
    cartan = sp.Matrix(d6_cartan())
    inverse = cartan.inv()
    bounds: list[int] = []
    for index in range(6):
        diagonal = Fraction(inverse[index, index])
        quotient = Fraction(max_square, 1) / diagonal
        bounds.append(isqrt(quotient.numerator // quotient.denominator) + 1)
    vectors: list[DominantD6Vector] = []
    for weights in product(*(range(bound + 1) for bound in bounds)):
        weight = sp.Matrix(weights)
        square = (weight.T * inverse * weight)[0]
        if square.q != 1 or square > max_square:
            continue
        coordinates = inverse * weight
        if any(entry.q != 1 for entry in coordinates):
            continue
        vector = DominantD6Vector(
            coordinates=tuple(int(entry) for entry in coordinates),
            weights=tuple(weights),
            square=int(square),
        )
        vectors.append(vector)
    return tuple(sorted(vectors, key=lambda vector: vector.coordinates))


@dataclass(frozen=True, slots=True)
class Candidate:
    root: Vector
    control: int
    norm_abs: int
    root_divisibility: int
    difference: int
    d6_square: int

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
    d6_by_square: dict[int, list[DominantD6Vector]] = {}
    for vector in dominant_d6_vectors(radius // 2):
        if vector.weights[5] < vector.weights[4]:
            continue
        d6_by_square.setdefault(vector.square, []).append(vector)

    candidates: list[Candidate] = []
    for difference in range(-isqrt(radius), 1):
        if (difference - control) % 2:
            continue
        a = (-control + difference) // 2
        b = (-control - difference) // 2
        z_bound = isqrt(radius // 4)
        for z1, z2, z3 in product(range(-z_bound, z_bound + 1), repeat=3):
            remainder = radius - difference * difference - 4 * (
                z1 * z1 + z2 * z2 + z3 * z3
            )
            if remainder < 0 or remainder % 2:
                continue
            d6_square = remainder // 2
            for vector in d6_by_square.get(d6_square, ()):
                root = (a, b, *vector.coordinates, z1, z2, z3)
                if not _primitive(root) or pair(root, root) != -norm_abs:
                    continue
                if any(pair(root, wall) < 0 for wall in height_zero_simple_roots()):
                    continue
                root_divisibility = divisibility(root)
                if (2 * root_divisibility) % norm_abs:
                    continue
                candidates.append(
                    Candidate(
                        root=root,
                        control=control,
                        norm_abs=norm_abs,
                        root_divisibility=root_divisibility,
                        difference=difference,
                        d6_square=d6_square,
                    )
                )

    result = tuple(sorted(candidates, key=lambda candidate: candidate.root))
    if any(
        pair(candidate.root, wall) < 0
        for candidate in result
        for wall in height_zero_simple_roots()
    ):
        raise AssertionError("dominant enumeration crossed a height-zero wall")
    return result


def ordered_level_keys(control_bound: int) -> tuple[tuple[int, int], ...]:
    if control_bound < 0:
        raise ValueError("control bound must be nonnegative")
    cutoff = Fraction(control_bound * control_bound, 2)
    keys = [
        (control, norm_abs)
        for norm_abs in (2, 4)
        for control in range(
            1,
            isqrt(
                (cutoff * norm_abs).numerator // (cutoff * norm_abs).denominator
            )
            + 2,
        )
        if Fraction(control * control, norm_abs) <= cutoff
    ]
    return tuple(sorted(keys, key=lambda item: (Fraction(item[0] ** 2, item[1]), item)))


@dataclass(frozen=True, slots=True)
class HeightLevel:
    control: int
    norm_abs: int
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
    coefficients = finite_gram.inv() * sp.Matrix([-1] * 10)
    root_matrix = sp.Matrix.hstack(*(sp.Matrix(root) for root in height_zero))
    tangent = root_matrix * coefficients
    perturbation = tuple(Fraction(entry) for entry in tangent)
    diagnostic_gram = tuple(
        tuple(Golden(-entry) for entry in row) for row in lattice_gram()
    )
    diagnostic_roots = tuple(tuple(Golden(-entry) for entry in root) for root in roots)
    for denominator in (10, 100, 1_000, 10_000):
        point = (
            Golden(1) + Golden(perturbation[0]) / denominator,
            Golden(1) + Golden(perturbation[1]) / denominator,
            *(Golden(entry) / denominator for entry in perturbation[2:]),
        )
        norm = sum(
            point[row] * diagnostic_gram[row][column] * point[column]
            for row in range(11)
            for column in range(11)
        )
        values = tuple(
            sum(
                point[row] * diagnostic_gram[row][column] * root[column]
                for row in range(11)
                for column in range(11)
            )
            for root in diagnostic_roots
        )
        if norm.sign() < 0 and all(value.sign() < 0 for value in values):
            return point
    raise AssertionError("failed to construct a strict chamber point")


def audit_stopping(roots: Sequence[Vector]) -> StoppingAudit:
    if len(roots) < 11:
        return StoppingAudit(False, False, None, None, None, None, None)
    gram = tuple(tuple(Golden(-entry) for entry in row) for row in lattice_gram())
    normals = tuple(tuple(Golden(-entry) for entry in root) for root in roots)
    geometry = diagnose_chamber_cone(
        normals, gram, _diagnostic_interior(roots), dimension=10
    )
    fired = not geometry.hyperideal_rays and not geometry.opposite_sheet_rays
    return StoppingAudit(
        attempted=True,
        fired=fired,
        finite_vertex_count=len(geometry.finite_vertices),
        cusp_count=len(geometry.cusps),
        hyperideal_ray_count=len(geometry.hyperideal_rays),
        opposite_sheet_ray_count=len(geometry.opposite_sheet_rays),
        geometry=geometry,
    )


@dataclass(frozen=True, slots=True)
class Search:
    control_bound: int
    levels: tuple[HeightLevel, ...]
    roots: tuple[Vector, ...]
    stopping: StoppingAudit


@lru_cache(maxsize=None)
def run_search(control_bound: int, *, attempt_stopping: bool = False) -> Search:
    roots = list(height_zero_simple_roots())
    levels: list[HeightLevel] = []
    for control, norm_abs in ordered_level_keys(control_bound):
        candidates = enumerate_level_candidates(control, norm_abs)
        accepted = tuple(
            candidate.root
            for candidate in candidates
            if all(pair(candidate.root, wall) >= 0 for wall in roots)
        )
        if any(
            pair(left, right) < 0
            for position, left in enumerate(accepted)
            for right in accepted[position + 1 :]
        ):
            raise AssertionError("same-height accepted roots are not pairwise acute")
        roots.extend(accepted)
        levels.append(HeightLevel(control, norm_abs, len(candidates), accepted))
    stopping = (
        audit_stopping(roots)
        if attempt_stopping
        else StoppingAudit(False, False, None, None, None, None, None)
    )
    return Search(control_bound, tuple(levels), tuple(roots), stopping)
