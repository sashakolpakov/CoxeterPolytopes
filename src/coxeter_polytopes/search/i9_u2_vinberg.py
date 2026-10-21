"""Exact Vinberg search on ``U(2) + <-1>^9``.

This rank-eleven Lorentz lattice is the negative of Scharlau--Kirschmer
``List10/Gram10`` entry 6.  It is also the globally half-scaled reflection
lattice occurring in the Alexeev ``(15,11,1)`` fixed section.  Global
scaling does not change any reflection.

Write a vector as ``r=(a,b;z_1,...,z_9)`` and use ``p=e+f`` as control
vector.  Put ``m=-(a+b)>0`` and ``d=a-b``.  If ``r^2=-q``, then

    d^2 + z_1^2 + ... + z_9^2 = m^2 + q,
    d == m (mod 2).

The height-zero system is ``A1+B9``.  A primitive integral reflection root
has ``q`` in ``{1,2,4}``: if some ``z_i`` is odd its divisibility is one,
while otherwise primitivity forces divisibility two.  The implementation
still checks primitivity, norm and reflection integrality for every vector.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt
from typing import Iterator, Sequence

import sympy as sp

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.verify.chamber import ChamberVerification, diagnose_chamber_cone


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


@lru_cache(maxsize=1)
def lattice_gram() -> Matrix:
    gram = [[0 for _ in range(11)] for _ in range(11)]
    gram[0][1] = gram[1][0] = 2
    for index in range(2, 11):
        gram[index][index] = -1
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
    roots: list[Vector] = [(1, -1, *(0,) * 9)]
    # B9: -e_i+e_{i+1} (i<9), followed by the short root -e_9.
    for index in range(8):
        tail = [0] * 9
        tail[index] = -1
        tail[index + 1] = 1
        roots.append((0, 0, *tail))
    roots.append((0, 0, *(0,) * 8, -1))
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


def _square_partitions(
    total: int, length: int, maximum: int
) -> Iterator[tuple[int, ...]]:
    """Yield nonincreasing nonnegative vectors with prescribed square sum."""

    if length == 0:
        if total == 0:
            yield ()
        return
    upper = min(maximum, isqrt(total))
    for value in range(upper, -1, -1):
        remainder = total - value * value
        if remainder < 0 or remainder > (length - 1) * value * value:
            continue
        for suffix in _square_partitions(remainder, length - 1, value):
            yield (value, *suffix)


@dataclass(frozen=True, slots=True)
class Candidate:
    root: Vector
    control: int
    norm_abs: int
    root_divisibility: int
    difference: int

    @property
    def height(self) -> Fraction:
        return Fraction(self.control * self.control, self.norm_abs)


@lru_cache(maxsize=None)
def enumerate_level_candidates(control: int, norm_abs: int) -> tuple[Candidate, ...]:
    if control <= 0:
        raise ValueError("control must be positive")
    if norm_abs not in (1, 2, 4):
        raise ValueError("supported root norms are -1, -2 and -4")

    radius = control * control + norm_abs
    candidates: list[Candidate] = []
    for difference in range(-isqrt(radius), 1):
        if (difference - control) % 2:
            continue
        a = (-control + difference) // 2
        b = (-control - difference) // 2
        remainder = radius - difference * difference
        for tail in _square_partitions(remainder, 9, isqrt(remainder)):
            root = (a, b, *tail)
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
                )
            )
    return tuple(sorted(candidates, key=lambda candidate: candidate.root))


def ordered_level_keys(height_bound: Fraction | int) -> tuple[tuple[int, int], ...]:
    bound = Fraction(height_bound)
    if bound < 0:
        raise ValueError("height bound must be nonnegative")
    keys: list[tuple[int, int]] = []
    for norm_abs in (1, 2, 4):
        control = 1
        while Fraction(control * control, norm_abs) <= bound:
            keys.append((control, norm_abs))
            control += 1
    return tuple(sorted(keys, key=lambda item: (Fraction(item[0] ** 2, item[1]), item)))


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
    positions = 0
    stopping = StoppingAudit(False, False, None, None, None, None, None)
    while positions < len(keys):
        height = Fraction(keys[positions][0] ** 2, keys[positions][1])
        same_height: list[tuple[int, int]] = []
        while positions < len(keys) and Fraction(
            keys[positions][0] ** 2, keys[positions][1]
        ) == height:
            same_height.append(keys[positions])
            positions += 1
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
            for position, left in enumerate(accepted)
            for right in accepted[position + 1 :]
        ):
            raise AssertionError("same-height accepted roots are not pairwise acute")
        roots.extend(accepted)
        levels.append(HeightLevel(height, tuple(same_height), len(candidates), accepted))
        if len(roots) >= 11:
            stopping = audit_stopping(roots)
            if stop_when_finite and stopping.fired:
                break
    return Search(bound, tuple(levels), tuple(roots), stopping)
