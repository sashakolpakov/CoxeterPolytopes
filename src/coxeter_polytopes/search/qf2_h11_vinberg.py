"""Exact Vinberg search for ``<2> + <-1>^11``.

This is the sign convention of the odd determinant-two Lorentzian lattice
``I_{11,1}([1^11 2^1]_2)``.  With control vector ``p=(1,0,...,0)``, write a
positive-height root as ``r=(-m;z_1,...,z_11)`` in the dominant ``B11``
chamber.  A root of square ``-q`` satisfies

    z_1^2 + ... + z_11^2 = 2*m^2 + q.

The only possible primitive integral reflection norms are ``-1,-2,-4``.
Every enumerated vector is nevertheless checked for norm, primitivity,
divisibility and reflection integrality.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt
from typing import Iterator, Sequence

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.verify.chamber import ChamberVerification, diagnose_chamber_cone


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


@lru_cache(maxsize=1)
def lattice_gram() -> Matrix:
    return tuple(
        tuple((2 if row == 0 else -1) if row == column else 0 for column in range(12))
        for row in range(12)
    )


def pair(left: Sequence[int], right: Sequence[int]) -> int:
    return 2 * left[0] * right[0] - sum(
        left[index] * right[index] for index in range(1, 12)
    )


def divisibility(root: Sequence[int]) -> int:
    result = abs(2 * root[0])
    for value in root[1:]:
        result = gcd(result, abs(value))
    return result


def _primitive(root: Sequence[int]) -> bool:
    result = 0
    for value in root:
        result = gcd(result, abs(value))
    return result == 1


@lru_cache(maxsize=1)
def height_zero_simple_roots() -> tuple[Vector, ...]:
    roots: list[Vector] = []
    for index in range(10):
        tail = [0] * 11
        tail[index] = -1
        tail[index + 1] = 1
        roots.append((0, *tail))
    roots.append((0, *(0,) * 10, -1))
    result = tuple(roots)
    if len(result) != 11 or any(
        pair(left, right) < 0
        for position, left in enumerate(result)
        for right in result[position + 1 :]
    ):
        raise AssertionError("invalid height-zero B11 chamber")
    return result


def _square_partitions(
    total: int, length: int, maximum: int
) -> Iterator[tuple[int, ...]]:
    """Yield nonincreasing nonnegative vectors with the given square sum."""

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

    @property
    def height(self) -> Fraction:
        return Fraction(4 * self.control * self.control, self.norm_abs)


@lru_cache(maxsize=None)
def enumerate_level_candidates(control: int, norm_abs: int) -> tuple[Candidate, ...]:
    if control <= 0:
        raise ValueError("control must be positive")
    if norm_abs not in (1, 2, 4):
        raise ValueError("supported root norms are -1, -2 and -4")
    total = 2 * control * control + norm_abs
    candidates = []
    for tail in _square_partitions(total, 11, isqrt(total)):
        root = (-control, *tail)
        if not _primitive(root) or pair(root, root) != -norm_abs:
            continue
        root_divisibility = divisibility(root)
        if (2 * root_divisibility) % norm_abs:
            continue
        if any(pair(root, wall) < 0 for wall in height_zero_simple_roots()):
            continue
        candidates.append(Candidate(root, control, norm_abs, root_divisibility))
    return tuple(sorted(candidates, key=lambda candidate: candidate.root))


def ordered_level_keys(height_bound: Fraction | int) -> tuple[tuple[int, int], ...]:
    bound = Fraction(height_bound)
    keys = []
    for norm_abs in (1, 2, 4):
        control = 1
        while Fraction(4 * control * control, norm_abs) <= bound:
            keys.append((control, norm_abs))
            control += 1
    return tuple(
        sorted(keys, key=lambda key: (Fraction(4 * key[0] ** 2, key[1]), key))
    )


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
    # With our inward wall orientations, use z_1<...<z_11<0 near p.
    for denominator in (100, 1_000, 10_000, 100_000):
        point = (Golden(1),) + tuple(
            -Golden(12 - index) / denominator for index in range(1, 12)
        )
        gram = tuple(tuple(Golden(-value) for value in row) for row in lattice_gram())
        normals = tuple(tuple(Golden(-value) for value in root) for root in roots)
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
    levels = []
    keys = ordered_level_keys(bound)
    position = 0
    stopping = StoppingAudit(False, False, None, None, None, None, None)
    while position < len(keys):
        height = Fraction(4 * keys[position][0] ** 2, keys[position][1])
        same_height = []
        while position < len(keys) and Fraction(
            4 * keys[position][0] ** 2, keys[position][1]
        ) == height:
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
            raise AssertionError("same-height Vinberg roots are not pairwise acute")
        roots.extend(accepted)
        levels.append(HeightLevel(height, tuple(same_height), len(candidates), accepted))
        stopping = audit_stopping(roots)
        if stop_when_finite and stopping.fired:
            break
    return Search(bound, tuple(levels), tuple(roots), stopping)
