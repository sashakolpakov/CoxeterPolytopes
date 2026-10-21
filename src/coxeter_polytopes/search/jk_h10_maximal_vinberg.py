"""Exact Vinberg search for the full reflection group of
``U + E8(-2) + <-4>``.

This is the intrinsic lattice of the 84-wall Jacobian--Kummer ``H10`` catch.
The induced divisor chamber need not be a chamber of the *maximal* integral
reflection group, so this module computes the latter independently.

Write a vector as ``r=(a,b;v,z)`` and use the control vector ``p=e+f``.  For
``n=-(p,r)=-(a+b)>0`` and ``d=a-b``, a root of square ``-q`` satisfies

    d^2 + 4*v_E8^2 + 8*z^2 = n^2 + 2*q,   d == n (mod 2).

The height-zero chamber is ``A1+E8+A1``.  Hence ``d<=0``, the E8 Dynkin
weights are nonnegative, and ``z>=0``.  Primitive integral reflection roots
can only have ``q`` in ``{2,4,8}``; the divisibility test is performed
directly rather than assumed.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt
from typing import Sequence

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.lattices import block_diagonal, e8_cartan_gram, hyperbolic_plane
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.s19_vinberg import dominant_e8_vectors
from coxeter_polytopes.verify.chamber import ChamberVerification, diagnose_chamber_cone


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


@lru_cache(maxsize=1)
def lattice_gram() -> Matrix:
    e8 = [[-2 * entry for entry in row] for row in e8_cartan_gram()]
    return tuple(
        tuple(row)
        for row in block_diagonal(hyperbolic_plane(), e8, [[-4]])
    )


def _mat_vec(matrix: Sequence[Sequence[int]], vector: Sequence[int]) -> tuple[int, ...]:
    return tuple(sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix)


def pair(left: Sequence[int], right: Sequence[int]) -> int:
    gram = lattice_gram()
    return sum(left[i] * gram[i][j] * right[j] for i in range(11) for j in range(11))


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
    roots.extend(
        (0, 0, *tuple(-1 if i == j else 0 for i in range(8)), 0)
        for j in range(8)
    )
    roots.append((0, 0, *(0,) * 8, -1))
    result = tuple(roots)
    if len(result) != 10:
        raise AssertionError("wrong height-zero rank")
    if any(pair(left, right) < 0 for i, left in enumerate(result) for right in result[i + 1 :]):
        raise AssertionError("height-zero roots are not pairwise acute")
    return result


@dataclass(frozen=True, slots=True)
class Candidate:
    root: Vector
    control: int
    norm_abs: int
    root_divisibility: int
    difference: int
    e8_square: int

    @property
    def height(self) -> Fraction:
        return Fraction(self.control * self.control, self.norm_abs)


@lru_cache(maxsize=None)
def enumerate_level_candidates(control: int, norm_abs: int) -> tuple[Candidate, ...]:
    if control <= 0:
        raise ValueError("control must be positive")
    if norm_abs not in (2, 4, 8):
        raise ValueError("supported root norms are -2, -4, and -8")

    radius = control * control + 2 * norm_abs
    e8_by_square: dict[int, list[object]] = {}
    for vector in dominant_e8_vectors(radius // 4):
        e8_by_square.setdefault(vector.square, []).append(vector)

    candidates: list[Candidate] = []
    for difference in range(-isqrt(radius), 1):
        if (difference - control) % 2:
            continue
        a = (-control + difference) // 2
        b = (-control - difference) // 2
        for z in range(isqrt(radius // 8) + 1):
            remainder = radius - difference * difference - 8 * z * z
            if remainder < 0 or remainder % 4:
                continue
            e8_square = remainder // 4
            for vector in e8_by_square.get(e8_square, ()):  # type: ignore[union-attr]
                root = (a, b, *vector.coordinates, z)  # type: ignore[union-attr]
                if not _primitive(root) or pair(root, root) != -norm_abs:
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
                        e8_square=e8_square,
                    )
                )

    result = tuple(sorted(candidates, key=lambda candidate: candidate.root))
    if any(pair(candidate.root, wall) < 0 for candidate in result for wall in height_zero_simple_roots()):
        raise AssertionError("dominant enumeration crossed a height-zero wall")
    return result


def ordered_level_keys(control_bound: int) -> tuple[tuple[int, int], ...]:
    if control_bound < 0:
        raise ValueError("control bound must be nonnegative")
    cutoff = Fraction(control_bound * control_bound, 2)
    keys = [
        (control, norm_abs)
        for norm_abs in (2, 4, 8)
        for control in range(1, isqrt((cutoff * norm_abs).numerator // (cutoff * norm_abs).denominator) + 2)
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
    inverse = _fraction_inverse(
        tuple(tuple(Fraction(entry) for entry in row) for row in e8_cartan_gram())
    )
    anti_weyl = tuple(-sum(row) for row in inverse)
    if any(entry.denominator != 1 for entry in anti_weyl):
        raise AssertionError("E8 anti-Weyl vector is not integral")
    perturbation = (1, -1, *(int(entry) for entry in anti_weyl), -1)
    diagnostic_gram = tuple(tuple(Golden(-entry) for entry in row) for row in lattice_gram())
    diagnostic_roots = tuple(tuple(Golden(-entry) for entry in root) for root in roots)
    for denominator in (10, 100, 1_000, 10_000):
        point = (
            Golden(1) + Golden(perturbation[0]) / denominator,
            Golden(1) + Golden(perturbation[1]) / denominator,
            *(Golden(entry) / denominator for entry in perturbation[2:]),
        )
        norm = sum(point[i] * diagnostic_gram[i][j] * point[j] for i in range(11) for j in range(11))
        values = tuple(
            sum(point[i] * diagnostic_gram[i][j] * root[j] for i in range(11) for j in range(11))
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
    geometry = diagnose_chamber_cone(normals, gram, _diagnostic_interior(roots), dimension=10)
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
        if any(pair(left, right) < 0 for i, left in enumerate(accepted) for right in accepted[i + 1 :]):
            raise AssertionError("same-height accepted roots are not pairwise acute")
        roots.extend(accepted)
        levels.append(HeightLevel(control, norm_abs, len(candidates), accepted))
    stopping = audit_stopping(roots) if attempt_stopping else StoppingAudit(False, False, None, None, None, None, None)
    return Search(control_bound, tuple(levels), tuple(roots), stopping)
