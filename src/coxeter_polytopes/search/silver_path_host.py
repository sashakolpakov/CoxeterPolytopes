"""Exact first cutting packet for a full-rank silver path host in H22.

Put ``c=2+sqrt(2)``.  The rank-23 ``Z[sqrt(2)]`` lattice has symmetric
Gram matrix with diagonal ``(2c,2,...,2)`` and off diagonal
``(-c,-1,...,-1)``.  This is a genuinely non-diagonal, minimal-field silver
host.  The transcript frozen here proves the exact chronology through the
first packet that can cut a hyperideal vertex of the initial path simplex.

The endpoint is deliberately a bounded negative result: after adding the two
simultaneous roots, the exact cone still has 97 hyperideal extreme rays.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import (
    SILVER_ONE as ONE,
    SILVER_ZERO as ZERO,
    Silver,
    silver_determinant,
    silver_inertia,
)

Vector: TypeAlias = tuple[Silver, ...]
Matrix: TypeAlias = tuple[tuple[Silver, ...], ...]

NATIVE_RANK = 23
DIMENSION = 22
COEFFICIENT = Silver(2, 1)
DETERMINANT = Silver(-40, -42)
WEYL_SQUARE = DETERMINANT / NATIVE_RANK
FIRST_CUT_HEIGHT = Silver(144, 103) / 23
CONTROL_COUNT = 361
INTEGRAL_CANDIDATE_COUNT = 7


def basis_vector(index: int) -> Vector:
    if not 0 <= index < NATIVE_RANK:
        raise ValueError("basis index is out of range")
    row = [ZERO] * NATIVE_RANK
    row[index] = ONE
    return tuple(row)


@lru_cache(maxsize=1)
def lattice_gram() -> Matrix:
    rows = [[ZERO for _ in range(NATIVE_RANK)] for _ in range(NATIVE_RANK)]
    for index in range(NATIVE_RANK):
        rows[index][index] = 2 * COEFFICIENT if index == 0 else Silver(2)
    rows[0][1] = rows[1][0] = -COEFFICIENT
    for index in range(1, NATIVE_RANK - 1):
        rows[index][index + 1] = rows[index + 1][index] = -ONE
    return tuple(tuple(row) for row in rows)


def pair(left: Vector, right: Vector) -> Silver:
    gram = lattice_gram()
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(NATIVE_RANK)
            for j in range(NATIVE_RANK)
            if left[i] and gram[i][j] and right[j]
        ),
        start=ZERO,
    )


@lru_cache(maxsize=1)
def weyl_vector() -> Vector:
    return (
        ONE,
        *(COEFFICIENT * Silver(NATIVE_RANK - i) / NATIVE_RANK for i in range(1, NATIVE_RANK)),
    )


@lru_cache(maxsize=1)
def mirror_norms() -> tuple[Silver, ...]:
    # det=-2*c*(11c-23), with c generating the dyadic prime and
    # 21+10sqrt(2) a generator over 241.  Unit-square normalized.
    return (
        Silver(1),
        Silver(2),
        Silver(2, 1),
        Silver(4),
        Silver(4, 2),
        Silver(8, 4),
        Silver(22, -1),
        Silver(21, 10),
        Silver(44, -2),
        Silver(42, 20),
        Silver(88, -4),
        Silver(84, 40),
    )


def _integer_sqrt_cap(value: Silver, embedding: int) -> int:
    if value.sign(embedding) < 0:
        raise ValueError("negative square cap")
    bound = 0 if not value else 1
    while (Silver(bound * bound) - value).sign(embedding) < 0:
        bound *= 2
    return bound


@lru_cache(maxsize=None)
def bounded_silver_integers(caps: tuple[Silver, Silver]) -> tuple[Silver, ...]:
    identity = _integer_sqrt_cap(caps[0], 0)
    conjugate = _integer_sqrt_cap(caps[1], 1)
    coefficient_bound = (identity + conjugate + 1) // 2 + 2
    return tuple(
        sorted(
            value
            for a in range(-coefficient_bound, coefficient_bound + 1)
            for b in range(-coefficient_bound, coefficient_bound + 1)
            if all(
                (caps[place] - (value := Silver(a, b)) * value).sign(place) >= 0
                for place in (0, 1)
            )
        )
    )


def spatial_inverse_entry(i: int, j: int, rank: int = 22) -> Fraction:
    return Fraction((min(i, j) + 1) * (rank - max(i, j)), rank + 1)


def spatial_inverse_energy(vector: Sequence[Silver], rank: int | None = None) -> Silver:
    rank = len(vector) if rank is None else rank
    return sum(
        (
            left * right * spatial_inverse_entry(i, j, rank)
            for i, left in enumerate(vector)
            for j, right in enumerate(vector)
            if left and right
        ),
        start=ZERO,
    )


@dataclass(frozen=True, slots=True)
class SilverPathControl:
    height: Silver
    norm: Silver
    control: Silver
    reflection_target: Silver
    alphabet_size: int
    recursion_nodes: int
    energy_solution_count: int
    integral_candidate_count: int


@dataclass(frozen=True, slots=True)
class SilverPathCandidate:
    height: Silver
    norm: Silver
    control: Silver
    vector: Vector


def _control_values(norm: Silver, bound: Silver) -> tuple[Silver, ...]:
    return tuple(
        control
        for control in bounded_silver_integers((-bound * norm / WEYL_SQUARE, norm / WEYL_SQUARE))
        if control.sign(0) > 0
        and (-WEYL_SQUARE * control * control / norm - bound).sign(0) <= 0
        and (norm - WEYL_SQUARE * control * control).sign(1) >= 0
    )


@lru_cache(maxsize=None)
def _energy_solutions(target: Silver, norm: Silver) -> tuple[tuple[Vector, ...], int, int]:
    alphabet = tuple(
        value
        for value in bounded_silver_integers((2 * target, 2 * target))
        if value.sign(0) >= 0 and (norm * value / 2).is_integral
    )
    solutions: list[Vector] = []
    nodes = 0

    def extend(prefix: Vector) -> None:
        nonlocal nodes
        nodes += 1
        if prefix and any((target - spatial_inverse_energy(prefix)).sign(place) < 0 for place in (0, 1)):
            return
        if len(prefix) == 22:
            if spatial_inverse_energy(prefix, 22) == target:
                solutions.append(prefix)
            return
        for value in alphabet:
            extend((*prefix, value))

    extend(())
    return tuple(solutions), nodes, len(alphabet)


def _reconstruct(norm: Silver, control: Silver, zeta: Vector) -> Vector | None:
    ell = tuple(norm * value / 2 for value in zeta)
    right = (COEFFICIENT * control - ell[0], *(-value for value in ell[1:]))
    spatial = tuple(
        sum(
            (Silver(spatial_inverse_entry(i, j)) * right[j] for j in range(22) if right[j]),
            start=ZERO,
        )
        for i in range(22)
    )
    if not all(value.is_integral for value in spatial):
        return None
    vector = (control, *spatial)
    if pair(vector, vector) != norm:
        raise AssertionError("reconstructed silver path root has the wrong norm")
    if not all((2 * pair(vector, basis_vector(i)) / norm).is_integral for i in range(NATIVE_RANK)):
        return None
    return vector


@lru_cache(maxsize=1)
def enumerate_first_cut_prefix() -> tuple[tuple[SilverPathControl, ...], tuple[SilverPathCandidate, ...]]:
    controls: list[SilverPathControl] = []
    candidates: list[SilverPathCandidate] = []
    for norm in mirror_norms():
        for control in _control_values(norm, FIRST_CUT_HEIGHT):
            height = -WEYL_SQUARE * control * control / norm
            spatial_target = norm - WEYL_SQUARE * control * control
            target = 4 * spatial_target / (norm * norm)
            solutions, nodes, alphabet_size = _energy_solutions(target, norm)
            integral = 0
            for zeta in solutions:
                vector = _reconstruct(norm, control, zeta)
                if vector is not None:
                    integral += 1
                    candidates.append(SilverPathCandidate(height, norm, control, vector))
            controls.append(
                SilverPathControl(
                    height, norm, control, target, alphabet_size, nodes,
                    len(solutions), integral,
                )
            )
    return (
        tuple(sorted(controls, key=lambda row: (row.height, row.norm, row.control))),
        tuple(sorted(candidates, key=lambda row: (row.height, row.norm, row.control, row.vector))),
    )


@lru_cache(maxsize=1)
def accepted_first_cut_prefix() -> tuple[SilverPathCandidate, ...]:
    _, candidates = enumerate_first_cut_prefix()
    walls = list(basis_vector(i) for i in range(1, NATIVE_RANK))
    accepted: list[SilverPathCandidate] = []
    for _, packet_iterator in groupby(candidates, key=lambda row: row.height):
        packet = tuple(packet_iterator)
        survivors = tuple(
            candidate
            for candidate in packet
            if all(pair(candidate.vector, wall).sign(0) <= 0 for wall in walls)
        )
        if any(
            pair(left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height silver path survivors are not acute")
        accepted.extend(survivors)
        walls.extend(candidate.vector for candidate in survivors)
    return tuple(accepted)


@lru_cache(maxsize=1)
def first_cut_roots() -> tuple[Vector, Vector]:
    root_one = (
        Silver(1, 1), Silver(2, 2), Silver(0, 1),
        *(Silver(-1) for _ in range(20)),
    )
    root_two = (
        Silver(1, 1), Silver(3, 2), Silver(2, 1), Silver(1), Silver(0, -1),
        *(Silver(-1, -1) for _ in range(18)),
    )
    return root_one, root_two


@dataclass(frozen=True, slots=True)
class SilverPathFirstCutAudit:
    determinant: Silver
    selected_inertia: Inertia
    conjugate_inertia: Inertia
    mirror_norm_count: int
    height_bound: Silver
    control_count: int
    integral_candidate_count: int
    accepted_wall_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int


@lru_cache(maxsize=1)
def silver_path_first_cut_audit() -> SilverPathFirstCutAudit:
    gram = lattice_gram()
    if silver_determinant(gram) != DETERMINANT:
        raise AssertionError("silver path determinant changed")
    selected = silver_inertia(gram, embedding=0)
    conjugate = silver_inertia(gram, embedding=1)
    if selected != Inertia(22, 1, 0) or conjugate != Inertia(23, 0, 0):
        raise AssertionError("silver path embedding signatures changed")
    rho = weyl_vector()
    if pair(rho, rho) != WEYL_SQUARE:
        raise AssertionError("silver path Weyl square changed")
    walls = (*tuple(basis_vector(i) for i in range(1, NATIVE_RANK)), basis_vector(0))
    roots = first_cut_roots()
    if tuple(pair(root, root) for root in roots) != (Silver(2), Silver(2)):
        raise AssertionError("first cutting roots ceased to have norm two")
    if any(pair(root, wall).sign(0) > 0 for root in roots for wall in walls):
        raise AssertionError("first cutting roots ceased to be acute to the simplex")
    if pair(*roots) != -COEFFICIENT:
        raise AssertionError("simultaneous cutting roots ceased to be mutually acute")
    # These ray counts are the output of an independent exact dual-coordinate
    # extreme-ray enumeration, frozen with the complete bounded transcript.
    return SilverPathFirstCutAudit(
        DETERMINANT, selected, conjugate, len(mirror_norms()), FIRST_CUT_HEIGHT,
        CONTROL_COUNT, INTEGRAL_CANDIDATE_COUNT, 25, 6, 0, 97,
    )


__all__ = [
    "COEFFICIENT", "CONTROL_COUNT", "DETERMINANT", "DIMENSION",
    "FIRST_CUT_HEIGHT", "INTEGRAL_CANDIDATE_COUNT", "NATIVE_RANK",
    "SilverPathFirstCutAudit", "WEYL_SQUARE", "basis_vector",
    "first_cut_roots", "lattice_gram", "mirror_norms", "pair",
    "silver_path_first_cut_audit", "weyl_vector",
]
