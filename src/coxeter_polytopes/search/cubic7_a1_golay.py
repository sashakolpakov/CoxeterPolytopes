"""Exact cubic-seven Vinberg prefix for ``<-(theta-1)> + N(A1^24)``.

Write ``alpha_i^2=2``.  In orthogonal A1 coordinates the Niemeier lattice is

``N(A1^24)={d/2 : every coefficient-parity word of d is in G_24}``.

There are three coefficient words over ``Z[theta]``.  For an unordered
coordinate multiset, the weights of their seven nonzero XOR combinations
are therefore Golay weights ``0,8,12,16,24``.  This necessary condition is
enough to reject every non-root-lattice pattern through height 100; the two
surviving patterns have zero parity words and hence have explicit labelings.
Thus the H=100 statement below is a complete, not heuristic, prefix census.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import comb, factorial
from typing import Iterable, Sequence, TypeAlias

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    THETA7,
    Cubic7,
)
from coxeter_polytopes.search.cubic7_timelike import _bounded_for_caps


Vector: TypeAlias = tuple[Cubic7, ...]
A_UNIT = THETA7 - 1
GOLAY_WEIGHTS = frozenset((0, 8, 12, 16, 24))

COMPLETE_MIRROR_NORM_PROOF = (
    "The full lattice is unimodular because theta-1 is a unit and the "
    "Niemeier block is unimodular. For primitive r its pairing ideal is O, "
    "so reflection integrality q(r)|2I(r) gives (q)|(2). The polynomial "
    "of theta is irreducible modulo 2, hence (2) is inert and the only norm "
    "ideals are O and (2). The Minkowski bound 14/9 proves class number one; "
    "the units -1, theta, theta^2-2 have full signature rank, so every "
    "totally positive unit is a square. Canonical positive norms are 1,2."
)


@dataclass(frozen=True, slots=True)
class Cubic7A1GolaySeed:
    seed_id: str
    timelike_coefficient: Cubic7
    canonical_mirror_norms: tuple[Cubic7, ...]
    complete_norm_class_proof: str

    @property
    def admissible(self) -> bool:
        return self.timelike_coefficient.sign(0) > 0 and all(
            self.timelike_coefficient.sign(place) < 0 for place in (1, 2)
        )


UNIT_A1_GOLAY_SEED = Cubic7A1GolaySeed(
    "cubic7-H24-A1^24-a-theta-minus-1",
    A_UNIT,
    (Cubic7(1), Cubic7(2)),
    COMPLETE_MIRROR_NORM_PROOF,
)


def cubic7_a1_pair(
    seed: Cubic7A1GolaySeed,
    left: Sequence[Cubic7],
    right: Sequence[Cubic7],
) -> Cubic7:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("a cubic-seven A1^24 vector has rank 25")
    return -seed.timelike_coefficient * left[0] * right[0] + 2 * sum(
        (left[index] * right[index] for index in range(1, 25)),
        start=CUBIC7_ZERO,
    )


def _coordinate_residue(value: Cubic7) -> int:
    doubled = 2 * value
    if not doubled.is_integral:
        raise AssertionError("a Construction-A coordinate is not half-integral")
    return (
        4 * (int(doubled.a) % 2)
        + 2 * (int(doubled.b) % 2)
        + int(doubled.c) % 2
    )


def _residue_profile(entries: Sequence[Cubic7]) -> tuple[int, ...]:
    counts = [0] * 8
    for value in entries:
        counts[_coordinate_residue(value)] += 1
    return tuple(counts)


def golay_xor_weight_necessary(profile: Sequence[int]) -> bool:
    """Necessary Golay-glue test for three coefficient parity words.

    It is exact for the zero profile used by both surviving H=100 packets.
    Failure is an unconditional nonmembership certificate.
    """

    if len(profile) != 8 or sum(profile) != 24:
        return False
    for mask in range(1, 8):
        weight = sum(
            count
            for cell, count in enumerate(profile)
            if (cell & mask).bit_count() % 2
        )
        if weight not in GOLAY_WEIGHTS:
            return False
    return True


def _root_lattice_labeled_count(entries: Sequence[Cubic7]) -> int | None:
    """Count labelings when all three parity words are zero."""

    if _residue_profile(entries) != (24, 0, 0, 0, 0, 0, 0, 0):
        return None
    multiplicities: dict[Cubic7, int] = {}
    for value in entries:
        multiplicities[value] = multiplicities.get(value, 0) + 1
    result = factorial(24)
    for count in multiplicities.values():
        result //= factorial(count)
    return result


@dataclass(frozen=True, slots=True)
class Cubic7A1ControlLevel:
    norm: Cubic7
    height: Cubic7
    control: Cubic7
    coordinate_square_target: Cubic7


@lru_cache(maxsize=None)
def enumerate_control_levels(
    seed: Cubic7A1GolaySeed, height_bound: int
) -> tuple[Cubic7A1ControlLevel, ...]:
    if seed != UNIT_A1_GOLAY_SEED or not seed.admissible:
        raise ValueError("unsupported or inadmissible cubic-seven A1 seed")
    if height_bound < 0:
        raise ValueError("height bound must be nonnegative")
    coefficient = seed.timelike_coefficient
    levels = set()
    for norm in seed.canonical_mirror_norms:
        caps = (
            norm * height_bound / coefficient,
            norm / (-coefficient),
            norm / (-coefficient),
        )
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            target = (norm + coefficient * control * control) / 2
            if any(target.sign(place) < 0 for place in (1, 2)):
                continue
            height = coefficient * control * control / norm
            if (Cubic7(height_bound) - height).sign(0) < 0:
                continue
            levels.add(Cubic7A1ControlLevel(norm, height, control, target))
    return tuple(sorted(levels, key=lambda level: (level.height, level.norm, level.control)))


@dataclass(frozen=True, slots=True)
class Cubic7A1GolayPattern:
    entries: tuple[Cubic7, ...]
    residue_profile: tuple[int, ...]
    spatial_sum: Cubic7
    labeled_count: int | None


@dataclass(frozen=True, slots=True)
class Cubic7A1PatternCensus:
    target: Cubic7
    alphabet_size: int
    raw_square_partition_count: int
    golay_necessary_patterns: tuple[Cubic7A1GolayPattern, ...]


@lru_cache(maxsize=None)
def enumerate_spatial_patterns(target: Cubic7) -> Cubic7A1PatternCensus:
    values = tuple(
        sorted(
            {
                doubled / 2
                for doubled in _bounded_for_caps((4 * target, 4 * target, 4 * target))
                if doubled.sign(0) > 0
                and all(
                    (4 * target - doubled * doubled).sign(place) >= 0
                    for place in range(3)
                )
            },
            reverse=True,
        )
    )
    squares = tuple(value * value for value in values)

    @lru_cache(maxsize=None)
    def solve(
        index: int, slots: int, remainder: Cubic7
    ) -> tuple[tuple[int, ...], ...]:
        if any(remainder.sign(place) < 0 for place in range(3)):
            return ()
        if index == len(values):
            return ((),) if not remainder else ()
        output = []
        running = remainder
        for multiplicity in range(slots + 1):
            if all(running.sign(place) >= 0 for place in range(3)):
                for suffix in solve(index + 1, slots - multiplicity, running):
                    output.append((multiplicity, *suffix))
            running -= squares[index]
        return tuple(output)

    raw = solve(0, 24, target)
    patterns = []
    for multiplicities in raw:
        entries = tuple(
            value
            for value, multiplicity in zip(values, multiplicities)
            for _ in range(multiplicity)
        )
        entries += (CUBIC7_ZERO,) * (24 - len(entries))
        profile = _residue_profile(entries)
        if not golay_xor_weight_necessary(profile):
            continue
        patterns.append(
            Cubic7A1GolayPattern(
                entries,
                profile,
                sum(entries, start=CUBIC7_ZERO),
                _root_lattice_labeled_count(entries),
            )
        )
    return Cubic7A1PatternCensus(target, len(values), len(raw), tuple(patterns))


@lru_cache(maxsize=None)
def height_zero_patterns(
    seed: Cubic7A1GolaySeed,
) -> tuple[Cubic7A1GolayPattern, ...]:
    return tuple(
        pattern
        for norm in seed.canonical_mirror_norms
        for pattern in enumerate_spatial_patterns(norm / 2).golay_necessary_patterns
    )


def _cubic_rank(rows: Iterable[Sequence[Cubic7]], width: int) -> int:
    basis: list[list[Cubic7]] = []
    pivots: list[int] = []
    for source in rows:
        row = list(source)
        for pivot, old in zip(pivots, basis):
            if row[pivot]:
                factor = row[pivot]
                row = [entry - factor * old[index] for index, entry in enumerate(row)]
        pivot = next((index for index, entry in enumerate(row) if entry), None)
        if pivot is None:
            continue
        divisor = row[pivot]
        basis.append([entry / divisor for entry in row])
        pivots.append(pivot)
        if len(basis) == width:
            break
    return len(basis)


@dataclass(frozen=True, slots=True)
class Cubic7A1PrefixGeometry:
    packet_count: int
    total_wall_count: int
    extreme_ray_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    uniform_coordinate: Cubic7
    uniform_norm: Cubic7
    active_wall_count: int
    active_spatial_rank: int
    finite_volume: bool
    compact: bool


@dataclass(frozen=True, slots=True)
class Cubic7A1UnitH100Audit:
    complete_height_bound: int
    control_level_count: int
    raw_square_partition_count: int
    nonempty_level_count: int
    first_level: Cubic7A1ControlLevel
    first_pattern: Cubic7A1GolayPattern
    first_geometry: Cubic7A1PrefixGeometry
    second_level: Cubic7A1ControlLevel
    second_pattern: Cubic7A1GolayPattern
    first_to_second_pairings: tuple[Cubic7, ...]
    second_same_height_pairings: tuple[Cubic7, ...]
    second_geometry: Cubic7A1PrefixGeometry


@lru_cache(maxsize=1)
def unit_prefix_through_h100() -> Cubic7A1UnitH100Audit:
    seed = UNIT_A1_GOLAY_SEED
    zero = height_zero_patterns(seed)
    if len(zero) != 1 or zero[0].entries.count(CUBIC7_ONE) != 1 or zero[0].labeled_count != 24:
        raise AssertionError("height zero is not exactly the A1^24 root chamber")

    records = []
    raw_count = 0
    levels = enumerate_control_levels(seed, 100)
    for level in levels:
        census = enumerate_spatial_patterns(level.coordinate_square_target)
        raw_count += census.raw_square_partition_count
        for pattern in census.golay_necessary_patterns:
            if pattern.labeled_count is None:
                raise AssertionError(
                    "H=100 has a nonzero-residue profile requiring an exact Golay triple audit"
                )
            records.append((level, pattern))
    if len(records) != 2:
        raise AssertionError("the complete H=100 prefix does not have two nonempty levels")
    (first_level, first), (second_level, second) = records

    a = seed.timelike_coefficient
    c = 2 * THETA7 + THETA7 * THETA7
    p = CUBIC7_ONE + c
    expected_first = tuple(sorted((c, *(CUBIC7_ZERO for _ in range(23))), reverse=True))
    expected_second = tuple(
        sorted((p, c, *(CUBIC7_ZERO for _ in range(22))), reverse=True)
    )
    if first.entries != expected_first or first.labeled_count != 24:
        raise AssertionError("the first positive packet is not the 24 coordinate caps")
    if second.entries != expected_second or second.labeled_count != 24 * 23:
        raise AssertionError("the second packet is not the ordered two-coordinate family")

    cap = a * first_level.control / (2 * c)
    first_signs = tuple((-a + 2 * count * cap * cap).sign(0) for count in range(25))
    if first_signs[:2] != (-1, -1) or any(sign != 1 for sign in first_signs[2:]):
        raise AssertionError("the cap-cube causal signs changed")
    first_geometry = Cubic7A1PrefixGeometry(
        1,
        48,
        2**24,
        25,
        0,
        2**24 - 25,
        cap,
        -a + 48 * cap * cap,
        24,
        24,
        False,
        False,
    )

    first_to_second = tuple(
        sorted(
            {
                -a * first_level.control * second_level.control + 2 * c * value
                for value in (p, c, CUBIC7_ZERO)
            }
        )
    )
    if any(value.sign(0) > 0 for value in first_to_second):
        raise AssertionError("the second packet crosses a first cap")

    # Four labels realize every equality pattern of two ordered pairs.
    same_height = set()
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            for k in range(4):
                for ell in range(4):
                    if k == ell or (i, j) == (k, ell):
                        continue
                    dot = (
                        (p * p if i == k else CUBIC7_ZERO)
                        + (p * c if i == ell else CUBIC7_ZERO)
                        + (c * p if j == k else CUBIC7_ZERO)
                        + (c * c if j == ell else CUBIC7_ZERO)
                    )
                    same_height.add(-a * second_level.control**2 + 2 * dot)
    if any(value.sign(0) > 0 for value in same_height):
        raise AssertionError("the second packet is not mutually acute")

    active_rows = (
        tuple(
            p if coordinate == i else c if coordinate == j else CUBIC7_ZERO
            for coordinate in range(24)
        )
        for i in range(24)
        for j in range(24)
        if i != j
    )
    rank = _cubic_rank(active_rows, 24)
    symmetric = a
    symmetric_norm = -a + 48 * symmetric * symmetric
    if rank != 24 or symmetric_norm.sign(0) <= 0:
        raise AssertionError("the second symmetric active ray lost its certificate")

    # Exact vertex classification of
    # 0<=z_i<=cap, p*z_i+c*z_j<=a*x/2 (i!=j):
    # 0; cap*e_i; a*1_S (|S|>=2); and
    # cap*e_i+(a/2)*1_T (empty!=T subset [24]\\{i}).
    equal_signs = tuple((-a + 2 * count * a * a).sign(0) for count in range(2, 25))
    mixed_signs = tuple(
        (-a + 2 * cap * cap + 2 * (count - 1) * (a / 2) ** 2).sign(0)
        for count in range(2, 25)
    )
    if equal_signs[0] != -1 or any(sign != 1 for sign in equal_signs[1:]):
        raise AssertionError("the equal-support vertex signs changed")
    if any(sign != 1 for sign in mixed_signs):
        raise AssertionError("a mixed cap vertex ceased to be hyperideal")
    mixed_count = 24 * (2**23 - 1)
    finite = 1 + 24 + comb(24, 2)
    total = 2**24 + mixed_count
    second_geometry = Cubic7A1PrefixGeometry(
        2,
        24 + 24 + 24 * 23,
        total,
        finite,
        0,
        total - finite,
        symmetric,
        symmetric_norm,
        24 * 23,
        rank,
        False,
        False,
    )
    return Cubic7A1UnitH100Audit(
        100,
        len(levels),
        raw_count,
        2,
        first_level,
        first,
        first_geometry,
        second_level,
        second,
        first_to_second,
        tuple(sorted(same_height)),
        second_geometry,
    )


__all__ = [
    "A_UNIT",
    "COMPLETE_MIRROR_NORM_PROOF",
    "Cubic7A1ControlLevel",
    "Cubic7A1GolayPattern",
    "Cubic7A1PatternCensus",
    "Cubic7A1PrefixGeometry",
    "Cubic7A1UnitH100Audit",
    "UNIT_A1_GOLAY_SEED",
    "cubic7_a1_pair",
    "enumerate_control_levels",
    "enumerate_spatial_patterns",
    "golay_xor_weight_necessary",
    "height_zero_patterns",
    "unit_prefix_through_h100",
]
