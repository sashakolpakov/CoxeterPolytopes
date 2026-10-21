"""Exact silver Vinberg shells for ``<-a> + N(A1^24)``.

The Niemeier lattice is written in orthogonal A1 coordinates, with
``alpha_i^2=2``::

    N(A1^24) = {d/2 : the two coefficient parity words of d lie in G_24}.

Here coefficients belong to ``Z[sqrt(2)]`` and ``G_24`` is the extended
binary Golay code.  Dominance for the height-zero A1 Weyl chamber is simply
coordinatewise nonnegativity at the distinguished real embedding.  Exact
Minkowski boxes, square-partition dynamic programming, and the joint Golay
residue table therefore give a complete *unlabelled* shell census.  A joint
table multiplicity supplies exact labelled counts without iterating millions
of coordinate permutations.

Three arithmetically admissible timelike coefficients are supported:
``1+sqrt(2)``, ``sqrt(2)``, and ``2sqrt(2)``.  The first two complete
positive packets for the unit coefficient are especially dense: 24
coordinate caps, followed by 6,072 octad/special-point walls.  Their common
symmetric ray is certified extreme and hyperideal.  A separate H=100 slope
audit deliberately ignores later-wall acuteness, hence gives a rigorous
*upper bound* for every accepted symmetric slope; even that relaxed bound
leaves an explicit hyperideal vector in the prefix cone.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby
from math import factorial
from typing import Iterable, Sequence, TypeAlias

from coxeter_polytopes.exact.silver import (
    SILVER_ONE,
    SILVER_ZERO,
    SQRT2,
    Silver,
)
from coxeter_polytopes.search.cubic7_leech import extended_binary_golay_codewords
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


Vector: TypeAlias = tuple[Silver, ...]
EPSILON = SILVER_ONE + SQRT2


UNIT_NORM_PROOF = (
    "Let pi=sqrt(2), so (2)=(pi)^2, and let a=epsilon=1+pi be a unit. "
    "For primitive r=(x,y), reflection integrality is (q(r))|2I(r), "
    "I(r)=(a*x,(y,N)). Away from (pi), unimodularity and primitivity make "
    "I(r)=O. Locally at (pi), if x is a unit then evenness gives "
    "q=-a*x^2 mod pi^2, so v_pi(q)=0. If x is divisible by pi, primitive "
    "y and unimodularity give I(r)=O; evenness gives v_pi(q)>=2 while "
    "q|2I(r) gives v_pi(q)<=2. Thus the norm ideals are O and (pi^2)=(2). "
    "The unit group is +/-epsilon^Z and epsilon has signs (+,-), so every "
    "totally positive unit is a square. Canonical positive norms are 1,2."
)

SQRT2_NORM_PROOF = (
    "Let pi=sqrt(2), so (2)=(pi)^2 and a=pi. For primitive r=(x,y), "
    "q(r)|2I(r). Away from (pi), unimodularity makes I(r)=O. If x is a "
    "pi-adic unit, q=-pi*x^2 mod pi^2 and v_pi(q)=1. Otherwise primitive "
    "y makes I(r)=O, while evenness gives v_pi(q)>=2 and q|2 gives "
    "v_pi(q)=2. Thus the ideals are (pi),(pi^2). Since every totally "
    "positive unit is a square, canonical positive generators are "
    "epsilon*pi=2+pi and pi^2=2."
)

TWO_SQRT2_NORM_PROOF = (
    "Let pi=sqrt(2), so a=2sqrt(2)=pi^3. Away from (pi), unimodularity "
    "and primitivity make the pairing ideal a unit. If y is locally "
    "primitive, evenness and q|2I force v_pi(q)=2. Otherwise y lies in "
    "pi*N and primitive x is a unit; then I(r)=(pi), the spatial norm is "
    "0 mod pi^4, and q=-pi^3*x^2 mod pi^4, so v_pi(q)=3. Hence the norm "
    "ideals are (pi^2),(pi^3). Modulo totally-positive unit squares their "
    "canonical positive generators are 2 and epsilon*pi^3=4+2pi. For the "
    "pi^3 class, reflection integrality is equivalently y in pi*N."
)


@dataclass(frozen=True, slots=True)
class SilverA1GolaySeed:
    seed_id: str
    timelike_coefficient: Silver
    canonical_mirror_norms: tuple[Silver, ...]
    complete_norm_class_proof: str

    @property
    def admissible(self) -> bool:
        return (
            self.timelike_coefficient.sign(0) > 0
            and self.timelike_coefficient.sign(1) < 0
        )


UNIT_A1_GOLAY_SEED = SilverA1GolaySeed(
    "silver-H24-A1^24-a-1+sqrt2",
    EPSILON,
    (SILVER_ONE, Silver(2)),
    UNIT_NORM_PROOF,
)
SQRT2_A1_GOLAY_SEED = SilverA1GolaySeed(
    "silver-H24-A1^24-a-sqrt2",
    SQRT2,
    (Silver(2), Silver(2, 1)),
    SQRT2_NORM_PROOF,
)
TWO_SQRT2_A1_GOLAY_SEED = SilverA1GolaySeed(
    "silver-H24-A1^24-a-2sqrt2",
    2 * SQRT2,
    (Silver(2), Silver(4, 2)),
    TWO_SQRT2_NORM_PROOF,
)
SUPPORTED_SILVER_A1_GOLAY_SEEDS = (
    UNIT_A1_GOLAY_SEED,
    SQRT2_A1_GOLAY_SEED,
    TWO_SQRT2_A1_GOLAY_SEED,
)


def silver_a1_pair(
    seed: SilverA1GolaySeed,
    left: Sequence[Silver],
    right: Sequence[Silver],
) -> Silver:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("a silver A1^24 vector has rank 25")
    return -seed.timelike_coefficient * left[0] * right[0] + 2 * sum(
        (left[index] * right[index] for index in range(1, 25)),
        start=SILVER_ZERO,
    )


def _coordinate_residue(value: Silver) -> tuple[int, int]:
    doubled = 2 * value
    if not doubled.is_integral:
        raise AssertionError("a Construction-A coordinate is not half-integral")
    return int(doubled.a) % 2, int(doubled.b) % 2


@lru_cache(maxsize=1)
def golay_joint_residue_table_multiplicities(
) -> dict[tuple[int, int, int, int], int]:
    """Count ordered Golay word pairs by their four joint residue counts."""

    words = extended_binary_golay_codewords()
    counts: Counter[tuple[int, int, int, int]] = Counter()
    for rational_word in words:
        rational_weight = rational_word.bit_count()
        for radical_word in words:
            radical_weight = radical_word.bit_count()
            both = (rational_word & radical_word).bit_count()
            rational_only = rational_weight - both
            radical_only = radical_weight - both
            neither = 24 - both - rational_only - radical_only
            counts[(neither, radical_only, rational_only, both)] += 1
    return dict(counts)


def _residue_counts(entries: Sequence[Silver]) -> tuple[int, int, int, int]:
    counts = [0, 0, 0, 0]
    for value in entries:
        rational, radical = _coordinate_residue(value)
        counts[2 * rational + radical] += 1
    return tuple(counts)


def _labeled_count(entries: Sequence[Silver]) -> int:
    residue_counts = _residue_counts(entries)
    word_pair_count = golay_joint_residue_table_multiplicities().get(
        residue_counts, 0
    )
    if not word_pair_count:
        return 0
    values_by_residue: dict[int, Counter[Silver]] = {}
    for value in entries:
        rational, radical = _coordinate_residue(value)
        values_by_residue.setdefault(2 * rational + radical, Counter())[value] += 1
    assignments = 1
    for residue, count in enumerate(residue_counts):
        assignments *= factorial(count)
        for multiplicity in values_by_residue.get(residue, {}).values():
            assignments //= factorial(multiplicity)
    return word_pair_count * assignments


def golay_membership(coordinates: Sequence[Silver]) -> bool:
    if len(coordinates) != 24:
        return False
    words = set(extended_binary_golay_codewords())
    coefficient_words = []
    for coefficient in ("a", "b"):
        doubled = [2 * getattr(value, coefficient) for value in coordinates]
        if any(value.denominator != 1 for value in doubled):
            return False
        coefficient_words.append(
            sum((int(value) % 2) << index for index, value in enumerate(doubled))
        )
    return all(word in words for word in coefficient_words)


@dataclass(frozen=True, slots=True)
class SilverA1ControlLevel:
    norm: Silver
    height: Silver
    control: Silver
    coordinate_square_target: Silver


@lru_cache(maxsize=None)
def enumerate_control_levels(
    seed: SilverA1GolaySeed, height_bound: int
) -> tuple[SilverA1ControlLevel, ...]:
    if seed not in SUPPORTED_SILVER_A1_GOLAY_SEEDS or not seed.admissible:
        raise ValueError("unsupported or arithmetically inadmissible silver A1 seed")
    levels = set()
    coefficient = seed.timelike_coefficient
    for norm in seed.canonical_mirror_norms:
        caps = (norm * height_bound / coefficient, norm / (-coefficient))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0 or not (2 * coefficient * control / norm).is_integral:
                continue
            coordinate_target = (norm + coefficient * control * control) / 2
            if coordinate_target.sign(1) < 0:
                continue
            height = coefficient * control * control / norm
            if (Silver(height_bound) - height).sign(0) < 0:
                continue
            levels.add(
                SilverA1ControlLevel(norm, height, control, coordinate_target)
            )
    return tuple(sorted(levels, key=lambda level: (level.height, level.norm, level.control)))


@dataclass(frozen=True, slots=True)
class SilverA1GolayPattern:
    entries: tuple[Silver, ...]
    residue_counts: tuple[int, int, int, int]
    spatial_sum: Silver
    labeled_count: int


@lru_cache(maxsize=None)
def enumerate_spatial_patterns(
    target: Silver,
) -> tuple[SilverA1GolayPattern, ...]:
    """Exhaust dominant Construction-A multisets with ``sum d_i^2=target``."""

    values = tuple(
        sorted(
            {
                doubled / 2
                for doubled in _bounded_for_caps((4 * target, 4 * target))
                if (doubled / 2).sign(0) >= 0
                and all(
                    (target - (doubled / 2) ** 2).sign(place) >= 0
                    for place in (0, 1)
                )
            },
            reverse=True,
        )
    )
    positive = tuple(value for value in values if value)
    squares = tuple(value * value for value in positive)

    @lru_cache(maxsize=None)
    def solve(
        index: int, slots: int, remainder: Silver
    ) -> tuple[tuple[int, ...], ...]:
        if any(remainder.sign(place) < 0 for place in (0, 1)):
            return ()
        if index == len(positive):
            return ((),) if not remainder else ()
        output = []
        running = remainder
        for multiplicity in range(slots + 1):
            if all(running.sign(place) >= 0 for place in (0, 1)):
                for suffix in solve(index + 1, slots - multiplicity, running):
                    output.append((multiplicity, *suffix))
            running -= squares[index]
        return tuple(output)

    patterns = []
    for multiplicities in solve(0, 24, target):
        entries = tuple(
            value
            for value, multiplicity in zip(positive, multiplicities)
            for _ in range(multiplicity)
        )
        entries += (SILVER_ZERO,) * (24 - len(entries))
        labeled_count = _labeled_count(entries)
        if not labeled_count:
            continue
        patterns.append(
            SilverA1GolayPattern(
                entries,
                _residue_counts(entries),
                sum(entries, start=SILVER_ZERO),
                labeled_count,
            )
        )
    return tuple(patterns)


def _reflective_pattern(
    seed: SilverA1GolaySeed,
    norm: Silver,
    pattern: SilverA1GolayPattern,
) -> SilverA1GolayPattern | None:
    if seed != TWO_SQRT2_A1_GOLAY_SEED or norm != Silver(4, 2):
        return pattern
    divided = tuple(value / SQRT2 for value in pattern.entries)
    if any(not (2 * value).is_integral for value in divided):
        return None
    labeled_count = _labeled_count(divided)
    if not labeled_count:
        return None
    return SilverA1GolayPattern(
        pattern.entries,
        pattern.residue_counts,
        pattern.spatial_sum,
        labeled_count,
    )


@lru_cache(maxsize=None)
def enumerate_level_patterns(
    seed: SilverA1GolaySeed,
    level: SilverA1ControlLevel,
) -> tuple[SilverA1GolayPattern, ...]:
    output = []
    for pattern in enumerate_spatial_patterns(level.coordinate_square_target):
        reflective = _reflective_pattern(seed, level.norm, pattern)
        if reflective is not None:
            output.append(reflective)
    return tuple(output)


@lru_cache(maxsize=None)
def height_zero_patterns(
    seed: SilverA1GolaySeed,
) -> tuple[SilverA1GolayPattern, ...]:
    output = []
    for norm in seed.canonical_mirror_norms:
        for pattern in enumerate_spatial_patterns(norm / 2):
            reflective = _reflective_pattern(seed, norm, pattern)
            if reflective is not None:
                output.append(reflective)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class SilverA1FirstPacketGeometry:
    height: Silver
    norm: Silver
    control: Silver
    cap_coordinate: Silver
    accepted_wall_count: int
    total_wall_count: int
    extreme_ray_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    uniform_ray_norm: Silver
    finite_volume: bool
    compact: bool


@lru_cache(maxsize=None)
def first_packet_cube_geometry(
    seed: SilverA1GolaySeed,
) -> SilverA1FirstPacketGeometry:
    """Find the first nonempty packet and classify its complete cube cone."""

    if len(height_zero_patterns(seed)) != 1:
        raise AssertionError("height zero is not exactly the A1^24 root chamber")
    selected = None
    for level in enumerate_control_levels(seed, 10):
        patterns = enumerate_level_patterns(seed, level)
        if patterns:
            selected = (level, patterns)
            break
    if selected is None:
        raise AssertionError("no silver A1 first packet through H=10")
    level, patterns = selected
    if len(patterns) != 1:
        raise AssertionError("the reflective first packet has more than one multiset")
    pattern = patterns[0]
    nonzero = tuple(value for value in pattern.entries if value)
    if len(nonzero) != 1 or pattern.labeled_count != 24:
        raise AssertionError("the first packet is not 24 coordinate caps")
    spatial = nonzero[0]
    cap = seed.timelike_coefficient * level.control / (2 * spatial)
    signs = tuple(
        (-seed.timelike_coefficient + 2 * Silver(count) * cap * cap).sign(0)
        for count in range(25)
    )
    binomials = tuple(
        factorial(24) // (factorial(count) * factorial(24 - count))
        for count in range(25)
    )
    finite = sum(weight for weight, sign in zip(binomials, signs) if sign < 0)
    cusps = sum(weight for weight, sign in zip(binomials, signs) if sign == 0)
    hyperideal = sum(weight for weight, sign in zip(binomials, signs) if sign > 0)
    uniform_norm = -seed.timelike_coefficient + 48 * cap * cap
    return SilverA1FirstPacketGeometry(
        level.height,
        level.norm,
        level.control,
        cap,
        pattern.labeled_count,
        24 + pattern.labeled_count,
        2**24,
        finite,
        cusps,
        hyperideal,
        uniform_norm,
        hyperideal == 0,
        hyperideal == 0 and cusps == 0,
    )


def _silver_rank(rows: Sequence[Sequence[Silver]], width: int) -> int:
    work = [list(row) for row in rows]
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for row in range(rank + 1, len(work)):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for index in range(column, width):
                work[row][index] -= factor * work[rank][index]
        rank += 1
        if rank == width:
            break
    return rank


def _independent_row_basis(
    rows: Iterable[Sequence[Silver]], width: int
) -> tuple[tuple[Silver, ...], ...]:
    """Stream a row echelon basis, stopping as soon as full rank is reached."""

    pivots: dict[int, list[Silver]] = {}
    for source in rows:
        row = list(source)
        for column, basis_row in sorted(pivots.items()):
            if not row[column]:
                continue
            factor = row[column]
            for index in range(column, width):
                row[index] -= factor * basis_row[index]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        value = row[pivot]
        row = [entry / value for entry in row]
        for column, basis_row in tuple(pivots.items()):
            if not basis_row[pivot]:
                continue
            factor = basis_row[pivot]
            pivots[column] = [
                basis_row[index] - factor * row[index] for index in range(width)
            ]
        pivots[pivot] = row
        if len(pivots) == width:
            break
    return tuple(tuple(row) for _pivot, row in sorted(pivots.items()))


@lru_cache(maxsize=1)
def _octads() -> tuple[int, ...]:
    words = tuple(
        word for word in extended_binary_golay_codewords() if word.bit_count() == 8
    )
    if len(words) != 759:
        raise AssertionError("the extended Golay code lost its 759 octads")
    return words


@dataclass(frozen=True, slots=True)
class SilverA1UnitDensePrefix:
    complete_height_bound: int
    first_height: Silver
    first_family_count: int
    second_height: Silver
    second_control: Silver
    second_pattern: tuple[tuple[Silver, int], ...]
    second_family_count: int
    second_family_rank: int
    second_same_height_margins: tuple[Silver, ...]
    symmetric_coordinate: Silver
    symmetric_norm: Silver
    symmetric_extreme: bool
    symmetric_hyperideal: bool


@lru_cache(maxsize=1)
def unit_dense_prefix_through_h20() -> SilverA1UnitDensePrefix:
    """Certify the first two complete nonempty unit-coefficient packets."""

    seed = UNIT_A1_GOLAY_SEED
    nonempty = []
    for height, iterator in groupby(
        enumerate_control_levels(seed, 20), key=lambda level: level.height
    ):
        levels = tuple(iterator)
        patterns = tuple(
            (level, pattern)
            for level in levels
            for pattern in enumerate_level_patterns(seed, level)
        )
        if patterns:
            nonempty.append((height, patterns))
    if len(nonempty) != 2:
        raise AssertionError("the H=20 unit prefix does not have two nonempty packets")
    first_height, first = nonempty[0]
    second_height, second = nonempty[1]
    if len(first) != 1 or first[0][1].labeled_count != 24 or len(second) != 1:
        raise AssertionError("unexpected unit packet decomposition")
    level, pattern = second[0]
    high = (Silver(3) + SQRT2) / 2
    low = EPSILON / 2
    expected = tuple(sorted((high, *(low for _ in range(7)), *(SILVER_ZERO for _ in range(16))), reverse=True))
    if pattern.entries != expected or pattern.labeled_count != 759 * 8:
        raise AssertionError("the second unit packet is not the octad/special family")
    cap = first_packet_cube_geometry(seed).cap_coordinate
    if any((cap * level.control - value).sign(0) < 0 for value in pattern.entries):
        raise AssertionError("the second packet crosses a coordinate cap")

    cases = (
        (0, 0, 0),
        *((intersection, indicators, equal) for intersection in (2, 4) for indicators, equal in ((0, 0), (1, 0), (2, 0), (2, 1))),
        (8, 2, 0),
    )
    threshold = seed.timelike_coefficient * level.control * level.control / 2
    margins = tuple(
        low * low * intersection + low * indicators + equal - threshold
        for intersection, indicators, equal in cases
    )
    if any(margin.sign(0) > 0 for margin in margins):
        raise AssertionError("the octad/special family is not mutually acute")

    rows = (
        tuple(
            high
            if index == special
            else low
            if octad & (1 << index)
            else SILVER_ZERO
            for index in range(24)
        )
        for octad in _octads()
        for special in range(24)
        if octad & (1 << special)
    )
    rank = len(_independent_row_basis(rows, 24))
    spatial_sum = pattern.spatial_sum
    coordinate = seed.timelike_coefficient * level.control / (2 * spatial_sum)
    norm = -seed.timelike_coefficient + 48 * coordinate * coordinate
    return SilverA1UnitDensePrefix(
        20,
        first_height,
        first[0][1].labeled_count,
        second_height,
        level.control,
        ((high, 1), (low, 7), (SILVER_ZERO, 16)),
        pattern.labeled_count,
        rank,
        margins,
        coordinate,
        norm,
        rank == 24,
        norm.sign(0) > 0,
    )


@dataclass(frozen=True, slots=True)
class SilverA1SlopeUpperBound:
    height_bound: int
    relaxed_pattern_count: int
    maximum_slope: Silver
    maximizing_height: Silver
    maximizing_control: Silver
    critical_slope_square: Silver
    slope_square_below_critical: bool
    feasible_uniform_coordinate: Silver
    feasible_uniform_norm: Silver
    hyperideal_prefix_witness: bool


@dataclass(frozen=True, slots=True)
class SilverA1UnitThirdPacket:
    complete_height_bound: int
    height: Silver
    control: Silver
    norm: Silver
    unlabelled_pattern_count: int
    labeled_wall_count: int
    prior_packet_rearrangement_margins: tuple[Silver, ...]
    possible_positive_same_height_pairings: tuple[Silver, ...]
    positive_pairing_difference_pattern_counts: tuple[int, ...]
    same_height_mutually_acute: bool
    accepted_wall_count: int
    total_positive_wall_count: int
    maximum_slope: Silver
    maximum_slope_family_count: int
    maximum_slope_family_rank: int
    symmetric_coordinate: Silver
    symmetric_norm: Silver
    symmetric_extreme: bool
    symmetric_hyperideal: bool


@lru_cache(maxsize=1)
def unit_third_packet_through_h30() -> SilverA1UnitThirdPacket:
    """Certify the complete million-wall packet at ``h=14+10sqrt(2)``.

    Rearrangement bounds prove acuteness to the earlier octad/special packet
    without expanding labels.  Mutual same-height acuteness follows from a
    short difference-shell argument.  Indeed, conjugate Cauchy bounds force
    a positive integral pairing of two norm-two roots to be 1 or sqrt(2).
    Their spatial difference has coordinate-square target 1 or 2-sqrt(2).
    The latter shell is empty; the former is one A1 coordinate root.  Two
    dominant equal-square vectors cannot differ by that root, since
    ``d_i-e_i=+/-1`` and ``d_i^2=e_i^2`` would give ``d_i+e_i=0``.
    """

    seed = UNIT_A1_GOLAY_SEED
    target_height = 14 + 10 * SQRT2
    levels = tuple(
        level
        for level in enumerate_control_levels(seed, 30)
        if level.height == target_height
    )
    rows = tuple(
        (level, pattern)
        for level in levels
        for pattern in enumerate_level_patterns(seed, level)
    )
    if len(rows) != 12 or {level.norm for level, _pattern in rows} != {Silver(2)}:
        raise AssertionError("the third unit packet lost its twelve norm-two patterns")
    controls = {level.control for level, _pattern in rows}
    if controls != {2 + 2 * SQRT2}:
        raise AssertionError("the third unit packet has an unexpected control")
    control = next(iter(controls))

    old = unit_dense_prefix_through_h20()
    old_pattern = tuple(
        value
        for value, multiplicity in old.second_pattern
        for _ in range(multiplicity)
    )
    old_descending = tuple(sorted(old_pattern, reverse=True))
    prior_margins = []
    for _level, pattern in rows:
        maximum_dot = sum(
            (
                left * right
                for left, right in zip(
                    tuple(sorted(pattern.entries, reverse=True)), old_descending
                )
            ),
            start=SILVER_ZERO,
        )
        prior_margins.append(
            maximum_dot
            - seed.timelike_coefficient * control * old.second_control / 2
        )
    if any(margin.sign(0) >= 0 for margin in prior_margins):
        raise AssertionError("a third-packet pattern crosses an earlier orbit")
    cap = first_packet_cube_geometry(seed).cap_coordinate
    if any(
        (cap * control - value).sign(0) < 0
        for _level, pattern in rows
        for value in pattern.entries
    ):
        raise AssertionError("a third-packet pattern crosses a coordinate cap")

    possible_pairings = tuple(
        value
        for value in _bounded_for_caps((Silver(4), Silver(4)))
        if value.is_integral
        and value.sign(0) > 0
        and (Silver(2) - value).sign(0) > 0
        and (Silver(2) - value.conjugate()).sign(0) >= 0
        and (Silver(2) + value.conjugate()).sign(0) >= 0
    )
    possible_pairings = tuple(sorted(set(possible_pairings)))
    if possible_pairings != (SQRT2, SILVER_ONE):
        # Silver's identity order puts sqrt(2) after 1, so use a set check
        # below and retain a deterministic mathematical order in the audit.
        if set(possible_pairings) != {SILVER_ONE, SQRT2}:
            raise AssertionError("positive same-height pairing interval changed")
        possible_pairings = (SILVER_ONE, SQRT2)
    difference_counts = tuple(
        len(enumerate_spatial_patterns(Silver(2) - pairing))
        for pairing in possible_pairings
    )
    if difference_counts != (1, 0):
        raise AssertionError("the two small difference shells changed")
    same_height_acute = True

    maximum_slope = max(
        (pattern.spatial_sum / level.control for level, pattern in rows)
    )
    maximum_rows = tuple(
        pattern
        for level, pattern in rows
        if pattern.spatial_sum / level.control == maximum_slope
    )
    if len(maximum_rows) != 1 or maximum_slope != Silver(9) / 2:
        raise AssertionError("the third-packet maximum slope changed")
    maximum_pattern = maximum_rows[0]
    high = EPSILON
    low = EPSILON / 2
    family_rows = (
        tuple(
            high
            if index == special
            else low
            if not octad & (1 << index)
            else SILVER_ZERO
            for index in range(24)
        )
        for octad in _octads()
        for special in range(24)
        if octad & (1 << special)
    )
    family_rank = len(_independent_row_basis(family_rows, 24))
    if maximum_pattern.labeled_count != 759 * 8 or family_rank != 24:
        raise AssertionError("the maximum-slope hexadecad family lost rank")
    coordinate = seed.timelike_coefficient / (2 * maximum_slope)
    norm = -seed.timelike_coefficient + 48 * coordinate * coordinate
    labeled = sum(pattern.labeled_count for _level, pattern in rows)
    return SilverA1UnitThirdPacket(
        30,
        target_height,
        control,
        Silver(2),
        len(rows),
        labeled,
        tuple(prior_margins),
        possible_pairings,
        difference_counts,
        same_height_acute,
        labeled,
        24 + old.second_family_count + labeled,
        maximum_slope,
        maximum_pattern.labeled_count,
        family_rank,
        coordinate,
        norm,
        family_rank == 24,
        norm.sign(0) > 0,
    )


@lru_cache(maxsize=None)
def unit_symmetric_slope_upper_bound(
    height_bound: int = 100,
) -> SilverA1SlopeUpperBound:
    """Bound every accepted unit-prefix slope by a relaxed exact census.

    Only height-zero dominance, Golay glue, reflection integrality, and the
    first 24 coordinate caps are imposed.  Omitting all later walls enlarges
    the candidate set, so the returned maximum is an upper bound for the
    actual Vinberg prefix.  Choosing ``c=a/(2*maximum)`` therefore gives an
    explicit vector ``(1;c^24)`` in the true prefix cone.
    """

    seed = UNIT_A1_GOLAY_SEED
    cap = first_packet_cube_geometry(seed).cap_coordinate
    maximum = None
    maximizing_level = None
    pattern_count = 0
    for level in enumerate_control_levels(seed, height_bound):
        for pattern in enumerate_level_patterns(seed, level):
            if any((cap * level.control - value).sign(0) < 0 for value in pattern.entries):
                continue
            pattern_count += 1
            slope = pattern.spatial_sum / level.control
            if maximum is None or (slope - maximum).sign(0) > 0:
                maximum = slope
                maximizing_level = level
    if maximum is None or maximizing_level is None:
        raise AssertionError("the relaxed slope census is empty")
    critical_square = 12 * seed.timelike_coefficient
    below = (critical_square - maximum * maximum).sign(0) > 0
    coordinate = seed.timelike_coefficient / (2 * maximum)
    norm = -seed.timelike_coefficient + 48 * coordinate * coordinate
    return SilverA1SlopeUpperBound(
        height_bound,
        pattern_count,
        maximum,
        maximizing_level.height,
        maximizing_level.control,
        critical_square,
        below,
        coordinate,
        norm,
        below and norm.sign(0) > 0,
    )


__all__ = [
    "EPSILON",
    "SQRT2_A1_GOLAY_SEED",
    "SUPPORTED_SILVER_A1_GOLAY_SEEDS",
    "SilverA1ControlLevel",
    "SilverA1FirstPacketGeometry",
    "SilverA1GolayPattern",
    "SilverA1GolaySeed",
    "SilverA1SlopeUpperBound",
    "SilverA1UnitDensePrefix",
    "SilverA1UnitThirdPacket",
    "TWO_SQRT2_A1_GOLAY_SEED",
    "UNIT_A1_GOLAY_SEED",
    "enumerate_control_levels",
    "enumerate_level_patterns",
    "enumerate_spatial_patterns",
    "first_packet_cube_geometry",
    "golay_membership",
    "height_zero_patterns",
    "silver_a1_pair",
    "unit_dense_prefix_through_h20",
    "unit_symmetric_slope_upper_bound",
    "unit_third_packet_through_h30",
]
