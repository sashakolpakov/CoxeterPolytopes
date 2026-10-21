"""Exact first two packets for ``<-2phi> + N(A1^24)``.

Let ``alpha_i^2=2``.  In alpha-coordinates the Niemeier lattice with root
system A1^24 is

``{ d/2 : d mod 2 is an extended binary Golay codeword }``.

The coefficient condition is imposed independently in both coefficients of
``Z[phi]``.  This is the glue-aware Construction-A model; in particular the
support of every half-integral coefficient vector is a Golay codeword.

At control x=1 the exact coefficient comparison leaves the 24 cap roots
``(1;phi e_i)``.  At x=phi the norm equation leaves the 276 pair roots
``(phi;phi(e_i+e_j))`` and the 759 octad roots
``(phi;(phi/2)1_O)``.  They form simultaneous, mutually acute Vinberg
packets.

The vector ``u=(1;(phi/4)^24)`` lies on every octad wall.  The rational rank
of the 759 octad incidence vectors is 24, so u is an exact extreme ray; its
norm is ``3+phi>0``.  A later root ``(x;d)`` contains u exactly when
``sum(d_i)<=4x`` because ``<u,(x;d)>=(phi/2)(sum d_i-4x)``.  The helper below
audits this persistence inequality on any exact root prefix.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from typing import Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.search.cubic7_leech import extended_binary_golay_codewords
from coxeter_polytopes.search.golden_dn import _coefficient_box


Vector = tuple[Golden, ...]


COMPLETE_MIRROR_NORM_PROOF = (
    "N(A1^24) is even unimodular, so a primitive reflective vector has norm "
    "ideal dividing (2). Evenness excludes the unit ideal. At the inert "
    "prime (2), primitivity and q|2I(r) force valuation exactly one. "
    "Q(sqrt(5)) has narrow class number one and every totally positive unit "
    "is a square, hence every mirror has canonical norm element 2."
)


def a1_golay_pair(left: Sequence[Golden], right: Sequence[Golden]) -> Golden:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("the golden A1^24 host has rank 25")
    return -2 * PHI * left[0] * right[0] + 2 * sum(
        (left[index] * right[index] for index in range(1, 25)), start=ZERO
    )


def golay_membership(coordinates: Sequence[Golden]) -> bool:
    """Construction-A membership in the alpha-coordinate model."""

    if len(coordinates) != 24:
        return False
    words = set(extended_binary_golay_codewords())
    coefficient_words = []
    for coefficient in ("a", "b"):
        doubled = [2 * getattr(value, coefficient) for value in coordinates]
        if any(value.denominator != 1 for value in doubled):
            return False
        word = sum((int(value) % 2) << index for index, value in enumerate(doubled))
        coefficient_words.append(word)
    return all(word in words for word in coefficient_words)


@lru_cache(maxsize=1)
def height_zero_walls() -> tuple[Vector, ...]:
    return tuple(
        (ZERO, *tuple(-ONE if coordinate == index else ZERO for coordinate in range(24)))
        for index in range(24)
    )


@lru_cache(maxsize=1)
def octads() -> tuple[int, ...]:
    result = tuple(word for word in extended_binary_golay_codewords() if word.bit_count() == 8)
    if len(result) != 759:
        raise AssertionError("the extended Golay code lost its 759 octads")
    return result


def _rational_rank(rows: Iterable[Sequence[int]], width: int) -> int:
    work = [[Fraction(value) for value in row] for row in rows]
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][column]
        work[rank] = [entry / value for entry in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                work[row][index] - multiplier * work[rank][index]
                for index in range(width)
            ]
        rank += 1
        if rank == width:
            break
    return rank


@lru_cache(maxsize=1)
def octad_incidence_rank() -> int:
    return _rational_rank(
        ([int(word & (1 << index) != 0) for index in range(24)] for word in octads()),
        24,
    )


@dataclass(frozen=True, slots=True)
class GoldenA1GolayRoot:
    orbit: str
    control: Golden
    vector: Vector


@lru_cache(maxsize=1)
def x1_packet() -> tuple[GoldenA1GolayRoot, ...]:
    roots = []
    for index in range(24):
        spatial = tuple(PHI if coordinate == index else ZERO for coordinate in range(24))
        if not golay_membership(spatial):
            raise AssertionError("root-lattice cap left Construction A")
        roots.append(GoldenA1GolayRoot("single", ONE, (ONE, *spatial)))
    return tuple(roots)


@lru_cache(maxsize=1)
def xphi_packet() -> tuple[GoldenA1GolayRoot, ...]:
    roots: list[GoldenA1GolayRoot] = []
    for first, second in combinations(range(24), 2):
        spatial = tuple(PHI if coordinate in (first, second) else ZERO for coordinate in range(24))
        roots.append(GoldenA1GolayRoot("pair", PHI, (PHI, *spatial)))
    for word in octads():
        spatial = tuple(PHI / 2 if word & (1 << coordinate) else ZERO for coordinate in range(24))
        if not golay_membership(spatial):
            raise AssertionError("Golay octad root failed Construction A")
        roots.append(GoldenA1GolayRoot("octad", PHI, (PHI, *spatial)))
    result = tuple(roots)
    if len(result) != 1035:
        raise AssertionError("x=phi packet does not have 276+759 roots")
    return result


@dataclass(frozen=True, slots=True)
class UniformRayAudit:
    vector: Vector
    norm: Golden
    octad_incidence_rank: int
    octad_wall_count: int
    all_octad_pairings_zero: bool
    x1_pairings_nonpositive: bool
    xphi_pairings_nonpositive: bool
    hyperideal: bool
    extreme: bool


@lru_cache(maxsize=1)
def uniform_ray_audit() -> UniformRayAudit:
    vector = (ONE, *(PHI / 4 for _ in range(24)))
    norm = a1_golay_pair(vector, vector)
    x1_ok = all(a1_golay_pair(vector, root.vector).sign(0) <= 0 for root in x1_packet())
    xphi_ok = all(a1_golay_pair(vector, root.vector).sign(0) <= 0 for root in xphi_packet())
    octad_zero = all(
        not a1_golay_pair(vector, root.vector)
        for root in xphi_packet()
        if root.orbit == "octad"
    )
    rank = octad_incidence_rank()
    return UniformRayAudit(
        vector,
        norm,
        rank,
        len(octads()),
        octad_zero,
        x1_ok,
        xphi_ok,
        norm.sign(0) > 0,
        rank == 24 and octad_zero,
    )


@dataclass(frozen=True, slots=True)
class PersistenceViolation:
    index: int
    control: Golden
    spatial_sum: Golden
    excess: Golden
    vector: Vector


@lru_cache(maxsize=1)
def golay_joint_residue_counts() -> frozenset[tuple[int, int, int, int]]:
    """Possible counts of the four joint coefficient parities.

    A coordinate value ``d=(A+B phi)/2`` has residue ``(A mod 2,B mod 2)``.
    Construction A asks the two 24-bit support words to lie independently in
    the extended Golay code.  Their contingency table is therefore the exact
    assignment criterion for an unordered coordinate multiset.
    """

    words = extended_binary_golay_codewords()
    tables = set()
    for rational_word in words:
        rational_weight = rational_word.bit_count()
        for radical_word in words:
            radical_weight = radical_word.bit_count()
            both = (rational_word & radical_word).bit_count()
            rational_only = rational_weight - both
            radical_only = radical_weight - both
            neither = 24 - both - rational_only - radical_only
            tables.add((neither, radical_only, rational_only, both))
    return frozenset(tables)


@dataclass(frozen=True, slots=True)
class UniformInequalityLevelAudit:
    control: Golden
    height: Golden
    target: Golden
    square_partition_count: int
    golay_admissible_pattern_count: int
    maximum_spatial_sum: Golden | None
    maximum_excess: Golden | None
    maximum_patterns: tuple[tuple[Golden, ...], ...]
    violating_pattern_count: int
    first_violating_pattern: tuple[Golden, ...] | None


@dataclass(frozen=True, slots=True)
class MaximumSlopeLevelAudit:
    """Maximum spatial sum without materializing every square partition."""

    control: Golden
    height: Golden
    target: Golden
    alphabet_size: int
    dynamic_state_count: int
    feasible_joint_profile_count: int
    maximum_spatial_sum: Golden | None
    slope: Golden | None
    maximum_patterns: tuple[tuple[Golden, ...], ...]
    beats_hyperbolic_threshold: bool


@dataclass(frozen=True, slots=True)
class ExplicitUniformCutterAudit:
    root: Vector
    control: Golden
    height: Golden
    norm: Golden
    spatial_sum: Golden
    uniform_pairing: Golden
    golay_member: bool
    height_zero_pairings_nonpositive: bool
    x1_pairings_nonpositive: bool
    xphi_pairings_nonpositive: bool
    earlier_relaxed_pairing_margins: tuple[Golden, ...]
    accepted_through_earlier_levels: bool
    labeled_family_count: int
    octad_intersection_sizes: tuple[int, ...]
    same_height_pairing_margins: tuple[tuple[tuple[int, int, int], Golden], ...]
    labeled_family_mutually_acute: bool


@dataclass(frozen=True, slots=True)
class ReplacementUniformRayAudit:
    """Exact symmetric ray after the first 6,072-wall cutter orbit."""

    vector: Vector
    spatial_coordinate: Golden
    norm: Golden
    active_control: Golden
    active_spatial_sum: Golden
    active_orbit_count: int
    active_spatial_rank: int
    all_active_pairings_zero: bool
    all_earlier_shells_feasible: bool
    full_active_shell_feasible: bool
    hyperideal: bool
    extreme: bool
    later_root_slope: Golden


def _golden_incremental_rank(rows: Iterable[Sequence[Golden]], width: int) -> int:
    """Rank over Q(phi), retaining only a small echelon basis."""

    basis: list[list[Golden]] = []
    pivots: list[int] = []
    for source in rows:
        if len(source) != width:
            raise ValueError("a rank row has the wrong width")
        row = list(source)
        for pivot, old in zip(pivots, basis):
            if not row[pivot]:
                continue
            factor = row[pivot]
            row = [entry - factor * old[index] for index, entry in enumerate(row)]
        pivot = next((index for index, entry in enumerate(row) if entry), None)
        if pivot is None:
            continue
        divisor = row[pivot]
        row = [entry / divisor for entry in row]
        basis.append(row)
        pivots.append(pivot)
        if len(basis) == width:
            break
    return len(basis)


def _first_cutter_spatial(octad: int, special: int) -> tuple[Golden, ...]:
    if not octad & (1 << special):
        raise ValueError("the distinguished point must lie in the complementary octad")
    hexadecad = ((1 << 24) - 1) ^ octad
    return tuple(
        PHI
        if index == special
        else (ONE + PHI) / 2
        if hexadecad & (1 << index)
        else ZERO
        for index in range(24)
    )


def _first_cutter_spatials() -> Iterable[tuple[Golden, ...]]:
    for octad in octads():
        for special in range(24):
            if octad & (1 << special):
                yield _first_cutter_spatial(octad, special)


@lru_cache(maxsize=1)
def first_cutter_spatial_rank() -> int:
    """Rank of the complete 759*8 active spatial incidence orbit."""

    return _golden_incremental_rank(_first_cutter_spatials(), 24)


def _coordinate_residue(value: Golden) -> tuple[int, int]:
    doubled = 2 * value
    if not doubled.is_integral:
        raise AssertionError("a Construction-A coordinate is not half-integral")
    return int(doubled.a) % 2, int(doubled.b) % 2


@lru_cache(maxsize=None)
def _complete_golay_patterns(
    control: Golden,
) -> tuple[int, tuple[tuple[Golden, ...], ...]]:
    """Return raw count and every Golay-admissible dominant multiset."""

    target = ONE + PHI * control * control
    positive = tuple(
        sorted(
            {
                value / 2
                for value in _coefficient_box(4 * target)
                if value.sign(0) > 0
                and all((4 * target - value * value).sign(place) >= 0 for place in (0, 1))
            },
            reverse=True,
        )
    )
    squares = tuple(value * value for value in positive)

    @lru_cache(maxsize=None)
    def solve(index: int, slots: int, remainder: Golden) -> tuple[tuple[int, ...], ...]:
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

    partitions = solve(0, 24, target)
    allowed = golay_joint_residue_counts()
    patterns = []
    for multiplicities in partitions:
        residue_counts = [0, 0, 0, 0]
        entries = []
        for value, multiplicity in zip(positive, multiplicities):
            rational, radical = _coordinate_residue(value)
            residue_counts[2 * rational + radical] += multiplicity
            entries.extend(value for _ in range(multiplicity))
        residue_counts[0] += 24 - sum(multiplicities)
        if tuple(residue_counts) in allowed:
            patterns.append(tuple(entries) + (ZERO,) * (24 - len(entries)))
    return len(partitions), tuple(patterns)


@lru_cache(maxsize=None)
def audit_uniform_inequality_level(control: Golden) -> UniformInequalityLevelAudit:
    """Exhaust one norm-two shell and test ``sum(d_i)<=4x``.

    Coordinates are enumerated as unordered nonnegative multisets.  The
    exact Golay contingency-table criterion proves whether some labeling of
    each multiset belongs to the Niemeier lattice; because the tested sum is
    symmetric, expanding those labelings is unnecessary.
    """

    if not control.is_integral or control.sign(0) <= 0:
        raise ValueError("a control must be a positive golden integer")
    target = ONE + PHI * control * control
    raw_count, patterns = _complete_golay_patterns(control)
    violating = 0
    maximum_sum = None
    maximum_excess = None
    maximum_patterns: list[tuple[Golden, ...]] = []
    first_violation = None
    for pattern in patterns:
        spatial_sum = sum(pattern, start=ZERO)
        excess = spatial_sum - 4 * control
        if maximum_sum is None or (spatial_sum - maximum_sum).sign(0) > 0:
            maximum_sum = spatial_sum
            maximum_excess = excess
            maximum_patterns = [pattern]
        elif spatial_sum == maximum_sum:
            maximum_patterns.append(pattern)
        if excess.sign(0) > 0:
            violating += 1
            if first_violation is None:
                first_violation = pattern
    return UniformInequalityLevelAudit(
        control,
        PHI * control * control,
        target,
        raw_count,
        len(patterns),
        maximum_sum,
        maximum_excess,
        tuple(maximum_patterns),
        violating,
        first_violation,
    )


def audit_uniform_inequality_prefix(
    controls: Sequence[Golden],
) -> tuple[UniformInequalityLevelAudit, ...]:
    return tuple(audit_uniform_inequality_level(control) for control in controls)


def maximum_rearranged_pairing_margin(
    later_control: Golden,
    later_spatial_pattern: Sequence[Golden],
    earlier_control: Golden,
) -> Golden:
    """Upper-bound pairing/2 with every root in an earlier complete shell.

    Rearrangement maximizes the Euclidean dot product over all labelings, so
    a nonpositive result is an exact, label-independent Vinberg acuteness
    certificate.  It may overestimate because not every permutation is a
    Golay labeling, which is harmless for this use.
    """

    if len(later_spatial_pattern) != 24:
        raise ValueError("a spatial pattern must have 24 coordinates")
    _raw_count, earlier_patterns = _complete_golay_patterns(earlier_control)
    later = tuple(sorted(later_spatial_pattern, reverse=True))
    return max(
        (
            sum(
                (left * right for left, right in zip(later, sorted(pattern, reverse=True))),
                start=ZERO,
            )
            - PHI * later_control * earlier_control
            for pattern in earlier_patterns
        ),
        default=Golden(-10**9),
    )


@lru_cache(maxsize=None)
def audit_maximum_slope_level(control: Golden) -> MaximumSlopeLevelAudit:
    """Compute ``max sum(d_i)/x`` by a profile-valued exact dynamic program.

    The state retains only the greatest coordinate sum for each of the 49
    exact joint Golay residue profiles.  This is substantially cheaper than
    retaining every unordered square partition at large controls.
    """

    if not control.is_integral or control.sign(0) <= 0:
        raise ValueError("a control must be a positive golden integer")
    target = ONE + PHI * control * control
    doubled_values = _coefficient_box(4 * target)
    positive = tuple(
        sorted(
            {
                value / 2
                for value in doubled_values
                if value.sign(0) > 0
                and all((4 * target - value * value).sign(place) >= 0 for place in (0, 1))
            },
            reverse=True,
        )
    )
    squares = tuple(value * value for value in positive)
    residues = tuple(2 * rational + radical for rational, radical in map(_coordinate_residue, positive))
    allowed = golay_joint_residue_counts()

    # Each record is (joint profile, maximum sum, one maximizing multiset).
    @lru_cache(maxsize=None)
    def solve(
        index: int, slots: int, remainder: Golden
    ) -> tuple[tuple[tuple[int, int, int, int], Golden, tuple[Golden, ...]], ...]:
        if any(remainder.sign(place) < 0 for place in (0, 1)):
            return ()
        if index == len(positive):
            if remainder:
                return ()
            return (((slots, 0, 0, 0), ZERO, (ZERO,) * slots),)
        best: dict[tuple[int, int, int, int], tuple[Golden, tuple[Golden, ...]]] = {}
        value = positive[index]
        square = squares[index]
        residue = residues[index]
        running = remainder
        for multiplicity in range(slots + 1):
            if all(running.sign(place) >= 0 for place in (0, 1)):
                for profile, suffix_sum, suffix_pattern in solve(
                    index + 1, slots - multiplicity, running
                ):
                    counts = list(profile)
                    counts[residue] += multiplicity
                    joint = tuple(counts)
                    if not any(
                        all(joint[cell] <= possible[cell] for cell in range(4))
                        for possible in allowed
                    ):
                        continue
                    spatial_sum = multiplicity * value + suffix_sum
                    old = best.get(joint)
                    if old is None or (spatial_sum - old[0]).sign(0) > 0:
                        best[joint] = (
                            spatial_sum,
                            (value,) * multiplicity + suffix_pattern,
                        )
            running -= square
        return tuple((profile, total, pattern) for profile, (total, pattern) in best.items())

    records = solve(0, 24, target)
    maximum = None
    patterns: list[tuple[Golden, ...]] = []
    for profile, spatial_sum, pattern in records:
        if profile not in allowed:
            raise AssertionError("the dynamic program returned a non-Golay profile")
        if maximum is None or (spatial_sum - maximum).sign(0) > 0:
            maximum = spatial_sum
            patterns = [pattern]
        elif spatial_sum == maximum:
            patterns.append(pattern)
    slope = None if maximum is None else maximum / control
    beats = maximum is not None and (maximum * maximum - 24 * PHI * control * control).sign(0) > 0
    return MaximumSlopeLevelAudit(
        control,
        PHI * control * control,
        target,
        len(positive),
        solve.cache_info().currsize,
        len(records),
        maximum,
        slope,
        tuple(patterns),
        beats,
    )


@lru_cache(maxsize=1)
def explicit_first_uniform_cutter() -> ExplicitUniformCutterAudit:
    """Construct an accepted H<30 root that cuts the uniform ray.

    Take a Golay hexadecad ``H`` (the complement of an octad), put
    ``phi^2/2`` on H, put ``phi`` at one point of the complementary octad,
    and zero at its other seven points.  Both coefficient residue words are
    H, so Construction A is automatic.  Rearrangement bounds against the
    three complete ``x=phi^2`` multisets and the single ``x=2phi`` multiset
    prove nonpositive pairing with every earlier root, independently of its
    Golay labeling.
    """

    octad = octads()[0]
    hexadecad = ((1 << 24) - 1) ^ octad
    special = next(index for index in range(24) if octad & (1 << index))
    spatial = tuple(
        PHI
        if index == special
        else (ONE + PHI) / 2
        if hexadecad & (1 << index)
        else ZERO
        for index in range(24)
    )
    control = ONE + 2 * PHI
    root = (control, *spatial)
    norm = a1_golay_pair(root, root)
    spatial_sum = sum(spatial, start=ZERO)
    uniform_pairing = a1_golay_pair(uniform_ray_audit().vector, root)
    height_zero_ok = all(
        a1_golay_pair(root, wall).sign(0) <= 0 for wall in height_zero_walls()
    )
    x1_ok = all(a1_golay_pair(root, wall.vector).sign(0) <= 0 for wall in x1_packet())
    xphi_ok = all(
        a1_golay_pair(root, wall.vector).sign(0) <= 0 for wall in xphi_packet()
    )

    earlier_patterns = (
        # The complete x=phi^2 multiset census.
        ((ONE + PHI) / 2,) * 7 + ((-ONE + PHI) / 2,),
        (ONE + PHI,) + (PHI / 2,) * 8,
        (ONE + PHI,) + (PHI,) * 2,
        # The complete x=2phi multiset census.
        (ONE + 2 * PHI,),
    )
    earlier_controls = (ONE + PHI, ONE + PHI, ONE + PHI, 2 * PHI)
    descending = tuple(sorted(spatial, reverse=True))
    margins = []
    for earlier_control, pattern in zip(earlier_controls, earlier_patterns):
        padded = tuple(sorted((*pattern, *(ZERO for _ in range(24 - len(pattern)))), reverse=True))
        maximum_dot = sum(
            (left * right for left, right in zip(descending, padded)), start=ZERO
        )
        # Pairing/2 is dot - phi*x*x_earlier.
        margins.append(maximum_dot - PHI * control * earlier_control)
    accepted = (
        norm == 2
        and golay_membership(spatial)
        and height_zero_ok
        and x1_ok
        and xphi_ok
        and all(margin.sign(0) <= 0 for margin in margins)
    )
    intersections = tuple(
        sorted(
            {
                (left & right).bit_count()
                for index, left in enumerate(octads())
                for right in octads()[index + 1 :]
            }
        )
    )
    if intersections != (0, 2, 4):
        raise AssertionError("Golay octads have an unexpected intersection size")
    # For two distinct labels (O,j),(P,k), write H=O^c,K=P^c.
    # Their dot product depends only on
    # (|H intersect K|, [j in K]+[k in H], [j=k]).
    cases = (
        (8, 2, 0),
        (10, 0, 0), (10, 0, 1), (10, 1, 0), (10, 2, 0),
        (12, 0, 0), (12, 0, 1), (12, 1, 0), (12, 2, 0),
        (16, 0, 0),
    )
    half_square = (ONE + PHI) / 2
    same_height_margins = tuple(
        (
            case,
            half_square * half_square * case[0]
            + half_square * PHI * case[1]
            + PHI * PHI * case[2]
            - PHI * control * control,
        )
        for case in cases
    )
    return ExplicitUniformCutterAudit(
        root,
        control,
        PHI * control * control,
        norm,
        spatial_sum,
        uniform_pairing,
        golay_membership(spatial),
        height_zero_ok,
        x1_ok,
        xphi_ok,
        tuple(margins),
        accepted,
        len(octads()) * 8,
        intersections,
        same_height_margins,
        all(margin.sign(0) <= 0 for _case, margin in same_height_margins),
    )


@lru_cache(maxsize=1)
def replacement_uniform_ray_audit() -> ReplacementUniformRayAudit:
    """Certify the new extreme hyperideal ray after ``x=1+2phi``.

    The complete unordered shell census has the unique maximum-sum pattern
    ``phi^1, ((1+phi)/2)^16, 0^7``.  Its Golay labelings are precisely a
    hexadecad and a distinguished point in the complementary octad, giving
    ``759*8`` active roots.  Their spatial rows have exact Q(phi)-rank 24.
    """

    cutter = explicit_first_uniform_cutter()
    shell = audit_uniform_inequality_level(cutter.control)
    expected_pattern = tuple(
        sorted(
            (PHI, *((ONE + PHI) / 2 for _ in range(16)), *(ZERO for _ in range(7))),
            reverse=True,
        )
    )
    if shell.maximum_spatial_sum != cutter.spatial_sum:
        raise AssertionError("the explicit cutter does not maximize its complete shell")
    if shell.maximum_patterns != (expected_pattern,):
        raise AssertionError("the first cutter shell has an unexpected maximizing multiset")

    coordinate = PHI * cutter.control / cutter.spatial_sum
    vector = (ONE, *(coordinate for _ in range(24)))
    norm = a1_golay_pair(vector, vector)
    active_pairings_zero = all(
        not a1_golay_pair(vector, (cutter.control, *spatial))
        for spatial in _first_cutter_spatials()
    )

    # At all earlier positive shells the old coordinate phi/4 was feasible.
    # Every dominant spatial coefficient is nonnegative and the replacement
    # coordinate is smaller, so those pairings only decrease.  We retain an
    # independent complete-shell check at the two intervening controls.
    intervening = tuple(
        audit_uniform_inequality_level(control) for control in (ONE + PHI, 2 * PHI)
    )
    earlier_feasible = (
        (coordinate - PHI / 4).sign(0) < 0
        and all(a1_golay_pair(vector, wall).sign(0) <= 0 for wall in height_zero_walls())
        and all(a1_golay_pair(vector, root.vector).sign(0) <= 0 for root in x1_packet())
        and all(a1_golay_pair(vector, root.vector).sign(0) <= 0 for root in xphi_packet())
        and all(
            audit.maximum_spatial_sum is not None
            and (
                -2 * PHI * audit.control
                + 2 * coordinate * audit.maximum_spatial_sum
            ).sign(0)
            <= 0
            for audit in intervening
        )
    )
    active_shell_feasible = shell.maximum_spatial_sum is not None and (
        -2 * PHI * shell.control + 2 * coordinate * shell.maximum_spatial_sum
    ).sign(0) <= 0
    rank = first_cutter_spatial_rank()
    return ReplacementUniformRayAudit(
        vector,
        coordinate,
        norm,
        cutter.control,
        cutter.spatial_sum,
        cutter.labeled_family_count,
        rank,
        active_pairings_zero,
        earlier_feasible,
        active_shell_feasible,
        norm.sign(0) > 0,
        rank == 24 and active_pairings_zero,
        cutter.spatial_sum / cutter.control,
    )


def first_uniform_ray_cutter(
    roots: Sequence[Vector],
) -> PersistenceViolation | None:
    for index, root in enumerate(roots):
        if len(root) != 25:
            raise ValueError("a candidate root must have rank 25")
        spatial_sum = sum(root[1:], start=ZERO)
        excess = spatial_sum - 4 * root[0]
        if excess.sign(0) > 0:
            return PersistenceViolation(index, root[0], spatial_sum, excess, root)
    return None


def build_first_two_packet_audit() -> dict[str, object]:
    uniform = uniform_ray_audit()
    cutter = explicit_first_uniform_cutter()
    cutter_shell = audit_uniform_inequality_level(cutter.control)
    replacement = replacement_uniform_ray_audit()
    roots = (*x1_packet(), *xphi_packet())
    if any(a1_golay_pair(root.vector, root.vector) != 2 for root in roots):
        raise AssertionError("first packets contain a non-norm-two vector")
    if any(
        a1_golay_pair(left.vector, right.vector).sign(0) > 0
        for index, left in enumerate(xphi_packet())
        for right in xphi_packet()[index + 1 :]
    ):
        raise AssertionError("x=phi packet is not mutually acute")
    return {
        "schema_version": "0.1.0",
        "lattice": "<-2phi> + N(A1^24, extended Golay glue)",
        "canonical_mirror_norms": ["2"],
        "complete_norm_class_proof": COMPLETE_MIRROR_NORM_PROOF,
        "height_zero_wall_count": len(height_zero_walls()),
        "x1_packet": {"single_orbit_count": len(x1_packet())},
        "xphi_packet": {
            "pair_orbit_count": sum(root.orbit == "pair" for root in xphi_packet()),
            "octad_orbit_count": sum(root.orbit == "octad" for root in xphi_packet()),
            "total_count": len(xphi_packet()),
        },
        "uniform_ray": {
            "vector": [str(value) for value in uniform.vector],
            "norm": str(uniform.norm),
            "octad_incidence_rank": uniform.octad_incidence_rank,
            "octad_wall_count": uniform.octad_wall_count,
            "all_octad_pairings_zero": uniform.all_octad_pairings_zero,
            "x1_pairings_nonpositive": uniform.x1_pairings_nonpositive,
            "xphi_pairings_nonpositive": uniform.xphi_pairings_nonpositive,
            "hyperideal": uniform.hyperideal,
            "extreme": uniform.extreme,
            "later_root_condition": "sum(d_i) <= 4*x",
        },
        "first_uniform_ray_cutter": {
            "control": str(cutter.control),
            "height": str(cutter.height),
            "spatial_pattern": [
                "phi^1",
                "((1+phi)/2)^16",
                "0^7",
            ],
            "spatial_sum": str(cutter.spatial_sum),
            "norm": str(cutter.norm),
            "complete_shell_pattern_count": cutter_shell.golay_admissible_pattern_count,
            "complete_shell_maximum_pattern_count": len(cutter_shell.maximum_patterns),
            "complete_shell_violating_pattern_count": cutter_shell.violating_pattern_count,
            "labeled_orbit_count": cutter.labeled_family_count,
            "accepted_through_all_earlier_complete_shells": cutter.accepted_through_earlier_levels,
            "labeled_orbit_mutually_acute": cutter.labeled_family_mutually_acute,
        },
        "replacement_uniform_ray": {
            "vector": [str(value) for value in replacement.vector],
            "spatial_coordinate": str(replacement.spatial_coordinate),
            "norm": str(replacement.norm),
            "active_orbit_count": replacement.active_orbit_count,
            "active_spatial_rank": replacement.active_spatial_rank,
            "all_active_pairings_zero": replacement.all_active_pairings_zero,
            "all_earlier_shells_feasible": replacement.all_earlier_shells_feasible,
            "full_active_shell_feasible": replacement.full_active_shell_feasible,
            "hyperideal": replacement.hyperideal,
            "extreme": replacement.extreme,
            "later_root_condition": (
                "sum(d_i) <= (" + str(replacement.later_root_slope) + ")*x"
            ),
        },
        "scope": (
            "The packets through x=1+2phi and the replacement hyperideal "
            "extreme ray are exact. Persistence against later roots is a "
            "separate inequality and is not claimed without a complete "
            "later-level audit."
        ),
    }


__all__ = [
    "COMPLETE_MIRROR_NORM_PROOF",
    "a1_golay_pair",
    "audit_maximum_slope_level",
    "audit_uniform_inequality_level",
    "audit_uniform_inequality_prefix",
    "build_first_two_packet_audit",
    "explicit_first_uniform_cutter",
    "first_cutter_spatial_rank",
    "first_uniform_ray_cutter",
    "golay_membership",
    "golay_joint_residue_counts",
    "height_zero_walls",
    "maximum_rearranged_pairing_margin",
    "octad_incidence_rank",
    "octads",
    "replacement_uniform_ray_audit",
    "uniform_ray_audit",
    "x1_packet",
    "xphi_packet",
]
