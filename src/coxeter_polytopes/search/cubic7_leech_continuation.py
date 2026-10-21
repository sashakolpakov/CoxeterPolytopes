"""Exact later prefixes for cubic-seven Leech-spatial hosts in ``H24``.

The original :mod:`coxeter_polytopes.search.cubic7_leech` audit stopped at
height ten, before the first nonempty packet of ``<-2*theta> + Lambda_24``.
This module continues that *actual* host.  It also puts the four small
coefficient directions supported at the inert prime two and the ramified
prime ``delta=2-theta`` under one local-arithmetic convention.

There are two deliberately different kinds of result here.

* For ``a=2*theta`` the first two nonempty packets are realised and accepted
  completely.  They are described by the minimal Leech shell and by ordered
  pairs of minimal vectors of inner product ``-1``.
* For every later control in the stated prefix, it is enough to enumerate
  necessary coefficient Gram matrices.  Adjoining one fixed norm-six Leech
  direction and checking the exact integral minimum of the resulting Gram
  quotient gives an upper bound for its pairing with *every possible wall*.
  A strict hyperideal ray surviving this larger necessary set certainly
  survives the actual Vinberg walls.

No Conway triple-orbit enumeration is assumed.  In particular, a necessary
Gram type is never silently promoted to a realised wall.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import isqrt
from typing import Sequence

from coxeter_polytopes.exact.cubic7 import CUBIC7_ONE, THETA7, Cubic7
from coxeter_polytopes.search.cubic7_lambda22 import (
    Lambda22GramType,
    lambda22_integer_gram_minimum_at_least,
    lambda22_rank3_gram_types,
)
from coxeter_polytopes.search.cubic7_leech import (
    LEECH_MINIMAL_VECTOR_COUNT,
    extended_binary_golay_codewords,
    iter_minimal_leech_vectors_scaled,
    minimal_leech_shell_audit,
)
from coxeter_polytopes.search.cubic7_root_blocks import _bounded_for_caps
from coxeter_polytopes.search.silver_leech import (
    P7_CUTTER_P_SCALED,
    P7_CUTTER_Q_SCALED,
    integer_determinant,
)


DELTA7 = 2 - THETA7
SCALAR_S = 2 * THETA7 + THETA7**2
LEECH_MINIMUM = 4
NORM_SIX_SHELL_COUNT = 16_773_120


LOCAL_MIRROR_NORM_PROOF = (
    "For primitive r=(x,y), q(r) divides 2I(r), where "
    "I(r)=(a*x,(y,Lambda_24)).  The Leech lattice is even unimodular. "
    "Thus no prime away from (2) and the squarefree support of a can divide "
    "q.  At an odd supported prime P, the valuation is zero or one; the "
    "positive case is equivalent to y in P*Lambda_24, and then "
    "q=-a*x^2 mod P^2 has valuation exactly one.  At inert (2), an odd "
    "coefficient permits valuations zero or one, while v_2(a)=1 forces "
    "valuation exactly one by evenness and primitivity.  Class number one "
    "and U^+=U^2 leave the displayed canonical totally-positive norms."
)


@dataclass(frozen=True, slots=True)
class Cubic7LeechSeed:
    seed_id: str
    timelike_coefficient: Cubic7
    canonical_mirror_norms: tuple[Cubic7, ...]
    ideal_support: str
    complete_mirror_norm_proof: str = LOCAL_MIRROR_NORM_PROOF

    @property
    def admissible(self) -> bool:
        return self.timelike_coefficient.sign(0) > 0 and all(
            self.timelike_coefficient.sign(place) < 0 for place in (1, 2)
        )


UNIT_LEECH_SEED = Cubic7LeechSeed(
    "cubic7-H24-Leech-a-theta",
    THETA7,
    (Cubic7(1), Cubic7(2)),
    "unit",
)
TWO_THETA_LEECH_SEED = Cubic7LeechSeed(
    "cubic7-H24-Leech-a-2theta",
    2 * THETA7,
    (Cubic7(2),),
    "(2)",
)
P7_LEECH_SEED = Cubic7LeechSeed(
    "cubic7-H24-Leech-a-theta-delta",
    THETA7 * DELTA7,
    (Cubic7(1), Cubic7(2), DELTA7, 2 * DELTA7),
    "(delta), delta=2-theta",
)
TWO_P7_LEECH_SEED = Cubic7LeechSeed(
    "cubic7-H24-Leech-a-2theta-delta",
    2 * THETA7 * DELTA7,
    (Cubic7(2), 2 * DELTA7),
    "(2)(delta), delta=2-theta",
)
SUPPORTED_LEECH_SEEDS = (
    UNIT_LEECH_SEED,
    TWO_THETA_LEECH_SEED,
    P7_LEECH_SEED,
    TWO_P7_LEECH_SEED,
)


@dataclass(frozen=True, slots=True, order=True)
class Cubic7LeechControlLevel:
    height: Cubic7
    norm: Cubic7
    control: Cubic7
    target_before_gate: Cubic7


@lru_cache(maxsize=None)
def cubic7_leech_control_levels(
    seed: Cubic7LeechSeed, height_bound: int
) -> tuple[Cubic7LeechControlLevel, ...]:
    """Exhaust every locally compatible control through exact height ``H``."""

    if seed not in SUPPORTED_LEECH_SEEDS or not seed.admissible:
        raise ValueError("unsupported or inadmissible cubic-seven Leech seed")
    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    coefficient = seed.timelike_coefficient
    levels: set[Cubic7LeechControlLevel] = set()
    for norm in seed.canonical_mirror_norms:
        caps = (
            norm * height_bound / coefficient,
            norm / (-coefficient),
            norm / (-coefficient),
        )
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if any(target.sign(place) < 0 for place in (1, 2)):
                continue
            height = coefficient * control * control / norm
            if (Cubic7(height_bound) - height).sign(0) < 0:
                continue
            levels.add(Cubic7LeechControlLevel(height, norm, control, target))
    return tuple(sorted(levels))


@dataclass(frozen=True, slots=True)
class Cubic7LeechSpatialGate:
    spatial_divisor: Cubic7
    target_after_gate: Cubic7 | None
    reason: str


def cubic7_leech_spatial_gate(
    level: Cubic7LeechControlLevel,
) -> Cubic7LeechSpatialGate:
    """Apply the exact ramified-prime spatial divisibility condition."""

    if level.norm not in (DELTA7, 2 * DELTA7):
        return Cubic7LeechSpatialGate(CUBIC7_ONE, level.target_before_gate, "none")
    reduced = level.target_before_gate / (DELTA7 * DELTA7)
    if not reduced.is_integral:
        return Cubic7LeechSpatialGate(
            DELTA7,
            None,
            "y in delta*Lambda_24, but y^2/delta^2 is nonintegral",
        )
    return Cubic7LeechSpatialGate(
        DELTA7,
        reduced,
        "y in delta*Lambda_24; enumerate z=y/delta",
    )


@lru_cache(maxsize=None)
def cubic7_leech_necessary_gram_types(
    level: Cubic7LeechControlLevel,
) -> tuple[Lambda22GramType, ...]:
    """Return every minimum-four coefficient Gram type for one control.

    ``lambda22_rank3_gram_types`` is lattice-independent at this stage: it
    solves the three coefficient equations and checks only evenness,
    positive semidefiniteness, and exact integral minimum at least four.
    Those conditions are necessary for a triple of Leech vectors.
    """

    gate = cubic7_leech_spatial_gate(level)
    if gate.target_after_gate is None:
        return ()
    return lambda22_rank3_gram_types(gate.target_after_gate, 4, 1)


@dataclass(frozen=True, slots=True)
class Cubic7LeechHeightZeroAudit:
    seed_id: str
    norm: Cubic7
    spatial_divisor: Cubic7
    target_after_gate: Cubic7 | None
    necessary_gram_type_count: int


def cubic7_leech_height_zero_audit(
    seed: Cubic7LeechSeed,
) -> tuple[Cubic7LeechHeightZeroAudit, ...]:
    """Prove the height-zero census empty in every complete norm class."""

    rows = []
    for norm in seed.canonical_mirror_norms:
        level = Cubic7LeechControlLevel(Cubic7(0), norm, Cubic7(0), norm)
        gate = cubic7_leech_spatial_gate(level)
        types = cubic7_leech_necessary_gram_types(level)
        if types:
            raise AssertionError("a small cubic-seven Leech seed gained a height-zero type")
        rows.append(
            Cubic7LeechHeightZeroAudit(
                seed.seed_id,
                norm,
                gate.spatial_divisor,
                gate.target_after_gate,
                0,
            )
        )
    return tuple(rows)


def _scaled_leech_dot(left: Sequence[int], right: Sequence[int]) -> int:
    if len(left) != 24 or len(right) != 24:
        raise ValueError("scaled Leech coordinates must have length 24")
    numerator = sum(a * b for a, b in zip(left, right))
    if numerator % 8:
        raise AssertionError("scaled Leech pairing ceased to be integral")
    return numerator // 8


NORM_SIX_DIRECTION_SCALED = tuple(
    left + right
    for left, right in zip(P7_CUTTER_P_SCALED, P7_CUTTER_Q_SCALED)
)


@lru_cache(maxsize=1)
def norm_six_direction_audit() -> dict[str, object]:
    """Check the explicit catalogue-compatible norm-six Leech direction."""

    support = sum(
        1 << index
        for index, coordinate in enumerate(P7_CUTTER_Q_SCALED)
        if coordinate
    )
    if support.bit_count() != 8 or support not in extended_binary_golay_codewords():
        raise AssertionError("the norm-six direction lost its Golay octad")
    if _scaled_leech_dot(P7_CUTTER_P_SCALED, P7_CUTTER_P_SCALED) != 4:
        raise AssertionError("the type-(4^2) summand is not minimal")
    if _scaled_leech_dot(P7_CUTTER_Q_SCALED, P7_CUTTER_Q_SCALED) != 4:
        raise AssertionError("the octad summand is not minimal")
    if _scaled_leech_dot(P7_CUTTER_P_SCALED, P7_CUTTER_Q_SCALED) != -1:
        raise AssertionError("the two minimal summands lost pairing -1")
    square = _scaled_leech_dot(NORM_SIX_DIRECTION_SCALED, NORM_SIX_DIRECTION_SCALED)
    if square != 6:
        raise AssertionError("the fixed direction lost norm six")
    return {
        "scaled_coordinates": NORM_SIX_DIRECTION_SCALED,
        "coordinate_scale": "sqrt(8)",
        "square": square,
        "construction": "sum of two minimal Leech vectors of inner product -1",
    }


def _rank_three_determinant(
    a: int, b: int, c: int, ab: int, ac: int, bc: int
) -> int:
    return (
        a * b * c
        + 2 * ab * ac * bc
        - a * bc * bc
        - b * ac * ac
        - c * ab * ab
    )


def _augmented_gram(
    gram_type: Lambda22GramType, pairings: tuple[int, int, int]
) -> tuple[tuple[int, ...], ...]:
    base = gram_type.matrix
    p, q, r = pairings
    return (
        (6, p, q, r),
        (p, base[0][0], base[0][1], base[0][2]),
        (q, base[1][0], base[1][1], base[1][2]),
        (r, base[2][0], base[2][1], base[2][2]),
    )


def _augmented_gram_is_positive_semidefinite(
    gram_type: Lambda22GramType, pairings: tuple[int, int, int]
) -> bool:
    """Check every new principal minor after adjoining the norm-six row."""

    base = gram_type.matrix
    p, q, r = pairings
    diagonals = (base[0][0], base[1][1], base[2][2])
    if any(value * value > 6 * diagonal for value, diagonal in zip(pairings, diagonals)):
        return False
    for i, j in combinations(range(3), 2):
        if _rank_three_determinant(
            6,
            base[i][i],
            base[j][j],
            pairings[i],
            pairings[j],
            base[i][j],
        ) < 0:
            return False
    return integer_determinant(_augmented_gram(gram_type, pairings)) >= 0


@dataclass(frozen=True, slots=True)
class NormSixPairingBound:
    maximum: Cubic7
    maximizing_pairings: tuple[int, int, int]
    positive_semidefinite_extension_count: int
    higher_extensions_rejected_by_minimum: int


@lru_cache(maxsize=None)
def norm_six_pairing_bound(gram_type: Lambda22GramType) -> NormSixPairingBound:
    """Maximise ``d.(u+theta*v+theta^2*w)`` exactly.

    All integral pairings allowed by Cauchy and positive semidefiniteness are
    sorted in the distinguished real embedding.  The first augmented Gram
    whose integral kernel quotient has minimum at least four is the exact
    maximum over this necessary superset.  This handles singular Gram types
    without choosing a noncanonical rational complement.
    """

    base = gram_type.matrix
    bounds = tuple(isqrt(6 * base[index][index]) for index in range(3))
    candidates: list[tuple[Cubic7, tuple[int, int, int]]] = []
    for pairings in product(*(range(-bound, bound + 1) for bound in bounds)):
        if not _augmented_gram_is_positive_semidefinite(gram_type, pairings):
            continue
        candidates.append((Cubic7(*pairings), pairings))
    candidates.sort(reverse=True, key=lambda item: item[0])
    for rejected, (value, pairings) in enumerate(candidates):
        if lambda22_integer_gram_minimum_at_least(
            _augmented_gram(gram_type, pairings), LEECH_MINIMUM
        ):
            return NormSixPairingBound(value, pairings, len(candidates), rejected)
    raise AssertionError("the zero-pairing norm-six extension disappeared")


@dataclass(frozen=True, slots=True)
class Cubic7LeechNecessaryLevel:
    level: Cubic7LeechControlLevel
    spatial_divisor: Cubic7
    target_after_gate: Cubic7
    gram_types: tuple[Lambda22GramType, ...]
    normalized_pairing_maximum: Cubic7
    actual_pairing_maximum: Cubic7
    maximizing_pairings: tuple[int, int, int]
    witness_margin: Cubic7


@dataclass(frozen=True, slots=True)
class Cubic7LeechPrefixObstruction:
    seed: Cubic7LeechSeed
    height_bound: int
    control_level_count: int
    empty_necessary_level_count: int
    nonempty_necessary_levels: tuple[Cubic7LeechNecessaryLevel, ...]
    witness_coordinate: Fraction
    witness_spatial_square: Fraction
    witness_lorentzian_norm: Cubic7
    minimum_wall_margin: Cubic7 | None

    @property
    def strict_hyperideal(self) -> bool:
        return self.witness_lorentzian_norm.sign(0) > 0 and (
            self.minimum_wall_margin is None
            or self.minimum_wall_margin.sign(0) > 0
        )


@lru_cache(maxsize=None)
def cubic7_leech_prefix_obstruction(
    seed: Cubic7LeechSeed,
    height_bound: int,
    witness_coordinate: Fraction,
) -> Cubic7LeechPrefixObstruction:
    """Bound a fixed norm-six ray against every possible prefix wall."""

    norm_six_direction_audit()
    levels = cubic7_leech_control_levels(seed, height_bound)
    rows = []
    empty = 0
    scale = Cubic7(witness_coordinate)
    for level in levels:
        gate = cubic7_leech_spatial_gate(level)
        gram_types = cubic7_leech_necessary_gram_types(level)
        if not gram_types:
            empty += 1
            continue
        bounds = tuple(norm_six_pairing_bound(gram) for gram in gram_types)
        maximum_bound = max(bounds, key=lambda bound: bound.maximum)
        actual_maximum = gate.spatial_divisor * maximum_bound.maximum
        margin = (
            seed.timelike_coefficient * level.control
            - scale * actual_maximum
        )
        if margin.sign(0) <= 0:
            raise AssertionError(
                f"the fixed ray fails at {seed.seed_id}, height {level.height}"
            )
        if gate.target_after_gate is None:  # impossible when gram_types is nonempty
            raise AssertionError("a nonempty level has no normalized target")
        rows.append(
            Cubic7LeechNecessaryLevel(
                level,
                gate.spatial_divisor,
                gate.target_after_gate,
                gram_types,
                maximum_bound.maximum,
                actual_maximum,
                maximum_bound.maximizing_pairings,
                margin,
            )
        )
    lorentzian_norm = Cubic7(6 * witness_coordinate * witness_coordinate) - seed.timelike_coefficient
    if lorentzian_norm.sign(0) <= 0:
        raise AssertionError("the fixed prefix witness is not hyperideal")
    minimum_margin = min((row.witness_margin for row in rows), default=None)
    result = Cubic7LeechPrefixObstruction(
        seed,
        height_bound,
        len(levels),
        empty,
        tuple(rows),
        witness_coordinate,
        6 * witness_coordinate * witness_coordinate,
        lorentzian_norm,
        minimum_margin,
    )
    if not result.strict_hyperideal:
        raise AssertionError("the prefix obstruction lost strictness")
    return result


@lru_cache(maxsize=1)
def minimal_leech_inner_product_distribution() -> dict[int, int]:
    """Inner-product distribution of the minimal shell about one vector."""

    fixed = (-4, -4) + (0,) * 22
    distribution: dict[int, int] = {}
    for _kind, vector in iter_minimal_leech_vectors_scaled():
        pairing = _scaled_leech_dot(fixed, vector)
        distribution[pairing] = distribution.get(pairing, 0) + 1
    expected = {-4: 1, -2: 4600, -1: 47104, 0: 93150, 1: 47104, 2: 4600, 4: 1}
    if distribution != expected:
        raise AssertionError("the minimal-shell inner-product distribution drifted")
    return dict(sorted(distribution.items()))


@lru_cache(maxsize=1)
def _second_packet_same_height_audit() -> tuple[int, Cubic7]:
    """Maximise the pairing of two distinct second-packet roots.

    A root is indexed by an ordered minimal pair ``(D,E)`` with ``D.E=-1``.
    Four cross pairings in ``[-4,4]`` determine the joint Gram.  Exhausting
    all ``9^4`` possibilities, exact PSD, and exact minimum four is a
    complete necessary comparison; its maximum is already negative.
    """

    s = SCALAR_S
    control = Cubic7(1, 4, 2)
    feasible = 0
    maximum: Cubic7 | None = None
    maximizing_cross_data: tuple[int, int, int, int] | None = None
    identical = (4, -1, -1, 4)
    for dd, de, ed, ee in product(range(-4, 5), repeat=4):
        if (dd, de, ed, ee) == identical:
            continue
        gram = (
            (4, -1, dd, de),
            (-1, 4, ed, ee),
            (dd, ed, 4, -1),
            (de, ee, -1, 4),
        )
        if any(
            integer_determinant(
                tuple(tuple(gram[i][j] for j in indices) for i in indices)
            ) < 0
            for size in range(1, 5)
            for indices in combinations(range(4), size)
        ):
            continue
        if not lambda22_integer_gram_minimum_at_least(gram, LEECH_MINIMUM):
            continue
        feasible += 1
        spatial_pairing = (
            (s + 1) * (s + 1) * dd
            - s * (s + 1) * (de + ed)
            + s * s * ee
        )
        root_pairing = -2 * THETA7 * control * control + spatial_pairing
        if maximum is None or root_pairing > maximum:
            maximum = root_pairing
            maximizing_cross_data = (dd, de, ed, ee)
    # The maximum is attained in the actual packet: for any (D,E), take the
    # second ordered pair (D',E')=(-E,-D).  Its cross data are exactly the
    # tuple below.
    if (
        feasible != 494
        or maximum != Cubic7(-1)
        or maximizing_cross_data != (1, -4, -4, 1)
    ):
        raise AssertionError("the second-packet same-height audit drifted")
    return feasible, maximum


@dataclass(frozen=True, slots=True)
class TwoThetaFirstTwoPackets:
    first_height: Cubic7
    first_control: Cubic7
    first_target: Cubic7
    first_gram_type: Lambda22GramType
    first_wall_count: int
    first_same_height_maximum: Cubic7
    second_height: Cubic7
    second_control: Cubic7
    second_target: Cubic7
    second_gram_type: Lambda22GramType
    minimal_pair_inner_product: int
    minimal_pairs_per_fixed_vector: int
    second_wall_count: int
    first_second_maximum: Cubic7
    second_same_height_cross_type_count: int
    second_same_height_maximum: Cubic7


@lru_cache(maxsize=1)
def two_theta_first_two_packets() -> TwoThetaFirstTwoPackets:
    """Certify the first two realised simultaneous Vinberg packets."""

    levels = cubic7_leech_control_levels(TWO_THETA_LEECH_SEED, 150)
    nonempty = tuple(
        (level, cubic7_leech_necessary_gram_types(level))
        for level in levels
        if cubic7_leech_necessary_gram_types(level)
    )
    if len(nonempty) < 2:
        raise AssertionError("the first two necessary packets disappeared")
    (first, first_types), (second, second_types) = nonempty[:2]
    if len(first_types) != 1 or len(second_types) != 1:
        raise AssertionError("an early packet gained another Gram type")
    first_type = first_types[0]
    second_type = second_types[0]
    if (
        first.height,
        first.control,
        first.target_before_gate,
        (
            first_type.u_square,
            first_type.v_square,
            first_type.w_square,
            first_type.u_dot_v,
            first_type.u_dot_w,
            first_type.v_dot_w,
        ),
    ) != (
        Cubic7(5, 14, 6),
        Cubic7(1, 2, 1),
        Cubic7(12, 28, 12),
        (0, 16, 4, 0, 0, 8),
    ):
        raise AssertionError("the first realised packet changed")
    if (
        second.height,
        second.control,
        second.target_before_gate,
        (
            second_type.u_square,
            second_type.v_square,
            second_type.w_square,
            second_type.u_dot_v,
            second_type.u_dot_w,
            second_type.v_dot_w,
        ),
    ) != (
        Cubic7(16, 45, 20),
        Cubic7(1, 4, 2),
        Cubic7(34, 90, 40),
        (4, 40, 10, 10, 5, 20),
    ):
        raise AssertionError("the second realised packet changed")

    # First packet: y=s*M with M^2=4.  The coefficient s is a unit and the
    # control is a unit, so all roots are primitive and reflective.
    if SCALAR_S.norm() != 1 or first.control.norm() != 1:
        raise AssertionError("the first-packet unit proof drifted")
    if (
        -TWO_THETA_LEECH_SEED.timelike_coefficient * first.control * first.control
        + 4 * SCALAR_S * SCALAR_S
        != Cubic7(2)
    ):
        raise AssertionError("the first-packet roots lost norm two")
    first_same_height_maximum = (
        -2 * THETA7 * first.control * first.control + 2 * SCALAR_S * SCALAR_S
    )
    if first_same_height_maximum != Cubic7(-4, -14, -6):
        raise AssertionError("the first packet lost acuteness")

    # Second packet: v=2w, put D=-u and E=D+w.  Then D and E are minimal
    # and D.E=-1.  Conversely every such ordered pair gives the Gram type.
    distribution = minimal_leech_inner_product_distribution()
    pairs_per_fixed = distribution[-1]
    second_count = LEECH_MINIMAL_VECTOR_COUNT * pairs_per_fixed
    if second_count != 9_258_762_240:
        raise AssertionError("the second-packet incidence count drifted")
    if second.control.norm() not in (-1, 1):
        raise AssertionError("the second-packet control ceased to be a unit")

    # Against the first packet, |B.M|<=5 and |D.M|<=4.  Equality is attained
    # by M=-D, so this is an exact maximum, not merely a bound.
    first_second_maximum = (
        -2 * THETA7 * first.control * second.control
        + SCALAR_S * (5 * SCALAR_S + 4)
    )
    if first_second_maximum != -(SCALAR_S * SCALAR_S):
        raise AssertionError("the second packet is not acute to the first")
    cross_count, second_maximum = _second_packet_same_height_audit()
    if second_maximum.sign(0) >= 0:
        raise AssertionError("the second equal-height packet is not acute")

    # The archived minimal shell is one Co0 orbit.  This is the only group
    # transitivity used in multiplying the fixed-vector incidence count.
    if minimal_leech_shell_audit()["co0_orbit_count"] != 1:
        raise AssertionError("the minimal shell lost its one-orbit certificate")
    return TwoThetaFirstTwoPackets(
        first.height,
        first.control,
        first.target_before_gate,
        first_type,
        LEECH_MINIMAL_VECTOR_COUNT,
        first_same_height_maximum,
        second.height,
        second.control,
        second.target_before_gate,
        second_type,
        -1,
        pairs_per_fixed,
        second_count,
        first_second_maximum,
        cross_count,
        second_maximum,
    )


SMALL_COEFFICIENT_PREFIXES = (
    (UNIT_LEECH_SEED, 300, Fraction(7, 15)),
    (P7_LEECH_SEED, 300, Fraction(1, 2)),
    (TWO_P7_LEECH_SEED, 300, Fraction(3, 5)),
)


def _cubic_json(value: Cubic7 | None) -> dict[str, str] | None:
    return None if value is None else value.to_json()


def _gram_json(gram: Lambda22GramType) -> dict[str, object]:
    return {
        "matrix": gram.matrix,
        "determinant": gram.determinant,
        "rank": gram.rank,
    }


def _prefix_json(prefix: Cubic7LeechPrefixObstruction) -> dict[str, object]:
    return {
        "seed_id": prefix.seed.seed_id,
        "timelike_coefficient": prefix.seed.timelike_coefficient.to_json(),
        "height_bound": prefix.height_bound,
        "control_level_count": prefix.control_level_count,
        "empty_necessary_level_count": prefix.empty_necessary_level_count,
        "nonempty_necessary_level_count": len(prefix.nonempty_necessary_levels),
        "total_necessary_gram_type_count": sum(
            len(row.gram_types) for row in prefix.nonempty_necessary_levels
        ),
        "witness": {
            "direction_square": 6,
            "coordinate": str(prefix.witness_coordinate),
            "spatial_square": str(prefix.witness_spatial_square),
            "lorentzian_norm": prefix.witness_lorentzian_norm.to_json(),
            "minimum_wall_margin": _cubic_json(prefix.minimum_wall_margin),
            "strict_hyperideal": prefix.strict_hyperideal,
        },
        "nonempty_levels": [
            {
                "height": row.level.height.to_json(),
                "norm": row.level.norm.to_json(),
                "control": row.level.control.to_json(),
                "target_before_gate": row.level.target_before_gate.to_json(),
                "spatial_divisor": row.spatial_divisor.to_json(),
                "target_after_gate": row.target_after_gate.to_json(),
                "necessary_gram_types": [_gram_json(gram) for gram in row.gram_types],
                "normalized_pairing_maximum": row.normalized_pairing_maximum.to_json(),
                "actual_pairing_maximum": row.actual_pairing_maximum.to_json(),
                "maximizing_pairings": row.maximizing_pairings,
                "witness_margin": row.witness_margin.to_json(),
            }
            for row in prefix.nonempty_necessary_levels
        ],
    }


def cubic7_leech_h24_continuation_payload() -> dict[str, object]:
    """Build the complete frozen H1000/H300 continuation record."""

    packets = two_theta_first_two_packets()
    deep = cubic7_leech_prefix_obstruction(
        TWO_THETA_LEECH_SEED, 1000, Fraction(2, 3)
    )
    comparisons = tuple(
        cubic7_leech_prefix_obstruction(seed, height, coordinate)
        for seed, height, coordinate in SMALL_COEFFICIENT_PREFIXES
    )
    return {
        "schema": "coxeter.cubic7_leech_h24_continuation.v1",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "hyperbolic_dimension": 24,
        "spatial_lattice": "Lambda_24",
        "local_arithmetic": {
            "delta": DELTA7.to_json(),
            "proof": LOCAL_MIRROR_NORM_PROOF,
            "seeds": [
                {
                    "seed_id": seed.seed_id,
                    "timelike_coefficient": seed.timelike_coefficient.to_json(),
                    "ideal_support": seed.ideal_support,
                    "canonical_mirror_norms": [
                        norm.to_json() for norm in seed.canonical_mirror_norms
                    ],
                    "admissible": seed.admissible,
                    "height_zero": [
                        {
                            "norm": row.norm.to_json(),
                            "spatial_divisor": row.spatial_divisor.to_json(),
                            "target_after_gate": _cubic_json(row.target_after_gate),
                            "necessary_gram_type_count": row.necessary_gram_type_count,
                        }
                        for row in cubic7_leech_height_zero_audit(seed)
                    ],
                }
                for seed in SUPPORTED_LEECH_SEEDS
            ],
        },
        "norm_six_direction": norm_six_direction_audit(),
        "two_theta_actual_packets": {
            "first": {
                "height": packets.first_height.to_json(),
                "control": packets.first_control.to_json(),
                "target": packets.first_target.to_json(),
                "gram_type": _gram_json(packets.first_gram_type),
                "description": "y=(2theta+theta^2)M, M^2=4",
                "wall_count": packets.first_wall_count,
                "same_height_pairing_maximum": packets.first_same_height_maximum.to_json(),
            },
            "second": {
                "height": packets.second_height.to_json(),
                "control": packets.second_control.to_json(),
                "target": packets.second_target.to_json(),
                "gram_type": _gram_json(packets.second_gram_type),
                "description": (
                    "y=(2theta+theta^2)(E-D)-D for ordered minimal "
                    "D,E with D.E=-1"
                ),
                "minimal_pair_inner_product": packets.minimal_pair_inner_product,
                "minimal_pairs_per_fixed_vector": packets.minimal_pairs_per_fixed_vector,
                "wall_count": packets.second_wall_count,
                "first_second_pairing_maximum": packets.first_second_maximum.to_json(),
                "same_height_cross_type_count": packets.second_same_height_cross_type_count,
                "same_height_pairing_maximum": packets.second_same_height_maximum.to_json(),
            },
            "minimal_shell_inner_product_distribution": minimal_leech_inner_product_distribution(),
            "complete_acceptance": True,
        },
        "two_theta_h1000": _prefix_json(deep),
        "small_coefficient_h300": [_prefix_json(prefix) for prefix in comparisons],
        "conclusion": {
            "new_compact_example": False,
            "new_finite_volume_example": False,
            "two_theta_status": (
                "two complete accepted packets; strict hyperideal point survives "
                "every possible wall through H=1000"
            ),
            "scope": (
                "Exact only for the four displayed Leech-spatial coefficient "
                "directions, complete local norm/gates, the stated H300/H1000 "
                "control boxes, and the first two accepted a=2theta packets. "
                "Later walls remain unresolved; this is not nonreflectivity "
                "or field-wide nonexistence."
            ),
        },
    }


__all__ = [
    "DELTA7",
    "LOCAL_MIRROR_NORM_PROOF",
    "NORM_SIX_DIRECTION_SCALED",
    "P7_LEECH_SEED",
    "SCALAR_S",
    "SMALL_COEFFICIENT_PREFIXES",
    "SUPPORTED_LEECH_SEEDS",
    "TWO_P7_LEECH_SEED",
    "TWO_THETA_LEECH_SEED",
    "UNIT_LEECH_SEED",
    "Cubic7LeechControlLevel",
    "Cubic7LeechHeightZeroAudit",
    "Cubic7LeechNecessaryLevel",
    "Cubic7LeechPrefixObstruction",
    "Cubic7LeechSeed",
    "NormSixPairingBound",
    "TwoThetaFirstTwoPackets",
    "cubic7_leech_control_levels",
    "cubic7_leech_h24_continuation_payload",
    "cubic7_leech_height_zero_audit",
    "cubic7_leech_necessary_gram_types",
    "cubic7_leech_prefix_obstruction",
    "cubic7_leech_spatial_gate",
    "minimal_leech_inner_product_distribution",
    "norm_six_direction_audit",
    "norm_six_pairing_bound",
    "two_theta_first_two_packets",
]
