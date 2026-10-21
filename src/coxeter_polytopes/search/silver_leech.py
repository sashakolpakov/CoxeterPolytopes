"""Exact early Vinberg packets for silver Leech-spatial lattices.

The primary seed is ``<-epsilon> + Lambda_24`` over ``Z[sqrt(2)]``, where
``epsilon=1+sqrt(2)`` and ``Lambda_24`` is the even unimodular Leech lattice.
There are no height-zero walls.  Its first nonempty positive packet is the
complete norm-ten Leech shell, without ever materialising its 4.6 billion
vectors.  The packet is bounded in the affine section but an explicit
minimal-Leech-vector direction leaves an open hyperideal patch, so boundedness
must not be confused with finite hyperbolic volume.

All finite shell exclusions below use ``y=u+sqrt(2)v`` with ``u,v`` in the
Leech lattice.  If ``y^2=A+B sqrt(2)``, then

``u^2+2v^2=A`` and ``u.v=B/2``.

Evenness, minimum four, Cauchy--Schwarz, and short integral combinations of
``u,v`` give a complete finite list of possible Gram triples.  This is enough
for the first packet and the next nonempty unit packet; no Conway-orbit shell
expansion is needed.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import isqrt
from typing import Sequence

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.silver_a1_golay import (
    SQRT2_NORM_PROOF,
    TWO_SQRT2_NORM_PROOF,
    UNIT_NORM_PROOF,
)
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


EPSILON = SILVER_ONE + SQRT2
LEECH_MINIMUM = 4
LEECH_SHELL_COUNTS = {
    4: 196_560,
    6: 16_773_120,
    8: 398_034_000,
    10: 4_629_381_120,
}


@dataclass(frozen=True, slots=True)
class SilverLeechSeed:
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


UNIT_LEECH_SEED = SilverLeechSeed(
    "silver-H24-Leech-a-1+sqrt2",
    EPSILON,
    (Silver(1), Silver(2)),
    UNIT_NORM_PROOF,
)
SQRT2_LEECH_SEED = SilverLeechSeed(
    "silver-H24-Leech-a-sqrt2",
    SQRT2,
    (Silver(2), Silver(2, 1)),
    SQRT2_NORM_PROOF,
)
TWO_SQRT2_LEECH_SEED = SilverLeechSeed(
    "silver-H24-Leech-a-2sqrt2",
    2 * SQRT2,
    (Silver(2), Silver(4, 2)),
    TWO_SQRT2_NORM_PROOF,
)
P7_COEFFICIENT = 2 * SQRT2 - SILVER_ONE
P7_GENERATOR = 3 + SQRT2
P7_NORM_PROOF = (
    "Put pi=sqrt(2) and rho=2*pi-1, whose field norm is -7. For a "
    "primitive r=(x,y), q(r)|2I(r), I(r)=(rho*x,(y,Lambda)). Away "
    "from (pi) and the prime p=(rho) above 7 singled out by the "
    "coefficient, unimodularity and "
    "primitivity make I(r)=O. At (pi), evenness and the same primitive "
    "split give v_pi(q)=0 or 2. At p, either y is p-primitive and "
    "v_p(q)=0, or y lies in p*Lambda, primitive x is a p-unit, and "
    "q=-rho*x^2 mod p^2 gives v_p(q)=1. Thus the norm ideals are "
    "O,(2),p,2p. Since every totally positive unit is a square, their "
    "canonical totally-positive generators are 1,2,3+pi,6+2*pi. "
    "For the p-divisible classes, reflection integrality additionally "
    "requires y in p*Lambda."
)
P7_LEECH_SEED = SilverLeechSeed(
    "silver-H24-Leech-a-2sqrt2-1",
    P7_COEFFICIENT,
    (Silver(1), Silver(2), P7_GENERATOR, 2 * P7_GENERATOR),
    P7_NORM_PROOF,
)

# Table 2a, row 1, of Dutour Sikiric--Schuermann--Vallentin gives this
# exceptional contact-polar vertex as one third of the following raw MOG
# row.  MOG coordinates are divided by sqrt(8), so the row has Euclidean
# square 192 and the displayed polar vertex has Leech square 192/(8*9)=8/3.
P7_EXCEPTIONAL_MOG_RAW: tuple[int, ...] = (
    2, 2, 2, 2, -2, -6, 2, -2, -2, 6, 2, 2,
    -2, 2, -2, 2, 2, 2, 2, 2, -2, 2, 2, 6,
)

# An explicit norm-six direction and an H=31 cutter in the same standard
# scaled-sqrt(8) model.  P and C are type-(4^2) minimal vectors.  Q is the
# all-positive vector on a Golay octad.  Put D=P+Q, V=D+C and U=D+V.
P7_CUTTER_P_SCALED: tuple[int, ...] = (-4, -4) + (0,) * 22
P7_CUTTER_Q_SCALED: tuple[int, ...] = (
    0, 2, 2, 0, 0, 0, 2, 2, 2, 0, 2, 0,
    2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
)
P7_CUTTER_C_SCALED: tuple[int, ...] = (-4, 0, 0, -4) + (0,) * 20
SUPPORTED_SILVER_LEECH_SEEDS = (
    UNIT_LEECH_SEED,
    SQRT2_LEECH_SEED,
    TWO_SQRT2_LEECH_SEED,
    P7_LEECH_SEED,
)


@dataclass(frozen=True, slots=True)
class SilverLeechControlLevel:
    norm: Silver
    height: Silver
    control: Silver
    spatial_square_target: Silver


@lru_cache(maxsize=None)
def enumerate_control_levels(
    seed: SilverLeechSeed, height_bound: int
) -> tuple[SilverLeechControlLevel, ...]:
    """Enumerate every arithmetically possible control below ``height_bound``."""

    if seed not in SUPPORTED_SILVER_LEECH_SEEDS or not seed.admissible:
        raise ValueError("unsupported or inadmissible silver Leech seed")
    coefficient = seed.timelike_coefficient
    levels: set[SilverLeechControlLevel] = set()
    for norm in seed.canonical_mirror_norms:
        for control in _bounded_for_caps(
            (norm * height_bound / coefficient, norm / (-coefficient))
        ):
            if control.sign(0) <= 0:
                continue
            # Vinberg crystallographic height requires 2*a*x/q in O_K.
            if not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if target.sign(1) < 0:
                continue
            height = coefficient * control * control / norm
            if (Silver(height_bound) - height).sign(0) < 0:
                continue
            levels.add(SilverLeechControlLevel(norm, height, control, target))
    return tuple(
        sorted(levels, key=lambda level: (level.height, level.norm, level.control))
    )


@dataclass(frozen=True, slots=True)
class LeechGramTriple:
    u_square: int
    v_square: int
    u_dot_v: int

    def combination_square(self, u_coefficient: int, v_coefficient: int) -> int:
        return (
            u_coefficient * u_coefficient * self.u_square
            + v_coefficient * v_coefficient * self.v_square
            + 2 * u_coefficient * v_coefficient * self.u_dot_v
        )


def candidate_leech_gram_triples(target: Silver) -> tuple[LeechGramTriple, ...]:
    """Return all triples surviving evenness, minimum four, and Cauchy.

    This is only a necessary census.  Callers must still exclude impossible
    triples using norms of integral combinations, or construct them.
    """

    if not target.is_integral or target.b.denominator != 1:
        return ()
    dot = target.b / 2
    if dot.denominator != 1:
        return ()
    rational_target = int(target.a)
    dot_int = int(dot)
    output: list[LeechGramTriple] = []
    for v_square in range(0, rational_target // 2 + 1, 2):
        u_square = rational_target - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if u_square not in (0,) and u_square < LEECH_MINIMUM:
            continue
        if v_square not in (0,) and v_square < LEECH_MINIMUM:
            continue
        if dot_int * dot_int > u_square * v_square:
            continue
        output.append(LeechGramTriple(u_square, v_square, dot_int))
    return tuple(output)


def reflective_candidate_leech_gram_triples(
    seed: SilverLeechSeed,
    norm: Silver,
    target: Silver,
) -> tuple[LeechGramTriple, ...]:
    """Apply the seed's local reflection-integrality divisibility gate.

    For ``a=2sqrt(2)=pi^3`` and the ``pi^3`` norm class, reflectivity
    requires ``y in pi*Lambda``.  Writing ``y=pi*z`` replaces the square
    target by ``target/2`` before the ordinary Leech Gram census.  Omitting
    this gate creates a spurious minimal-shell packet.
    """

    if seed == TWO_SQRT2_LEECH_SEED and norm == Silver(4, 2):
        target = target / 2
    if seed == P7_LEECH_SEED and norm in (
        P7_GENERATOR,
        2 * P7_GENERATOR,
    ):
        # Here p=(rho) and p^2=(7).  Writing y=rho*z replaces y^2 by
        # rho^2*z^2, not merely by the rational field norm seven.
        target = target / (P7_COEFFICIENT * P7_COEFFICIENT)
    return candidate_leech_gram_triples(target)


def has_forbidden_norm_two_combination(triple: LeechGramTriple) -> bool:
    """Detect the short combinations sufficient for all archived exclusions."""

    return any(
        triple.combination_square(a, b) == 2
        for a, b in ((1, -1), (1, -2), (2, -1), (1, 1), (1, 2), (2, 1))
    )


@dataclass(frozen=True, slots=True)
class UnitLeechFirstPacketAudit:
    complete_through_height: Silver
    earlier_control_count: int
    earlier_candidate_triples: tuple[tuple[Silver, tuple[LeechGramTriple, ...]], ...]
    first_norm: Silver
    first_height: Silver
    first_control: Silver
    first_target: Silver
    surviving_triple: LeechGramTriple
    shell_norm: int
    wall_count: int
    same_height_pairing_bound: Silver
    span_witness_determinant: int
    affine_section_bounded: bool
    strict_witness_coordinate: Silver
    strict_witness_spatial_square: Silver
    boundary_witness_coordinate: Silver
    boundary_witness_spatial_square: Silver
    timelike_threshold: Silver
    hyperideal_open_patch: bool
    finite_volume: bool
    compact: bool


# Coordinates in the exact integral MOG Leech basis used by GAP/ALCO.  Every
# row has Gram square ten; the 24x24 coordinate determinant is two.
MOG_NORM10_SPAN_WITNESS: tuple[tuple[int, ...], ...] = (
    (0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
    (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
    (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
    (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
    (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
    (0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1),
    (0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
    (0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
    (0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
    (0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
    (0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
    (0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
    (0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
    (0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0),
    (0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0),
    (0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0),
    (0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
    (0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1),
    (0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
    (0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
    (0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0),
    (-1,1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0),
)


def integer_determinant(matrix: Sequence[Sequence[int]]) -> int:
    """Fraction-free Bareiss determinant for a small integral certificate."""

    work = [list(map(int, row)) for row in matrix]
    size = len(work)
    if any(len(row) != size for row in work):
        raise ValueError("matrix must be square")
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column]
        for row in range(column + 1, size):
            for index in range(column + 1, size):
                numerator = (
                    work[row][index] * pivot_value
                    - work[row][column] * work[column][index]
                )
                if numerator % previous:
                    raise AssertionError("Bareiss division ceased to be exact")
                work[row][index] = numerator // previous
        previous = pivot_value
    return sign * work[-1][-1]


def unit_leech_first_packet_audit() -> UnitLeechFirstPacketAudit:
    levels = enumerate_control_levels(UNIT_LEECH_SEED, 29)
    first_height = Silver(14, 10)
    first = next(
        level
        for level in levels
        if level.height == first_height and level.norm == 2
    )
    earlier = tuple(level for level in levels if level.height < first_height)
    earlier_triples = tuple(
        (level.height, candidate_leech_gram_triples(level.spatial_square_target))
        for level in earlier
    )
    for _height, triples in earlier_triples:
        if any(not has_forbidden_norm_two_combination(triple) for triple in triples):
            raise AssertionError("an unexplained earlier Leech Gram triple survived")

    triples = candidate_leech_gram_triples(first.spatial_square_target)
    survivors = tuple(
        triple for triple in triples if not has_forbidden_norm_two_combination(triple)
    )
    expected = LeechGramTriple(10, 10, 10)
    if survivors != (expected,):
        raise AssertionError("the unit-Leech first shell ceased to be scalar")

    determinant = integer_determinant(MOG_NORM10_SPAN_WITNESS)
    if abs(determinant) != 2:
        raise AssertionError("the archived MOG norm-ten span witness drifted")

    strict_coordinate = Silver(Fraction(4, 5))
    boundary_coordinate = 2 * EPSILON / 5
    strict_square = 4 * strict_coordinate * strict_coordinate
    boundary_square = 4 * boundary_coordinate * boundary_coordinate
    return UnitLeechFirstPacketAudit(
        complete_through_height=first_height,
        earlier_control_count=len(earlier),
        earlier_candidate_triples=earlier_triples,
        first_norm=first.norm,
        first_height=first.height,
        first_control=first.control,
        first_target=first.spatial_square_target,
        surviving_triple=expected,
        shell_norm=10,
        wall_count=LEECH_SHELL_COUNTS[10],
        # For distinct b,c, (b-c)^2>=4 gives b.c<=8.
        same_height_pairing_bound=EPSILON * EPSILON * (8 - 4 * EPSILON),
        span_witness_determinant=determinant,
        affine_section_bounded=True,
        # If u^2=4 and b^2=10, then |u.b|<=5.  The strict witness
        # z=4u/5 has |z.b|<=4<2epsilon and is already spacelike.
        strict_witness_coordinate=strict_coordinate,
        strict_witness_spatial_square=strict_square,
        boundary_witness_coordinate=boundary_coordinate,
        boundary_witness_spatial_square=boundary_square,
        timelike_threshold=EPSILON,
        hyperideal_open_patch=strict_square > EPSILON,
        finite_volume=False,
        compact=False,
    )


@dataclass(frozen=True, slots=True)
class UnitLeechSecondPacketAudit:
    height: Silver
    norm: Silver
    control: Silver
    target: Silver
    candidate_triples: tuple[LeechGramTriple, ...]
    surviving_triple: LeechGramTriple
    decomposition: str
    first_packet_pairing_bound_status: str
    old_boundary_witness_is_cut: bool
    replacement_witness_coordinate: Silver
    replacement_witness_spatial_square: Silver
    replacement_first_packet_margin: Silver
    replacement_hyperideal: bool


def unit_leech_second_packet_audit() -> UnitLeechSecondPacketAudit:
    """Audit the next nonempty shell after the norm-ten packet.

    The q=1 level at the same height has no Gram triple.  For q=2 the only
    realizable triple is ``(30,14,20)``.  Writing ``w=u-v`` turns it into
    ``w^2=4, v^2=14, w.v=6`` and ``y=w+epsilon*v``.
    """

    height = Silver(28, 20)
    levels = tuple(
        level
        for level in enumerate_control_levels(UNIT_LEECH_SEED, 57)
        if level.height == height
    )
    if tuple(level.norm for level in levels) != (Silver(1), Silver(2)):
        raise AssertionError("the second unit-Leech height packet drifted")
    if candidate_leech_gram_triples(levels[0].spatial_square_target):
        raise AssertionError("the q=1 second-height shell should be empty")
    level = levels[1]
    candidates = candidate_leech_gram_triples(level.spatial_square_target)
    survivors = tuple(
        triple for triple in candidates if not has_forbidden_norm_two_combination(triple)
    )
    expected = LeechGramTriple(30, 14, 20)
    if survivors != (expected,):
        raise AssertionError("the unit-Leech second shell Gram census drifted")

    control = level.control
    old_cut_margin = 2 * EPSILON * (16 * EPSILON / 5 - control)
    replacement = EPSILON * control / (4 + 7 * EPSILON)
    replacement_square = 4 * replacement * replacement
    return UnitLeechSecondPacketAudit(
        height=height,
        norm=level.norm,
        control=control,
        target=level.spatial_square_target,
        candidate_triples=candidates,
        surviving_triple=expected,
        decomposition="y=w+epsilon*v, w^2=4, v^2=14, w.v=6",
        first_packet_pairing_bound_status=(
            "survivors are the Gram-(30,14,20) roots satisfying "
            "y.b <= 2*epsilon*x for every norm-ten b"
        ),
        old_boundary_witness_is_cut=old_cut_margin.sign(0) > 0,
        replacement_witness_coordinate=replacement,
        replacement_witness_spatial_square=replacement_square,
        replacement_first_packet_margin=2 * EPSILON - 5 * replacement,
        replacement_hyperideal=replacement_square > EPSILON,
    )


@dataclass(frozen=True, slots=True)
class MinimalDirectionBound:
    height: Silver
    control: Silver
    necessary_triple_count: int
    maximum_pairing: Silver
    wall_margin: Silver


def _minimum_direction_pairing_bound(triple: LeechGramTriple) -> Silver:
    """Maximise ``u0.y`` under exact elementary Leech constraints.

    Here ``u0^2=4`` is fixed and ``y=u+sqrt(2)v``.  We choose
    ``w=u-kv`` with minimal square, then enumerate the integral pairings
    ``a=u0.w`` and ``b=u0.v``.  The exact 3x3 Gram determinant must be
    nonnegative, and each of ``u0+-w,u0+-v`` is either zero or has norm at
    least four.  These are necessary constraints for an actual Leech pair,
    so their maximum is a rigorous upper bound even if some enumerated
    pairing data are not realised.
    """

    transformations = []
    for coefficient in range(-8, 9):
        w_square = triple.combination_square(1, -coefficient)
        if w_square < 0:
            continue
        w_dot_v = triple.u_dot_v - coefficient * triple.v_square
        transformations.append(
            (w_square, abs(w_dot_v), abs(coefficient), coefficient, w_dot_v)
        )
    if not transformations:
        raise AssertionError("no nonnegative transformed square")
    w_square, _abs_dot, _abs_coefficient, coefficient, w_dot_v = min(
        transformations
    )
    a_bound = isqrt(4 * w_square)
    b_bound = isqrt(4 * triple.v_square)
    best: Silver | None = None
    for a in range(-a_bound, a_bound + 1):
        for b in range(-b_bound, b_bound + 1):
            if a * a > 4 * w_square or b * b > 4 * triple.v_square:
                continue
            determinant = (
                4 * (w_square * triple.v_square - w_dot_v * w_dot_v)
                - a * (a * triple.v_square - b * w_dot_v)
                + b * (a * w_dot_v - b * w_square)
            )
            if determinant < 0:
                continue
            pair_norms = (
                4 + w_square + 2 * a,
                4 + w_square - 2 * a,
                4 + triple.v_square + 2 * b,
                4 + triple.v_square - 2 * b,
            )
            if any(value != 0 and value < 4 for value in pair_norms):
                continue
            value = Silver(a + coefficient * b, b)
            if best is None or value > best:
                best = value
    if best is None:
        raise AssertionError("no elementary pairing data survived")
    return best


@lru_cache(maxsize=1)
def unit_leech_minimal_direction_prefix_h82() -> tuple[MinimalDirectionBound, ...]:
    """Certify one hyperideal ray through every possible wall up to H=82."""

    witness = Silver(4, 10) / 23
    output: list[MinimalDirectionBound] = []
    # The endpoint 82+58sqrt(2) is about 164.02 at the identity embedding.
    for level in enumerate_control_levels(UNIT_LEECH_SEED, 165):
        if level.height <= Silver(14, 10):
            continue
        triples = tuple(
            triple
            for triple in reflective_candidate_leech_gram_triples(
                UNIT_LEECH_SEED, level.norm, level.spatial_square_target
            )
            if not has_forbidden_norm_two_combination(triple)
        )
        if not triples:
            continue
        maximum = max(_minimum_direction_pairing_bound(triple) for triple in triples)
        margin = EPSILON * level.control - witness * maximum
        if margin.sign(0) < 0:
            raise AssertionError("the H82 minimal-direction certificate failed")
        output.append(
            MinimalDirectionBound(
                height=level.height,
                control=level.control,
                necessary_triple_count=len(triples),
                maximum_pairing=maximum,
                wall_margin=margin,
            )
        )
    expected = (
        (Silver(28, 20), Silver(10, 6)),
        (Silver(41, 29), Silver(11, 8)),
        (Silver(56, 40), Silver(14, 10)),
        (Silver(63, 45), Silver(15, 11)),
        (Silver(82, 58), Silver(17, 12)),
    )
    if tuple((row.height, row.maximum_pairing) for row in output) != expected:
        raise AssertionError("the complete H82 nonempty-level census drifted")
    return tuple(output)


@dataclass(frozen=True, slots=True)
class NonunitLeechFirstPacketAudit:
    seed_id: str
    earlier_control_count: int
    first_height: Silver
    first_norm: Silver
    first_control: Silver
    first_target: Silver
    scalar_shell_norm: int
    wall_count: int
    slab_bound: Silver
    witness_description: str
    witness_spatial_square: Silver
    timelike_threshold: Silver
    strict_hyperideal_witness: bool


@dataclass(frozen=True, slots=True)
class NonunitLeechSecondPacketAudit:
    seed_id: str
    complete_height: Silver
    norm: Silver
    control: Silver
    target: Silver
    surviving_triple: LeechGramTriple
    decomposition: str
    witness_coordinate: Silver
    first_packet_margin: Silver
    second_packet_margin: Silver
    witness_spatial_square: Silver
    timelike_threshold: Silver
    hyperideal_prefix_witness: bool


@dataclass(frozen=True, slots=True)
class P7LeechFirstPacketAudit:
    complete_norm_class_proof: str
    complete_through_height: Silver
    earlier_control_count: int
    earlier_surviving_triple_count: int
    first_norm: Silver
    first_height: Silver
    first_control: Silver
    first_target: Silver
    surviving_triple: LeechGramTriple
    scalar_shell_norm: int
    wall_count: int
    normalized_slab_bound: Silver
    contact_polar_reference_bound: int
    contact_polar_vertex_orbit_count: int
    exceptional_mog_raw: tuple[int, ...]
    exceptional_mog_raw_square: int
    exceptional_mog_divisor: int
    exceptional_mog_coordinate_scale_square: int
    exceptional_contact_polar_square: Fraction
    exceptional_incidence: int
    exceptional_stabilizer: str
    scaled_exceptional_square: Silver
    timelike_threshold: Silver
    exceptional_lorentzian_norm: Silver
    hyperideal_open_patch: bool
    finite_volume: bool
    compact: bool


@dataclass(frozen=True, slots=True)
class P7LeechPrefixBound:
    height: Silver
    norm: Silver
    control: Silver
    necessary_triples: tuple[LeechGramTriple, ...]
    maximum_pairing: Silver
    wall_margin: Silver


@dataclass(frozen=True, slots=True)
class P7LeechSecondPacketAudit:
    complete_through_height: Silver
    norm: Silver
    control: Silver
    target: Silver
    candidate_triples: tuple[LeechGramTriple, ...]
    surviving_triples: tuple[LeechGramTriple, ...]
    explicit_scaled_sqrt8_w: tuple[int, ...]
    explicit_scaled_sqrt8_v: tuple[int, ...]
    explicit_scaled_sqrt8_u: tuple[int, ...]
    explicit_realized_triple: LeechGramTriple
    old_witness_is_cut: bool
    replacement_direction_square: int
    replacement_witness_coordinate: Silver
    first_packet_margin: Silver
    replacement_spatial_square: Silver
    timelike_threshold: Silver
    replacement_hyperideal: bool
    complete_prefix_bounds: tuple[P7LeechPrefixBound, ...]


@lru_cache(maxsize=1)
def p7_leech_first_packet_audit() -> P7LeechFirstPacketAudit:
    """Exact first-packet geometry for ``a=2sqrt(2)-1``.

    The first packet is the complete minimal Leech shell, so its affine
    section is a scaled polar of the Leech contact polytope.  The exact
    facet classification of Dutour Sikirić--Schürmann--Vallentin has 232
    polar-vertex orbits; its unique maximum-square orbit has square ``8/3``
    and 552 incident minimal-vector facets.  Scaling their reference bound
    two to the wall bound ``4-sqrt(2)`` gives an exact hyperideal orbit.
    """

    endpoint = Silver(5, 4)
    levels = enumerate_control_levels(P7_LEECH_SEED, 11)
    first = next(
        level
        for level in levels
        if level.height == endpoint and level.norm == Silver(2)
    )
    earlier = tuple(level for level in levels if level.height < endpoint)
    earlier_survivors = []
    for level in earlier:
        for triple in reflective_candidate_leech_gram_triples(
            P7_LEECH_SEED, level.norm, level.spatial_square_target
        ):
            if not has_forbidden_norm_two_combination(triple):
                earlier_survivors.append((level, triple))
    if earlier_survivors:
        raise AssertionError("a lower p7-Leech Gram packet survived")

    triples = tuple(
        triple
        for triple in reflective_candidate_leech_gram_triples(
            P7_LEECH_SEED, first.norm, first.spatial_square_target
        )
        if not has_forbidden_norm_two_combination(triple)
    )
    expected = LeechGramTriple(4, 4, 4)
    if triples != (expected,):
        raise AssertionError("the p7-Leech first packet ceased to be minimal")

    # y=epsilon*u and x/epsilon=sqrt(2).  Hence the two oriented roots
    # +/-u give |z.u| <= a*sqrt(2)=4-sqrt(2) in the t=1 section.
    bound = 4 - SQRT2
    raw_square = sum(coordinate * coordinate for coordinate in P7_EXCEPTIONAL_MOG_RAW)
    if len(P7_EXCEPTIONAL_MOG_RAW) != 24 or raw_square != 192:
        raise AssertionError("the archived exceptional MOG row drifted")
    exceptional_reference_square = Fraction(8, 3)
    if Fraction(raw_square, 8 * 3 * 3) != exceptional_reference_square:
        raise AssertionError("the exceptional MOG normalization drifted")
    scaled_square = exceptional_reference_square * (bound / 2) ** 2
    lorentzian_norm = scaled_square - P7_COEFFICIENT
    if lorentzian_norm.sign(0) <= 0:
        raise AssertionError("the exceptional contact-polar orbit is not hyperideal")
    return P7LeechFirstPacketAudit(
        complete_norm_class_proof=P7_NORM_PROOF,
        complete_through_height=endpoint,
        earlier_control_count=len(earlier),
        earlier_surviving_triple_count=0,
        first_norm=first.norm,
        first_height=first.height,
        first_control=first.control,
        first_target=first.spatial_square_target,
        surviving_triple=expected,
        scalar_shell_norm=4,
        wall_count=LEECH_SHELL_COUNTS[4],
        normalized_slab_bound=bound,
        contact_polar_reference_bound=2,
        contact_polar_vertex_orbit_count=232,
        exceptional_mog_raw=P7_EXCEPTIONAL_MOG_RAW,
        exceptional_mog_raw_square=raw_square,
        exceptional_mog_divisor=3,
        exceptional_mog_coordinate_scale_square=8,
        exceptional_contact_polar_square=exceptional_reference_square,
        exceptional_incidence=552,
        exceptional_stabilizer="Co3",
        scaled_exceptional_square=scaled_square,
        timelike_threshold=P7_COEFFICIENT,
        exceptional_lorentzian_norm=lorentzian_norm,
        hyperideal_open_patch=True,
        finite_volume=False,
        compact=False,
    )


def _fixed_norm_six_pairing_bound(triple: LeechGramTriple) -> Silver:
    """Bound ``d.(u+sqrt(2)v)`` for a fixed Leech vector ``d^2=6``.

    Choose ``w=u-kv`` of minimum square.  For the integral pairings
    ``p=d.w`` and ``s=d.v``, positivity of the exact three-vector Gram
    matrix and ``min(Lambda)=4`` for ``d+-w`` and ``d+-v`` are necessary
    conditions.  Maximising over this finite superset therefore gives a
    rigorous upper bound for every actual Leech pair of the given type.
    """

    if triple.v_square:
        quotient = triple.u_dot_v // triple.v_square
        coefficients = range(quotient - 2, quotient + 3)
    else:
        coefficients = (0,)
    transformations = []
    for coefficient in coefficients:
        w_square = triple.combination_square(1, -coefficient)
        if w_square < 0:
            continue
        w_dot_v = triple.u_dot_v - coefficient * triple.v_square
        transformations.append(
            (w_square, abs(w_dot_v), abs(coefficient), coefficient, w_dot_v)
        )
    if not transformations:
        raise AssertionError("no nonnegative transformed square")
    w_square, _abs_dot, _abs_coefficient, coefficient, w_dot_v = min(
        transformations
    )

    direction_square = 6
    p_bound = isqrt(direction_square * w_square)
    s_bound = isqrt(direction_square * triple.v_square)
    best: Silver | None = None
    for p in range(-p_bound, p_bound + 1):
        for s in range(-s_bound, s_bound + 1):
            if p * p > direction_square * w_square:
                continue
            if s * s > direction_square * triple.v_square:
                continue
            determinant = (
                direction_square
                * (w_square * triple.v_square - w_dot_v * w_dot_v)
                - p * (p * triple.v_square - s * w_dot_v)
                + s * (p * w_dot_v - s * w_square)
            )
            if determinant < 0:
                continue
            pair_norms = (
                direction_square + w_square + 2 * p,
                direction_square + w_square - 2 * p,
                direction_square + triple.v_square + 2 * s,
                direction_square + triple.v_square - 2 * s,
            )
            if any(value != 0 and value < LEECH_MINIMUM for value in pair_norms):
                continue
            value = Silver(p + coefficient * s, s)
            if best is None or value > best:
                best = value
    if best is None:
        raise AssertionError("no norm-six pairing data survived")
    return best


def _scaled_mog_dot(left: Sequence[int], right: Sequence[int]) -> int:
    """Leech inner product for integer coordinates divided by ``sqrt(8)``."""

    if len(left) != 24 or len(right) != 24:
        raise ValueError("expected 24 scaled MOG coordinates")
    numerator = sum(a * b for a, b in zip(left, right))
    if numerator % 8:
        raise AssertionError("the archived MOG pairing is not integral")
    return numerator // 8


def _scaled_mog_add(
    left: Sequence[int], right: Sequence[int]
) -> tuple[int, ...]:
    if len(left) != 24 or len(right) != 24:
        raise ValueError("expected 24 scaled MOG coordinates")
    return tuple(a + b for a, b in zip(left, right))


@dataclass(frozen=True, slots=True)
class P7LeechPrefixLevelAudit:
    height: Silver
    norm: Silver
    control: Silver
    target: Silver
    spatial_normalizer: str
    candidate_triples: tuple[LeechGramTriple, ...]
    maximum_norm_six_pairing: Silver
    normalized_wall_bound: Silver
    replacement_margin: Silver


@dataclass(frozen=True, slots=True)
class P7LeechPostFirstAudit:
    complete_height: Silver
    control_level_count: int
    empty_necessary_control_count: int
    nonempty_necessary_levels: tuple[P7LeechPrefixLevelAudit, ...]
    earliest_later_height: Silver
    earliest_later_norm: Silver
    earliest_later_control: Silver
    earliest_later_target: Silver
    earliest_later_triples: tuple[LeechGramTriple, ...]
    cutter_p_scaled: tuple[int, ...]
    cutter_q_scaled: tuple[int, ...]
    cutter_c_scaled: tuple[int, ...]
    cutter_direction_square: int
    cutter_u_v_triple: LeechGramTriple
    cutter_direction_pairing: Silver
    cutter_first_packet_pairing_bound: Silver
    cutter_first_packet_margin: Silver
    exceptional_ray_coordinate: Silver
    exceptional_ray_cut_excess: Silver
    replacement_coordinate: Fraction
    replacement_first_packet_margin: Silver
    replacement_spatial_square: Fraction
    replacement_lorentzian_norm: Silver
    strict_hyperideal_prefix_witness: bool
    finite_volume_through_prefix: bool


@lru_cache(maxsize=1)
def p7_leech_post_first_prefix_h82() -> P7LeechPostFirstAudit:
    """Cut the exceptional vertex, then certify a hyperideal H82 prefix.

    The four recorded rows are every nonempty *necessary* Gram level after
    the first packet through ``82+58sqrt(2)``; all other exact controls are
    empty.  Bounding the larger necessary sets is sufficient for every
    actual or Vinberg-accepted root.

    For a p-divisible norm class the Gram triple describes ``z=y/rho``.
    Dividing the root by ``rho`` makes its section inequality
    ``point.z <= x``.  This normalization is essential: using ``a*x`` here
    would give an invalid, overly generous margin.
    """

    endpoint = Silver(82, 58)
    levels = tuple(
        level
        for level in enumerate_control_levels(P7_LEECH_SEED, 165)
        if Silver(5, 4) < level.height <= endpoint
    )
    replacement = Fraction(37, 67)
    replacement_silver = Silver(replacement)
    rows: list[P7LeechPrefixLevelAudit] = []
    empty = 0
    p_norms = (P7_GENERATOR, 2 * P7_GENERATOR)
    for level in levels:
        triples = tuple(
            triple
            for triple in reflective_candidate_leech_gram_triples(
                P7_LEECH_SEED, level.norm, level.spatial_square_target
            )
            if not has_forbidden_norm_two_combination(triple)
        )
        if not triples:
            empty += 1
            continue
        maximum = max(_fixed_norm_six_pairing_bound(triple) for triple in triples)
        if level.norm in p_norms:
            normalizer = "rho=2*sqrt(2)-1"
            wall_bound = level.control
        else:
            normalizer = "1"
            wall_bound = P7_COEFFICIENT * level.control
        margin = wall_bound - replacement_silver * maximum
        if margin.sign(0) <= 0:
            raise AssertionError("the corrected p7 H82 witness lost feasibility")
        rows.append(
            P7LeechPrefixLevelAudit(
                height=level.height,
                norm=level.norm,
                control=level.control,
                target=level.spatial_square_target,
                spatial_normalizer=normalizer,
                candidate_triples=triples,
                maximum_norm_six_pairing=maximum,
                normalized_wall_bound=wall_bound,
                replacement_margin=margin,
            )
        )

    expected_rows = (
        (Silver(31, 22), Silver(14, 8)),
        (Silver(56, 40), Silver(21, 15)),
        (Silver(62, 44), Silver(17, 13)),
        (Silver(82, 58), Silver(21, 15)),
    )
    if tuple((row.height, row.maximum_norm_six_pairing) for row in rows) != expected_rows:
        raise AssertionError("the complete p7 H82 necessary-level census drifted")
    if len(levels) != 27 or empty != 23:
        raise AssertionError("the p7 H82 control count drifted")

    # Certify the explicit H=31 root.  P,C are type-(4^2) vectors and Q is
    # an even-sign octad vector, hence all lie in Lambda.
    from coxeter_polytopes.search.cubic7_leech import extended_binary_golay_codewords

    q_support = sum(
        1 << index
        for index, coordinate in enumerate(P7_CUTTER_Q_SCALED)
        if coordinate
    )
    if q_support.bit_count() != 8 or q_support not in extended_binary_golay_codewords():
        raise AssertionError("the archived p7 cutter support ceased to be a Golay octad")
    if any(
        _scaled_mog_dot(vector, vector) != 4
        for vector in (P7_CUTTER_P_SCALED, P7_CUTTER_Q_SCALED, P7_CUTTER_C_SCALED)
    ):
        raise AssertionError("an archived p7 cutter generator ceased to be minimal")
    direction = _scaled_mog_add(P7_CUTTER_P_SCALED, P7_CUTTER_Q_SCALED)
    v = _scaled_mog_add(direction, P7_CUTTER_C_SCALED)
    u = _scaled_mog_add(direction, v)
    cutter_triple = LeechGramTriple(
        _scaled_mog_dot(u, u),
        _scaled_mog_dot(v, v),
        _scaled_mog_dot(u, v),
    )
    if cutter_triple != LeechGramTriple(36, 14, 22):
        raise AssertionError("the explicit p7 H31 cutter Gram data drifted")
    direction_square = _scaled_mog_dot(direction, direction)
    direction_pairing = Silver(
        _scaled_mog_dot(direction, u), _scaled_mog_dot(direction, v)
    )
    if direction_square != 6 or direction_pairing != Silver(14, 8):
        raise AssertionError("the explicit p7 cutter direction drifted")

    first_wall_bound_for_cutter = (4 - SQRT2) * Silver(4, 3)
    first_pairing_bound = _minimum_direction_pairing_bound(cutter_triple)
    first_margin = first_wall_bound_for_cutter - first_pairing_bound
    if first_margin.sign(0) <= 0:
        raise AssertionError("the explicit H31 cutter is rejected by the first packet")

    exceptional_coordinate = (4 - SQRT2) / 3
    cutter_wall_bound = P7_COEFFICIENT * Silver(4, 3)
    cut_excess = exceptional_coordinate * direction_pairing - cutter_wall_bound
    if cut_excess.sign(0) <= 0:
        raise AssertionError("the H31 root ceased to cut the exceptional orbit")

    first_margin_replacement = (4 - SQRT2) - 3 * replacement_silver
    replacement_square = 6 * replacement * replacement
    lorentzian_norm = Silver(replacement_square) - P7_COEFFICIENT
    if first_margin_replacement.sign(0) <= 0 or lorentzian_norm.sign(0) <= 0:
        raise AssertionError("the replacement p7 direction is not strict hyperideal")
    return P7LeechPostFirstAudit(
        complete_height=endpoint,
        control_level_count=len(levels),
        empty_necessary_control_count=empty,
        nonempty_necessary_levels=tuple(rows),
        earliest_later_height=rows[0].height,
        earliest_later_norm=rows[0].norm,
        earliest_later_control=rows[0].control,
        earliest_later_target=rows[0].target,
        earliest_later_triples=rows[0].candidate_triples,
        cutter_p_scaled=P7_CUTTER_P_SCALED,
        cutter_q_scaled=P7_CUTTER_Q_SCALED,
        cutter_c_scaled=P7_CUTTER_C_SCALED,
        cutter_direction_square=direction_square,
        cutter_u_v_triple=cutter_triple,
        cutter_direction_pairing=direction_pairing,
        cutter_first_packet_pairing_bound=first_pairing_bound,
        cutter_first_packet_margin=first_margin,
        exceptional_ray_coordinate=exceptional_coordinate,
        exceptional_ray_cut_excess=cut_excess,
        replacement_coordinate=replacement,
        replacement_first_packet_margin=first_margin_replacement,
        replacement_spatial_square=replacement_square,
        replacement_lorentzian_norm=lorentzian_norm,
        strict_hyperideal_prefix_witness=True,
        finite_volume_through_prefix=False,
    )


def _direction_pairing_bound(
    direction_square: int, triple: LeechGramTriple
) -> Silver:
    """Rigorous pairing bound with a fixed Leech direction.

    This is the square-``direction_square`` analogue of
    :func:`_minimum_direction_pairing_bound`.  The enumeration imposes
    positive semidefiniteness of the exact three-vector Gram matrix and the
    Leech minimum on the four sums and differences with a reduced basis
    ``w=u-kv``.  These are necessary conditions for every actual triple, so
    the maximum is an upper bound even when some pairing data are unrealised.
    """

    transformations = []
    for coefficient in range(-8, 9):
        w_square = triple.combination_square(1, -coefficient)
        if w_square < 0:
            continue
        w_dot_v = triple.u_dot_v - coefficient * triple.v_square
        transformations.append(
            (w_square, abs(w_dot_v), abs(coefficient), coefficient, w_dot_v)
        )
    if not transformations:
        raise AssertionError("no nonnegative transformed square")
    w_square, _abs_dot, _abs_coefficient, coefficient, w_dot_v = min(
        transformations
    )
    a_bound = isqrt(direction_square * w_square)
    b_bound = isqrt(direction_square * triple.v_square)
    best: Silver | None = None
    for a in range(-a_bound, a_bound + 1):
        for b in range(-b_bound, b_bound + 1):
            if (
                a * a > direction_square * w_square
                or b * b > direction_square * triple.v_square
            ):
                continue
            determinant = (
                direction_square
                * (w_square * triple.v_square - w_dot_v * w_dot_v)
                - a * (a * triple.v_square - b * w_dot_v)
                + b * (a * w_dot_v - b * w_square)
            )
            if determinant < 0:
                continue
            pair_norms = (
                direction_square + w_square + 2 * a,
                direction_square + w_square - 2 * a,
                direction_square + triple.v_square + 2 * b,
                direction_square + triple.v_square - 2 * b,
            )
            if any(value != 0 and value < LEECH_MINIMUM for value in pair_norms):
                continue
            value = Silver(a + coefficient * b, b)
            if best is None or value > best:
                best = value
    if best is None:
        raise AssertionError("no elementary pairing data survived")
    return best


@lru_cache(maxsize=1)
def p7_leech_minimal_direction_prefix_h82() -> tuple[P7LeechPrefixBound, ...]:
    """Keep an explicit norm-six direction hyperideal through exact H82."""

    endpoint = Silver(82, 58)
    witness = Silver(37) / 67
    output: list[P7LeechPrefixBound] = []
    # The endpoint is about 164.02 at the identity embedding, so the rational
    # cap 165 includes it without silently truncating the radical endpoint.
    for level in enumerate_control_levels(P7_LEECH_SEED, 165):
        if level.height <= Silver(5, 4) or level.height > endpoint:
            continue
        triples = tuple(
            triple
            for triple in reflective_candidate_leech_gram_triples(
                P7_LEECH_SEED, level.norm, level.spatial_square_target
            )
            if not has_forbidden_norm_two_combination(triple)
        )
        if not triples:
            continue
        maximum = max(_direction_pairing_bound(6, triple) for triple in triples)
        # In the p-divisible classes y=rho*z and q is represented by
        # delta=epsilon*rho (or 2delta).  Cancelling rho in the section wall
        # leaves the threshold x, not rho*x.
        threshold = (
            level.control
            if level.norm in (P7_GENERATOR, 2 * P7_GENERATOR)
            else P7_COEFFICIENT * level.control
        )
        margin = threshold - witness * maximum
        if margin.sign(0) < 0:
            raise AssertionError("the p7-Leech H82 direction certificate failed")
        output.append(
            P7LeechPrefixBound(
                height=level.height,
                norm=level.norm,
                control=level.control,
                necessary_triples=triples,
                maximum_pairing=maximum,
                wall_margin=margin,
            )
        )
    expected = (
        (Silver(31, 22), Silver(2), Silver(14, 8)),
        (Silver(56, 40), 2 * P7_GENERATOR, Silver(21, 15)),
        (Silver(62, 44), Silver(2), Silver(17, 13)),
        (Silver(82, 58), 2 * P7_GENERATOR, Silver(21, 15)),
    )
    if tuple((row.height, row.norm, row.maximum_pairing) for row in output) != expected:
        raise AssertionError("the complete p7-Leech H82 census drifted")
    return tuple(output)


@lru_cache(maxsize=1)
def p7_leech_second_packet_audit() -> P7LeechSecondPacketAudit:
    """Audit the next nonempty packet and replace the first ray witness."""

    height = Silver(31, 22)
    levels = tuple(
        level
        for level in enumerate_control_levels(P7_LEECH_SEED, 63)
        if level.height == height
    )
    if tuple(level.norm for level in levels) != (Silver(1), Silver(2)):
        raise AssertionError("the p7-Leech second height packet drifted")
    q1_candidates = candidate_leech_gram_triples(levels[0].spatial_square_target)
    if tuple(
        triple
        for triple in q1_candidates
        if not has_forbidden_norm_two_combination(triple)
    ):
        raise AssertionError("the p7-Leech second-height q=1 shell survived")
    level = levels[1]
    candidates = candidate_leech_gram_triples(level.spatial_square_target)
    survivors = tuple(
        triple for triple in candidates if not has_forbidden_norm_two_combination(triple)
    )
    expected = (
        LeechGramTriple(36, 14, 22),
        LeechGramTriple(32, 16, 22),
    )
    if survivors != expected:
        raise AssertionError("the p7-Leech second shell Gram census drifted")

    # Integral coordinates scaled by sqrt(8).  The three vectors have true
    # norms w^2=6, v^2=14, u^2=36 and u.v=22, proving the first Gram type is
    # realised.  The construction below explicitly writes them as sums of
    # three standard minimal Leech coordinate vectors p,q,c.
    p = (-4,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    q = (0,2,2,0,0,0,2,2,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,2)
    c = (-4,0,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    if any(
        sum(value * value for value in vector) != 32
        for vector in (p, q, c)
    ):
        raise AssertionError("the archived p7-Leech minimal summands drifted")
    w = tuple(p[index] + q[index] for index in range(24))
    v = tuple(w[index] + c[index] for index in range(24))
    u = tuple(w[index] + v[index] for index in range(24))
    if (
        sum(value * value for value in w),
        sum(value * value for value in v),
        sum(value * value for value in u),
        sum(u[index] * v[index] for index in range(24)),
    ) != (48, 112, 288, 176):
        raise AssertionError("the explicit p7-Leech second-packet vectors drifted")

    old_witness = Silver(3) / 5
    old_margin = P7_COEFFICIENT * level.control - old_witness * Silver(14, 8)
    witness = Silver(37) / 67
    square = Silver(6) * witness * witness
    first_margin = (4 - SQRT2) - Silver(3) * witness
    prefix = p7_leech_minimal_direction_prefix_h82()
    if prefix[0].height != height:
        raise AssertionError("the p7-Leech prefix lost its second packet")
    return P7LeechSecondPacketAudit(
        complete_through_height=height,
        norm=level.norm,
        control=level.control,
        target=level.spatial_square_target,
        candidate_triples=candidates,
        surviving_triples=survivors,
        explicit_scaled_sqrt8_w=w,
        explicit_scaled_sqrt8_v=v,
        explicit_scaled_sqrt8_u=u,
        explicit_realized_triple=expected[0],
        old_witness_is_cut=old_margin.sign(0) < 0,
        replacement_direction_square=6,
        replacement_witness_coordinate=witness,
        first_packet_margin=first_margin,
        replacement_spatial_square=square,
        timelike_threshold=P7_COEFFICIENT,
        replacement_hyperideal=square > P7_COEFFICIENT,
        complete_prefix_bounds=prefix,
    )


def _first_reflective_scalar_level(
    seed: SilverLeechSeed, height_bound: int
) -> tuple[SilverLeechControlLevel, int, int]:
    earlier = 0
    for level in enumerate_control_levels(seed, height_bound):
        triples = tuple(
            triple
            for triple in reflective_candidate_leech_gram_triples(
                seed, level.norm, level.spatial_square_target
            )
            if not has_forbidden_norm_two_combination(triple)
        )
        if not triples:
            earlier += 1
            continue
        if len(triples) != 1:
            raise AssertionError("first nonunit packet is not a scalar shell")
        triple = triples[0]
        if not (triple.u_square == triple.v_square == triple.u_dot_v):
            raise AssertionError("first nonunit packet ceased to be scalar")
        return level, triple.u_square, earlier
    raise AssertionError("no scalar packet in requested range")


def sqrt2_leech_first_packet_audit() -> NonunitLeechFirstPacketAudit:
    level, shell_norm, earlier = _first_reflective_scalar_level(
        SQRT2_LEECH_SEED, 17
    )
    if shell_norm != 6:
        raise AssertionError("sqrt2 Leech first shell drifted")
    # Divide y=epsilon*b out of the root: the section slabs are
    # |z.b|<=2sqrt(2).  For u^2=4,b^2=6, minimum four gives |u.b|<=3.
    coefficient = 2 * SQRT2 / 3
    square = 4 * coefficient * coefficient
    return NonunitLeechFirstPacketAudit(
        seed_id=SQRT2_LEECH_SEED.seed_id,
        earlier_control_count=earlier,
        first_height=level.height,
        first_norm=level.norm,
        first_control=level.control,
        first_target=level.spatial_square_target,
        scalar_shell_norm=shell_norm,
        wall_count=LEECH_SHELL_COUNTS[shell_norm],
        slab_bound=2 * SQRT2,
        witness_description="z=(2sqrt(2)/3)u for any minimal u^2=4",
        witness_spatial_square=square,
        timelike_threshold=SQRT2,
        strict_hyperideal_witness=square > SQRT2,
    )


def two_sqrt2_leech_first_packet_audit() -> NonunitLeechFirstPacketAudit:
    level, shell_norm, earlier = _first_reflective_scalar_level(
        TWO_SQRT2_LEECH_SEED, 17
    )
    if shell_norm != 6 or level.norm != 2:
        raise AssertionError("2sqrt2 Leech first reflective shell drifted")
    # The pi^3-class candidate at h=epsilon is killed by y in pi*Lambda.
    # Dividing y=epsilon*b from the q=2 root gives normalized control sqrt2
    # and slabs |z.b|<=4.  For minimal u, |u.b|<=3, hence z=6u/5 is strict.
    coefficient = Silver(Fraction(6, 5))
    square = 4 * coefficient * coefficient
    return NonunitLeechFirstPacketAudit(
        seed_id=TWO_SQRT2_LEECH_SEED.seed_id,
        earlier_control_count=earlier,
        first_height=level.height,
        first_norm=level.norm,
        first_control=level.control,
        first_target=level.spatial_square_target,
        scalar_shell_norm=shell_norm,
        wall_count=LEECH_SHELL_COUNTS[shell_norm],
        slab_bound=Silver(4),
        witness_description="z=(6/5)u for any minimal u^2=4",
        witness_spatial_square=square,
        timelike_threshold=2 * SQRT2,
        strict_hyperideal_witness=square > 2 * SQRT2,
    )


def nonunit_leech_second_packet_audit(
    seed: SilverLeechSeed,
) -> NonunitLeechSecondPacketAudit:
    """Certify a minimal-vector hyperideal ray through the second packet."""

    if seed not in (SQRT2_LEECH_SEED, TWO_SQRT2_LEECH_SEED):
        raise ValueError("expected a nonunit silver Leech seed")
    height = Silver(24, 17)
    levels = tuple(
        level
        for level in enumerate_control_levels(seed, 49)
        if level.height == height
    )
    if len(levels) != 1 or levels[0].norm != 2:
        raise AssertionError("the nonunit second control level drifted")
    level = levels[0]
    triples = tuple(
        triple
        for triple in reflective_candidate_leech_gram_triples(
            seed, level.norm, level.spatial_square_target
        )
        if not has_forbidden_norm_two_combination(triple)
    )
    expected = LeechGramTriple(26, 12, 17)
    if triples != (expected,):
        raise AssertionError("the nonunit second Leech Gram triple drifted")
    # w=u-v has w^2=4, v^2=12, w.v=5.  For fixed minimal u0,
    # |u0.w|<=4 and |u0.v|<=6, so |u0.y|<=4+6epsilon=10+6sqrt2.
    maximum = Silver(10, 6)
    if seed == SQRT2_LEECH_SEED:
        witness = Silver(57, 19) / 140
        first_bound = 2 * SQRT2
        coefficient = SQRT2
    else:
        witness = Silver(38, 57) / 140
        first_bound = Silver(4)
        coefficient = 2 * SQRT2
    square = 4 * witness * witness
    return NonunitLeechSecondPacketAudit(
        seed_id=seed.seed_id,
        complete_height=height,
        norm=level.norm,
        control=level.control,
        target=level.spatial_square_target,
        surviving_triple=expected,
        decomposition="y=w+epsilon*v, w^2=4, v^2=12, w.v=5",
        witness_coordinate=witness,
        first_packet_margin=first_bound - 3 * witness,
        second_packet_margin=coefficient * level.control - witness * maximum,
        witness_spatial_square=square,
        timelike_threshold=coefficient,
        hyperideal_prefix_witness=square > coefficient,
    )
