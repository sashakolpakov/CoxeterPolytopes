"""Exact silver screens for the odd-Leech section ``S23`` in ``H23``.

Let ``O24`` be the catalogue odd Leech lattice and put

``S23 = O24 cap (1^24)^perp``.

The exact model, determinant-three Smith form, minimum-three shell, and the
norm-four test vector are certified in :mod:`cubic7_odd_leech_section`.
Here the Lorentzian forms are ``<-a> + S23`` over ``Z[sqrt(2)]`` for the two
admissible coefficients ``a=1+sqrt(2)`` and ``a=sqrt(2)``.

Both complete prefixes through their first nonempty packet have a strict
hyperideal point.  The result is packet-local: later mirrors can still cut
that point.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import comb

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.cubic7_odd_leech_section import (
    RADICAL_NORM_EXCLUSION,
    SECTION_WITNESS_SCALED,
    OddLeechSectionModelAudit,
    odd_leech_section_model_audit,
)
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps
from coxeter_polytopes.search.silver_odd_leech import (
    OddLeechGramTriple,
    _binary_minimum,
    catalog_golay_codewords,
    necessary_odd_leech_gram_triples,
)


EPSILON = SILVER_ONE + SQRT2

# S23 is locally unimodular away from 3 and has one elementary divisor 3.
# At the dyadic prime pi=sqrt(2), the possible exponents are 0,1,2.  The
# totally-positive representatives and their products with 3 are below.
CANONICAL_MIRROR_NORMS: tuple[Silver, ...] = (
    Silver(1),
    Silver(2, 1),
    Silver(2),
    Silver(3),
    Silver(6, 3),
    Silver(6),
)
THREE_GATED_NORMS = frozenset(CANONICAL_MIRROR_NORMS[3:])

COMPLETE_MIRROR_NORM_PROOF = (
    "Put pi=sqrt(2).  The section S23 has Smith invariants 1^22,3.  For a "
    "primitive r=(x,y), reflection integrality gives q(r)|2I(r), where "
    "I(r)=(a*x,(y,S23)).  Local unimodularity and primitivity exclude every "
    "prime away from (pi) and (3).  At (3), positive q-valuation puts x in "
    "(3) and every coefficient of y in the one-dimensional pairing radical "
    "modulo 3.  That radical is primitive and the discriminant has exponent "
    "one, so I(r) has valuation exactly one and v_3(q)<=1.  At (pi), if a "
    "is a unit then I(r) is a unit and v_pi(q)<=v_pi(2)=2.  If a=pi and y "
    "has a unit pairing, the same bound holds.  Otherwise nondegeneracy of "
    "S23 modulo pi gives y in pi*S23; primitivity makes x a unit and "
    "q=-pi*x^2 modulo pi^2, so v_pi(q)=1.  Thus v_pi(q) is 0,1,or 2 and "
    "v_3(q) is 0 or 1.  Z[sqrt(2)] has class number one and its totally "
    "positive units are squares.  The complete canonical list is therefore "
    "1,2+pi,2,3,6+3pi,6.  In the final three classes both silver "
    "coefficients of y lie in the mod-three pairing radical."
)

HEIGHT_ZERO_PROOF = (
    "Write y=u+sqrt(2)*v.  Then Tr(y^2)=2(u^2+2v^2).  Since S23 has "
    "minimum three, norms 1,2+sqrt(2),2 cannot occur at height zero.  For "
    "norms divisible by 3, u and v lie in the mod-three pairing radical.  "
    "The certified radical exclusion gives minimum at least seven, whereas "
    "the three targets have trace at most twelve.  Hence no height-zero "
    "wall survives for either coefficient."
)


def _filtered_gram_triples(
    target: Silver, coefficient_minimum: int
) -> tuple[OddLeechGramTriple, ...]:
    """Apply the exact minimum bound to all necessary binary Gram triples."""

    return tuple(
        triple
        for triple in necessary_odd_leech_gram_triples(target)
        if _binary_minimum(triple) >= coefficient_minimum
    )


@dataclass(frozen=True, slots=True)
class NormSixSectionShellAudit:
    """Complete norm-six census without storing 3.7 million vectors."""

    odd_coset_count: int
    even_weight_zero_count: int
    even_weight_eight_count: int
    even_weight_twelve_count: int
    total_count: int
    witness_pairing_histogram: tuple[tuple[int, int], ...]
    witness_maximum_pairing: int


NORM_SIX_SHELL_PROOF = (
    "Use raw coordinates scaled by sqrt(8).  A norm-six vector has raw "
    "square 48 and coordinate sum zero.  In the odd coset all coordinates "
    "are odd, so exactly three have absolute value 3.  Relative to its "
    "Golay sign word those three signs reverse, and the zero-sum equation "
    "would equate one of 24,8,0,-8,-24 with four times a sum of three signs, "
    "which is impossible.  In the even coset write the raw vector as 2z.  "
    "Then z^2=12 and z mod 2 is a Golay word, hence has weight 0,8,or 12.  "
    "Weight zero would require three signed unit entries after dividing by "
    "two and cannot have sum zero.  At weight twelve, z has twelve entries "
    "+/-1 with six of each: 2576*C(12,6)=2380224 vectors.  At weight eight, "
    "there are eight entries +/-1 and one outside entry +/-2; the zero-sum "
    "condition gives 759*16*2*C(8,3)=1360128 vectors.  The Construction-B "
    "residual sum is even in both cases.  This is the complete norm-six "
    "shell of 3740352 vectors."
)


@lru_cache(maxsize=1)
def norm_six_section_shell_audit() -> NormSixSectionShellAudit:
    """Count the shell and its pairings with ``(4,-4,0^22)/sqrt(8)``."""

    histogram: Counter[int] = Counter()
    weight_eight_count = 0
    weight_twelve_count = 0

    for word in catalog_golay_codewords():
        weight = word.bit_count()
        in_zero = bool(word & 1)
        in_one = bool(word & 2)

        if weight == 12:
            weight_twelve_count += comb(12, 6)
            if not in_zero and not in_one:
                histogram[0] += comb(12, 6)
            elif in_zero != in_one:
                # The supported distinguished coordinate is +1 or -1.
                histogram[-1] += comb(11, 5)
                histogram[1] += comb(11, 6)
            else:
                # Equal signs give pairing zero; opposite signs give +/-2.
                histogram[0] += 2 * comb(10, 4)
                histogram[-2] += comb(10, 5)
                histogram[2] += comb(10, 5)
            continue

        if weight != 8:
            continue

        support = tuple(index for index in range(24) if word & (1 << index))
        outside = tuple(index for index in range(24) if not word & (1 << index))
        distinguished_support = tuple(index for index in (0, 1) if index in support)
        remaining = 8 - len(distinguished_support)
        for outside_index in outside:
            for outside_sign, plus_count in ((1, 3), (-1, 5)):
                for mask in range(1 << len(distinguished_support)):
                    plus_distinguished = {
                        index
                        for bit, index in enumerate(distinguished_support)
                        if mask & (1 << bit)
                    }
                    needed = plus_count - len(plus_distinguished)
                    if not 0 <= needed <= remaining:
                        continue
                    multiplicity = comb(remaining, needed)

                    def coordinate(index: int) -> int:
                        if index == outside_index:
                            return 2 * outside_sign
                        if index in support:
                            return 1 if index in plus_distinguished else -1
                        return 0

                    histogram[coordinate(0) - coordinate(1)] += multiplicity
                    weight_eight_count += multiplicity

    expected_histogram = Counter(
        {-3: 12320, -2: 229152, -1: 924000, 0: 1409408,
         1: 924000, 2: 229152, 3: 12320}
    )
    if (
        weight_twelve_count != 2_380_224
        or weight_eight_count != 1_360_128
        or histogram != expected_histogram
    ):
        raise AssertionError("the exact norm-six section census changed")
    return NormSixSectionShellAudit(
        odd_coset_count=0,
        even_weight_zero_count=0,
        even_weight_eight_count=weight_eight_count,
        even_weight_twelve_count=weight_twelve_count,
        total_count=weight_eight_count + weight_twelve_count,
        witness_pairing_histogram=tuple(sorted(histogram.items())),
        witness_maximum_pairing=max(histogram),
    )


@dataclass(frozen=True, slots=True, order=True)
class SilverSectionControlLevel:
    height: Silver
    norm: Silver
    control: Silver
    target: Silver
    coefficient_minimum: int
    necessary_gram_triples: tuple[OddLeechGramTriple, ...]


UNIT_COEFFICIENT = EPSILON
UNIT_FIRST_HEIGHT = Silver(4, 3) / 2
UNIT_FIRST_NORM = Silver(2, 1)
UNIT_FIRST_CONTROL = EPSILON
UNIT_FIRST_TARGET = Silver(9, 6)

SQRT2_COEFFICIENT = SQRT2
SQRT2_FIRST_HEIGHT = Silver(8, 6)
SQRT2_FIRST_Q1_CONTROL = Silver(2, 1)
SQRT2_FIRST_Q2_CONTROL = Silver(2, 2)


def _controls_through(
    coefficient: Silver, endpoint: Silver, integer_cap: int
) -> tuple[SilverSectionControlLevel, ...]:
    levels = set()
    for norm in CANONICAL_MIRROR_NORMS:
        caps = (norm * integer_cap / coefficient, norm / (-coefficient))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if target.sign(1) < 0:
                continue
            height = coefficient * control * control / norm
            if (endpoint - height).sign(0) < 0:
                continue
            minimum = 7 if norm in THREE_GATED_NORMS else 3
            levels.add(
                SilverSectionControlLevel(
                    height=height,
                    norm=norm,
                    control=control,
                    target=target,
                    coefficient_minimum=minimum,
                    necessary_gram_triples=_filtered_gram_triples(target, minimum),
                )
            )
    return tuple(sorted(levels))


@lru_cache(maxsize=1)
def unit_controls_through_first_packet() -> tuple[SilverSectionControlLevel, ...]:
    """Complete control census below the exact unit endpoint (<5)."""

    return _controls_through(UNIT_COEFFICIENT, UNIT_FIRST_HEIGHT, 5)


@lru_cache(maxsize=1)
def sqrt2_controls_through_first_packet() -> tuple[SilverSectionControlLevel, ...]:
    """Complete control census below the exact sqrt(2) endpoint (<17)."""

    return _controls_through(SQRT2_COEFFICIENT, SQRT2_FIRST_HEIGHT, 17)


@dataclass(frozen=True, slots=True)
class SilverSectionPrefixAudit:
    coefficient: Silver
    model: OddLeechSectionModelAudit
    canonical_mirror_norms: tuple[Silver, ...]
    height_zero_ambient_types: tuple[tuple[Silver, tuple[OddLeechGramTriple, ...]], ...]
    height_zero_gated_types: tuple[tuple[Silver, tuple[OddLeechGramTriple, ...]], ...]
    controls: tuple[SilverSectionControlLevel, ...]
    first_height: Silver
    first_packet_norm_counts: tuple[tuple[Silver, int], ...]
    imprimitive_candidate_count: int
    accepted_wall_count: int
    packet_pairing_upper_bounds: tuple[tuple[str, Silver], ...]
    norm_six_shell: NormSixSectionShellAudit | None
    witness_spatial_scale: Fraction
    witness_square: Silver
    closest_wall_pairing: Silver
    strict_hyperideal: bool
    finite_volume: bool
    compact: bool


def _height_zero_censuses() -> tuple[
    tuple[tuple[Silver, tuple[OddLeechGramTriple, ...]], ...],
    tuple[tuple[Silver, tuple[OddLeechGramTriple, ...]], ...],
]:
    ambient = tuple(
        (norm, necessary_odd_leech_gram_triples(norm))
        for norm in CANONICAL_MIRROR_NORMS
    )
    gated = tuple(
        (
            norm,
            _filtered_gram_triples(norm, 7 if norm in THREE_GATED_NORMS else 3),
        )
        for norm in CANONICAL_MIRROR_NORMS
    )
    if tuple(len(types) for _norm, types in ambient) != (0, 0, 0, 1, 0, 2):
        raise AssertionError("the height-zero ambient census changed")
    if any(types for _norm, types in gated):
        raise AssertionError("a height-zero section wall survived its local gate")
    return ambient, gated


@lru_cache(maxsize=1)
def unit_section_first_packet_audit() -> SilverSectionPrefixAudit:
    model = odd_leech_section_model_audit()
    ambient, gated = _height_zero_censuses()
    controls = unit_controls_through_first_packet()
    nonempty = tuple(level for level in controls if level.necessary_gram_triples)
    expected = OddLeechGramTriple(3, 3, 3)
    if (
        len(controls) != 6
        or len(nonempty) != 1
        or nonempty[0].height != UNIT_FIRST_HEIGHT
        or nonempty[0].norm != UNIT_FIRST_NORM
        or nonempty[0].control != UNIT_FIRST_CONTROL
        or nonempty[0].target != UNIT_FIRST_TARGET
        or nonempty[0].necessary_gram_triples != (expected,)
    ):
        raise AssertionError("the complete silver-unit first prefix changed")

    # The Gram triple forces u=v=b, b^2=3.  These are all 2576 shell vectors.
    if (
        -EPSILON * EPSILON * EPSILON
        + 3 * EPSILON * EPSILON
        != UNIT_FIRST_NORM
        or not (2 * UNIT_COEFFICIENT * UNIT_FIRST_CONTROL / UNIT_FIRST_NORM).is_integral
        or not (2 * EPSILON / UNIT_FIRST_NORM).is_integral
    ):
        raise AssertionError("the displayed unit packet root lost norm or integrality")
    # For distinct minimum vectors b,c, b.c is at most one.  Thus every
    # off-diagonal root pairing is <= eps^2*(1-eps)<0.
    closest = EPSILON * EPSILON * (Silver(1) - EPSILON)
    witness_square = Silver(model.witness_square) - EPSILON
    if closest.sign(0) >= 0 or witness_square.sign(0) <= 0:
        raise AssertionError("the unit strict hyperideal witness failed")
    return SilverSectionPrefixAudit(
        coefficient=UNIT_COEFFICIENT,
        model=model,
        canonical_mirror_norms=CANONICAL_MIRROR_NORMS,
        height_zero_ambient_types=ambient,
        height_zero_gated_types=gated,
        controls=controls,
        first_height=UNIT_FIRST_HEIGHT,
        first_packet_norm_counts=((UNIT_FIRST_NORM, model.minimal_shell_count),),
        imprimitive_candidate_count=0,
        accepted_wall_count=model.minimal_shell_count,
        packet_pairing_upper_bounds=(("q=(2+sqrt2) versus itself", closest),),
        norm_six_shell=None,
        witness_spatial_scale=Fraction(1),
        witness_square=witness_square,
        closest_wall_pairing=EPSILON * (Silver(1) - EPSILON),
        strict_hyperideal=True,
        finite_volume=False,
        compact=False,
    )


@lru_cache(maxsize=1)
def sqrt2_section_first_packet_audit() -> SilverSectionPrefixAudit:
    model = odd_leech_section_model_audit()
    shell6 = norm_six_section_shell_audit()
    ambient, gated = _height_zero_censuses()
    controls = sqrt2_controls_through_first_packet()
    nonempty = tuple(level for level in controls if level.necessary_gram_triples)
    q1_triple = OddLeechGramTriple(3, 3, 3)
    q2_imprimitive = OddLeechGramTriple(12, 3, 6)
    q2_primitive = OddLeechGramTriple(6, 6, 6)
    if (
        len(controls) != 10
        or len(nonempty) != 2
        or any(level.height != SQRT2_FIRST_HEIGHT for level in nonempty)
        or nonempty[0].norm != Silver(1)
        or nonempty[0].control != SQRT2_FIRST_Q1_CONTROL
        or nonempty[0].necessary_gram_triples != (q1_triple,)
        or nonempty[1].norm != Silver(2)
        or nonempty[1].control != SQRT2_FIRST_Q2_CONTROL
        or nonempty[1].necessary_gram_triples != (q2_imprimitive, q2_primitive)
    ):
        raise AssertionError("the complete sqrt(2) first prefix changed")

    # q=1: (x,y)=(pi*eps,eps*b), b^2=3.
    q1_square = -SQRT2 * (SQRT2 * EPSILON) ** 2 + 3 * EPSILON ** 2
    # q=2 primitive: (x,y)=(2*eps,eps*c), c^2=6.
    q2_square = -SQRT2 * (2 * EPSILON) ** 2 + 6 * EPSILON ** 2
    if q1_square != 1 or q2_square != 2:
        raise AssertionError("a displayed sqrt(2) packet root lost its norm")
    if (
        not (2 * SQRT2 * SQRT2_FIRST_Q1_CONTROL).is_integral
        or not (2 * EPSILON).is_integral
        or not (2 * SQRT2 * SQRT2_FIRST_Q2_CONTROL / 2).is_integral
        or not EPSILON.is_integral
    ):
        raise AssertionError("a displayed sqrt(2) packet root lost reflection integrality")

    # The other q=2 Gram type is exactly pi times a q=1 shell root.
    # Acuteness follows from minimum three: b.c<=3 for norms 3 and 6,
    # while distinct norm-six vectors have c.d<=4.
    q1_q1_closest = EPSILON ** 2 * (Silver(1) - 2 * SQRT2)
    q1_q2_closest = EPSILON ** 2 * (Silver(3) - Silver(4))
    q2_q2_closest = EPSILON ** 2 * (Silver(4) - 4 * SQRT2)
    if any(value.sign(0) >= 0 for value in (q1_q1_closest, q1_q2_closest, q2_q2_closest)):
        raise AssertionError("the simultaneous sqrt(2) packet ceased to be acute")

    # z=(1,(2/3)e), e^2=4.  The shell maxima are 1 and 3.
    scale = Fraction(2, 3)
    q1_closest = EPSILON * (Silver(scale) - Silver(2))
    q2_closest = EPSILON * (Silver(3 * scale) - 2 * SQRT2)
    closest = max(q1_closest, q2_closest)
    witness_square = Silver(Fraction(4) * scale * scale) - SQRT2
    if closest.sign(0) >= 0 or witness_square.sign(0) <= 0:
        raise AssertionError("the sqrt(2) strict hyperideal witness failed")
    accepted = model.minimal_shell_count + shell6.total_count
    return SilverSectionPrefixAudit(
        coefficient=SQRT2_COEFFICIENT,
        model=model,
        canonical_mirror_norms=CANONICAL_MIRROR_NORMS,
        height_zero_ambient_types=ambient,
        height_zero_gated_types=gated,
        controls=controls,
        first_height=SQRT2_FIRST_HEIGHT,
        first_packet_norm_counts=(
            (Silver(1), model.minimal_shell_count),
            (Silver(2), shell6.total_count),
        ),
        imprimitive_candidate_count=model.minimal_shell_count,
        accepted_wall_count=accepted,
        packet_pairing_upper_bounds=(
            ("q=1 versus q=1", q1_q1_closest),
            ("q=1 versus q=2", q1_q2_closest),
            ("q=2 versus q=2", q2_q2_closest),
        ),
        norm_six_shell=shell6,
        witness_spatial_scale=scale,
        witness_square=witness_square,
        closest_wall_pairing=closest,
        strict_hyperideal=True,
        finite_volume=False,
        compact=False,
    )


__all__ = [name for name in globals() if not name.startswith("_")]
