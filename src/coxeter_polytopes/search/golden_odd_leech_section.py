"""Exact golden first packets for the odd-Leech section in ``H23``.

Let ``S23=(1^24)^perp cap O24`` be the determinant-three, minimum-three
section certified in :mod:`cubic7_odd_leech_section`.  This module studies

``<-a> + S23`` over ``Z[phi]``

for the two elementary admissible timelike coefficients ``a=phi`` and
``a=2*phi``.  The embeddings of either coefficient have signs ``(+,-)``,
so the conjugate quadratic form is positive definite.

The calculation is deliberately prefix-local.  For each coefficient it
proves the complete local mirror-norm alphabet, eliminates height zero,
enumerates every control through the first nonempty positive packet, and
exhibits a strict spacelike point in that finite prefix cone.  Later mirrors
can still cut the displayed point.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import comb

from coxeter_polytopes.exact.golden import PHI, Golden
from coxeter_polytopes.search.cubic7_odd_leech_section import (
    RADICAL_NORM_EXCLUSION,
    OddLeechSectionModelAudit,
    odd_leech_section_model_audit,
)
from coxeter_polytopes.search.golden_leech_manifest import _bounded_golden_for_caps
from coxeter_polytopes.search.silver_odd_leech import (
    OddLeechGramTriple,
    _binary_minimum,
    catalog_golay_codewords,
)


CANONICAL_MIRROR_NORMS: tuple[Golden, ...] = tuple(
    Golden(value) for value in (1, 2, 3, 6)
)
THREE_GATED_NORMS = frozenset((Golden(3), Golden(6)))

UNIT_COEFFICIENT = PHI
TWICE_UNIT_COEFFICIENT = 2 * PHI

COMPLETE_MIRROR_NORM_PROOF = (
    "The primes 2 and 3 are inert in Z[phi], the class number is one, and "
    "the totally-positive units are precisely the even powers of phi, hence "
    "unit squares.  For a primitive root r=(x,y), reflection integrality is "
    "(q(r))|2I(r), I(r)=(a*x,(y,S23)).  Since S23 has Smith invariants "
    "1^22,3, every prime away from 2, 3, and the coefficient would divide "
    "both x and y, contradicting primitivity.  For a=phi the coefficient is "
    "a unit, so v_2(q)<=v_2(2)=1.  For a=2*phi, either y has a unit pairing "
    "and the same bound holds, or y lies in 2*S23; then primitivity makes x "
    "a 2-adic unit and q=-2*phi*x^2 mod 4, so v_2(q)=1.  At 3, positive "
    "q-valuation forces x into (3) and both coefficients of y into the "
    "one-dimensional pairing radical modulo 3.  Its primitive pairing ideal "
    "has valuation exactly one, so v_3(q)<=1.  Thus, modulo totally-positive "
    "unit-square scaling, the complete norm alphabet for both coefficients "
    "is 1,2,3,6.  In the last two classes both spatial coefficients lie in "
    "the mod-three pairing radical."
)

HEIGHT_ZERO_PROOF = (
    "Write y=u+phi*v.  If y^2 is the rational canonical norm q, then "
    "u^2+v^2=q and 2(u.v)+v^2=0.  Minimum three excludes q=1,2.  Before "
    "the 3-adic gate the only surviving binary Gram types for q=3,6 are "
    "(3,0,0) and (6,0,0).  The gate puts u and v in the mod-three pairing "
    "radical, whose certified minimum is at least seven, so both types are "
    "excluded.  Hence height zero is empty for both timelike coefficients."
)


def necessary_section_gram_triples(
    target: Golden, *, pairing_modulus: int = 1, coefficient_minimum: int = 3
) -> tuple[OddLeechGramTriple, ...]:
    """All binary Gram types allowed by the exact target and local gates."""

    if not target.is_integral or target.a < 0:
        return ()
    rational_target = int(target.a)
    phi_target = int(target.b)
    output = []
    for v_square in range(rational_target + 1):
        u_square = rational_target - v_square
        dot_numerator = phi_target - v_square
        if dot_numerator % 2:
            continue
        u_dot_v = dot_numerator // 2
        if u_square not in (0,) and u_square < 3:
            continue
        if v_square not in (0,) and v_square < 3:
            continue
        if u_dot_v * u_dot_v > u_square * v_square:
            continue
        if any(
            value % pairing_modulus
            for value in (u_square, v_square, u_dot_v)
        ):
            continue
        triple = OddLeechGramTriple(u_square, v_square, u_dot_v)
        if _binary_minimum(triple) < coefficient_minimum:
            continue
        output.append(triple)
    return tuple(output)


@dataclass(frozen=True, slots=True)
class GoldenSectionShellAudit:
    norm: int
    even_coset_count: int
    odd_coset_count: int
    component_counts: tuple[tuple[str, int], ...]
    total_count: int


NORM_FOUR_SHELL_PROOF = (
    "On the raw sqrt(8) scale, norm four means square 32 and coordinate sum "
    "zero.  In the even coset write the vector as 2t.  The equation t^2=8 "
    "has only two Golay-compatible parity supports: weight eight, giving "
    "759*C(8,4)=53130 sign-balanced octad vectors, and weight zero, giving "
    "24*23=552 ordered choices of one +4 and one -4 raw coordinate.  In the "
    "odd coset square 32 would mean one coordinate of absolute value three "
    "and 23 of absolute value one.  Flipping that coordinate changes a "
    "Golay sign-word sum by plus or minus four, which would require codeword "
    "weight 10 or 14; neither exists.  Thus the complete section shell has "
    "53682 vectors."
)

NORM_FIVE_SHELL_PROOF = (
    "On the raw sqrt(8) scale, norm five means square 40 and coordinate sum "
    "zero.  In the even coset, after division by two, every partition of ten "
    "into squares has odd-coordinate support of weight 2, 6, or 10, none a "
    "Golay weight.  In the odd coset exactly two coordinates have absolute "
    "value three and the other 22 have absolute value one.  Relative to a "
    "Golay sign word, the two large coordinates reverse sign.  Zero sum "
    "allows: two zero positions for weight 8, one zero and one one position "
    "for weight 12, or two one positions for weight 16.  Hence the complete "
    "count is 759*C(16,2)+2576*12*12+759*C(16,2)=553104."
)


@lru_cache(maxsize=1)
def section_norm_four_shell_audit() -> GoldenSectionShellAudit:
    distribution = Counter(word.bit_count() for word in catalog_golay_codewords())
    weight_eight = distribution[8] * comb(8, 4)
    weight_zero = 24 * 23
    components = (
        ("even Golay weight 0", weight_zero),
        ("even Golay weight 8", weight_eight),
        ("odd coset", 0),
    )
    total = sum(count for _label, count in components)
    if weight_zero != 552 or weight_eight != 53_130 or total != 53_682:
        raise AssertionError("the exact norm-four section census changed")
    return GoldenSectionShellAudit(4, total, 0, components, total)


@lru_cache(maxsize=1)
def section_norm_five_shell_audit() -> GoldenSectionShellAudit:
    distribution = Counter(word.bit_count() for word in catalog_golay_codewords())
    weight_eight = distribution[8] * comb(16, 2)
    weight_twelve = distribution[12] * 12 * 12
    weight_sixteen = distribution[16] * comb(16, 2)
    components = (
        ("even coset", 0),
        ("odd Golay weight 8", weight_eight),
        ("odd Golay weight 12", weight_twelve),
        ("odd Golay weight 16", weight_sixteen),
    )
    total = sum(count for _label, count in components)
    if (
        weight_eight != 91_080
        or weight_twelve != 370_944
        or weight_sixteen != 91_080
        or total != 553_104
    ):
        raise AssertionError("the exact norm-five section census changed")
    return GoldenSectionShellAudit(5, 0, total, components, total)


@dataclass(frozen=True, slots=True, order=True)
class GoldenSectionControlLevel:
    height: Golden
    norm: Golden
    control: Golden
    target: Golden
    pairing_modulus: int
    coefficient_minimum: int
    necessary_gram_triples: tuple[OddLeechGramTriple, ...]


UNIT_FIRST_HEIGHT = Golden(Fraction(3, 2), Fraction(5, 2))
UNIT_FIRST_NORM = Golden(2)
UNIT_FIRST_CONTROL = 1 + PHI
UNIT_FIRST_TARGET = Golden(5, 5)
UNIT_FIRST_SCALE = PHI

TWICE_UNIT_FIRST_HEIGHT = Golden(1, 2)
TWICE_UNIT_FIRST_NORM = Golden(2)
TWICE_UNIT_FIRST_CONTROL = PHI
TWICE_UNIT_FIRST_TARGET = Golden(4, 4)
TWICE_UNIT_FIRST_SCALE = PHI


def _controls_through(
    coefficient: Golden, endpoint: Golden, integer_height_cap: int
) -> tuple[GoldenSectionControlLevel, ...]:
    levels = set()
    for norm in CANONICAL_MIRROR_NORMS:
        controls = _bounded_golden_for_caps(
            (norm * integer_height_cap / coefficient, norm / (-coefficient))
        )
        for control in controls:
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
            gated = norm in THREE_GATED_NORMS
            modulus = 3 if gated else 1
            minimum = 7 if gated else 3
            levels.add(
                GoldenSectionControlLevel(
                    height,
                    norm,
                    control,
                    target,
                    modulus,
                    minimum,
                    necessary_section_gram_triples(
                        target,
                        pairing_modulus=modulus,
                        coefficient_minimum=minimum,
                    ),
                )
            )
    return tuple(sorted(levels))


@lru_cache(maxsize=1)
def unit_controls_through_first_packet() -> tuple[GoldenSectionControlLevel, ...]:
    """Complete controls through the first ``a=phi`` packet (endpoint < 6)."""

    return _controls_through(UNIT_COEFFICIENT, UNIT_FIRST_HEIGHT, 6)


@lru_cache(maxsize=1)
def twice_unit_controls_through_first_packet() -> tuple[GoldenSectionControlLevel, ...]:
    """Complete controls through the first ``a=2phi`` packet (endpoint < 5)."""

    return _controls_through(
        TWICE_UNIT_COEFFICIENT, TWICE_UNIT_FIRST_HEIGHT, 5
    )


@dataclass(frozen=True, slots=True)
class GoldenSectionPrefixAudit:
    coefficient: Golden
    model: OddLeechSectionModelAudit
    canonical_mirror_norms: tuple[Golden, ...]
    height_zero_ambient_types: tuple[
        tuple[Golden, tuple[OddLeechGramTriple, ...]], ...
    ]
    height_zero_gated_types: tuple[
        tuple[Golden, tuple[OddLeechGramTriple, ...]], ...
    ]
    controls: tuple[GoldenSectionControlLevel, ...]
    first_height: Golden
    first_norm: Golden
    first_control: Golden
    first_target: Golden
    first_scale: Golden
    first_gram_type: OddLeechGramTriple
    first_spatial_norm: int
    first_packet_root_count: int
    packet_pairing_upper_bound: Golden
    witness_label: str
    witness_spatial_scale: Fraction
    witness_square: Golden
    closest_wall_pairing: Golden
    strict_hyperideal: bool
    finite_volume: bool
    compact: bool


def _height_zero_censuses() -> tuple[
    tuple[tuple[Golden, tuple[OddLeechGramTriple, ...]], ...],
    tuple[tuple[Golden, tuple[OddLeechGramTriple, ...]], ...],
]:
    ambient = tuple(
        (norm, necessary_section_gram_triples(norm))
        for norm in CANONICAL_MIRROR_NORMS
    )
    gated = tuple(
        (
            norm,
            necessary_section_gram_triples(
                norm,
                pairing_modulus=3 if norm in THREE_GATED_NORMS else 1,
                coefficient_minimum=7 if norm in THREE_GATED_NORMS else 3,
            ),
        )
        for norm in CANONICAL_MIRROR_NORMS
    )
    if tuple(len(types) for _norm, types in ambient) != (0, 0, 1, 1):
        raise AssertionError("the golden height-zero ambient census changed")
    if any(types for _norm, types in gated):
        raise AssertionError("a golden height-zero wall survived its local gate")
    return ambient, gated


@lru_cache(maxsize=1)
def unit_section_first_packet_audit() -> GoldenSectionPrefixAudit:
    model = odd_leech_section_model_audit()
    shell = section_norm_five_shell_audit()
    ambient, gated = _height_zero_censuses()
    controls = unit_controls_through_first_packet()
    nonempty = tuple(level for level in controls if level.necessary_gram_triples)
    expected = OddLeechGramTriple(0, 5, 0)
    if (
        len(controls) != 7
        or len(nonempty) != 1
        or nonempty[0].height != UNIT_FIRST_HEIGHT
        or nonempty[0].norm != UNIT_FIRST_NORM
        or nonempty[0].control != UNIT_FIRST_CONTROL
        or nonempty[0].target != UNIT_FIRST_TARGET
        or nonempty[0].necessary_gram_triples != (expected,)
    ):
        raise AssertionError("the complete golden-unit first prefix changed")

    root_square = (
        -UNIT_COEFFICIENT * UNIT_FIRST_CONTROL**2
        + Golden(5) * UNIT_FIRST_SCALE**2
    )
    integrality = (
        2 * UNIT_COEFFICIENT * UNIT_FIRST_CONTROL / UNIT_FIRST_NORM,
        2 * UNIT_FIRST_SCALE / UNIT_FIRST_NORM,
    )
    if root_square != UNIT_FIRST_NORM or not all(
        value.is_integral for value in integrality
    ):
        raise AssertionError("the displayed norm-five packet roots lost validity")

    # Distinct norm-five vectors have inner product at most three because
    # the section minimum is three.
    packet_bound = (
        -UNIT_COEFFICIENT * UNIT_FIRST_CONTROL**2
        + Golden(3) * UNIT_FIRST_SCALE**2
    )
    witness_scale = Fraction(2, 3)
    witness_square = Golden(4 * witness_scale * witness_scale) - UNIT_COEFFICIENT
    # The norm-four witness pairs at most three with a norm-five vector.
    closest = (
        -UNIT_COEFFICIENT * UNIT_FIRST_CONTROL
        + Golden(3 * witness_scale) * UNIT_FIRST_SCALE
    )
    if (
        packet_bound != -2 * PHI
        or closest != -1
        or packet_bound.sign(0) >= 0
        or closest.sign(0) >= 0
        or witness_square.sign(0) <= 0
    ):
        raise AssertionError("the golden-unit hyperideal certificate failed")
    return GoldenSectionPrefixAudit(
        UNIT_COEFFICIENT,
        model,
        CANONICAL_MIRROR_NORMS,
        ambient,
        gated,
        controls,
        UNIT_FIRST_HEIGHT,
        UNIT_FIRST_NORM,
        UNIT_FIRST_CONTROL,
        UNIT_FIRST_TARGET,
        UNIT_FIRST_SCALE,
        expected,
        5,
        shell.total_count,
        packet_bound,
        "norm-four section vector (4,-4,0^22)/sqrt(8)",
        witness_scale,
        witness_square,
        closest,
        True,
        False,
        False,
    )


@lru_cache(maxsize=1)
def twice_unit_section_first_packet_audit() -> GoldenSectionPrefixAudit:
    model = odd_leech_section_model_audit()
    shell = section_norm_four_shell_audit()
    ambient, gated = _height_zero_censuses()
    controls = twice_unit_controls_through_first_packet()
    nonempty = tuple(level for level in controls if level.necessary_gram_triples)
    expected = OddLeechGramTriple(0, 4, 0)
    if (
        len(controls) != 2
        or len(nonempty) != 1
        or nonempty[0].height != TWICE_UNIT_FIRST_HEIGHT
        or nonempty[0].norm != TWICE_UNIT_FIRST_NORM
        or nonempty[0].control != TWICE_UNIT_FIRST_CONTROL
        or nonempty[0].target != TWICE_UNIT_FIRST_TARGET
        or nonempty[0].necessary_gram_triples != (expected,)
    ):
        raise AssertionError("the complete twice-golden-unit first prefix changed")

    root_square = (
        -TWICE_UNIT_COEFFICIENT * TWICE_UNIT_FIRST_CONTROL**2
        + Golden(4) * TWICE_UNIT_FIRST_SCALE**2
    )
    integrality = (
        2
        * TWICE_UNIT_COEFFICIENT
        * TWICE_UNIT_FIRST_CONTROL
        / TWICE_UNIT_FIRST_NORM,
        2 * TWICE_UNIT_FIRST_SCALE / TWICE_UNIT_FIRST_NORM,
    )
    if root_square != TWICE_UNIT_FIRST_NORM or not all(
        value.is_integral for value in integrality
    ):
        raise AssertionError("the displayed norm-four packet roots lost validity")

    # Distinct norm-four vectors have inner product at most two.
    packet_bound = (
        -TWICE_UNIT_COEFFICIENT * TWICE_UNIT_FIRST_CONTROL**2
        + Golden(2) * TWICE_UNIT_FIRST_SCALE**2
    )
    witness_scale = Fraction(4, 3)
    witness_square = (
        Golden(3 * witness_scale * witness_scale) - TWICE_UNIT_COEFFICIENT
    )
    # A norm-three vector pairs at most two with a norm-four vector.
    closest = (
        -TWICE_UNIT_COEFFICIENT * TWICE_UNIT_FIRST_CONTROL
        + Golden(2 * witness_scale) * TWICE_UNIT_FIRST_SCALE
    )
    if (
        packet_bound != -2 * PHI
        or packet_bound.sign(0) >= 0
        or closest.sign(0) >= 0
        or witness_square.sign(0) <= 0
    ):
        raise AssertionError("the twice-golden-unit hyperideal certificate failed")
    return GoldenSectionPrefixAudit(
        TWICE_UNIT_COEFFICIENT,
        model,
        CANONICAL_MIRROR_NORMS,
        ambient,
        gated,
        controls,
        TWICE_UNIT_FIRST_HEIGHT,
        TWICE_UNIT_FIRST_NORM,
        TWICE_UNIT_FIRST_CONTROL,
        TWICE_UNIT_FIRST_TARGET,
        TWICE_UNIT_FIRST_SCALE,
        expected,
        4,
        shell.total_count,
        packet_bound,
        "any norm-three vector in the certified 2576-vector section shell",
        witness_scale,
        witness_square,
        closest,
        True,
        False,
        False,
    )


__all__ = [name for name in globals() if not name.startswith("_")]
