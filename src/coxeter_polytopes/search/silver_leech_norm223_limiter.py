"""Realization and Vinberg acceptance of the norm-223 H=114998 limiter.

The H=250000 norm-223 audit deliberately treated every binary Leech Gram
type as merely necessary.  This module resolves its unique limiting component.
Exactly two binary types attain the norm-four direction bound.  One would
create a norm-two Leech vector and is impossible; the other has an explicit
Construction-A/MOG realization and gives a primitive, reflective, accepted
norm-two root.

Acceptance is stronger than a sample test.  An exact Schur/Cauchy inequality
makes the displayed root strictly obtuse to every root in every earlier
necessary component, and a difference-shell argument handles its complete
same-height cohort.  The result is still bounded at H=250000 and does not
assert Vinberg termination or finite volume.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import isqrt

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.cubic7_leech_continuation import (
    NORM_SIX_DIRECTION_SCALED,
    _scaled_leech_dot,
    norm_six_direction_audit,
)
from coxeter_polytopes.search.silver_leech import LeechGramTriple
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)
from coxeter_polytopes.search.silver_leech_norm223 import (
    DIRECTION_SQUARE,
    NORM223_COEFFICIENT,
    _candidate_triples,
    d4_pairing_bound,
    norm223_control_skeleton,
    norm223_leech_prefix_audit,
)


LIMITING_HEIGHT = Silver(114_998, 81_316)
LIMITING_NORM = Silver(2)
LIMITING_CONTROL = Silver(130, 92)
LIMITING_TARGET = Silver(229_998, 162_632)
LIMITING_MAXIMUM = Silver(679, 479)
REALIZED_TRIPLE = LeechGramTriple(115_266, 57_366, 81_316)
IMPOSSIBLE_TRIPLE = LeechGramTriple(115_270, 57_364, 81_316)

# Standard MOG coordinates are divided by sqrt(8).  D is a type-(4^2)
# minimal vector.  E is the previously audited sum of two minimal vectors of
# pairing -1, hence belongs to the Leech lattice and has norm six.
LIMITING_D_SCALED: tuple[int, ...] = (4, -4) + (0,) * 22
LIMITING_E_SCALED: tuple[int, ...] = NORM_SIX_DIRECTION_SCALED
LIMITING_U_SCALED: tuple[int, ...] = tuple(
    170 * d + e for d, e in zip(LIMITING_D_SCALED, LIMITING_E_SCALED)
)
LIMITING_V_SCALED: tuple[int, ...] = tuple(
    120 * d + e for d, e in zip(LIMITING_D_SCALED, LIMITING_E_SCALED)
)


@dataclass(frozen=True, slots=True)
class LimitingMaximizer:
    triple: LeechGramTriple
    d_dot_u: int
    d_dot_v: int
    augmented_determinant: int
    shifted_e_square: int
    shifted_f_square: int
    shifted_e_dot_f: int
    shifted_difference_square: int
    realized: bool
    obstruction: str | None


@dataclass(frozen=True, slots=True)
class EarlierSchurCertificate:
    control_ordinal: int
    height: Silver
    norm: Silver
    control: Silver
    spatial_divisor: Silver
    actual_spatial_square: Silver
    actual_direction_bound: Silver
    perpendicular_residual: Silver
    threshold_gap: Silver
    monotonicity_discriminant: Silver
    rejection_discriminant: Silver


@dataclass(frozen=True, slots=True)
class Norm223LimiterAudit:
    limiting_height: Silver
    limiting_norm: Silver
    limiting_control: Silver
    limiting_target: Silver
    surviving_triple_count: int
    maximizing_types: tuple[LimitingMaximizer, ...]
    explicit_d_scaled: tuple[int, ...]
    explicit_e_scaled: tuple[int, ...]
    explicit_u_scaled: tuple[int, ...]
    explicit_v_scaled: tuple[int, ...]
    d_e_gram: tuple[tuple[int, int], tuple[int, int]]
    realized_u_v_gram: LeechGramTriple
    direction_pairing: Silver
    root_square: Silver
    primitive: bool
    reflection_integral: bool
    height_zero_empty: bool
    earlier_component_count: int
    earlier_schur_certificates: tuple[EarlierSchurCertificate, ...]
    tight_earlier_control_ordinal: int
    tight_earlier_rejection_discriminant: Silver
    same_height_q1_empty: bool
    positive_distinct_same_height_pairings: tuple[Silver, ...]
    same_height_difference_targets: tuple[Silver, ...]
    same_height_acute: bool
    accepted: bool


def _d4_pairings_at_value(
    triple: LeechGramTriple, target: Silver
) -> tuple[tuple[int, int, int], ...]:
    """Return every local-Gram extension attaining ``target``.

    Each row is ``(d.u,d.v,det Gram(d,u,v))``.  This independently replays
    the reduced-basis enumeration used by :func:`d4_pairing_bound`.
    """

    quotient = triple.u_dot_v // triple.v_square if triple.v_square else 0
    transforms = []
    for coefficient in range(quotient - 2, quotient + 3):
        w_square = triple.combination_square(1, -coefficient)
        if w_square >= 0:
            transforms.append(
                (
                    w_square,
                    abs(triple.u_dot_v - coefficient * triple.v_square),
                    abs(coefficient),
                    coefficient,
                    triple.u_dot_v - coefficient * triple.v_square,
                )
            )
    w_square, _abs_dot, _abs_coefficient, coefficient, w_dot_v = min(transforms)
    determinant = w_square * triple.v_square - w_dot_v * w_dot_v
    output: set[tuple[int, int, int]] = set()
    second_bound = isqrt(DIRECTION_SQUARE * triple.v_square)
    for second in range(-second_bound, second_bound + 1):
        if any(
            value != 0 and value < 4
            for value in (
                DIRECTION_SQUARE + triple.v_square + 2 * second,
                DIRECTION_SQUARE + triple.v_square - 2 * second,
            )
        ):
            continue
        if determinant == 0:
            first_values = (
                (w_dot_v * second // triple.v_square,)
                if (w_dot_v * second) % triple.v_square == 0
                else ()
            )
        else:
            radius = isqrt(
                determinant
                * (DIRECTION_SQUARE * triple.v_square - second * second)
            )
            lower_numerator = w_dot_v * second - radius
            upper_numerator = w_dot_v * second + radius
            lower = -((-lower_numerator) // triple.v_square)
            upper = upper_numerator // triple.v_square
            first_values = range(lower, upper + 1)
        for first in first_values:
            augmented_determinant = (
                DIRECTION_SQUARE * determinant
                - first * (first * triple.v_square - second * w_dot_v)
                + second * (first * w_dot_v - second * w_square)
            )
            if augmented_determinant < 0:
                continue
            if any(
                value != 0 and value < 4
                for value in (
                    DIRECTION_SQUARE + w_square + 2 * first,
                    DIRECTION_SQUARE + w_square - 2 * first,
                    DIRECTION_SQUARE + triple.v_square + 2 * second,
                    DIRECTION_SQUARE + triple.v_square - 2 * second,
                )
            ):
                continue
            d_dot_u = first + coefficient * second
            if Silver(d_dot_u, second) == target:
                output.add((d_dot_u, second, augmented_determinant))
    return tuple(sorted(output))


def _maximizing_types() -> tuple[LimitingMaximizer, ...]:
    survivors = tuple(
        triple
        for triple in _candidate_triples(LIMITING_TARGET)
        if not has_short_binary_combination(triple)
    )
    hits = tuple(
        triple for triple in survivors if d4_pairing_bound(triple) == LIMITING_MAXIMUM
    )
    if hits != (IMPOSSIBLE_TRIPLE, REALIZED_TRIPLE):
        raise AssertionError("the norm-223 maximizing Gram locus drifted")
    output = []
    for triple in hits:
        pairings = _d4_pairings_at_value(triple, LIMITING_MAXIMUM)
        if len(pairings) != 1 or pairings[0][:2] != (679, 479):
            raise AssertionError("a limiting augmented pairing ceased to be unique")
        d_dot_u, d_dot_v, determinant = pairings[0]
        e_square = triple.u_square - 2 * 170 * d_dot_u + 170 * 170 * 4
        f_square = triple.v_square - 2 * 120 * d_dot_v + 120 * 120 * 4
        e_dot_f = (
            triple.u_dot_v
            - 170 * d_dot_v
            - 120 * d_dot_u
            + 170 * 120 * 4
        )
        difference = e_square + f_square - 2 * e_dot_f
        if triple == IMPOSSIBLE_TRIPLE:
            realized = False
            obstruction = "(u-170d)-(v-120d) would have Leech norm 2"
        else:
            realized = True
            obstruction = None
        output.append(
            LimitingMaximizer(
                triple,
                d_dot_u,
                d_dot_v,
                determinant,
                e_square,
                f_square,
                e_dot_f,
                difference,
                realized,
                obstruction,
            )
        )
    if output[0].shifted_difference_square != 2:
        raise AssertionError("the positive-definite maximizing type lost its obstruction")
    if (
        output[1].shifted_e_square,
        output[1].shifted_f_square,
        output[1].shifted_e_dot_f,
        output[1].shifted_difference_square,
    ) != (6, 6, 6, 0):
        raise AssertionError("the singular maximizing type lost its equality relation")
    return tuple(output)


def _explicit_root_audit() -> tuple[bool, bool]:
    direction = norm_six_direction_audit()
    if direction["scaled_coordinates"] != LIMITING_E_SCALED:
        raise AssertionError("the norm-six Leech direction changed")
    if _scaled_leech_dot(LIMITING_D_SCALED, LIMITING_D_SCALED) != 4:
        raise AssertionError("the displayed type-(4^2) vector is not minimal")
    if _scaled_leech_dot(LIMITING_E_SCALED, LIMITING_E_SCALED) != 6:
        raise AssertionError("the displayed Leech vector does not have norm six")
    if _scaled_leech_dot(LIMITING_D_SCALED, LIMITING_E_SCALED) != -1:
        raise AssertionError("the displayed rank-two Leech Gram changed")
    realized = LeechGramTriple(
        _scaled_leech_dot(LIMITING_U_SCALED, LIMITING_U_SCALED),
        _scaled_leech_dot(LIMITING_V_SCALED, LIMITING_V_SCALED),
        _scaled_leech_dot(LIMITING_U_SCALED, LIMITING_V_SCALED),
    )
    if realized != REALIZED_TRIPLE:
        raise AssertionError("the explicit limiting realization changed")
    pairing = Silver(
        _scaled_leech_dot(LIMITING_D_SCALED, LIMITING_U_SCALED),
        _scaled_leech_dot(LIMITING_D_SCALED, LIMITING_V_SCALED),
    )
    if pairing != LIMITING_MAXIMUM:
        raise AssertionError("the explicit direction pairing changed")
    spatial_square = Silver(
        realized.u_square + 2 * realized.v_square,
        2 * realized.u_dot_v,
    )
    root_square = spatial_square - NORM223_COEFFICIENT * LIMITING_CONTROL**2
    if spatial_square != LIMITING_TARGET or root_square != LIMITING_NORM:
        raise AssertionError("the explicit vector ceased to be a norm-two root")

    # If the root were divisible by a nonunit, its square two permits only
    # the prime sqrt(2), hence u would lie in 2*Lambda.  But d.u=679 is odd.
    primitive = _scaled_leech_dot(LIMITING_D_SCALED, LIMITING_U_SCALED) % 2 == 1
    # For q=2, 2(z.r)/q=z.r belongs to O_K for every lattice vector z.
    reflection_integral = root_square == 2
    return primitive, reflection_integral


def _height_zero_empty() -> bool:
    prefix = norm223_leech_prefix_audit()
    return all(
        not _candidate_triples(
            row.canonical_norm / (row.odd_spatial_divisor**2)
        )
        for row in prefix.local_classes
    )


def _earlier_schur_certificates() -> tuple[EarlierSchurCertificate, ...]:
    """Bound the displayed root against every earlier necessary component."""

    prefix = norm223_leech_prefix_audit()
    c = LIMITING_MAXIMUM / 4
    beta = Silver(1, 1)
    perpendicular_square = Fraction(23, 4)  # (e+d/4)^2
    output = []
    for row in prefix.components:
        if row.height >= LIMITING_HEIGHT:
            continue
        divisor = row.odd_spatial_divisor
        spatial_square = row.target_after_gate * divisor**2
        direction_bound = row.maximum_pairing * divisor
        residual = spatial_square - direction_bound**2 / 4
        threshold_gap = (
            NORM223_COEFFICIENT * LIMITING_CONTROL * row.control
            - c * direction_bound
        )
        # f(D)=cD+beta*sqrt(K(N-D^2/4)) is increasing through D=M if
        # 16c^2N-M^2(4c^2+beta^2K)>0.  Its value at M is below the time
        # threshold precisely when the positive gap has square larger than
        # beta^2*K*residual.
        monotonicity = (
            16 * c**2 * spatial_square
            - direction_bound**2
            * (4 * c**2 + beta**2 * perpendicular_square)
        )
        rejection = (
            threshold_gap**2
            - beta**2 * perpendicular_square * residual
        )
        if any(
            value.sign(0) <= 0
            for value in (residual, threshold_gap, monotonicity, rejection)
        ):
            raise AssertionError(
                f"an earlier Schur certificate failed at {row.control_ordinal}"
            )
        output.append(
            EarlierSchurCertificate(
                row.control_ordinal,
                row.height,
                row.norm,
                row.control,
                divisor,
                spatial_square,
                direction_bound,
                residual,
                threshold_gap,
                monotonicity,
                rejection,
            )
        )
    if len(output) != 34:
        raise AssertionError("the earlier nonempty component count changed")
    tight = min(output, key=lambda row: row.rejection_discriminant)
    if (
        tight.control_ordinal != 565
        or tight.threshold_gap != Silver(Fraction(63, 4), Fraction(45, 4))
        or tight.perpendicular_residual != Silver(Fraction(55, 4), 9)
        or tight.rejection_discriminant != Silver(57, 41)
    ):
        raise AssertionError("the tight earlier-wall certificate changed")
    return tuple(output)


def _same_height_audit() -> tuple[bool, tuple[Silver, ...], tuple[Silver, ...]]:
    rows = tuple(
        row for row in norm223_control_skeleton() if row[0] == LIMITING_HEIGHT
    )
    if tuple((row[1], row[2]) for row in rows) != (
        (Silver(1), Silver(92, 65)),
        (LIMITING_NORM, LIMITING_CONTROL),
    ):
        raise AssertionError("the limiting height packet changed")
    q1_empty = not _candidate_triples(rows[0][3])
    if not q1_empty:
        raise AssertionError("the simultaneous q=1 component became nonempty")

    conjugate_target = LIMITING_TARGET.conjugate()
    if (Silver(7) - 4 * conjugate_target).sign(0) <= 0:
        raise AssertionError("the conjugate same-height Cauchy cap changed")
    possible = set()
    for rational in range(-4, 5):
        for radical in range(-4, 5):
            pairing = Silver(rational, radical)
            if (
                pairing.sign(0) > 0
                and (Silver(2) - pairing).sign(0) >= 0
                and (Silver(2) - pairing.conjugate()).sign(0) >= 0
                and (Silver(2) + pairing.conjugate()).sign(0) >= 0
            ):
                possible.add(pairing)
    if possible != {Silver(1), Silver(0, 1), Silver(2)}:
        raise AssertionError("the same-height pairing interval changed")
    positive_distinct = (Silver(1), Silver(0, 1))
    difference_targets = tuple(2 * (Silver(2) - value) for value in positive_distinct)
    if any(_candidate_triples(target) for target in difference_targets):
        raise AssertionError("a positive same-height difference shell became nonempty")
    # Pairing two is equality in the positive-definite spatial Cauchy bound,
    # so it represents the same oriented root, not a distinct wall.
    return q1_empty, positive_distinct, difference_targets


@lru_cache(maxsize=1)
def norm223_limiting_root_audit() -> Norm223LimiterAudit:
    prefix = norm223_leech_prefix_audit()
    limiting_component = next(
        row for row in prefix.components if row.control_ordinal == 1084
    )
    if (
        limiting_component.height != LIMITING_HEIGHT
        or limiting_component.norm != LIMITING_NORM
        or limiting_component.control != LIMITING_CONTROL
        or limiting_component.target_after_gate != LIMITING_TARGET
        or limiting_component.maximum_pairing != LIMITING_MAXIMUM
    ):
        raise AssertionError("the frozen necessary limiter changed")
    maximizers = _maximizing_types()
    primitive, reflection_integral = _explicit_root_audit()
    earlier = _earlier_schur_certificates()
    q1_empty, same_pairings, difference_targets = _same_height_audit()
    height_zero_empty = _height_zero_empty()
    accepted = (
        primitive
        and reflection_integral
        and height_zero_empty
        and all(row.rejection_discriminant.sign(0) > 0 for row in earlier)
        and q1_empty
        and not any(_candidate_triples(target) for target in difference_targets)
    )
    if not accepted:
        raise AssertionError("the explicit norm-223 limiting root was not accepted")
    return Norm223LimiterAudit(
        LIMITING_HEIGHT,
        LIMITING_NORM,
        LIMITING_CONTROL,
        LIMITING_TARGET,
        limiting_component.surviving_triple_count,
        maximizers,
        LIMITING_D_SCALED,
        LIMITING_E_SCALED,
        LIMITING_U_SCALED,
        LIMITING_V_SCALED,
        ((4, -1), (-1, 6)),
        REALIZED_TRIPLE,
        LIMITING_MAXIMUM,
        LIMITING_NORM,
        primitive,
        reflection_integral,
        height_zero_empty,
        len(earlier),
        earlier,
        565,
        Silver(57, 41),
        q1_empty,
        same_pairings,
        difference_targets,
        True,
        accepted,
    )


__all__ = [
    "IMPOSSIBLE_TRIPLE",
    "LIMITING_CONTROL",
    "LIMITING_D_SCALED",
    "LIMITING_E_SCALED",
    "LIMITING_HEIGHT",
    "LIMITING_MAXIMUM",
    "LIMITING_NORM",
    "LIMITING_TARGET",
    "LIMITING_U_SCALED",
    "LIMITING_V_SCALED",
    "REALIZED_TRIPLE",
    "EarlierSchurCertificate",
    "LimitingMaximizer",
    "Norm223LimiterAudit",
    "norm223_limiting_root_audit",
]
