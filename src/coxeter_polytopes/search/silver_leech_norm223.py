"""Exact necessary-wall prefix for the norm-223 silver Leech seed.

The lattice under study is ``<-a> + Lambda_24`` over ``Z[sqrt(2)]``, where
``a=-13+14*sqrt(2)=(15+sqrt(2))*(sqrt(2)-1)``.  The routines below enumerate
the complete arithmetic control skeleton through a rational height cutoff and
apply the inert/split-prime spatial gate before testing binary Leech Grams.

This is deliberately a necessary *super-set* of Vinberg walls.  A surviving
Gram triple need not be realized by a Leech pair, primitive, or accepted by
Vinberg's algorithm.  Consequently a point satisfying every inequality here
is a rigorous bounded-prefix obstruction, while a new limiting component is
only a threat until those additional realization and acceptance tests pass.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import isqrt

from coxeter_polytopes.exact.silver import SILVER_ONE, Silver
from coxeter_polytopes.search.silver_leech import LEECH_MINIMUM, LeechGramTriple
from coxeter_polytopes.search.silver_leech_manifest import (
    MirrorNormClass,
    has_short_binary_combination,
    mirror_norm_classes,
)


NORM223_COEFFICIENT = Silver(-13, 14)
P223_DIVISOR = Silver(15, 1)
SIGN_CHANGING_UNIT = Silver(-1, 1)
DEFAULT_HEIGHT_BOUND = 250_000
DIRECTION_SQUARE = 4



@dataclass(frozen=True, slots=True)
class Norm223ComponentAudit:
    control_ordinal: int
    packet_ordinal: int
    height: Silver
    norm: Silver
    control: Silver
    target_before_gate: Silver
    odd_spatial_divisor: Silver
    target_after_gate: Silver
    raw_triple_count: int
    short_combination_exclusion_count: int
    surviving_triple_count: int
    maximum_pairing: Silver
    normalized_wall_bound: Silver
    slab_ratio: Silver


@dataclass(frozen=True, slots=True)
class Norm223PrefixAudit:
    complete_height_bound: int
    local_classes: tuple[MirrorNormClass, ...]
    control_count: int
    control_skeleton_sha256: str
    empty_control_count: int
    nonempty_component_count: int
    nonempty_height_count: int
    raw_candidate_triple_count: int
    short_combination_exclusion_count: int
    surviving_triple_count: int
    components: tuple[Norm223ComponentAudit, ...]
    limiting_control_ordinal: int
    limiting_height: Silver
    limiting_norm: Silver
    limiting_control: Silver
    limiting_coordinate: Silver
    limiting_hyperideal_excess: Silver
    shrink_denominator: int
    witness_coordinate: Silver
    minimum_wall_margin: Silver
    witness_hyperideal_excess: Silver


def norm223_local_classes() -> tuple[MirrorNormClass, ...]:
    """Return and verify the four normalized local mirror-norm classes."""

    if NORM223_COEFFICIENT != P223_DIVISOR * SIGN_CHANGING_UNIT:
        raise AssertionError("the norm-223 coefficient factorization drifted")
    if P223_DIVISOR.norm() != 223 or NORM223_COEFFICIENT.norm() != -223:
        raise AssertionError("the norm-223 field norms drifted")
    rows = mirror_norm_classes(NORM223_COEFFICIENT)
    expected = (
        MirrorNormClass(Silver(1), SILVER_ONE),
        MirrorNormClass(Silver(2), SILVER_ONE),
        MirrorNormClass(P223_DIVISOR, P223_DIVISOR),
        MirrorNormClass(2 * P223_DIVISOR, P223_DIVISOR),
    )
    if rows != expected:
        raise AssertionError("the norm-223 local classes drifted")
    return rows


def _embedding_square_root_upper(value: Silver, embedding: int) -> int:
    """Least integer ``n`` with ``n^2 >= value`` at one real place."""

    if value.sign(embedding) < 0:
        raise ValueError("a Minkowski square cap is negative")
    lower = 0
    upper = 1
    while (Silver(upper * upper) - value).sign(embedding) < 0:
        upper *= 2
    while upper - lower > 1:
        midpoint = (lower + upper) // 2
        if (Silver(midpoint * midpoint) - value).sign(embedding) >= 0:
            upper = midpoint
        else:
            lower = midpoint
    return upper


def _floor_n_sqrt2(value: int) -> int:
    """Return ``floor(value*sqrt(2))`` without floating-point arithmetic."""

    if value >= 0:
        return isqrt(2 * value * value)
    # sqrt(2) is irrational, so the nonzero negative case is never integral.
    return -isqrt(2 * value * value) - 1


@lru_cache(maxsize=None)
def norm223_control_skeleton(
    height_bound: int = DEFAULT_HEIGHT_BOUND,
) -> tuple[tuple[Silver, Silver, Silver, Silver, Silver], ...]:
    """Enumerate every locally normalized arithmetic control below ``H``.

    The thin Minkowski strip is exact.  If ``x=m+n*sqrt(2)``, the conjugate
    square cap bounds ``m-n*sqrt(2)``; ``_floor_n_sqrt2`` then supplies the
    finite integral interval for ``m`` without a floating-point decision.
    """

    if height_bound <= 0:
        raise ValueError("height_bound must be positive")
    coefficient = NORM223_COEFFICIENT
    output: set[tuple[Silver, Silver, Silver, Silver, Silver]] = set()
    for norm_class in norm223_local_classes():
        norm = norm_class.canonical_norm
        identity_cap = norm * height_bound / coefficient
        conjugate_cap = norm / (-coefficient)
        identity_bound = _embedding_square_root_upper(identity_cap, 0)
        conjugate_bound = _embedding_square_root_upper(conjugate_cap, 1)
        radical_bound = (identity_bound + conjugate_bound + 1) // 2
        for radical in range(-radical_bound, radical_bound + 1):
            center = _floor_n_sqrt2(radical)
            # Two extra integers on either side make the interval generous;
            # the exact cap tests below make the generosity harmless.
            for rational in range(
                center - conjugate_bound - 2,
                center + conjugate_bound + 3,
            ):
                control = Silver(rational, radical)
                control_square = control * control
                if control.sign(0) <= 0:
                    continue
                if (identity_cap - control_square).sign(0) < 0:
                    continue
                if (conjugate_cap - control_square).sign(1) < 0:
                    continue
                if not (2 * coefficient * control / norm).is_integral:
                    continue
                target = norm + coefficient * control_square
                if target.sign(1) < 0:
                    continue
                height = coefficient * control_square / norm
                if height > Silver(height_bound):
                    continue
                output.add(
                    (
                        height,
                        norm,
                        control,
                        target,
                        norm_class.odd_spatial_divisor,
                    )
                )
    return tuple(sorted(output, key=lambda row: (row[0], row[1], row[2], row[4])))


def _control_skeleton_sha256(
    skeleton: tuple[tuple[Silver, Silver, Silver, Silver, Silver], ...],
) -> str:
    def encode(value: Silver) -> tuple[str, str]:
        return str(value.a), str(value.b)

    payload = [
        [encode(height), encode(norm), encode(control), encode(target), encode(divisor)]
        for height, norm, control, target, divisor in skeleton
    ]
    return sha256(
        (json.dumps(payload, separators=(",", ":")) + "\n").encode()
    ).hexdigest()


@lru_cache(maxsize=None)
def _candidate_triples(target: Silver) -> tuple[LeechGramTriple, ...]:
    """Fast exact equivalent of the ordinary binary Leech Gram census."""

    if not target.is_integral:
        return ()
    dot = target.b / 2
    if dot.denominator != 1:
        return ()
    rational_target = int(target.a)
    dot_int = int(dot)
    if rational_target < 0 or rational_target % 2:
        return ()
    discriminant = rational_target * rational_target - 8 * dot_int * dot_int
    if discriminant < 0:
        return ()
    radical = isqrt(discriminant)
    # Cauchy's inequality confines V to the interval between the roots of
    # 2V^2-A*V+C^2.  The two-unit padding is followed by exact tests.
    lower = max(0, (rational_target - radical) // 4 - 2)
    upper = min(rational_target // 2, (rational_target + radical) // 4 + 2)
    output = []
    for v_square in range(lower, upper + 1):
        if v_square % 2:
            continue
        u_square = rational_target - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if (u_square and u_square < LEECH_MINIMUM) or (
            v_square and v_square < LEECH_MINIMUM
        ):
            continue
        if dot_int * dot_int > u_square * v_square:
            continue
        output.append(LeechGramTriple(u_square, v_square, dot_int))
    return tuple(output)


def _silver_integer_sign(rational: int, radical: int) -> int:
    """Sign of ``rational+radical*sqrt(2)`` using integer arithmetic."""

    if rational == 0:
        return (radical > 0) - (radical < 0)
    if radical == 0:
        return (rational > 0) - (rational < 0)
    if (rational > 0) == (radical > 0):
        return (rational > 0) - (rational < 0)
    norm = rational * rational - 2 * radical * radical
    norm_sign = (norm > 0) - (norm < 0)
    return norm_sign if rational > 0 else -norm_sign


@lru_cache(maxsize=None)
def d4_pairing_bound(triple: LeechGramTriple) -> Silver:
    """Exact minimum-four/3x3-Gram maximum for a norm-four direction.

    This is the streamed, integer-only specialization of the general pairing
    oracle.  For each second coordinate, the pairing increases with the first,
    so only the largest admissible first coordinate in its exact PSD interval
    can be optimal.
    """

    direction_square = DIRECTION_SQUARE
    v_square = triple.v_square
    dot = triple.u_dot_v
    quotient = dot // v_square if v_square else 0
    transformations = []
    for coefficient in range(quotient - 2, quotient + 3):
        w_square = triple.combination_square(1, -coefficient)
        if w_square >= 0:
            transformations.append(
                (
                    w_square,
                    abs(dot - coefficient * v_square),
                    abs(coefficient),
                    coefficient,
                    dot - coefficient * v_square,
                )
            )
    if not transformations:
        raise AssertionError("no reduced binary Gram basis survived")
    w_square, _abs_dot, _abs_coefficient, coefficient, w_dot_v = min(
        transformations
    )
    determinant = w_square * v_square - w_dot_v * w_dot_v
    if determinant < 0:
        raise AssertionError("the binary Gram matrix is not semidefinite")

    best: tuple[int, int] | None = None

    def consider(first: int, second: int) -> bool:
        nonlocal best
        gram_determinant = (
            direction_square * determinant
            - first * (first * v_square - second * w_dot_v)
            + second * (first * w_dot_v - second * w_square)
        )
        if gram_determinant < 0:
            return False
        adjacent_norms = (
            direction_square + w_square + 2 * first,
            direction_square + w_square - 2 * first,
            direction_square + v_square + 2 * second,
            direction_square + v_square - 2 * second,
        )
        if any(value != 0 and value < LEECH_MINIMUM for value in adjacent_norms):
            return False
        candidate = first + coefficient * second, second
        if best is None or _silver_integer_sign(
            candidate[0] - best[0], candidate[1] - best[1]
        ) > 0:
            best = candidate
        return True

    if v_square == 0:
        if w_dot_v:
            raise AssertionError("a zero-norm Leech vector has nonzero pairing")
        first_bound = isqrt(direction_square * w_square)
        for first in range(-first_bound, first_bound + 1):
            consider(first, 0)
    elif determinant == 0:
        second_bound = isqrt(direction_square * v_square)
        for second in range(-second_bound, second_bound + 1):
            if (w_dot_v * second) % v_square == 0:
                consider(w_dot_v * second // v_square, second)
    else:
        second_bound = isqrt(direction_square * v_square)
        for second in range(-second_bound, second_bound + 1):
            if any(
                value != 0 and value < LEECH_MINIMUM
                for value in (
                    direction_square + v_square + 2 * second,
                    direction_square + v_square - 2 * second,
                )
            ):
                continue
            radius = isqrt(
                determinant * (direction_square * v_square - second * second)
            )
            lower_numerator = w_dot_v * second - radius
            upper_numerator = w_dot_v * second + radius
            lower = -((-lower_numerator) // v_square)
            upper = upper_numerator // v_square
            # Only finitely many first coordinates can fail the minimum-four
            # adjacency gate, so this normally accepts on the first attempt.
            for first in range(upper, lower - 1, -1):
                if consider(first, second):
                    break
    if best is None:
        raise AssertionError("no norm-four fixed-direction pairing survived")
    return Silver(best[0], best[1])


def _smallest_strict_shrink_denominator(limiting: Silver) -> int:
    def is_hyperideal(denominator: int) -> bool:
        witness = limiting * Fraction(denominator - 1, denominator)
        return (
            DIRECTION_SQUARE * witness * witness - NORM223_COEFFICIENT
        ).sign(0) > 0

    lower = 1
    upper = 2
    while not is_hyperideal(upper):
        lower = upper
        upper *= 2
    while upper - lower > 1:
        midpoint = (lower + upper) // 2
        if is_hyperideal(midpoint):
            upper = midpoint
        else:
            lower = midpoint
    return upper


@lru_cache(maxsize=None)
def norm223_leech_prefix_audit(
    height_bound: int = DEFAULT_HEIGHT_BOUND,
) -> Norm223PrefixAudit:
    """Certify the complete norm-223 necessary prefix exactly."""

    skeleton = norm223_control_skeleton(height_bound)
    provisional = []
    nonempty_heights: list[Silver] = []
    raw_total = 0
    surviving_total = 0
    for ordinal, (height, norm, control, target, divisor) in enumerate(skeleton):
        reduced_target = target / (divisor * divisor)
        raw = _candidate_triples(reduced_target)
        surviving = tuple(
            triple for triple in raw if not has_short_binary_combination(triple)
        )
        raw_total += len(raw)
        surviving_total += len(surviving)
        if not surviving:
            continue
        if not nonempty_heights or height != nonempty_heights[-1]:
            nonempty_heights.append(height)
        maximum = max(d4_pairing_bound(triple) for triple in surviving)
        wall_bound = NORM223_COEFFICIENT * control / divisor
        provisional.append(
            (
                ordinal,
                height,
                norm,
                control,
                target,
                divisor,
                reduced_target,
                len(raw),
                len(surviving),
                maximum,
                wall_bound,
                wall_bound / maximum,
            )
        )
    if not provisional:
        raise AssertionError("the norm-223 prefix has no necessary components")
    packet_ordinals = {height: index for index, height in enumerate(nonempty_heights)}
    components = tuple(
        Norm223ComponentAudit(
            control_ordinal=ordinal,
            packet_ordinal=packet_ordinals[height],
            height=height,
            norm=norm,
            control=control,
            target_before_gate=target,
            odd_spatial_divisor=divisor,
            target_after_gate=reduced_target,
            raw_triple_count=raw_count,
            short_combination_exclusion_count=raw_count - surviving_count,
            surviving_triple_count=surviving_count,
            maximum_pairing=maximum,
            normalized_wall_bound=wall_bound,
            slab_ratio=ratio,
        )
        for (
            ordinal,
            height,
            norm,
            control,
            target,
            divisor,
            reduced_target,
            raw_count,
            surviving_count,
            maximum,
            wall_bound,
            ratio,
        ) in provisional
    )
    limiting = min(components, key=lambda row: row.slab_ratio)
    limiting_coordinate = limiting.slab_ratio
    limiting_excess = (
        DIRECTION_SQUARE * limiting_coordinate * limiting_coordinate
        - NORM223_COEFFICIENT
    )
    if limiting_excess.sign(0) <= 0:
        raise AssertionError("the norm-223 necessary prefix reached the light cone")
    shrink_denominator = _smallest_strict_shrink_denominator(limiting_coordinate)
    witness = limiting_coordinate * Fraction(
        shrink_denominator - 1, shrink_denominator
    )
    margins = tuple(
        row.normalized_wall_bound - witness * row.maximum_pairing
        for row in components
    )
    witness_excess = (
        DIRECTION_SQUARE * witness * witness - NORM223_COEFFICIENT
    )
    if any(margin.sign(0) <= 0 for margin in margins):
        raise AssertionError("a strict norm-223 wall margin was lost")
    if witness_excess.sign(0) <= 0:
        raise AssertionError("the norm-223 witness ceased to be hyperideal")
    return Norm223PrefixAudit(
        complete_height_bound=height_bound,
        local_classes=norm223_local_classes(),
        control_count=len(skeleton),
        control_skeleton_sha256=_control_skeleton_sha256(skeleton),
        empty_control_count=len(skeleton) - len(components),
        nonempty_component_count=len(components),
        nonempty_height_count=len(nonempty_heights),
        raw_candidate_triple_count=raw_total,
        short_combination_exclusion_count=raw_total - surviving_total,
        surviving_triple_count=surviving_total,
        components=components,
        limiting_control_ordinal=limiting.control_ordinal,
        limiting_height=limiting.height,
        limiting_norm=limiting.norm,
        limiting_control=limiting.control,
        limiting_coordinate=limiting_coordinate,
        limiting_hyperideal_excess=limiting_excess,
        shrink_denominator=shrink_denominator,
        witness_coordinate=witness,
        minimum_wall_margin=min(margins),
        witness_hyperideal_excess=witness_excess,
    )
