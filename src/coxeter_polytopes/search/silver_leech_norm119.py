"""Exact necessary-wall prefix for the norm-119 silver Leech seed.

The lattice under study is

``L = <-a> + Lambda_24`` over ``Z[sqrt(2)]``, with ``a=3+8*sqrt(2)``.

This module deliberately enumerates a *necessary super-set* of Vinberg
walls.  It does not test Leech realization, primitivity, or Vinberg
acceptance.  Consequently a hyperideal point satisfying every inequality
enumerated here is a valid obstruction for the bounded Vinberg prefix, while
the surviving Gram triples themselves must not be advertised as roots.

The odd-prime normalization is kept explicit.  If ``y=d*w`` for the spatial
divisor belonging to a mirror-norm class, the affine root inequality is

``point.w <= a*x/d``

and the binary Leech target is ``w^2=(q+a*x^2)/d^2``.  This simultaneous
division is important when the canonical totally-positive generator ``q``
differs from ``d`` by a totally-positive unit square.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import isqrt

from coxeter_polytopes.exact.silver import SILVER_ONE, Silver
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps
from coxeter_polytopes.search.silver_leech import (
    LEECH_MINIMUM,
    LEECH_SHELL_COUNTS,
    LeechGramTriple,
    candidate_leech_gram_triples,
)
from coxeter_polytopes.search.silver_leech_manifest import (
    MirrorNormClass,
    has_short_binary_combination,
    mirror_norm_classes,
)


NORM119_COEFFICIENT = Silver(3, 8)
NORM119_ABSOLUTE_NORM = 119
P7_DIVISOR = Silver(3, 1)
P17_DIVISOR = Silver(5, 2)
P119_DIVISOR = P7_DIVISOR * P17_DIVISOR
SIGN_CHANGING_UNIT = Silver(-1, 1)
INVERSE_UNIT_SQUARE = Silver(3, -2)
DEFAULT_HEIGHT_BOUND = 10_000
SUPPORTED_DIRECTION_SQUARES = (4, 6, 8, 10)


@dataclass(frozen=True, slots=True)
class Norm119LocalClassAudit:
    canonical_norm: Silver
    dyadic_factor: int
    odd_spatial_divisor: Silver
    unit_square_multiplier: Silver


@dataclass(frozen=True, slots=True)
class Norm119DirectionComponentBound:
    direction_square: int
    maximum_pairing: Silver
    slab_ratio: Silver


@dataclass(frozen=True, slots=True)
class Norm119ControlAudit:
    ordinal: int
    height: Silver
    norm: Silver
    control: Silver
    target_before_gate: Silver
    odd_spatial_divisor: Silver
    target_after_gate: Silver
    normalized_wall_bound: Silver
    raw_candidate_triples: tuple[LeechGramTriple, ...]
    surviving_triples: tuple[LeechGramTriple, ...]
    direction_bounds: tuple[Norm119DirectionComponentBound, ...]
    nonempty_packet_ordinal: int | None


@dataclass(frozen=True, slots=True)
class Norm119DirectionPrefixAudit:
    direction_square: int
    limiting_control_ordinal: int
    limiting_height: Silver
    limiting_norm: Silver
    limiting_control: Silver
    limiting_coordinate: Silver
    limiting_hyperideal_excess: Silver
    strict_hyperideal_witness_exists: bool
    shrink_denominator: int | None
    witness_coordinate: Silver | None
    component_wall_margins: tuple[Silver, ...]
    minimum_wall_margin: Silver | None
    witness_hyperideal_excess: Silver | None


@dataclass(frozen=True, slots=True)
class Norm119LeechPrefixAudit:
    coefficient: Silver
    absolute_norm: int
    complete_height_bound: int
    coefficient_factorization_verified: bool
    local_classes: tuple[Norm119LocalClassAudit, ...]
    local_normalization_proof: str
    control_count: int
    empty_control_count: int
    nonempty_component_count: int
    nonempty_height_count: int
    raw_candidate_triple_count: int
    short_combination_exclusion_count: int
    surviving_triple_count: int
    controls: tuple[Norm119ControlAudit, ...]
    direction_audits: tuple[Norm119DirectionPrefixAudit, ...]


EXPECTED_LOCAL_CLASSES = (
    Norm119LocalClassAudit(Silver(1), 1, Silver(1), Silver(1)),
    Norm119LocalClassAudit(Silver(2), 2, Silver(1), Silver(1)),
    Norm119LocalClassAudit(P7_DIVISOR, 1, P7_DIVISOR, Silver(1)),
    Norm119LocalClassAudit(
        Silver(13, -5), 1, P119_DIVISOR, INVERSE_UNIT_SQUARE
    ),
    Norm119LocalClassAudit(P17_DIVISOR, 1, P17_DIVISOR, Silver(1)),
    Norm119LocalClassAudit(2 * P7_DIVISOR, 2, P7_DIVISOR, Silver(1)),
    Norm119LocalClassAudit(
        Silver(26, -10), 2, P119_DIVISOR, INVERSE_UNIT_SQUARE
    ),
    Norm119LocalClassAudit(2 * P17_DIVISOR, 2, P17_DIVISOR, Silver(1)),
)


LOCAL_NORMALIZATION_PROOF = (
    "a=(3+sqrt(2))(5+2sqrt(2))(-1+sqrt(2)); the two totally-positive "
    "odd factors have norms 7 and 17, while -1+sqrt(2) is a unit of norm "
    "-1. At either odd prime dividing (a), reflection integrality and local "
    "primitivity force the mirror-norm valuation to be either zero, with a "
    "spatially primitive vector, or one, with y in that prime times Lambda. "
    "The choices are independent. At the ramified dyadic prime, even "
    "unimodularity gives valuations zero or two, represented by factors 1 "
    "and 2. For the product divisor d=(3+sqrt(2))(5+2sqrt(2))="
    "19+11sqrt(2), its Pell-domain norm generator is "
    "d*(1+sqrt(2))^(-2)=13-5sqrt(2). Writing y=d*w therefore gives exactly "
    "w^2=(q+a*x^2)/d^2 and point.w<=a*x/d."
)


def norm119_local_classes() -> tuple[Norm119LocalClassAudit, ...]:
    """Recheck the factorization and all eight normalized mirror classes."""

    if NORM119_COEFFICIENT.norm() != -NORM119_ABSOLUTE_NORM:
        raise AssertionError("the norm-119 coefficient norm drifted")
    if P7_DIVISOR.norm() != 7 or P17_DIVISOR.norm() != 17:
        raise AssertionError("an odd prime divisor norm drifted")
    if P119_DIVISOR != Silver(19, 11) or P119_DIVISOR.norm() != 119:
        raise AssertionError("the combined odd divisor drifted")
    if NORM119_COEFFICIENT != P119_DIVISOR * SIGN_CHANGING_UNIT:
        raise AssertionError("the coefficient factorization drifted")
    if INVERSE_UNIT_SQUARE != Silver(1, 1) ** -2:
        raise AssertionError("the canonical unit-square normalization drifted")

    manifest_classes = mirror_norm_classes(NORM119_COEFFICIENT)
    expected_manifest_classes = tuple(
        MirrorNormClass(row.canonical_norm, row.odd_spatial_divisor)
        for row in EXPECTED_LOCAL_CLASSES
    )
    if manifest_classes != expected_manifest_classes:
        raise AssertionError("the independently expected local classes drifted")
    for row in EXPECTED_LOCAL_CLASSES:
        reconstructed = (
            row.dyadic_factor
            * row.odd_spatial_divisor
            * row.unit_square_multiplier
        )
        if reconstructed != row.canonical_norm:
            raise AssertionError("a local norm/divisor normalization is inconsistent")
        if not row.unit_square_multiplier.is_integral:
            raise AssertionError("a unit-square multiplier ceased to be integral")
        if row.unit_square_multiplier.norm() != 1:
            raise AssertionError("a unit-square multiplier ceased to have norm one")
        if any(row.unit_square_multiplier.sign(place) <= 0 for place in (0, 1)):
            raise AssertionError("a unit-square multiplier ceased to be positive")
    return EXPECTED_LOCAL_CLASSES


def _fast_fixed_direction_pairing_bound(
    triple: LeechGramTriple, direction_square: int
) -> Silver:
    """Exact 3x3-Gram/minimum-four bound, using an interval enumeration.

    For a reduced basis ``(w,v)=(u-kv,v)``, write ``f=d.w`` and ``s=d.v``.
    Positive semidefiniteness of their Gram matrix is equivalent to

    ``(V*f-C*s)^2 <= det(w,v)*(D*V-s^2)``.

    This gives a complete short integer interval for ``f`` at each ``s`` and
    is algebraically identical to the rectangular oracle in
    ``silver_leech_manifest``, but is much faster on the later near-rank-one
    packets.
    """

    if direction_square < LEECH_MINIMUM or direction_square % 2:
        raise ValueError("expected an even Leech direction norm at least four")
    v_square = triple.v_square
    dot = triple.u_dot_v
    quotient = dot // v_square if v_square else 0
    transformations = []
    for coefficient in range(quotient - 2, quotient + 3):
        w_square = triple.combination_square(1, -coefficient)
        if w_square < 0:
            continue
        w_dot_v = dot - coefficient * v_square
        transformations.append(
            (
                w_square,
                abs(w_dot_v),
                abs(coefficient),
                coefficient,
                w_dot_v,
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

    pairs: list[tuple[int, int]] = []
    if v_square == 0:
        if w_dot_v:
            raise AssertionError("a zero-norm Leech vector has nonzero pairing")
        bound = isqrt(direction_square * w_square)
        pairs.extend((first, 0) for first in range(-bound, bound + 1))
    elif determinant == 0:
        second_bound = isqrt(direction_square * v_square)
        pairs.extend(
            (w_dot_v * second // v_square, second)
            for second in range(-second_bound, second_bound + 1)
            if (w_dot_v * second) % v_square == 0
        )
    else:
        second_bound = isqrt(direction_square * v_square)
        for second in range(-second_bound, second_bound + 1):
            radicand = determinant * (
                direction_square * v_square - second * second
            )
            radius = isqrt(radicand)
            lower_numerator = w_dot_v * second - radius
            upper_numerator = w_dot_v * second + radius
            lower = -((-lower_numerator) // v_square)
            upper = upper_numerator // v_square
            pairs.extend((first, second) for first in range(lower, upper + 1))

    best: Silver | None = None
    for first, second in pairs:
        gram_determinant = (
            direction_square * determinant
            - first * (first * v_square - second * w_dot_v)
            + second * (first * w_dot_v - second * w_square)
        )
        if gram_determinant < 0:
            raise AssertionError("the interval enumeration admitted a negative Gram")
        adjacent_norms = (
            direction_square + w_square + 2 * first,
            direction_square + w_square - 2 * first,
            direction_square + v_square + 2 * second,
            direction_square + v_square - 2 * second,
        )
        if any(
            value != 0 and value < LEECH_MINIMUM for value in adjacent_norms
        ):
            continue
        pairing = Silver(first + coefficient * second, second)
        if best is None or pairing > best:
            best = pairing
    if best is None:
        raise AssertionError("no fixed-direction pairing data survived")
    return best


def _enumerate_control_skeleton(
    height_bound: int,
) -> tuple[tuple[Silver, Silver, Silver, Silver, Silver], ...]:
    """Independently enumerate every arithmetic control below the cutoff."""

    if height_bound <= 0:
        raise ValueError("height_bound must be positive")
    coefficient = NORM119_COEFFICIENT
    levels: set[tuple[Silver, Silver, Silver, Silver, Silver]] = set()
    for norm_class in norm119_local_classes():
        norm = norm_class.canonical_norm
        for control in _bounded_for_caps(
            (norm * height_bound / coefficient, norm / (-coefficient))
        ):
            if control.sign(0) <= 0:
                continue
            if not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if target.sign(1) < 0:
                continue
            height = coefficient * control * control / norm
            if (Silver(height_bound) - height).sign(0) < 0:
                continue
            levels.add(
                (
                    height,
                    norm,
                    control,
                    target,
                    norm_class.odd_spatial_divisor,
                )
            )
    return tuple(sorted(levels, key=lambda row: (row[0], row[1], row[2], row[4])))


def _smallest_strict_shrink_denominator(
    direction_square: int, limiting: Silver
) -> int:
    def is_hyperideal(denominator: int) -> bool:
        witness = limiting * (denominator - 1) / denominator
        excess = (
            direction_square * witness * witness - NORM119_COEFFICIENT
        )
        return excess.sign(0) > 0

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
def norm119_leech_prefix_audit(
    height_bound: int = DEFAULT_HEIGHT_BOUND,
) -> Norm119LeechPrefixAudit:
    """Enumerate and certify the complete necessary prefix exactly."""

    local_classes = norm119_local_classes()
    skeleton = _enumerate_control_skeleton(height_bound)
    provisional: list[
        tuple[
            int,
            Silver,
            Silver,
            Silver,
            Silver,
            Silver,
            Silver,
            Silver,
            tuple[LeechGramTriple, ...],
            tuple[LeechGramTriple, ...],
            tuple[Norm119DirectionComponentBound, ...],
        ]
    ] = []
    nonempty_heights: list[Silver] = []
    for ordinal, (height, norm, control, target, divisor) in enumerate(skeleton):
        normalized_target = target / (divisor * divisor)
        raw_triples = candidate_leech_gram_triples(normalized_target)
        surviving_triples = tuple(
            triple
            for triple in raw_triples
            if not has_short_binary_combination(triple)
        )
        wall_bound = NORM119_COEFFICIENT * control / divisor
        direction_bounds = ()
        if surviving_triples:
            if not nonempty_heights or nonempty_heights[-1] != height:
                nonempty_heights.append(height)
            direction_bounds = tuple(
                Norm119DirectionComponentBound(
                    direction_square=direction_square,
                    maximum_pairing=(
                        max(
                            _fast_fixed_direction_pairing_bound(
                                triple, direction_square
                            )
                            for triple in surviving_triples
                        )
                    ),
                    slab_ratio=SILVER_ONE,
                )
                for direction_square in SUPPORTED_DIRECTION_SQUARES
                if LEECH_SHELL_COUNTS.get(direction_square, 0) > 0
            )
            if len(direction_bounds) != len(SUPPORTED_DIRECTION_SQUARES):
                raise AssertionError("a requested Leech direction shell is empty")
            direction_bounds = tuple(
                Norm119DirectionComponentBound(
                    direction_square=row.direction_square,
                    maximum_pairing=row.maximum_pairing,
                    slab_ratio=wall_bound / row.maximum_pairing,
                )
                for row in direction_bounds
            )
        provisional.append(
            (
                ordinal,
                height,
                norm,
                control,
                target,
                divisor,
                normalized_target,
                wall_bound,
                raw_triples,
                surviving_triples,
                direction_bounds,
            )
        )

    packet_ordinals = {height: index for index, height in enumerate(nonempty_heights)}
    controls = tuple(
        Norm119ControlAudit(
            ordinal=ordinal,
            height=height,
            norm=norm,
            control=control,
            target_before_gate=target,
            odd_spatial_divisor=divisor,
            target_after_gate=normalized_target,
            normalized_wall_bound=wall_bound,
            raw_candidate_triples=raw_triples,
            surviving_triples=surviving_triples,
            direction_bounds=direction_bounds,
            nonempty_packet_ordinal=(
                packet_ordinals[height] if surviving_triples else None
            ),
        )
        for (
            ordinal,
            height,
            norm,
            control,
            target,
            divisor,
            normalized_target,
            wall_bound,
            raw_triples,
            surviving_triples,
            direction_bounds,
        ) in provisional
    )
    nonempty_controls = tuple(row for row in controls if row.surviving_triples)
    if not nonempty_controls:
        raise AssertionError("the norm-119 prefix unexpectedly has no necessary walls")

    direction_audits = []
    for direction_square in SUPPORTED_DIRECTION_SQUARES:
        rows = tuple(
            (
                next(
                    bound
                    for bound in row.direction_bounds
                    if bound.direction_square == direction_square
                ).slab_ratio,
                row,
                next(
                    bound
                    for bound in row.direction_bounds
                    if bound.direction_square == direction_square
                ).maximum_pairing,
            )
            for row in nonempty_controls
        )
        limiting, limiting_row, _limiting_maximum = min(
            rows, key=lambda item: item[0]
        )
        limiting_excess = (
            direction_square * limiting * limiting - NORM119_COEFFICIENT
        )
        if limiting_excess.sign(0) > 0:
            shrink_denominator = _smallest_strict_shrink_denominator(
                direction_square, limiting
            )
            witness = limiting * (shrink_denominator - 1) / shrink_denominator
            margins = tuple(
                row.normalized_wall_bound - witness * maximum
                for _ratio, row, maximum in rows
            )
            witness_excess = (
                direction_square * witness * witness - NORM119_COEFFICIENT
            )
            if any(margin.sign(0) <= 0 for margin in margins):
                raise AssertionError("a strict norm-119 wall margin was lost")
            if witness_excess.sign(0) <= 0:
                raise AssertionError("the norm-119 witness ceased to be hyperideal")
            minimum_margin: Silver | None = min(margins)
        else:
            shrink_denominator = None
            witness = None
            margins = ()
            minimum_margin = None
            witness_excess = None
        direction_audits.append(
            Norm119DirectionPrefixAudit(
                direction_square=direction_square,
                limiting_control_ordinal=limiting_row.ordinal,
                limiting_height=limiting_row.height,
                limiting_norm=limiting_row.norm,
                limiting_control=limiting_row.control,
                limiting_coordinate=limiting,
                limiting_hyperideal_excess=limiting_excess,
                strict_hyperideal_witness_exists=limiting_excess.sign(0) > 0,
                shrink_denominator=shrink_denominator,
                witness_coordinate=witness,
                component_wall_margins=margins,
                minimum_wall_margin=minimum_margin,
                witness_hyperideal_excess=witness_excess,
            )
        )

    raw_count = sum(len(row.raw_candidate_triples) for row in controls)
    surviving_count = sum(len(row.surviving_triples) for row in controls)
    return Norm119LeechPrefixAudit(
        coefficient=NORM119_COEFFICIENT,
        absolute_norm=NORM119_ABSOLUTE_NORM,
        complete_height_bound=height_bound,
        coefficient_factorization_verified=True,
        local_classes=local_classes,
        local_normalization_proof=LOCAL_NORMALIZATION_PROOF,
        control_count=len(controls),
        empty_control_count=sum(not row.surviving_triples for row in controls),
        nonempty_component_count=len(nonempty_controls),
        nonempty_height_count=len(nonempty_heights),
        raw_candidate_triple_count=raw_count,
        short_combination_exclusion_count=raw_count - surviving_count,
        surviving_triple_count=surviving_count,
        controls=controls,
        direction_audits=tuple(direction_audits),
    )


def strongest_norm119_prefix_witness(
    height_bound: int = DEFAULT_HEIGHT_BOUND,
) -> Norm119DirectionPrefixAudit:
    """Return the surviving direction with greatest limiting clearance."""

    candidates = tuple(
        row
        for row in norm119_leech_prefix_audit(height_bound).direction_audits
        if row.strict_hyperideal_witness_exists
    )
    if not candidates:
        raise AssertionError("no fixed-direction witness survives this prefix")
    return max(candidates, key=lambda row: row.limiting_hyperideal_excess)
