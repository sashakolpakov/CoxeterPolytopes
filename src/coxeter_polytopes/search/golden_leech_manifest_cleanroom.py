"""Clean-room replay of the norm-128 golden Leech first-prefix manifest.

This module deliberately does not import :mod:`golden_leech_manifest`.  It
reconstructs the prime-ideal directions, primitive local mirror classes,
first necessary binary-Leech packets, and strict fixed-direction witnesses
directly in ``Z[phi]``.  The witnesses certify only the finite prefix through
the first necessary packet; later Vinberg walls may cut the displayed open
hyperideal patch.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gcd, isqrt

from coxeter_polytopes.exact.golden import ONE, PHI, Golden


UNIT_SQUARE = PHI * PHI
LEECH_MINIMUM = 4
LEECH_DIRECTION_NORMS = (4, 6, 8, 10)
NORM_BOUND = 128
HEIGHT_BOUND = 5001

# Standard MOG coordinates are divided by sqrt(8).  P and R are type-(4^2)
# minimal vectors and Q is supported on the standard Golay octad
# {1,2,6,7,8,10,12,23}.  These rows explicitly realize norms 4,6,8,10.
_MOG_P = (-4, -4) + (0,) * 22
_MOG_R = (0, 0, -4, -4) + (0,) * 20
_MOG_Q = (
    0, 2, 2, 0, 0, 0, 2, 2, 2, 0, 2, 0,
    2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
)
EXPLICIT_DIRECTION_VECTORS_SCALED: tuple[tuple[int, ...], ...] = (
    _MOG_P,
    tuple(left + right for left, right in zip(_MOG_P, _MOG_Q)),
    tuple(left + right for left, right in zip(_MOG_P, _MOG_R)),
    tuple(left - right for left, right in zip(_MOG_P, _MOG_Q)),
)


def mog_scaled_square(vector: tuple[int, ...]) -> int:
    numerator = sum(coordinate * coordinate for coordinate in vector)
    if numerator % 8:
        raise ValueError("nonintegral MOG square")
    return numerator // 8


@dataclass(frozen=True, slots=True)
class CleanGoldenAtom:
    label: str
    ideal_norm: int
    generator: Golden
    rational_prime: int
    splitting_type: str


@dataclass(frozen=True, slots=True)
class CleanGoldenDirection:
    coefficient: Golden
    ideal_norm: int
    atom_labels: tuple[str, ...]
    atom_generators: tuple[Golden, ...]


@dataclass(frozen=True, slots=True)
class CleanGoldenMirrorClass:
    canonical_norm: Golden
    odd_spatial_divisor: Golden


@dataclass(frozen=True, slots=True)
class CleanGoldenTriple:
    u_square: int
    v_square: int
    u_dot_v: int

    def square(self, first: int, second: int) -> int:
        return (
            self.u_square * first * first
            + 2 * self.u_dot_v * first * second
            + self.v_square * second * second
        )


@dataclass(frozen=True, slots=True)
class CleanGoldenComponent:
    norm: Golden
    control: Golden
    target_before_gate: Golden
    odd_spatial_divisor: Golden
    target_after_gate: Golden
    surviving_triples: tuple[CleanGoldenTriple, ...]


@dataclass(frozen=True, slots=True)
class CleanGoldenRow:
    coefficient: Golden
    ideal_norm: int
    atom_labels: tuple[str, ...]
    mirror_norm_classes: tuple[CleanGoldenMirrorClass, ...]
    first_necessary_height: Golden
    first_packet_components: tuple[CleanGoldenComponent, ...]


@dataclass(frozen=True, slots=True)
class CleanGoldenPrefixWitness:
    coefficient: Golden
    ideal_norm: int
    atom_labels: tuple[str, ...]
    first_height: Golden
    direction_square: int
    component_pairing_bounds: tuple[Golden, ...]
    component_wall_bounds: tuple[Golden, ...]
    limiting_coordinate: Golden
    limiting_hyperideal_excess: Golden
    shrink_denominator: int
    witness_coordinate: Golden
    minimum_wall_margin: Golden
    witness_hyperideal_excess: Golden


def _is_prime(value: int) -> bool:
    return value >= 2 and all(
        value % divisor for divisor in range(2, isqrt(value) + 1)
    )


def _canonical_tp(value: Golden) -> Golden:
    """Normalize a totally positive generator modulo ``phi^(2 Z)``."""

    if value.sign(0) <= 0 or value.sign(1) <= 0:
        raise ValueError("expected a totally positive generator")
    norm = Golden(value.norm())
    lower = norm / UNIT_SQUARE
    upper = norm * UNIT_SQUARE
    while (value * value - lower).sign(0) <= 0:
        value *= UNIT_SQUARE
    while (value * value - upper).sign(0) > 0:
        value /= UNIT_SQUARE
    if not value.is_integral:
        raise AssertionError("unit normalization escaped Z[phi]")
    return value


def _canonical_admissible(value: Golden) -> Golden:
    """Normalize a sign-(+,-) coefficient modulo unit squares."""

    if value.sign(0) <= 0 or value.sign(1) >= 0:
        raise ValueError("expected sign pattern (+,-)")
    absolute_norm = Golden(-value.norm())
    lower = absolute_norm / UNIT_SQUARE
    upper = absolute_norm * UNIT_SQUARE
    while (value * value - lower).sign(0) < 0:
        value *= UNIT_SQUARE
    while (value * value - upper).sign(0) >= 0:
        value /= UNIT_SQUARE
    if not value.is_integral:
        raise AssertionError("unit normalization escaped Z[phi]")
    return value


def _split_generators(prime: int) -> tuple[Golden, Golden]:
    """Find the two narrow-principal generators over a split prime."""

    found: set[Golden] = set()
    radius = prime + 6
    for rational in range(-radius, radius + 1):
        for phi_coefficient in range(-radius, radius + 1):
            value = Golden(rational, phi_coefficient)
            if value.norm() != prime:
                continue
            if value.sign(0) < 0 and value.sign(1) < 0:
                value = -value
            if value.sign(0) > 0 and value.sign(1) > 0:
                found.add(_canonical_tp(value))
    if len(found) != 2:
        raise AssertionError(f"split prime {prime} yielded {len(found)} generators")
    first, second = tuple(sorted(found))
    return first, second


@lru_cache(maxsize=None)
def clean_golden_atoms(
    norm_bound: int = NORM_BOUND,
) -> tuple[CleanGoldenAtom, ...]:
    """Prime ideals through the absolute ideal-norm cutoff."""

    output: list[CleanGoldenAtom] = []
    for prime in range(2, norm_bound + 1):
        if not _is_prime(prime):
            continue
        if prime == 5:
            output.append(
                CleanGoldenAtom("p5", 5, _canonical_tp(2 + PHI), 5, "ramified")
            )
            continue
        # Splitting is detected independently by the defining polynomial
        # X^2-X-1 modulo p.  An irreducible prime contributes norm p^2.
        split = any(
            (residue * residue - residue - 1) % prime == 0
            for residue in range(prime)
        )
        if split:
            first, second = _split_generators(prime)
            output.extend(
                (
                    CleanGoldenAtom(f"p{prime}+", prime, first, prime, "split"),
                    CleanGoldenAtom(f"p{prime}-", prime, second, prime, "split"),
                )
            )
        elif prime * prime <= norm_bound:
            output.append(
                CleanGoldenAtom(f"p{prime}", prime * prime, Golden(prime), prime, "inert")
            )
    return tuple(sorted(output, key=lambda atom: (atom.ideal_norm, atom.label)))


@lru_cache(maxsize=None)
def clean_golden_directions(
    norm_bound: int = NORM_BOUND,
) -> tuple[CleanGoldenDirection, ...]:
    """All squarefree products of the independently reconstructed atoms."""

    products: list[
        tuple[int, tuple[CleanGoldenAtom, ...], Golden]
    ] = [(1, (), ONE)]
    for atom in clean_golden_atoms(norm_bound):
        snapshot = tuple(products)
        for norm, lineage, generator in snapshot:
            new_norm = norm * atom.ideal_norm
            if new_norm <= norm_bound:
                products.append((new_norm, lineage + (atom,), generator * atom.generator))
    rows = tuple(
        CleanGoldenDirection(
            _canonical_admissible(generator * PHI),
            norm,
            tuple(atom.label for atom in lineage),
            tuple(atom.generator for atom in lineage),
        )
        for norm, lineage, generator in products
    )
    if len({row.atom_labels for row in rows}) != len(rows):
        raise AssertionError("duplicate ideal lineages")
    if len({(row.ideal_norm, row.coefficient) for row in rows}) != len(rows):
        raise AssertionError("two ideal lineages collapsed after normalization")
    return tuple(sorted(rows, key=lambda row: (row.ideal_norm, row.atom_labels)))


def clean_golden_mirror_classes(
    direction: CleanGoldenDirection,
) -> tuple[CleanGoldenMirrorClass, ...]:
    """Local primitive valuation product, with inert dyadic rule explicit."""

    contains_two = "p2" in direction.atom_labels
    dyadic_factors = (Golden(2),) if contains_two else (ONE, Golden(2))
    odd_generators = tuple(
        generator
        for label, generator in zip(direction.atom_labels, direction.atom_generators)
        if label != "p2"
    )
    output: set[CleanGoldenMirrorClass] = set()
    for dyadic in dyadic_factors:
        for mask in range(1 << len(odd_generators)):
            divisor = ONE
            for index, generator in enumerate(odd_generators):
                if mask & (1 << index):
                    divisor *= generator
            output.add(
                CleanGoldenMirrorClass(_canonical_tp(dyadic * divisor), divisor)
            )
    return tuple(
        sorted(output, key=lambda item: (item.canonical_norm, item.odd_spatial_divisor))
    )


def _sqrt_ceiling(value: Golden, embedding: int) -> int:
    if value.sign(embedding) < 0:
        raise ValueError("negative square cap")
    upper = 1
    while (Golden(upper * upper) - value).sign(embedding) < 0:
        upper *= 2
    lower = 0
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if (Golden(middle * middle) - value).sign(embedding) >= 0:
            upper = middle
        else:
            lower = middle
    return upper


@lru_cache(maxsize=None)
def _clean_controls(identity_cap: Golden, conjugate_cap: Golden) -> tuple[Golden, ...]:
    """Enumerate ``m+n phi`` from certified separate embedding bounds."""

    identity = _sqrt_ceiling(identity_cap, 0)
    conjugate = _sqrt_ceiling(conjugate_cap, 1)
    # n=(x-x')/sqrt(5) and
    # m=(phi*x'-(1-phi)*x)/sqrt(5).  The displayed integral boxes are
    # generous because 1/sqrt(5)<1/2, phi/sqrt(5)<1, and
    # |1-phi|/sqrt(5)<1/2.
    n_bound = (identity + conjugate + 1) // 2 + 1
    m_bound = conjugate + (identity + 1) // 2 + 2
    output = []
    for rational in range(-m_bound, m_bound + 1):
        for phi_coefficient in range(-n_bound, n_bound + 1):
            value = Golden(rational, phi_coefficient)
            square = value * value
            if (identity_cap - square).sign(0) < 0:
                continue
            if (conjugate_cap - square).sign(1) < 0:
                continue
            output.append(value)
    return tuple(sorted(output))


def clean_golden_candidate_triples(
    target: Golden,
) -> tuple[CleanGoldenTriple, ...]:
    """Solve the exact binary Gram equations for ``u+phi*v``."""

    if not target.is_integral or target.a < 0:
        return ()
    rational_target = int(target.a)
    phi_target = int(target.b)
    output = []
    for v_square in range(0, rational_target + 1, 2):
        u_square = rational_target - v_square
        dot_numerator = phi_target - v_square
        if u_square % 2 or dot_numerator % 2:
            continue
        dot = dot_numerator // 2
        if u_square and u_square < LEECH_MINIMUM:
            continue
        if v_square and v_square < LEECH_MINIMUM:
            continue
        if dot * dot > u_square * v_square:
            continue
        output.append(CleanGoldenTriple(u_square, v_square, dot))
    return tuple(output)


def clean_golden_has_short_combination(triple: CleanGoldenTriple) -> bool:
    determinant = triple.u_square * triple.v_square - triple.u_dot_v**2
    if determinant < 0:
        return True
    if determinant == 0:
        return gcd(
            gcd(triple.u_square, triple.v_square), abs(triple.u_dot_v)
        ) < LEECH_MINIMUM
    first_bound = isqrt((3 * triple.v_square) // determinant) + 1
    second_bound = isqrt((3 * triple.u_square) // determinant) + 1
    return any(
        0 < triple.square(first, second) < LEECH_MINIMUM
        for first in range(-first_bound, first_bound + 1)
        for second in range(-second_bound, second_bound + 1)
        if first or second
    )


def _clean_levels(
    direction: CleanGoldenDirection,
    height_bound: int = HEIGHT_BOUND,
) -> tuple[tuple[Golden, Golden, Golden, Golden, Golden], ...]:
    coefficient = direction.coefficient
    levels: set[tuple[Golden, Golden, Golden, Golden, Golden]] = set()
    for mirror in clean_golden_mirror_classes(direction):
        norm = mirror.canonical_norm
        controls = _clean_controls(
            norm * height_bound / coefficient,
            norm / (-coefficient),
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
            if (Golden(height_bound) - height).sign(0) < 0:
                continue
            levels.add(
                (height, norm, control, target, mirror.odd_spatial_divisor)
            )
    return tuple(sorted(levels, key=lambda item: (item[0], item[1], item[2], item[4])))


@lru_cache(maxsize=None)
def clean_golden_first_packet(
    direction: CleanGoldenDirection,
) -> tuple[Golden, tuple[CleanGoldenComponent, ...]]:
    levels = _clean_levels(direction)
    index = 0
    while index < len(levels):
        height = levels[index][0]
        same_height = []
        while index < len(levels) and levels[index][0] == height:
            same_height.append(levels[index])
            index += 1
        components = []
        for _height, norm, control, target, divisor in same_height:
            gated_target = target / (divisor * divisor)
            triples = tuple(
                triple
                for triple in clean_golden_candidate_triples(gated_target)
                if not clean_golden_has_short_combination(triple)
            )
            if triples:
                components.append(
                    CleanGoldenComponent(
                        norm, control, target, divisor, gated_target, triples
                    )
                )
        if components:
            return height, tuple(components)
    raise AssertionError(f"no first necessary packet through height {HEIGHT_BOUND}")


@lru_cache(maxsize=None)
def clean_golden_manifest(norm_bound: int = NORM_BOUND) -> tuple[CleanGoldenRow, ...]:
    rows = []
    for direction in clean_golden_directions(norm_bound):
        height, components = clean_golden_first_packet(direction)
        rows.append(
            CleanGoldenRow(
                direction.coefficient,
                direction.ideal_norm,
                direction.atom_labels,
                clean_golden_mirror_classes(direction),
                height,
                components,
            )
        )
    return tuple(rows)


def clean_golden_pairing_bound(
    triple: CleanGoldenTriple, direction_square: int
) -> Golden:
    """Necessary 3x3-Gram/minimum bound on ``d.(u+phi*v)``."""

    if direction_square not in LEECH_DIRECTION_NORMS:
        raise ValueError("unsupported explicit Leech direction norm")
    quotient = triple.u_dot_v // triple.v_square if triple.v_square else 0
    reductions = []
    for coefficient in range(quotient - 3, quotient + 4):
        w_square = triple.square(1, -coefficient)
        if w_square < 0:
            continue
        w_dot_v = triple.u_dot_v - coefficient * triple.v_square
        reductions.append(
            (w_square, abs(w_dot_v), abs(coefficient), coefficient, w_dot_v)
        )
    w_square, _abs_dot, _abs_coefficient, coefficient, w_dot_v = min(reductions)
    first_bound = isqrt(direction_square * w_square)
    second_bound = isqrt(direction_square * triple.v_square)
    best: Golden | None = None
    for first_pairing in range(-first_bound, first_bound + 1):
        for second_pairing in range(-second_bound, second_bound + 1):
            determinant = (
                direction_square
                * (w_square * triple.v_square - w_dot_v * w_dot_v)
                - first_pairing
                * (first_pairing * triple.v_square - second_pairing * w_dot_v)
                + second_pairing
                * (first_pairing * w_dot_v - second_pairing * w_square)
            )
            if determinant < 0:
                continue
            adjacent = (
                direction_square + w_square + 2 * first_pairing,
                direction_square + w_square - 2 * first_pairing,
                direction_square + triple.v_square + 2 * second_pairing,
                direction_square + triple.v_square - 2 * second_pairing,
            )
            if any(value and value < LEECH_MINIMUM for value in adjacent):
                continue
            pairing = Golden(
                first_pairing + coefficient * second_pairing,
                second_pairing,
            )
            if best is None or pairing > best:
                best = pairing
    if best is None:
        raise AssertionError("no pairing survived necessary Gram constraints")
    return best


@lru_cache(maxsize=None)
def clean_golden_prefix_witnesses(
    norm_bound: int = NORM_BOUND,
) -> tuple[CleanGoldenPrefixWitness, ...]:
    output = []
    for row in clean_golden_manifest(norm_bound):
        best: tuple[
            Golden, int, Golden, tuple[Golden, ...], tuple[Golden, ...]
        ] | None = None
        for direction_square in LEECH_DIRECTION_NORMS:
            pairing_bounds = tuple(
                max(
                    clean_golden_pairing_bound(triple, direction_square)
                    for triple in component.surviving_triples
                )
                for component in row.first_packet_components
            )
            wall_bounds = tuple(
                row.coefficient * component.control / component.odd_spatial_divisor
                for component in row.first_packet_components
            )
            limiting = min(
                wall / pairing
                for wall, pairing in zip(wall_bounds, pairing_bounds)
            )
            clearance = direction_square * limiting * limiting - row.coefficient
            candidate = (
                clearance,
                direction_square,
                limiting,
                pairing_bounds,
                wall_bounds,
            )
            if best is None or clearance > best[0]:
                best = candidate
        if best is None or best[0].sign(0) <= 0:
            raise AssertionError(f"no strict limiting witness for {row.coefficient}")
        clearance, direction_square, limiting, pairing_bounds, wall_bounds = best
        denominator = 2
        while True:
            witness = limiting * (denominator - 1) / denominator
            excess = direction_square * witness * witness - row.coefficient
            if excess.sign(0) > 0:
                break
            denominator += 1
            if denominator > 100_000:
                raise AssertionError("failed to shrink to a strict witness")
        margins = tuple(
            wall - witness * pairing
            for wall, pairing in zip(wall_bounds, pairing_bounds)
        )
        if any(margin.sign(0) <= 0 for margin in margins):
            raise AssertionError("non-strict wall margin")
        output.append(
            CleanGoldenPrefixWitness(
                row.coefficient,
                row.ideal_norm,
                row.atom_labels,
                row.first_necessary_height,
                direction_square,
                pairing_bounds,
                wall_bounds,
                limiting,
                clearance,
                denominator,
                witness,
                min(margins),
                excess,
            )
        )
    return tuple(output)
