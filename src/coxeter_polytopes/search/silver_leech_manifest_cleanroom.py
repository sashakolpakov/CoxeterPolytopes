"""Independent exact replay of the low-norm silver Leech manifest.

This module intentionally does not import ``silver_leech_manifest`` or any of
its enumeration helpers.  It reconstructs the prime-ideal directions, local
mirror classes, first necessary Leech packets, and fixed-direction
first-prefix witnesses directly from arithmetic in ``Z[sqrt(2)]``.

The final witnesses concern a finite prefix only.  A later Vinberg mirror can
cut the open hyperideal patch, so nothing here is a verdict on the volume of
the final reflection chamber.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gcd, isqrt

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver


EPSILON = SILVER_ONE + SQRT2
EPSILON_SQUARE = EPSILON * EPSILON
LEECH_MINIMUM = 4
LEECH_DIRECTION_NORMS = (4, 6, 8, 10)
NORM_BOUND = 128
HEIGHT_BOUND = 701

# Standard MOG coordinates are divided by sqrt(8).  P and R are type-(4^2)
# minimal vectors.  Q is the all-positive vector on the standard Golay octad
# {1,2,6,7,8,10,12,23}.  Hence P, P+Q, P+R, and P-Q explicitly witness the
# four direction norms used by the uniform obstruction oracle.
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
        raise ValueError("the MOG square is not integral")
    return numerator // 8


@dataclass(frozen=True, slots=True)
class CleanAtom:
    label: str
    ideal_norm: int
    generator: Silver
    rational_prime: int
    splitting_type: str


@dataclass(frozen=True, slots=True)
class CleanDirection:
    coefficient: Silver
    ideal_norm: int
    atom_labels: tuple[str, ...]
    atom_generators: tuple[Silver, ...]


@dataclass(frozen=True, slots=True)
class CleanMirrorClass:
    canonical_norm: Silver
    odd_spatial_divisor: Silver


@dataclass(frozen=True, slots=True)
class CleanTriple:
    u_square: int
    v_square: int
    u_dot_v: int

    def square(self, u_coefficient: int, v_coefficient: int) -> int:
        return (
            self.u_square * u_coefficient * u_coefficient
            + 2 * self.u_dot_v * u_coefficient * v_coefficient
            + self.v_square * v_coefficient * v_coefficient
        )


@dataclass(frozen=True, slots=True)
class CleanComponent:
    norm: Silver
    control: Silver
    target_before_gate: Silver
    odd_spatial_divisor: Silver
    target_after_gate: Silver
    surviving_triples: tuple[CleanTriple, ...]


@dataclass(frozen=True, slots=True)
class CleanRow:
    coefficient: Silver
    ideal_norm: int
    atom_labels: tuple[str, ...]
    mirror_norm_classes: tuple[CleanMirrorClass, ...]
    first_necessary_height: Silver
    first_packet_components: tuple[CleanComponent, ...]


@dataclass(frozen=True, slots=True)
class CleanPrefixWitness:
    coefficient: Silver
    ideal_norm: int
    first_height: Silver
    direction_square: int
    component_pairing_bounds: tuple[Silver, ...]
    component_wall_bounds: tuple[Silver, ...]
    limiting_coordinate: Silver
    limiting_hyperideal_excess: Silver
    shrink_denominator: int
    witness_coordinate: Silver
    minimum_wall_margin: Silver
    witness_hyperideal_excess: Silver


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def _canonical_tp(value: Silver) -> Silver:
    """Canonical totally-positive associate in a half-open Pell domain."""

    if value.sign(0) <= 0 or value.sign(1) <= 0:
        raise ValueError("expected a totally-positive element")
    norm = value.norm()
    lower = Silver(norm) / EPSILON_SQUARE
    upper = Silver(norm) * EPSILON_SQUARE
    while (value * value - lower).sign(0) <= 0:
        value *= EPSILON_SQUARE
    while (value * value - upper).sign(0) > 0:
        value /= EPSILON_SQUARE
    if not value.is_integral:
        raise AssertionError("unit normalization left Z[sqrt(2)]")
    return value


def _canonical_admissible(value: Silver) -> Silver:
    """Canonical sign-(+,-) associate modulo totally-positive unit squares."""

    if value.sign(0) <= 0 or value.sign(1) >= 0:
        raise ValueError("expected sign pattern (+,-)")
    absolute_norm = -value.norm()
    lower = Silver(absolute_norm) / EPSILON_SQUARE
    upper = Silver(absolute_norm) * EPSILON_SQUARE
    while (value * value - lower).sign(0) < 0:
        value *= EPSILON_SQUARE
    while (value * value - upper).sign(0) >= 0:
        value /= EPSILON_SQUARE
    if not value.is_integral:
        raise AssertionError("unit normalization left Z[sqrt(2)]")
    return value


def _split_generators(prime: int) -> tuple[Silver, Silver]:
    """Find the two TP generators above a split odd rational prime."""

    found: set[Silver] = set()
    for rational in range(1, prime + 5):
        for radical in range(-prime - 4, prime + 5):
            value = Silver(rational, radical)
            if value.norm() != prime:
                continue
            if value.sign(0) > 0 and value.sign(1) > 0:
                found.add(_canonical_tp(value))
    if len(found) != 2:
        raise AssertionError(f"prime {prime} did not yield two split generators")
    return tuple(sorted(found))  # type: ignore[return-value]


@lru_cache(maxsize=None)
def clean_atoms(norm_bound: int = NORM_BOUND) -> tuple[CleanAtom, ...]:
    atoms = [CleanAtom("p2", 2, Silver(2, 1), 2, "ramified")]
    for prime in range(3, norm_bound + 1):
        if not _is_prime(prime):
            continue
        split = any(value * value % prime == 2 for value in range(prime))
        if split:
            first, second = _split_generators(prime)
            atoms.extend(
                (
                    CleanAtom(f"p{prime}+", prime, first, prime, "split"),
                    CleanAtom(f"p{prime}-", prime, second, prime, "split"),
                )
            )
        elif prime * prime <= norm_bound:
            atoms.append(CleanAtom(f"p{prime}", prime * prime, Silver(prime), prime, "inert"))
    return tuple(sorted(atoms, key=lambda atom: (atom.ideal_norm, atom.label)))


@lru_cache(maxsize=None)
def clean_directions(norm_bound: int = NORM_BOUND) -> tuple[CleanDirection, ...]:
    products: list[tuple[int, tuple[CleanAtom, ...], Silver]] = [(1, (), Silver(1))]
    for atom in clean_atoms(norm_bound):
        products.extend(
            (norm * atom.ideal_norm, lineage + (atom,), generator * atom.generator)
            for norm, lineage, generator in tuple(products)
            if norm * atom.ideal_norm <= norm_bound
        )
    rows = tuple(
        CleanDirection(
            coefficient=_canonical_admissible(generator * EPSILON),
            ideal_norm=norm,
            atom_labels=tuple(atom.label for atom in lineage),
            atom_generators=tuple(atom.generator for atom in lineage),
        )
        for norm, lineage, generator in products
    )
    if len({row.atom_labels for row in rows}) != len(rows):
        raise AssertionError("duplicate squarefree ideal lineages")
    return tuple(sorted(rows, key=lambda row: (row.ideal_norm, row.atom_labels)))


def clean_mirror_classes(direction: CleanDirection) -> tuple[CleanMirrorClass, ...]:
    """Independent local valuation product, keeping the two p=7 ideals apart."""

    has_dyadic_atom = "p2" in direction.atom_labels
    dyadic = (Silver(2, 1), Silver(2)) if has_dyadic_atom else (Silver(1), Silver(2))
    odd_generators = tuple(
        generator
        for label, generator in zip(direction.atom_labels, direction.atom_generators)
        if label != "p2"
    )
    output: set[CleanMirrorClass] = set()
    for dyadic_factor in dyadic:
        for mask in range(1 << len(odd_generators)):
            divisor = Silver(1)
            for index, generator in enumerate(odd_generators):
                if mask & (1 << index):
                    divisor *= generator
            output.add(CleanMirrorClass(_canonical_tp(dyadic_factor * divisor), divisor))
    return tuple(sorted(output, key=lambda item: (item.canonical_norm, item.odd_spatial_divisor)))


def _embedding_sqrt_ceiling(value: Silver, embedding: int) -> int:
    if value.sign(embedding) < 0:
        raise ValueError("negative square cap")
    upper = 1
    while (Silver(upper * upper) - value).sign(embedding) < 0:
        upper *= 2
    lower = 0
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if (Silver(middle * middle) - value).sign(embedding) >= 0:
            upper = middle
        else:
            lower = middle
    return upper


@lru_cache(maxsize=None)
def _bounded_controls(identity_cap: Silver, conjugate_cap: Silver) -> tuple[Silver, ...]:
    first = _embedding_sqrt_ceiling(identity_cap, 0)
    second = _embedding_sqrt_ceiling(conjugate_cap, 1)
    coefficient_bound = (first + second + 1) // 2
    output = []
    for rational in range(-coefficient_bound, coefficient_bound + 1):
        for radical in range(-coefficient_bound, coefficient_bound + 1):
            value = Silver(rational, radical)
            square = value * value
            if (identity_cap - square).sign(0) < 0:
                continue
            if (conjugate_cap - square).sign(1) < 0:
                continue
            output.append(value)
    return tuple(sorted(output))


def clean_candidate_triples(target: Silver) -> tuple[CleanTriple, ...]:
    if not target.is_integral:
        return ()
    dot = target.b / 2
    if dot.denominator != 1:
        return ()
    rational_target = int(target.a)
    dot_int = int(dot)
    output = []
    for v_square in range(0, rational_target // 2 + 1, 2):
        u_square = rational_target - 2 * v_square
        if u_square < 0 or u_square % 2:
            continue
        if u_square and u_square < LEECH_MINIMUM:
            continue
        if v_square and v_square < LEECH_MINIMUM:
            continue
        if dot_int * dot_int > u_square * v_square:
            continue
        output.append(CleanTriple(u_square, v_square, dot_int))
    return tuple(output)


def clean_has_short_combination(triple: CleanTriple) -> bool:
    """Exhaust every binary combination that could have square 1, 2, or 3."""

    determinant = triple.u_square * triple.v_square - triple.u_dot_v**2
    if determinant < 0:
        return True
    if determinant == 0:
        return gcd(gcd(triple.u_square, triple.v_square), abs(triple.u_dot_v)) < 4
    # q(i,j)>=det*i^2/V and >=det*j^2/U.  Since q is integral, q<4
    # implies both displayed bounds with numerator 3.
    i_bound = isqrt((3 * triple.v_square) // determinant) + 1
    j_bound = isqrt((3 * triple.u_square) // determinant) + 1
    return any(
        0 < triple.square(i, j) < 4
        for i in range(-i_bound, i_bound + 1)
        for j in range(-j_bound, j_bound + 1)
        if i or j
    )


def _clean_levels(
    direction: CleanDirection, height_bound: int = HEIGHT_BOUND
) -> tuple[tuple[Silver, Silver, Silver, Silver, Silver], ...]:
    coefficient = direction.coefficient
    levels: set[tuple[Silver, Silver, Silver, Silver, Silver]] = set()
    for mirror in clean_mirror_classes(direction):
        norm = mirror.canonical_norm
        controls = _bounded_controls(
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
            if (Silver(height_bound) - height).sign(0) < 0:
                continue
            levels.add((height, norm, control, target, mirror.odd_spatial_divisor))
    return tuple(sorted(levels, key=lambda item: (item[0], item[1], item[2], item[4])))


@lru_cache(maxsize=None)
def clean_first_packet(direction: CleanDirection) -> tuple[Silver, tuple[CleanComponent, ...]]:
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
                for triple in clean_candidate_triples(gated_target)
                if not clean_has_short_combination(triple)
            )
            if triples:
                components.append(
                    CleanComponent(norm, control, target, divisor, gated_target, triples)
                )
        if components:
            return height, tuple(components)
    raise AssertionError(f"no necessary packet through height {HEIGHT_BOUND}")


@lru_cache(maxsize=None)
def clean_manifest(norm_bound: int = NORM_BOUND) -> tuple[CleanRow, ...]:
    rows = []
    for direction in clean_directions(norm_bound):
        height, components = clean_first_packet(direction)
        rows.append(
            CleanRow(
                direction.coefficient,
                direction.ideal_norm,
                direction.atom_labels,
                clean_mirror_classes(direction),
                height,
                components,
            )
        )
    return tuple(rows)


def clean_pairing_bound(triple: CleanTriple, direction_square: int) -> Silver:
    """A necessary 3x3 Gram/minimum bound on d.(u+sqrt(2)v)."""

    if direction_square not in LEECH_DIRECTION_NORMS:
        raise ValueError("unsupported Leech direction norm")
    quotient = triple.u_dot_v // triple.v_square if triple.v_square else 0
    choices = []
    for coefficient in range(quotient - 3, quotient + 4):
        w_square = triple.square(1, -coefficient)
        if w_square < 0:
            continue
        w_dot_v = triple.u_dot_v - coefficient * triple.v_square
        choices.append((w_square, abs(w_dot_v), abs(coefficient), coefficient, w_dot_v))
    w_square, _absolute_dot, _absolute_coefficient, coefficient, w_dot_v = min(choices)
    first_bound = isqrt(direction_square * w_square)
    second_bound = isqrt(direction_square * triple.v_square)
    best: Silver | None = None
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
            # Every nonzero Leech vector has norm at least four.  These four
            # combinations are enough for a rigorous upper bound; omitting
            # other minimum constraints can only make the bound weaker.
            adjacent_norms = (
                direction_square + w_square + 2 * first_pairing,
                direction_square + w_square - 2 * first_pairing,
                direction_square + triple.v_square + 2 * second_pairing,
                direction_square + triple.v_square - 2 * second_pairing,
            )
            if any(value and value < 4 for value in adjacent_norms):
                continue
            pairing = Silver(first_pairing + coefficient * second_pairing, second_pairing)
            if best is None or pairing > best:
                best = pairing
    if best is None:
        raise AssertionError("no pairing candidate survived the necessary Gram tests")
    return best


@lru_cache(maxsize=None)
def clean_prefix_witnesses(norm_bound: int = NORM_BOUND) -> tuple[CleanPrefixWitness, ...]:
    output = []
    for row in clean_manifest(norm_bound):
        best: tuple[Silver, int, Silver, tuple[Silver, ...], tuple[Silver, ...]] | None = None
        for direction_square in LEECH_DIRECTION_NORMS:
            pairing_bounds = tuple(
                max(clean_pairing_bound(triple, direction_square) for triple in component.surviving_triples)
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
            candidate = (clearance, direction_square, limiting, pairing_bounds, wall_bounds)
            if best is None or clearance > best[0]:
                best = candidate
        if best is None or best[0].sign(0) <= 0:
            raise AssertionError(f"no first-prefix witness for {row.coefficient}")
        clearance, direction_square, limiting, pairing_bounds, wall_bounds = best
        denominator = 2
        while True:
            witness = limiting * (denominator - 1) / denominator
            excess = direction_square * witness * witness - row.coefficient
            if excess.sign(0) > 0:
                break
            denominator += 1
            if denominator > 100_000:
                raise AssertionError("failed to make the limiting witness strict")
        margins = tuple(
            wall - witness * pairing
            for wall, pairing in zip(wall_bounds, pairing_bounds)
        )
        if any(margin.sign(0) <= 0 for margin in margins):
            raise AssertionError("non-strict first-prefix wall margin")
        output.append(
            CleanPrefixWitness(
                row.coefficient,
                row.ideal_norm,
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
