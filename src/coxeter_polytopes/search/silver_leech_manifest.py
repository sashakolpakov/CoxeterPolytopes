"""Finite low-norm coefficient-direction manifest for silver Leech H24.

This module is intentionally a *necessary-shell ranker*, not a Vinberg
termination engine.  Its primary census enumerates every squarefree principal
ideal of norm at most 128, including ramified, inert, and both split-prime
choices.  A sign-changing unit and a totally-positive unit square give a
canonical admissible coefficient in the Pell domain.  A retained helper also
enumerates the narrower ``gcd(A,B)=1`` slice for regression comparisons.

Odd prime-power factors of ``(a)`` give spatial divisibility gates.  The
local mirror-norm ideal support is the Cartesian product of those factors
with the dyadic choices; for each odd divisor ``d``, the shell census uses
``y=d*z`` and therefore divides the square target by ``d^2``.

Multiplication of ``a`` by a nonunit field square is used only to identify a
diagonal quadratic-space direction for search prioritization.  It is not an
assertion that the corresponding integral ``O_K``-lattices are isometric.
Likewise, a hyperideal patch surviving the first necessary packet is only a
finite-prefix obstruction: later Vinberg mirrors can cut that patch.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import gcd, isqrt

from coxeter_polytopes.exact.silver import SILVER_ONE, SQRT2, Silver
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps
from coxeter_polytopes.search.silver_leech import (
    EPSILON,
    LEECH_MINIMUM,
    LEECH_SHELL_COUNTS,
    LeechGramTriple,
    candidate_leech_gram_triples,
)


EPSILON_SQUARE = EPSILON * EPSILON
DEFAULT_NORM_BOUND = 128
DEFAULT_HEIGHT_BOUND = 701


@dataclass(frozen=True, slots=True)
class PrimitiveSilverCoefficient:
    coefficient: Silver
    absolute_norm: int


@dataclass(frozen=True, slots=True)
class SilverIdealAtom:
    label: str
    ideal_norm: int
    generator: Silver
    rational_prime: int
    splitting_type: str


@dataclass(frozen=True, slots=True)
class SilverIdealSquareclass:
    coefficient: Silver
    ideal_norm: int
    atom_labels: tuple[str, ...]
    atom_generators: tuple[Silver, ...]


@dataclass(frozen=True, slots=True)
class MirrorNormClass:
    canonical_norm: Silver
    odd_spatial_divisor: Silver


@dataclass(frozen=True, slots=True)
class ManifestShellComponent:
    norm: Silver
    control: Silver
    target_before_gate: Silver
    odd_spatial_divisor: Silver
    target_after_gate: Silver
    surviving_triples: tuple[LeechGramTriple, ...]


@dataclass(frozen=True, slots=True)
class SilverLeechManifestRow:
    coefficient: Silver
    absolute_norm: int
    mirror_norm_classes: tuple[MirrorNormClass, ...]
    first_necessary_height: Silver
    first_packet_components: tuple[ManifestShellComponent, ...]
    first_packet_triple_count: int
    all_first_packet_components_rank_one: bool
    first_packet_shell_norms: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class FixedDirectionObstruction:
    coefficient: Silver
    first_height: Silver
    direction_square: int
    maximum_pairing: Silver
    wall_bound: Silver
    witness_coordinate: Silver
    wall_margin: Silver
    hyperideal_excess: Silver


@dataclass(frozen=True, slots=True)
class UniformFirstPacketObstruction:
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


def _prime_power_factorization(value: int) -> tuple[tuple[int, int], ...]:
    output: list[tuple[int, int]] = []
    candidate = 2
    while candidate * candidate <= value:
        if value % candidate:
            candidate += 1
            continue
        exponent = 0
        while value % candidate == 0:
            value //= candidate
            exponent += 1
        output.append((candidate, exponent))
        candidate += 1
    if value > 1:
        output.append((value, 1))
    return tuple(output)


def _is_rational_prime(value: int) -> bool:
    if value < 2:
        return False
    return not any(value % divisor == 0 for divisor in range(2, isqrt(value) + 1))


def _canonical_admissible_unit_associate(value: Silver) -> Silver:
    """Give a sign-(+,-) element its exact Pell-domain representative."""

    if value.sign(0) <= 0 or value.sign(1) >= 0:
        raise ValueError("expected sign pattern (+,-)")
    norm = -value.norm()
    lower_squared = Silver(norm) / EPSILON_SQUARE
    upper_squared = Silver(norm) * EPSILON_SQUARE
    while (value * value - lower_squared).sign(0) < 0:
        value *= EPSILON_SQUARE
    while (value * value - upper_squared).sign(0) >= 0:
        value /= EPSILON_SQUARE
    if not value.is_integral:
        raise AssertionError("an admissible unit associate ceased to be integral")
    return value


@lru_cache(maxsize=None)
def silver_ideal_atoms(norm_bound: int = DEFAULT_NORM_BOUND) -> tuple[SilverIdealAtom, ...]:
    """Prime ideals of ``Z[sqrt(2)]`` whose norm is at most the cutoff."""

    if norm_bound < 2:
        return ()
    # (sqrt(2))=(2+sqrt(2)) and the latter generator is totally positive.
    atoms = [SilverIdealAtom("p2", 2, 2 + SQRT2, 2, "ramified")]
    for prime in range(3, norm_bound + 1):
        if not _is_rational_prime(prime):
            continue
        roots = tuple(value for value in range(prime) if value * value % prime == 2)
        if roots:
            # Each split prime has two conjugate principal ideals.  Find
            # generators of *positive* norm p, canonicalise by TP units, and
            # retain the conjugate choice separately.  Positive norm ensures
            # every atom is totally positive and one final epsilon supplies
            # the desired (+,-) coefficient sign.
            candidates: set[Silver] = set()
            search = prime + 4
            for rational in range(-search, search + 1):
                for radical in range(-search, search + 1):
                    element = Silver(rational, radical)
                    if element.norm() != prime or not element:
                        continue
                    # Convert to a positive-at-both-places generator of this
                    # ideal; this is only lineage, not the final coefficient.
                    for unit in (Silver(1), EPSILON, -Silver(1), -EPSILON):
                        adjusted = element * unit
                        if adjusted.sign(0) > 0 and adjusted.sign(1) > 0:
                            candidates.add(_canonical_totally_positive_associate(adjusted))
                            break
            if len(candidates) != 2:
                raise AssertionError(f"split prime {prime} lost its two ideal generators")
            ordered = tuple(sorted(candidates))
            atoms.extend(
                (
                    SilverIdealAtom(f"p{prime}+", prime, ordered[0], prime, "split"),
                    SilverIdealAtom(f"p{prime}-", prime, ordered[1], prime, "split"),
                )
            )
        elif prime * prime <= norm_bound:
            atoms.append(
                SilverIdealAtom(f"p{prime}", prime * prime, Silver(prime), prime, "inert")
            )
    return tuple(sorted(atoms, key=lambda atom: (atom.ideal_norm, atom.label)))


@lru_cache(maxsize=None)
def silver_ideal_squareclasses(
    norm_bound: int = DEFAULT_NORM_BOUND,
) -> tuple[SilverIdealSquareclass, ...]:
    """Enumerate all squarefree principal ideal classes through the cutoff.

    The ring has class number one.  A product of prime-ideal atoms gives a
    principal generator; multiplying by ``epsilon`` exactly when necessary
    changes (++/--) to (+,-), and a TP unit square then chooses the Pell
    representative.  The unit ideal itself contributes ``epsilon``.
    """

    atoms = silver_ideal_atoms(norm_bound)
    products: list[tuple[int, tuple[SilverIdealAtom, ...], Silver]] = [
        (1, (), Silver(1))
    ]
    for atom in atoms:
        products.extend(
            (norm * atom.ideal_norm, lineage + (atom,), generator * atom.generator)
            for norm, lineage, generator in tuple(products)
            if norm * atom.ideal_norm <= norm_bound
        )
    output = []
    for ideal_norm, lineage, generator in products:
        # A totally-positive ideal generator needs the sign-changing unit;
        # a sign-(+,-) generator does not.  The archived atom generators are
        # TP, so this uniform multiplication also handles every product.
        coefficient = _canonical_admissible_unit_associate(generator * EPSILON)
        output.append(
            SilverIdealSquareclass(
                coefficient=coefficient,
                ideal_norm=ideal_norm,
                atom_labels=tuple(atom.label for atom in lineage),
                atom_generators=tuple(atom.generator for atom in lineage),
            )
        )
    # Distinct squarefree ideals must stay distinct, including conjugate
    # primes over a split rational prime.
    if len({row.atom_labels for row in output}) != len(output):
        raise AssertionError("the ideal-lineage census acquired duplicates")
    return tuple(sorted(output, key=lambda row: (row.ideal_norm, row.atom_labels)))


def _canonical_totally_positive_associate(value: Silver) -> Silver:
    """Choose the associate with ``sqrt(N)/epsilon <= value < sqrt(N)*epsilon``."""

    if value.sign(0) <= 0 or value.sign(1) <= 0:
        raise ValueError("expected a totally positive element")
    norm = value.norm()
    lower_squared = Silver(norm) / EPSILON_SQUARE
    upper_squared = Silver(norm) * EPSILON_SQUARE
    while (value * value - lower_squared).sign(0) <= 0:
        value *= EPSILON_SQUARE
    while (value * value - upper_squared).sign(0) > 0:
        value /= EPSILON_SQUARE
    if not value.is_integral:
        raise AssertionError("a unit associate ceased to be integral")
    return value


def _canonical_positive_divisor(coefficient: Silver, rational_norm: int) -> Silver:
    """Find the totally-positive divisor of ``coefficient`` of given norm."""

    search_bound = rational_norm + 16
    found: set[Silver] = set()
    for rational in range(1, search_bound + 1):
        for radical in range(-search_bound, search_bound + 1):
            divisor = Silver(rational, radical)
            if divisor.norm() != rational_norm:
                continue
            if divisor.sign(0) <= 0 or divisor.sign(1) <= 0:
                continue
            if not (coefficient / divisor).is_integral:
                continue
            found.add(_canonical_totally_positive_associate(divisor))
    if len(found) != 1:
        raise AssertionError(
            f"expected one divisor of norm {rational_norm} for {coefficient}; "
            f"found {sorted(found)!r}"
        )
    return next(iter(found))


@lru_cache(maxsize=None)
def primitive_admissible_coefficients(
    norm_bound: int = DEFAULT_NORM_BOUND,
) -> tuple[PrimitiveSilverCoefficient, ...]:
    """Enumerate squareclasses in a Pell fundamental domain exactly."""

    if norm_bound <= 0:
        raise ValueError("norm bound must be positive")
    # If a is in the fundamental domain, a0^2 < N*epsilon^2.  Since
    # B*sqrt(2)=(a0-a1)/2 and |a1|=N/a0, this finite cap is generous.
    coefficient_cap = 2 * norm_bound + 4
    output: list[PrimitiveSilverCoefficient] = []
    for radical in range(1, coefficient_cap + 1):
        for rational in range(-coefficient_cap, coefficient_cap + 1):
            if gcd(abs(rational), radical) != 1:
                continue
            coefficient = Silver(rational, radical)
            norm = -coefficient.norm()
            if norm.denominator != 1 or not 1 <= norm <= norm_bound:
                continue
            if coefficient.sign(0) <= 0 or coefficient.sign(1) >= 0:
                continue
            # Exact squared version of sqrt(N)/epsilon <= a0 < sqrt(N)*epsilon.
            if (coefficient * coefficient - Silver(norm) / EPSILON_SQUARE).sign(0) < 0:
                continue
            if (coefficient * coefficient - Silver(norm) * EPSILON_SQUARE).sign(0) >= 0:
                continue
            output.append(PrimitiveSilverCoefficient(coefficient, int(norm)))
    return tuple(sorted(output, key=lambda row: (row.absolute_norm, row.coefficient)))


@lru_cache(maxsize=None)
def mirror_norm_classes(coefficient: Silver) -> tuple[MirrorNormClass, ...]:
    """Return the local norm classes and their odd spatial gates.

    For each odd prime-power ideal dividing ``(a)``, reflectivity gives
    valuation either zero or its full coefficient valuation; in the latter
    case ``y`` lies in that ideal times the Leech lattice.  At the ramified
    dyadic prime, the primitive-even alternatives have canonical factors
    ``1,2`` when ``N(a)`` is odd and ``sqrt(2),2`` when it is even.
    """

    norm = -coefficient.norm()
    if norm.denominator != 1 or norm <= 0:
        raise ValueError("coefficient must be primitive admissible integral")
    factors = _prime_power_factorization(int(norm))
    dyadic_exponent = dict(factors).get(2, 0)
    # Factor by prime *ideals*, not just rational primes.  This distinction
    # is essential when both conjugate primes over (say) 7 divide a: their
    # reflection valuations and spatial gates vary independently.
    atoms_list: list[Silver] = []
    for ideal_atom in silver_ideal_atoms(int(norm)):
        if ideal_atom.rational_prime == 2:
            continue
        quotient = coefficient
        exponent = 0
        while (quotient / ideal_atom.generator).is_integral:
            quotient /= ideal_atom.generator
            exponent += 1
        if exponent:
            atoms_list.append(ideal_atom.generator**exponent)
    atoms = tuple(atoms_list)
    dyadic_factors = (Silver(1), Silver(2)) if not dyadic_exponent else (Silver(2, 1), Silver(2))
    output: set[MirrorNormClass] = set()
    for dyadic in dyadic_factors:
        for mask in range(1 << len(atoms)):
            divisor = Silver(1)
            for index, atom in enumerate(atoms):
                if mask & (1 << index):
                    divisor *= atom
            output.add(
                MirrorNormClass(
                    _canonical_totally_positive_associate(dyadic * divisor),
                    divisor,
                )
            )
    return tuple(sorted(output, key=lambda row: (row.canonical_norm, row.odd_spatial_divisor)))


def has_short_binary_combination(triple: LeechGramTriple) -> bool:
    """Test *all* binary combinations for forbidden norm 1--3."""

    u_square = triple.u_square
    v_square = triple.v_square
    dot = triple.u_dot_v
    determinant = u_square * v_square - dot * dot
    if determinant < 0:
        return True
    if determinant == 0:
        # Rank one: the primitive direction has square gcd(U,V,M).
        return gcd(gcd(u_square, v_square), abs(dot)) < LEECH_MINIMUM
    # q(i,j)>=D*j^2/U and >=D*i^2/V after completing squares.  Hence any
    # q<4 lies inside these exact finite coefficient bounds.
    u_bound = isqrt((3 * v_square) // determinant) + 1
    v_bound = isqrt((3 * u_square) // determinant) + 1
    return any(
        0 < triple.combination_square(u_coefficient, v_coefficient) < LEECH_MINIMUM
        for u_coefficient in range(-u_bound, u_bound + 1)
        for v_coefficient in range(-v_bound, v_bound + 1)
        if u_coefficient or v_coefficient
    )


def _control_levels(
    coefficient: Silver,
    classes: tuple[MirrorNormClass, ...],
    height_bound: int,
) -> tuple[tuple[Silver, Silver, Silver, Silver, Silver], ...]:
    levels: set[tuple[Silver, Silver, Silver, Silver, Silver]] = set()
    for norm_class in classes:
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
            if height > height_bound:
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


@lru_cache(maxsize=None)
def first_necessary_packet(
    coefficient: Silver,
    height_bound: int = DEFAULT_HEIGHT_BOUND,
) -> tuple[Silver, tuple[ManifestShellComponent, ...]]:
    classes = mirror_norm_classes(coefficient)
    levels = _control_levels(coefficient, classes, height_bound)
    index = 0
    while index < len(levels):
        height = levels[index][0]
        same_height = []
        while index < len(levels) and levels[index][0] == height:
            same_height.append(levels[index])
            index += 1
        components: list[ManifestShellComponent] = []
        for _height, norm, control, target, divisor in same_height:
            gated_target = target / (divisor * divisor)
            triples = tuple(
                triple
                for triple in candidate_leech_gram_triples(gated_target)
                if not has_short_binary_combination(triple)
            )
            if triples:
                components.append(
                    ManifestShellComponent(
                        norm,
                        control,
                        target,
                        divisor,
                        gated_target,
                        triples,
                    )
                )
        if components:
            return height, tuple(components)
    raise AssertionError(f"no necessary Leech packet below height {height_bound}")


def _rank_one_shell_norm(triple: LeechGramTriple) -> int | None:
    if triple.u_square * triple.v_square != triple.u_dot_v * triple.u_dot_v:
        return None
    return gcd(gcd(triple.u_square, triple.v_square), abs(triple.u_dot_v))


def fixed_direction_pairing_bound(
    triple: LeechGramTriple, direction_square: int
) -> Silver:
    """Rigorous 3x3-Gram/minimum bound for ``d.(u+sqrt(2)v)``."""

    if direction_square < LEECH_MINIMUM or direction_square % 2:
        raise ValueError("expected an even Leech direction norm at least four")
    quotient = triple.u_dot_v // triple.v_square if triple.v_square else 0
    transformations = []
    for coefficient in range(quotient - 2, quotient + 3):
        w_square = triple.combination_square(1, -coefficient)
        if w_square < 0:
            continue
        w_dot_v = triple.u_dot_v - coefficient * triple.v_square
        transformations.append(
            (w_square, abs(w_dot_v), abs(coefficient), coefficient, w_dot_v)
        )
    w_square, _abs_dot, _abs_coefficient, coefficient, w_dot_v = min(transformations)
    first_bound = isqrt(direction_square * w_square)
    second_bound = isqrt(direction_square * triple.v_square)
    best: Silver | None = None
    for first in range(-first_bound, first_bound + 1):
        for second in range(-second_bound, second_bound + 1):
            determinant = (
                direction_square
                * (w_square * triple.v_square - w_dot_v * w_dot_v)
                - first * (first * triple.v_square - second * w_dot_v)
                + second * (first * w_dot_v - second * w_square)
            )
            if determinant < 0:
                continue
            norms = (
                direction_square + w_square + 2 * first,
                direction_square + w_square - 2 * first,
                direction_square + triple.v_square + 2 * second,
                direction_square + triple.v_square - 2 * second,
            )
            if any(value != 0 and value < LEECH_MINIMUM for value in norms):
                continue
            pairing = Silver(first + coefficient * second, second)
            if best is None or pairing > best:
                best = pairing
    if best is None:
        raise AssertionError("no fixed-direction pairing data survived")
    return best


def certify_fixed_direction_obstruction(
    coefficient: Silver,
    direction_square: int,
    witness_coordinate: Silver,
) -> FixedDirectionObstruction:
    """Certify a strict hyperideal direction through the first packet.

    Every first-packet component must have the same normalized section bound;
    this holds for the priority examples archived below.  Odd-divisible
    components use control/divisor after cancelling their spatial divisor.
    """

    height, components = first_necessary_packet(coefficient)
    maxima = []
    bounds = []
    for component in components:
        maxima.extend(
            fixed_direction_pairing_bound(triple, direction_square)
            for triple in component.surviving_triples
        )
        if component.odd_spatial_divisor == 1:
            bounds.append(coefficient * component.control)
        else:
            bounds.append(
                coefficient * component.control / component.odd_spatial_divisor
            )
    wall_bound = min(bounds)
    maximum = max(maxima)
    margin = wall_bound - witness_coordinate * maximum
    excess = direction_square * witness_coordinate * witness_coordinate - coefficient
    if margin.sign(0) <= 0 or excess.sign(0) <= 0:
        raise AssertionError("the requested fixed direction is not a strict obstruction")
    return FixedDirectionObstruction(
        coefficient,
        height,
        direction_square,
        maximum,
        wall_bound,
        witness_coordinate,
        margin,
        excess,
    )


@lru_cache(maxsize=1)
def norm89_priority_obstructions() -> tuple[FixedDirectionObstruction, ...]:
    """The two norm-89 squareclasses, currently closest to shell closure."""

    return (
        certify_fixed_direction_obstruction(Silver(-3, 7), 4, Silver(4) / 3),
        certify_fixed_direction_obstruction(Silver(3, 7), 4, Silver(9) / 5),
    )


@lru_cache(maxsize=1)
def scalar_shell_priority_obstructions() -> tuple[FixedDirectionObstruction, ...]:
    """Strict first-shell obstructions for all four new single-shell rows."""

    norm62 = Silver(-6, 7)
    norm63 = Silver(3, 6)
    norm82 = Silver(-4, 7)
    norm98 = Silver(0, 7)
    return (
        # First packet is the full minimal shell.  A fixed norm-six direction
        # pairs by at most three with every minimal vector.
        certify_fixed_direction_obstruction(norm62, 6, Silver(9) * norm62 / 28),
        # After dividing the scalar spatial multiplier, the slab bound is
        # a*(2-sqrt(2)); again a norm-six direction pairs by at most three.
        certify_fixed_direction_obstruction(
            norm63, 6, Silver(9) * norm63 * Silver(2, -1) / 28
        ),
        # First packet is the full norm-six shell.  A fixed minimal direction
        # pairs by at most three with every norm-six vector.
        certify_fixed_direction_obstruction(norm82, 4, Silver(9) * norm82 / 28),
        # The p7-divisible class is first divided by 3-sqrt(2), giving slab
        # bound a/(3-sqrt(2)) against the minimal shell.
        certify_fixed_direction_obstruction(
            norm98, 6, Silver(9) * norm98 / (28 * Silver(3, -1))
        ),
    )


@lru_cache(maxsize=None)
def silver_leech_low_norm_manifest(
    norm_bound: int = DEFAULT_NORM_BOUND,
) -> tuple[SilverLeechManifestRow, ...]:
    output = []
    for item in primitive_admissible_coefficients(norm_bound):
        classes = mirror_norm_classes(item.coefficient)
        height, components = first_necessary_packet(item.coefficient)
        shell_norms = tuple(
            sorted(
                {
                    shell_norm
                    for component in components
                    for triple in component.surviving_triples
                    if (shell_norm := _rank_one_shell_norm(triple)) is not None
                }
            )
        )
        triple_count = sum(len(component.surviving_triples) for component in components)
        output.append(
            SilverLeechManifestRow(
                coefficient=item.coefficient,
                absolute_norm=item.absolute_norm,
                mirror_norm_classes=classes,
                first_necessary_height=height,
                first_packet_components=components,
                first_packet_triple_count=triple_count,
                all_first_packet_components_rank_one=(
                    all(
                        _rank_one_shell_norm(triple) is not None
                        for component in components
                        for triple in component.surviving_triples
                    )
                ),
                first_packet_shell_norms=shell_norms,
            )
        )
    return tuple(output)


@lru_cache(maxsize=None)
def silver_leech_ideal_squareclass_manifest(
    norm_bound: int = DEFAULT_NORM_BOUND,
) -> tuple[SilverLeechManifestRow, ...]:
    """The full squarefree-principal-ideal manifest through ``norm_bound``."""

    output = []
    for item in silver_ideal_squareclasses(norm_bound):
        classes = mirror_norm_classes(item.coefficient)
        height, components = first_necessary_packet(item.coefficient, 2001)
        shell_norms = tuple(
            sorted(
                {
                    shell_norm
                    for component in components
                    for triple in component.surviving_triples
                    if (shell_norm := _rank_one_shell_norm(triple)) is not None
                }
            )
        )
        triple_count = sum(len(component.surviving_triples) for component in components)
        output.append(
            SilverLeechManifestRow(
                coefficient=item.coefficient,
                absolute_norm=item.ideal_norm,
                mirror_norm_classes=classes,
                first_necessary_height=height,
                first_packet_components=components,
                first_packet_triple_count=triple_count,
                all_first_packet_components_rank_one=all(
                    _rank_one_shell_norm(triple) is not None
                    for component in components
                    for triple in component.surviving_triples
                ),
                first_packet_shell_norms=shell_norms,
            )
        )
    return tuple(output)


@lru_cache(maxsize=None)
def full_manifest_first_packet_obstructions(
    norm_bound: int = DEFAULT_NORM_BOUND,
) -> tuple[UniformFirstPacketObstruction, ...]:
    """Give a strict hyperideal first-prefix witness for every manifest row.

    For each available Leech direction norm ``D in {4,6,8,10}``, the exact
    3x3-Gram oracle bounds its pairing with every necessary Gram component.
    The limiting scale is ``min(B_i/M_i)``.  We choose the direction with
    largest exact light-cone clearance, then shrink that scale by the explicit
    rational factor ``(n-1)/n`` to make every wall inequality strict.

    These witnesses prove that the finite prefix through the first necessary
    packet is not yet finite-volume.  They do not decide the volume of the
    final chamber, because later-height mirrors may remove the patch.
    """

    output = []
    for row in silver_leech_ideal_squareclass_manifest(norm_bound):
        best: tuple[
            Silver, int, Silver, tuple[Silver, ...], tuple[Silver, ...]
        ] | None = None
        for direction_square in (4, 6, 8, 10):
            if LEECH_SHELL_COUNTS.get(direction_square, 0) <= 0:
                raise AssertionError("the requested Leech direction shell is empty")
            pairing_bounds = tuple(
                max(
                    fixed_direction_pairing_bound(triple, direction_square)
                    for triple in component.surviving_triples
                )
                for component in row.first_packet_components
            )
            wall_bounds = tuple(
                row.coefficient
                * component.control
                / component.odd_spatial_divisor
                for component in row.first_packet_components
            )
            limiting = min(
                wall_bound / pairing_bound
                for wall_bound, pairing_bound in zip(wall_bounds, pairing_bounds)
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
            raise AssertionError(
                f"the first-packet fixed-direction filter failed at {row.coefficient}"
            )
        limiting_clearance, direction_square, limiting, pairing_bounds, wall_bounds = best
        shrink_denominator = 2
        while True:
            witness = limiting * (shrink_denominator - 1) / shrink_denominator
            witness_excess = direction_square * witness * witness - row.coefficient
            if witness_excess.sign(0) > 0:
                break
            shrink_denominator += 1
            if shrink_denominator > 100_000:
                raise AssertionError("failed to find a strict explicit shrink factor")
        margins = tuple(
            wall_bound - witness * pairing_bound
            for wall_bound, pairing_bound in zip(wall_bounds, pairing_bounds)
        )
        if any(margin.sign(0) <= 0 for margin in margins):
            raise AssertionError("the shrunken first-packet witness is not strict")
        output.append(
            UniformFirstPacketObstruction(
                coefficient=row.coefficient,
                ideal_norm=row.absolute_norm,
                first_height=row.first_necessary_height,
                direction_square=direction_square,
                component_pairing_bounds=pairing_bounds,
                component_wall_bounds=wall_bounds,
                limiting_coordinate=limiting,
                limiting_hyperideal_excess=limiting_clearance,
                shrink_denominator=shrink_denominator,
                witness_coordinate=witness,
                minimum_wall_margin=min(margins),
                witness_hyperideal_excess=witness_excess,
            )
        )
    return tuple(output)
