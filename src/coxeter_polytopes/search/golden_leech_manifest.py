"""Low-norm ideal-squareclass screen for golden Leech hosts in H24.

The hosts are ``<-a> + Lambda_24`` over ``O=Z[phi]``, where
``phi^2=phi+1`` and ``a`` has signs ``(+,-)``.  This module enumerates every
squarefree prime-ideal product of norm at most 128, not merely coefficients
whose two displayed coordinates are coprime.  Both primes over a split
rational prime are retained with their exact generator lineage.

This is a necessary first-packet screen, not a Vinberg termination proof.
For every complete local mirror-norm class it enumerates the first height at
which a Leech Gram triple survives the minimum-four test.  The triples need
not be realised by Leech vectors and need not be accepted chamber walls.

The dyadic rule is specific to the golden field.  The rational prime 2 is
inert.  If ``(2)`` does not divide ``(a)``, primitive reflective vectors have
dyadic norm valuation zero or one.  If ``(2)`` divides squarefree ``(a)``,
evenness forces valuation at least one and primitivity plus reflection
integrality forces it to be exactly one.  In either case there is no spatial
divisibility gate at 2.  Odd prime ideals dividing ``(a)`` independently
give valuation zero or one; valuation one forces the spatial vector into
that prime ideal times the Leech lattice.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import isqrt

from coxeter_polytopes.exact.golden import ONE, PHI, Golden
from coxeter_polytopes.search.silver_leech import (
    LEECH_MINIMUM,
    LEECH_SHELL_COUNTS,
    LeechGramTriple,
)
from coxeter_polytopes.search.silver_leech_manifest import has_short_binary_combination


FUNDAMENTAL_UNIT = PHI
FUNDAMENTAL_UNIT_SQUARE = PHI * PHI
DEFAULT_NORM_BOUND = 128
DEFAULT_HEIGHT_BOUND = 5001
HEIGHT_ZERO_MINIMUM_TRACE = 8

CLASS_AND_UNIT_PROOF = (
    "Z[phi] is norm-Euclidean, hence has class number one. Its units are "
    "+/-phi^n. Since phi has norm -1 and signs (+,-), the narrow and "
    "ordinary class groups agree. The totally positive units are exactly "
    "phi^(2n), and each is the square of the unit phi^n. Thus a squarefree "
    "ideal product has a (+,-) generator unique up to the chosen "
    "phi^(2Z) fundamental domain."
)

MIRROR_CLASS_PROOF = (
    "For a primitive r=(x,y), q(r) divides 2(a*x,(y,Lambda)) because the "
    "Leech lattice is even unimodular. Away from primes over 2a, local "
    "unimodularity and primitivity exclude positive q-valuation. At an odd "
    "prime P dividing squarefree (a), positive q-valuation forces y in "
    "P*Lambda; then x is a P-unit and q=-a*x^2 mod P^2, so the valuation is "
    "exactly one. These odd choices are independent. The prime (2) is "
    "inert. If 2 does not divide a, a primitive x gives dyadic valuation "
    "zero, while nonunit x and primitive y allow only valuation one. If 2 "
    "divides squarefree a, evenness makes q divisible by 2 and q|2I(r) "
    "bounds its valuation by one whenever y is primitive; if y is not "
    "primitive then x is a unit and q=-a*x^2 mod 4 again has valuation one. "
    "Thus the dyadic exponents are {0,1} for 2 not dividing a and {1} for "
    "2 dividing a, with no dyadic spatial gate."
)

HEIGHT_ZERO_PROOF = (
    "A primitive height-zero mirror would have q=y^2 generating (2): "
    "evenness gives (2)|q and q|2(y,Lambda) plus primitive unimodularity "
    "gives the reverse valuation bound. After unit rescaling, y^2=2 and "
    "Tr(y^2)=4. But writing y=u+phi*v gives Tr(y^2)="
    "u^2+(u+v)^2+2v^2. The Leech minimum is four, so every nonzero y has "
    "trace at least eight (if v=0, the first two terms coincide; if v is "
    "nonzero, 2v^2 is already at least eight). Therefore there are no "
    "height-zero mirrors."
)


@dataclass(frozen=True, slots=True)
class GoldenIdealAtom:
    label: str
    ideal_norm: int
    generator: Golden
    rational_prime: int
    splitting_type: str


@dataclass(frozen=True, slots=True)
class GoldenIdealSquareclass:
    coefficient: Golden
    ideal_norm: int
    atom_labels: tuple[str, ...]
    atom_generators: tuple[Golden, ...]


@dataclass(frozen=True, slots=True)
class GoldenMirrorNormClass:
    canonical_norm: Golden
    odd_spatial_divisor: Golden


@dataclass(frozen=True, slots=True)
class GoldenManifestShellComponent:
    norm: Golden
    control: Golden
    target_before_gate: Golden
    odd_spatial_divisor: Golden
    target_after_gate: Golden
    surviving_triples: tuple[LeechGramTriple, ...]


@dataclass(frozen=True, slots=True)
class GoldenLeechManifestRow:
    coefficient: Golden
    ideal_norm: int
    atom_labels: tuple[str, ...]
    mirror_norm_classes: tuple[GoldenMirrorNormClass, ...]
    first_necessary_height: Golden
    first_packet_components: tuple[GoldenManifestShellComponent, ...]
    first_packet_triple_count: int
    all_first_packet_components_rank_one: bool
    first_packet_shell_norms: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class GoldenUniformFirstPacketObstruction:
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


def _is_rational_prime(value: int) -> bool:
    if value < 2:
        return False
    return not any(value % divisor == 0 for divisor in range(2, isqrt(value) + 1))


def _canonical_totally_positive_associate(value: Golden) -> Golden:
    """Canonicalise a totally positive generator modulo ``phi^(2Z)``."""

    if value.sign(0) <= 0 or value.sign(1) <= 0:
        raise ValueError("expected a totally positive element")
    norm = value.norm()
    lower_squared = Golden(norm) / FUNDAMENTAL_UNIT_SQUARE
    upper_squared = Golden(norm) * FUNDAMENTAL_UNIT_SQUARE
    while (value * value - lower_squared).sign(0) <= 0:
        value *= FUNDAMENTAL_UNIT_SQUARE
    while (value * value - upper_squared).sign(0) > 0:
        value /= FUNDAMENTAL_UNIT_SQUARE
    if not value.is_integral:
        raise AssertionError("a totally-positive unit associate ceased to be integral")
    return value


def _canonical_admissible_unit_associate(value: Golden) -> Golden:
    """Canonicalise a sign-(+,-) coefficient modulo unit squares."""

    if value.sign(0) <= 0 or value.sign(1) >= 0:
        raise ValueError("expected sign pattern (+,-)")
    absolute_norm = -value.norm()
    lower_squared = Golden(absolute_norm) / FUNDAMENTAL_UNIT_SQUARE
    upper_squared = Golden(absolute_norm) * FUNDAMENTAL_UNIT_SQUARE
    while (value * value - lower_squared).sign(0) < 0:
        value *= FUNDAMENTAL_UNIT_SQUARE
    while (value * value - upper_squared).sign(0) >= 0:
        value /= FUNDAMENTAL_UNIT_SQUARE
    if not value.is_integral:
        raise AssertionError("an admissible unit associate ceased to be integral")
    return value


@lru_cache(maxsize=None)
def golden_ideal_atoms(norm_bound: int = DEFAULT_NORM_BOUND) -> tuple[GoldenIdealAtom, ...]:
    """Prime ideals of ``Z[phi]`` whose absolute norm is at most the cutoff."""

    if norm_bound < 4:
        return ()
    atoms: list[GoldenIdealAtom] = []
    for prime in range(2, norm_bound + 1):
        if not _is_rational_prime(prime):
            continue
        if prime == 5:
            if prime <= norm_bound:
                atoms.append(
                    GoldenIdealAtom(
                        "p5", 5, _canonical_totally_positive_associate(2 + PHI),
                        5, "ramified",
                    )
                )
            continue
        roots = tuple(
            residue
            for residue in range(prime)
            if (residue * residue - residue - 1) % prime == 0
        )
        if roots:
            candidates: set[Golden] = set()
            search = prime + 4
            for rational in range(-search, search + 1):
                for phi_coefficient in range(-search, search + 1):
                    element = Golden(rational, phi_coefficient)
                    if element.norm() != prime or not element:
                        continue
                    if element.sign(0) < 0 and element.sign(1) < 0:
                        element = -element
                    if element.sign(0) > 0 and element.sign(1) > 0:
                        candidates.add(_canonical_totally_positive_associate(element))
            if len(candidates) != 2:
                raise AssertionError(
                    f"split prime {prime} lost its two ideal generators: {sorted(candidates)!r}"
                )
            first, second = tuple(sorted(candidates))
            atoms.extend(
                (
                    GoldenIdealAtom(f"p{prime}+", prime, first, prime, "split"),
                    GoldenIdealAtom(f"p{prime}-", prime, second, prime, "split"),
                )
            )
        elif prime * prime <= norm_bound:
            atoms.append(
                GoldenIdealAtom(f"p{prime}", prime * prime, Golden(prime), prime, "inert")
            )
    return tuple(sorted(atoms, key=lambda atom: (atom.ideal_norm, atom.label)))


@lru_cache(maxsize=None)
def golden_ideal_squareclasses(
    norm_bound: int = DEFAULT_NORM_BOUND,
) -> tuple[GoldenIdealSquareclass, ...]:
    """Every squarefree prime-ideal product through ``norm_bound``."""

    if norm_bound <= 0:
        raise ValueError("norm bound must be positive")
    products: list[tuple[int, tuple[GoldenIdealAtom, ...], Golden]] = [(1, (), ONE)]
    for atom in golden_ideal_atoms(norm_bound):
        products.extend(
            (norm * atom.ideal_norm, lineage + (atom,), generator * atom.generator)
            for norm, lineage, generator in tuple(products)
            if norm * atom.ideal_norm <= norm_bound
        )
    output = tuple(
        GoldenIdealSquareclass(
            coefficient=_canonical_admissible_unit_associate(generator * PHI),
            ideal_norm=norm,
            atom_labels=tuple(atom.label for atom in lineage),
            atom_generators=tuple(atom.generator for atom in lineage),
        )
        for norm, lineage, generator in products
    )
    if len({row.atom_labels for row in output}) != len(output):
        raise AssertionError("the ideal-lineage census acquired duplicates")
    if len({(row.ideal_norm, row.coefficient) for row in output}) != len(output):
        raise AssertionError("distinct ideal products collapsed as principal generators")
    return tuple(sorted(output, key=lambda row: (row.ideal_norm, row.atom_labels)))


@lru_cache(maxsize=None)
def mirror_norm_classes(coefficient: Golden) -> tuple[GoldenMirrorNormClass, ...]:
    """Complete local mirror-norm classes for a squarefree coefficient."""

    absolute_norm = -coefficient.norm()
    if not coefficient.is_integral or absolute_norm.denominator != 1 or absolute_norm <= 0:
        raise ValueError("coefficient must be admissible and integral")
    odd_factors: list[Golden] = []
    dyadic_divides_coefficient = (coefficient / 2).is_integral
    # Reuse the complete default atom table for the bounded manifest.  Calling
    # the generator separately at every smaller norm gives the same factors
    # but needlessly repeats the split-prime generator search.
    atom_bound = max(DEFAULT_NORM_BOUND, int(absolute_norm))
    for atom in golden_ideal_atoms(atom_bound):
        if atom.rational_prime == 2:
            continue
        quotient = coefficient
        exponent = 0
        while (quotient / atom.generator).is_integral:
            quotient /= atom.generator
            exponent += 1
        if exponent:
            odd_factors.append(atom.generator**exponent)
    dyadic_factors = (Golden(2),) if dyadic_divides_coefficient else (ONE, Golden(2))
    output: set[GoldenMirrorNormClass] = set()
    for dyadic in dyadic_factors:
        for mask in range(1 << len(odd_factors)):
            divisor = ONE
            for index, factor in enumerate(odd_factors):
                if mask & (1 << index):
                    divisor *= factor
            output.add(
                GoldenMirrorNormClass(
                    canonical_norm=_canonical_totally_positive_associate(dyadic * divisor),
                    odd_spatial_divisor=divisor,
                )
            )
    return tuple(
        sorted(output, key=lambda row: (row.canonical_norm, row.odd_spatial_divisor))
    )


def _embedding_integer_upper(value: Golden, embedding: int) -> int:
    if value.sign(embedding) < 0:
        raise ValueError("negative embedding cap")
    upper = 0 if not value else 1
    while (Golden(upper * upper) - value).sign(embedding) < 0:
        upper *= 2
    return upper


@lru_cache(maxsize=None)
def _bounded_golden_for_caps(caps: tuple[Golden, Golden]) -> tuple[Golden, ...]:
    """Golden integers bounded by separate squared embedding caps."""

    identity_bound = _embedding_integer_upper(caps[0], 0)
    conjugate_bound = _embedding_integer_upper(caps[1], 1)
    # If s0=m+n*phi and s1=m+n*(1-phi), then
    # n=(s0-s1)/sqrt(5) and m=(phi*s1-(1-phi)*s0)/sqrt(5).
    # sqrt(5)>2, phi<2, and |1-phi|<1 give these certified caps.
    n_bound = (identity_bound + conjugate_bound + 1) // 2 + 1
    m_bound = conjugate_bound + (identity_bound + 1) // 2 + 2
    return tuple(
        sorted(
            value
            for m in range(-m_bound, m_bound + 1)
            for n in range(-n_bound, n_bound + 1)
            if all(
                (caps[place] - (value := Golden(m, n)) * value).sign(place) >= 0
                for place in (0, 1)
            )
        )
    )


def candidate_golden_leech_gram_triples(target: Golden) -> tuple[LeechGramTriple, ...]:
    """Necessary Leech Gram triples for ``y=u+phi*v`` and ``y^2=target``."""

    if not target.is_integral:
        return ()
    rational_target = int(target.a)
    phi_target = int(target.b)
    if rational_target < 0:
        return ()
    output = []
    for v_square in range(0, rational_target + 1, 2):
        u_square = rational_target - v_square
        dot_numerator = phi_target - v_square
        if dot_numerator % 2:
            continue
        dot = dot_numerator // 2
        if u_square % 2:
            continue
        if u_square not in (0,) and u_square < LEECH_MINIMUM:
            continue
        if v_square not in (0,) and v_square < LEECH_MINIMUM:
            continue
        if dot * dot > u_square * v_square:
            continue
        output.append(LeechGramTriple(u_square, v_square, dot))
    return tuple(output)


def _control_levels(
    coefficient: Golden,
    classes: tuple[GoldenMirrorNormClass, ...],
    height_bound: int,
) -> tuple[tuple[Golden, Golden, Golden, Golden, Golden], ...]:
    levels: set[tuple[Golden, Golden, Golden, Golden, Golden]] = set()
    for norm_class in classes:
        norm = norm_class.canonical_norm
        for control in _bounded_golden_for_caps(
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
            if (Golden(height_bound) - height).sign(0) < 0:
                continue
            levels.add(
                (height, norm, control, target, norm_class.odd_spatial_divisor)
            )
    return tuple(sorted(levels, key=lambda row: (row[0], row[1], row[2], row[4])))


@lru_cache(maxsize=None)
def first_necessary_packet(
    coefficient: Golden,
    height_bound: int = DEFAULT_HEIGHT_BOUND,
) -> tuple[Golden, tuple[GoldenManifestShellComponent, ...]]:
    levels = _control_levels(coefficient, mirror_norm_classes(coefficient), height_bound)
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
                for triple in candidate_golden_leech_gram_triples(gated_target)
                if not has_short_binary_combination(triple)
            )
            if triples:
                components.append(
                    GoldenManifestShellComponent(
                        norm=norm,
                        control=control,
                        target_before_gate=target,
                        odd_spatial_divisor=divisor,
                        target_after_gate=gated_target,
                        surviving_triples=triples,
                    )
                )
        if components:
            return height, tuple(components)
    raise AssertionError(f"no necessary Leech packet below height {height_bound}")


def _rank_one_shell_norm(triple: LeechGramTriple) -> int | None:
    if triple.u_square * triple.v_square != triple.u_dot_v * triple.u_dot_v:
        return None
    from math import gcd

    return gcd(gcd(triple.u_square, triple.v_square), abs(triple.u_dot_v))


def fixed_direction_pairing_bound(
    triple: LeechGramTriple, direction_square: int
) -> Golden:
    """A rigorous 3-by-3 Gram/minimum bound for ``d.(u+phi*v)``."""

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
    best: Golden | None = None
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
            pairing = Golden(first + coefficient * second, second)
            if best is None or pairing > best:
                best = pairing
    if best is None:
        raise AssertionError("no fixed-direction pairing data survived")
    return best


@lru_cache(maxsize=None)
def golden_leech_ideal_squareclass_manifest(
    norm_bound: int = DEFAULT_NORM_BOUND,
) -> tuple[GoldenLeechManifestRow, ...]:
    output = []
    for item in golden_ideal_squareclasses(norm_bound):
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
        output.append(
            GoldenLeechManifestRow(
                coefficient=item.coefficient,
                ideal_norm=item.ideal_norm,
                atom_labels=item.atom_labels,
                mirror_norm_classes=classes,
                first_necessary_height=height,
                first_packet_components=components,
                first_packet_triple_count=sum(
                    len(component.surviving_triples) for component in components
                ),
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
) -> tuple[GoldenUniformFirstPacketObstruction, ...]:
    """Seek a strict fixed-direction hyperideal ray for every first packet."""

    output = []
    for row in golden_leech_ideal_squareclass_manifest(norm_bound):
        best: tuple[
            Golden, int, Golden, tuple[Golden, ...], tuple[Golden, ...]
        ] | None = None
        for direction_square in (4, 6, 8, 10):
            if not LEECH_SHELL_COUNTS.get(direction_square):
                raise AssertionError("the requested Leech direction shell is empty")
            pairing_bounds = tuple(
                max(
                    fixed_direction_pairing_bound(triple, direction_square)
                    for triple in component.surviving_triples
                )
                for component in row.first_packet_components
            )
            wall_bounds = tuple(
                row.coefficient * component.control / component.odd_spatial_divisor
                for component in row.first_packet_components
            )
            limiting = min(
                wall_bound / pairing_bound
                for wall_bound, pairing_bound in zip(wall_bounds, pairing_bounds)
            )
            clearance = direction_square * limiting * limiting - row.coefficient
            candidate = (clearance, direction_square, limiting, pairing_bounds, wall_bounds)
            if best is None or clearance > best[0]:
                best = candidate
        if best is None or best[0].sign(0) <= 0:
            raise AssertionError(
                f"first packet has no fixed-direction obstruction at {row.coefficient}"
            )
        clearance, direction_square, limiting, pairing_bounds, wall_bounds = best
        shrink_denominator = 2
        while True:
            witness = limiting * (shrink_denominator - 1) / shrink_denominator
            witness_excess = direction_square * witness * witness - row.coefficient
            if witness_excess.sign(0) > 0:
                break
            shrink_denominator += 1
            if shrink_denominator > 100_000:
                raise AssertionError("failed to make the limiting witness strict")
        margins = tuple(
            wall_bound - witness * pairing_bound
            for wall_bound, pairing_bound in zip(wall_bounds, pairing_bounds)
        )
        if any(margin.sign(0) <= 0 for margin in margins):
            raise AssertionError("the shrunken witness is not strict")
        output.append(
            GoldenUniformFirstPacketObstruction(
                coefficient=row.coefficient,
                ideal_norm=row.ideal_norm,
                atom_labels=row.atom_labels,
                first_height=row.first_necessary_height,
                direction_square=direction_square,
                component_pairing_bounds=pairing_bounds,
                component_wall_bounds=wall_bounds,
                limiting_coordinate=limiting,
                limiting_hyperideal_excess=clearance,
                shrink_denominator=shrink_denominator,
                witness_coordinate=witness,
                minimum_wall_margin=min(margins),
                witness_hyperideal_excess=witness_excess,
            )
        )
    return tuple(output)
