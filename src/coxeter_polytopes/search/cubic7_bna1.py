"""Exact Vinberg prefixes for ``B_{n-1} x A_1`` cubic-seven lattices.

For a hyperbolic dimension ``n`` this module studies

``diag(-theta, 1, ..., 1, t)``

with ``n-1`` ordinary positive coordinates and one distinguished coordinate,
where ``t`` is either ``2`` or ``2-theta``.  The implementation is dimension
parametric; in particular the rank-six cases are the first ``H^5`` stretch of
the non-unimodular rank-five search.

The mirror-norm list is complete, not a norm cutoff.  If ``r`` is primitive,
reflection integrality says ``q(r) | 2 I(r)``, where

``I(r)=(-theta*r0, r1, ..., r_{n-1}, t*rn)``.

Away from primes dividing ``t``, primitivity makes ``I(r)`` locally the unit
ideal.  At the (squarefree) prime dividing ``t``, either this is still true or
all other coordinates lie in that prime and the last coordinate is a unit;
then the diagonal congruence makes the norm valuation exactly one.  Since 2
is inert, the possible norm ideals have canonical generators ``1,2`` for
``t=2`` and ``1,2,2-theta,2(2-theta)`` for ``t=2-theta``.  There is no hidden
unit family: by Dirichlet's unit theorem ``O^*/O^{*2}`` has order eight,
while ``-1``, ``theta``, and ``1+theta`` have independent sign vectors
``---``, ``+--``, and ``++-``.  Thus the signature map is an isomorphism and
every totally positive unit is a square.

Ordinary spatial coordinates are enumerated as ordered multisets.  This is
essential in dimensions 22--29: the height-zero Weyl group is ``B_{n-1}``,
so expanding all coordinate permutations during the positive-height search
would be pure duplication.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import factorial
from itertools import groupby
from time import monotonic
from typing import Sequence

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    THETA7,
    Cubic7,
)
from coxeter_polytopes.search.cubic7_diagonal import bounded_cubic7_integers


Vector = tuple[Cubic7, ...]
Matrix = tuple[tuple[Cubic7, ...], ...]
DELTA7 = 2 - THETA7


@dataclass(frozen=True, slots=True)
class Cubic7BnA1Seed:
    seed_id: str
    last_coefficient: Cubic7
    canonical_mirror_norms: tuple[Cubic7, ...]
    complete_norm_class_proof: str

    def coefficients(self, dimension: int) -> tuple[Cubic7, ...]:
        if dimension < 3:
            raise ValueError("hyperbolic dimension must be at least three")
        return (
            -THETA7,
            *(CUBIC7_ONE for _ in range(dimension - 1)),
            self.last_coefficient,
        )


L2_BNA1_SEED = Cubic7BnA1Seed(
    seed_id="L2-BnA1",
    last_coefficient=Cubic7(2),
    canonical_mirror_norms=(Cubic7(1), Cubic7(2)),
    complete_norm_class_proof=(
        "For primitive r, (q(r)) divides 2I(r). Away from the inert prime "
        "(2), all diagonal coefficients are local units and primitivity makes "
        "I(r)=O, hence no other prime divides q(r). At (2), either I(r)=O, "
        "giving v_2(q)<=1, or every non-final coordinate lies in (2) and the "
        "final coordinate is a unit; then q(r)=2*r_final^2 mod (4), so "
        "v_2(q)=1. Thus (q) is (1) or (2). Finally O*/O*^2 has order eight "
        "by Dirichlet, and the units -1, theta, 1+theta have independent "
        "signatures ---, +--, ++-; the signature map is therefore bijective, "
        "so every totally positive unit is a square. The canonical element "
        "norms are exactly 1 and 2 modulo totally-positive unit squares."
    ),
)

L7_BNA1_SEED = Cubic7BnA1Seed(
    seed_id="L7-BnA1",
    last_coefficient=DELTA7,
    canonical_mirror_norms=(Cubic7(1), Cubic7(2), DELTA7, 2 * DELTA7),
    complete_norm_class_proof=(
        "For primitive r, (q(r)) divides 2I(r). Away from the inert prime "
        "(2) and p=(2-theta), all diagonal coefficients are local units and "
        "primitivity makes I(r)=O, so no other prime divides q(r). At (2), "
        "2-theta is a unit, hence I(r)=O and v_2(q)<=1. At p, either I(r)=O "
        "and v_p(q)=0, or every non-final coordinate lies in p while the "
        "final coordinate is a unit; then q(r)=(2-theta)*r_final^2 mod p^2, "
        "so v_p(q)=1. Thus (q) is (1),(2),p,or (2)p. Finally O*/O*^2 has "
        "order eight by Dirichlet, and the units -1, theta, 1+theta have "
        "independent signatures ---, +--, ++-; therefore every totally "
        "positive unit is a square. The canonical element norms are exactly "
        "1, 2, 2-theta, and 2(2-theta) modulo such squares."
    ),
)


def bna1_diagonal_pair(
    seed: Cubic7BnA1Seed,
    dimension: int,
    left: Sequence[Cubic7],
    right: Sequence[Cubic7],
) -> Cubic7:
    coefficients = seed.coefficients(dimension)
    if len(left) != dimension + 1 or len(right) != dimension + 1:
        raise ValueError("a vector must have hyperbolic dimension plus one entries")
    return sum(
        (coefficient * x * y for coefficient, x, y in zip(coefficients, left, right)),
        start=CUBIC7_ZERO,
    )


def _embedding_integer_upper(value: Cubic7, embedding: int) -> int:
    if value.sign(embedding) < 0:
        raise ValueError("a square bound must be nonnegative at its selected embedding")
    if not value:
        return 0
    upper = 1
    while (Cubic7(upper) - value).sign(embedding) < 0:
        upper *= 2
    return upper


@lru_cache(maxsize=None)
def _bounded_for_caps(caps: tuple[Cubic7, Cubic7, Cubic7]) -> tuple[Cubic7, ...]:
    rational_caps = tuple(
        _embedding_integer_upper(cap, embedding) for embedding, cap in enumerate(caps)
    )
    return bounded_cubic7_integers(rational_caps)


def _reflection_is_integral(
    seed: Cubic7BnA1Seed,
    dimension: int,
    vector: Vector,
    norm: Cubic7,
) -> bool:
    return all(
        (2 * coefficient * coordinate / norm).is_integral
        for coefficient, coordinate in zip(seed.coefficients(dimension), vector)
    )


def _nonnegative_values(caps: tuple[Cubic7, Cubic7, Cubic7]) -> tuple[Cubic7, ...]:
    return tuple(value for value in _bounded_for_caps(caps) if value.sign(0) >= 0)


def _square_partition_solutions(
    values: Sequence[Cubic7], coordinate_count: int, target: Cubic7
) -> tuple[tuple[Cubic7, ...], ...]:
    """Solve a sum of squares as dominant ordered coordinate multisets.

    ``values`` must contain exactly the nonnegative-at-place-zero choices.
    Nonzero values have distinct squares, because ``x^2=y^2`` implies
    ``x=+/-y`` and only one sign is nonnegative at the distinguished place.
    The recursion stores multiplicities and appends zero coordinates once;
    it never expands Weyl-equivalent permutations.
    """

    if coordinate_count < 0:
        return ()
    positive = tuple(
        sorted((value for value in values if value), reverse=True)
    )
    squares = tuple(value * value for value in positive)

    @lru_cache(maxsize=None)
    def solve(index: int, slots: int, remainder: Cubic7) -> tuple[tuple[int, ...], ...]:
        if any(remainder.sign(embedding) < 0 for embedding in range(3)):
            return ()
        if index == len(positive):
            return ((),) if not remainder else ()
        square = squares[index]
        output: list[tuple[int, ...]] = []
        running = remainder
        for multiplicity in range(slots + 1):
            if all(running.sign(embedding) >= 0 for embedding in range(3)):
                for suffix in solve(index + 1, slots - multiplicity, running):
                    output.append((multiplicity, *suffix))
            else:
                break
            running -= square
        return tuple(output)

    solutions = []
    for multiplicities in solve(0, coordinate_count, target):
        used = sum(multiplicities)
        coordinates = tuple(
            value
            for value, multiplicity in zip(positive, multiplicities)
            for _ in range(multiplicity)
        ) + (CUBIC7_ZERO,) * (coordinate_count - used)
        solutions.append(coordinates)
    return tuple(solutions)


@dataclass(frozen=True, slots=True)
class HeightZeroPattern:
    norm: Cubic7
    last_coordinate: Cubic7
    ordinary_coordinates: tuple[Cubic7, ...]
    orbit_size: int


def _signed_permutation_orbit_size(coordinates: Sequence[Cubic7]) -> int:
    nonzero = [value for value in coordinates if value]
    multiplicities: dict[Cubic7, int] = {}
    for value in nonzero:
        square = value * value
        multiplicities[square] = multiplicities.get(square, 0) + 1
    permutation_count = factorial(len(coordinates))
    permutation_count //= factorial(len(coordinates) - len(nonzero))
    for count in multiplicities.values():
        permutation_count //= factorial(count)
    return permutation_count * 2 ** len(nonzero)


@lru_cache(maxsize=None)
def enumerate_height_zero_patterns(
    seed: Cubic7BnA1Seed, dimension: int
) -> tuple[HeightZeroPattern, ...]:
    """Exhaust height-zero roots modulo the ordinary ``B_{n-1}`` action."""

    ordinary_count = dimension - 1
    patterns: list[HeightZeroPattern] = []
    for norm in seed.canonical_mirror_norms:
        ordinary = _nonnegative_values((norm, norm, norm))
        weighted_caps = tuple(norm / seed.last_coefficient for _ in range(3))
        # Keep one representative of the +/- last-coordinate pair.  Its
        # sign contributes a factor two to the orbit size below.
        last_values = _nonnegative_values(weighted_caps)  # type: ignore[arg-type]
        for last in last_values:
            remainder = norm - seed.last_coefficient * last * last
            if any(remainder.sign(embedding) < 0 for embedding in range(3)):
                continue
            for spatial in _square_partition_solutions(ordinary, ordinary_count, remainder):
                vector = (CUBIC7_ZERO, *spatial, last)
                if not _reflection_is_integral(seed, dimension, vector, norm):
                    continue
                orbit = _signed_permutation_orbit_size(spatial) * (2 if last else 1)
                patterns.append(HeightZeroPattern(norm, last, spatial, orbit))
    return tuple(sorted(patterns, key=lambda item: (item.norm, item.last_coordinate, item.ordinary_coordinates)))


def height_zero_is_bn_minus_1_times_a1(seed: Cubic7BnA1Seed, dimension: int) -> bool:
    """Certify that the exact pattern list is precisely ``B_{n-1} x A1``."""

    ordinary_count = dimension - 1
    one = CUBIC7_ONE
    zero = CUBIC7_ZERO
    expected = {
        (Cubic7(1), zero, (one,) + (zero,) * (ordinary_count - 1)),
        (Cubic7(2), zero, (one, one) + (zero,) * (ordinary_count - 2)),
        (seed.last_coefficient, one, (zero,) * ordinary_count),
    }
    actual = {
        (pattern.norm, pattern.last_coordinate, pattern.ordinary_coordinates)
        for pattern in enumerate_height_zero_patterns(seed, dimension)
    }
    return actual == expected and sum(pattern.orbit_size for pattern in enumerate_height_zero_patterns(seed, dimension)) == (
        2 * ordinary_count * ordinary_count + 2
    )


def _cone_roots(dimension: int) -> tuple[Vector, ...]:
    ordinary_count = dimension - 1
    roots: list[Vector] = []
    for index in range(ordinary_count - 1):
        root = [CUBIC7_ZERO for _ in range(dimension + 1)]
        root[1 + index] = -CUBIC7_ONE
        root[2 + index] = CUBIC7_ONE
        roots.append(tuple(root))
    ordinary_end = [CUBIC7_ZERO for _ in range(dimension + 1)]
    ordinary_end[ordinary_count] = -CUBIC7_ONE
    roots.append(tuple(ordinary_end))
    final = [CUBIC7_ZERO for _ in range(dimension + 1)]
    final[-1] = -CUBIC7_ONE
    roots.append(tuple(final))
    return tuple(roots)


def dominant_bna1_representative(
    seed: Cubic7BnA1Seed,
    dimension: int,
    vector: Sequence[Cubic7],
) -> Vector:
    """Return the unique positive-root representative in the Weyl chamber.

    The complete height-zero reflection group is ``B_(d-1) x A1``.  Its
    ``B`` factor independently changes signs of, and permutes, the ordinary
    coordinates; the ``A1`` factor changes the sign of the final coordinate.
    Root orientation first makes a nonzero control coordinate positive.
    Sorting the ordinary absolute values at the distinguished real embedding
    therefore gives one representative of every positive-height mirror
    orbit, not merely a convenient sector of the root equation.

    Equality at the distinguished embedding is equality in the field, so
    repeated coordinates introduce stabilizers but no ambiguous orbit
    representatives.
    """

    entries = tuple(Cubic7.coerce(value) for value in vector)
    if len(entries) != dimension + 1:
        raise ValueError("a vector must have hyperbolic dimension plus one entries")
    # This also validates that the seed supports the requested dimension.
    seed.coefficients(dimension)
    if entries[0].sign(0) < 0:
        entries = tuple(-entry for entry in entries)
    ordinary = tuple(
        sorted(
            (
                entry if entry.sign(0) >= 0 else -entry
                for entry in entries[1:-1]
            ),
            reverse=True,
        )
    )
    final = entries[-1] if entries[-1].sign(0) >= 0 else -entries[-1]
    return (entries[0], *ordinary, final)


@dataclass(frozen=True, slots=True)
class Cubic7BnA1Root:
    norm: Cubic7
    height: Cubic7
    vector: Vector

    @property
    def control(self) -> Cubic7:
        return self.vector[0]


@dataclass(frozen=True, slots=True)
class Cubic7BnA1RootPool:
    seed: Cubic7BnA1Seed
    dimension: int
    height_bound: int
    cone_roots: tuple[Vector, ...]
    roots: tuple[Cubic7BnA1Root, ...]
    control_counts_by_norm: tuple[tuple[Cubic7, int], ...]


def _controls(norm: Cubic7, height_bound: int) -> tuple[Cubic7, ...]:
    caps = (
        norm * height_bound / THETA7,
        norm / (-THETA7),
        norm / (-THETA7),
    )
    return tuple(
        control
        for control in _bounded_for_caps(caps)
        if control.sign(0) > 0
        and (norm * height_bound - THETA7 * control * control).sign(0) >= 0
        and all((norm + THETA7 * control * control).sign(embedding) >= 0 for embedding in (1, 2))
    )


def enumerate_bna1_root_pool(
    seed: Cubic7BnA1Seed, dimension: int, height_bound: int
) -> Cubic7BnA1RootPool:
    """Exhaust normalized positive-height roots through Vinberg height ``H``."""

    if seed not in (L2_BNA1_SEED, L7_BNA1_SEED):
        raise ValueError("the complete local proof currently covers L2 and L7")
    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    if not height_zero_is_bn_minus_1_times_a1(seed, dimension):
        raise AssertionError("the asserted B_(n-1) x A1 height-zero chamber is incomplete")

    ordinary_count = dimension - 1
    roots: set[Cubic7BnA1Root] = set()
    control_counts: list[tuple[Cubic7, int]] = []
    for norm in seed.canonical_mirror_norms:
        controls = _controls(norm, height_bound)
        control_counts.append((norm, len(controls)))
        target_caps = (norm * (height_bound + 1), norm, norm)
        ordinary_pool = _nonnegative_values(target_caps)
        weighted_caps = tuple(cap / seed.last_coefficient for cap in target_caps)
        last_pool = _nonnegative_values(weighted_caps)  # type: ignore[arg-type]
        for control in controls:
            target = norm + THETA7 * control * control
            ordinary = tuple(
                value
                for value in ordinary_pool
                if all((target - value * value).sign(embedding) >= 0 for embedding in range(3))
            )
            for last in last_pool:
                remainder = target - seed.last_coefficient * last * last
                if any(remainder.sign(embedding) < 0 for embedding in range(3)):
                    continue
                for spatial in _square_partition_solutions(ordinary, ordinary_count, remainder):
                    vector = (control, *spatial, last)
                    if bna1_diagonal_pair(seed, dimension, vector, vector) != norm:
                        raise AssertionError("partition enumerator returned a vector of wrong norm")
                    if not _reflection_is_integral(seed, dimension, vector, norm):
                        continue
                    roots.add(
                        Cubic7BnA1Root(
                            norm=norm,
                            height=THETA7 * control * control / norm,
                            vector=vector,
                        )
                    )

    return Cubic7BnA1RootPool(
        seed=seed,
        dimension=dimension,
        height_bound=height_bound,
        cone_roots=_cone_roots(dimension),
        roots=tuple(sorted(roots, key=lambda root: (root.height, root.norm, root.vector))),
        control_counts_by_norm=tuple(control_counts),
    )


def _kernel(rows: Sequence[Sequence[Cubic7]]) -> Vector | None:
    if not rows:
        return None
    width = len(rows[0])
    work = [list(row) for row in rows]
    if any(len(row) != width for row in work):
        raise ValueError("kernel input has inconsistent row widths")
    pivots: list[int] = []
    pivot_row = 0
    for column in range(width):
        pivot = next((row for row in range(pivot_row, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                work[row][index] - multiplier * work[pivot_row][index]
                for index in range(width)
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(work):
            break
    if width - pivot_row != 1:
        return None
    free = next(column for column in range(width) if column not in pivots)
    vector = [CUBIC7_ZERO for _ in range(width)]
    vector[free] = CUBIC7_ONE
    for row, column in enumerate(pivots):
        vector[column] = -work[row][free]
    return tuple(vector)


def _row_rank(rows: Sequence[Sequence[Cubic7]], width: int) -> int:
    """Exact row rank over the conductor-seven field."""

    work = [list(row) for row in rows]
    if any(len(row) != width for row in work):
        raise ValueError("row-rank input has inconsistent widths")
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for row in range(rank + 1, len(work)):
            if not work[row][column]:
                continue
            multiplier = work[row][column] / pivot_value
            for index in range(column, width):
                work[row][index] -= multiplier * work[rank][index]
        rank += 1
        if rank == len(work):
            break
    return rank


def _bna1_support_kernel(
    dimension: int,
    functionals: Sequence[Vector],
    support: Sequence[int],
) -> Vector | None:
    """Exploit the ``B_(n-1) x A1`` cone before doing field elimination.

    The first ``dimension`` walls are the height-zero walls.  Their active
    equalities merely identify adjacent ordinary coordinates, anchor the
    final ordinary block at zero, or anchor the A1 coordinate at zero.  We
    first write an integral basis for their common kernel.  If ``p`` positive
    walls are active, the remaining field calculation is only a ``p`` by
    ``p+1`` kernel, independently of the ambient dimension.
    """

    ambient_rank = dimension + 1
    ordinary_count = dimension - 1
    active = set(support)
    active_positive = tuple(index for index in support if index >= dimension)

    basis: list[Vector] = []
    control = [CUBIC7_ZERO for _ in range(ambient_rank)]
    control[0] = CUBIC7_ONE
    basis.append(tuple(control))

    # Active chain walls join adjacent ordinary coordinates.  A missing chain
    # wall starts a new component.  The final ordinary wall anchors the last
    # component at zero; every earlier component stays free.
    components: list[list[int]] = [[1]]
    for chain in range(ordinary_count - 1):
        coordinate = chain + 2
        if chain in active:
            components[-1].append(coordinate)
        else:
            components.append([coordinate])
    ordinary_end_wall = ordinary_count - 1
    anchored_component = len(components) - 1 if ordinary_end_wall in active else None
    for component_index, component in enumerate(components):
        if component_index == anchored_component:
            continue
        vector = [CUBIC7_ZERO for _ in range(ambient_rank)]
        for coordinate in component:
            vector[coordinate] = CUBIC7_ONE
        basis.append(tuple(vector))

    last_wall = dimension - 1
    if last_wall not in active:
        vector = [CUBIC7_ZERO for _ in range(ambient_rank)]
        vector[-1] = CUBIC7_ONE
        basis.append(tuple(vector))

    reduced_rows = tuple(
        tuple(
            sum(
                (functionals[wall][index] * basis_vector[index] for index in range(ambient_rank)),
                start=CUBIC7_ZERO,
            )
            for basis_vector in basis
        )
        for wall in active_positive
    )
    if not reduced_rows:
        return basis[0] if len(basis) == 1 else None
    reduced = _kernel(reduced_rows)
    if reduced is None:
        return None
    return tuple(
        sum(
            (coefficient * basis_vector[index] for coefficient, basis_vector in zip(reduced, basis)),
            start=CUBIC7_ZERO,
        )
        for index in range(ambient_rank)
    )


def _projective_key(vector: Sequence[Cubic7]) -> tuple[Cubic7, ...]:
    pivot = next((entry for entry in vector if entry), None)
    if pivot is None:
        raise ValueError("zero vector is not a projective ray")
    return tuple(entry / pivot for entry in vector)


def _functional_value(functional: Sequence[Cubic7], vector: Sequence[Cubic7]) -> Cubic7:
    """Sparse exact functional evaluation (the cone rows have <=2 entries)."""

    return sum(
        (
            coefficient * coordinate
            for coefficient, coordinate in zip(functional, vector)
            if coefficient and coordinate
        ),
        start=CUBIC7_ZERO,
    )


@dataclass(frozen=True, slots=True)
class Cubic7BnA1ExtremeRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Cubic7
    kind: str
    finite_link_verified: bool


@dataclass(frozen=True, slots=True)
class Cubic7BnA1Geometry:
    extreme_rays: tuple[Cubic7BnA1ExtremeRay, ...]
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    irredundant: bool
    consistent_time_sheet: bool
    finite_volume: bool
    compact: bool


def _classify_bna1_ray_pairs(
    seed: Cubic7BnA1Seed,
    dimension: int,
    functionals: Sequence[Vector],
    ray_pairs: Sequence[tuple[Vector, Sequence[int]]],
) -> Cubic7BnA1Geometry:
    """Classify an already exhaustive projective-ray set exactly."""

    rays: dict[tuple[Cubic7, ...], Cubic7BnA1ExtremeRay] = {}
    for vector, active_facets in ray_pairs:
        norm = bna1_diagonal_pair(seed, dimension, vector, vector)
        sign = norm.sign(0)
        kind = "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal"
        ray = Cubic7BnA1ExtremeRay(
            vector=vector,
            active_facets=tuple(sorted(active_facets)),
            norm=norm,
            kind=kind,
            # A rank-d active normal space perpendicular to a nonspacelike
            # ray has the required elliptic/parabolic link by Lorentz
            # signature.  Every producer supplies such a rank-d support.
            finite_link_verified=True,
        )
        rays.setdefault(_projective_key(vector), ray)

    ordered = tuple(
        sorted(
            rays.values(),
            key=lambda ray: (_projective_key(ray.vector), ray.active_facets),
        )
    )
    finite = sum(ray.kind == "finite" for ray in ordered)
    cusps = sum(ray.kind == "cusp" for ray in ordered)
    hyperideal = sum(ray.kind == "hyperideal" for ray in ordered)
    irredundant = bool(ordered) and all(
        any(
            _functional_value(functionals[facet], ray.vector).sign(0) < 0
            for ray in ordered
        )
        for facet in range(len(functionals))
    )
    nonspacelike = tuple(ray for ray in ordered if ray.kind != "hyperideal")
    sheet_ok = all(
        bna1_diagonal_pair(seed, dimension, left.vector, right.vector).sign(0) <= 0
        for index, left in enumerate(nonspacelike)
        for right in nonspacelike[index + 1 :]
    )
    full_rank = _row_rank(functionals, dimension + 1) == dimension + 1
    links_ok = all(
        ray.finite_link_verified for ray in ordered if ray.kind != "hyperideal"
    )
    finite_volume = (
        bool(ordered)
        and full_rank
        and irredundant
        and sheet_ok
        and links_ok
        and hyperideal == 0
    )
    return Cubic7BnA1Geometry(
        extreme_rays=ordered,
        finite_vertex_count=finite,
        cusp_count=cusps,
        hyperideal_ray_count=hyperideal,
        irredundant=irredundant,
        consistent_time_sheet=sheet_ok,
        finite_volume=finite_volume,
        compact=finite_volume and cusps == 0,
    )


def analyze_bna1_cone(
    seed: Cubic7BnA1Seed,
    dimension: int,
    walls: Sequence[Vector],
    *,
    combination_limit: int = 2_000_000,
) -> Cubic7BnA1Geometry:
    """Enumerate every exact extreme ray of the current wall cone."""

    ambient_rank = dimension + 1
    facet_count = len(walls)
    combination_count = 0
    if facet_count >= dimension:
        # Avoid importing floating point just to compute a binomial integer.
        numerator = 1
        denominator = 1
        for index in range(1, dimension + 1):
            numerator *= facet_count - dimension + index
            denominator *= index
        combination_count = numerator // denominator
    if combination_count > combination_limit:
        raise RuntimeError(
            f"exact cone enumeration needs {combination_count} supports, above limit {combination_limit}"
        )
    coefficients = seed.coefficients(dimension)
    functionals = tuple(
        tuple(coefficients[index] * wall[index] for index in range(ambient_rank))
        for wall in walls
    )
    rays: dict[tuple[Cubic7, ...], tuple[Vector, tuple[int, ...]]] = {}
    for support in combinations(range(facet_count), dimension):
        vector = _bna1_support_kernel(dimension, functionals, support)
        if vector is None:
            continue
        pairings = tuple(_functional_value(row, vector) for row in functionals)
        if all(value.sign(0) <= 0 for value in pairings):
            oriented = vector
            oriented_pairings = pairings
        elif all(value.sign(0) >= 0 for value in pairings):
            oriented = tuple(-entry for entry in vector)
            oriented_pairings = tuple(-entry for entry in pairings)
        else:
            continue
        active = tuple(index for index, value in enumerate(oriented_pairings) if not value)
        # The selected support already has rank ``dimension``: the reduced
        # system above had nullity exactly one.  Extra zero pairings can only
        # enlarge ``active`` without changing that fact.
        rays.setdefault(_projective_key(oriented), (oriented, active))

    return _classify_bna1_ray_pairs(
        seed,
        dimension,
        functionals,
        tuple(rays.values()),
    )


def analyze_bna1_cone_incremental(
    seed: Cubic7BnA1Seed,
    dimension: int,
    walls: Sequence[Vector],
    *,
    ray_limit: int | None = None,
    time_limit_seconds: float | None = None,
) -> Cubic7BnA1Geometry:
    """Enumerate all exact rays by incremental double description.

    The support enumerator above is a useful independent oracle, but its
    ``binomial(f,d)`` support count is unsuitable once an H22 prefix has more
    than a handful of positive-height walls.  Here the ``d`` height-zero
    simple roots and the first positive-height root form a canonical
    independent seed cone.  Each remaining half-space is inserted exactly;
    two rays are adjacent precisely when their common active functionals have
    rank ``d-1``.  This is the standard double-description adjacency test and
    remains valid for nonsimple intermediate cones.
    """

    ambient_rank = dimension + 1
    if ray_limit is not None and ray_limit < ambient_rank:
        raise ValueError("ray_limit must accommodate the initial simplicial cone")
    if time_limit_seconds is not None and time_limit_seconds <= 0:
        raise ValueError("time_limit_seconds must be positive")
    started = monotonic()

    def check_limits(current_count: int, facet_index: int) -> None:
        if ray_limit is not None and current_count > ray_limit:
            raise RuntimeError(
                f"incremental ray limit {ray_limit} exceeded while inserting facet "
                f"{facet_index + 1}: at least {current_count} rays"
            )
        if (
            time_limit_seconds is not None
            and monotonic() - started > time_limit_seconds
        ):
            raise RuntimeError(
                f"incremental time limit {time_limit_seconds:g}s exceeded while "
                f"inserting facet {facet_index + 1}: at least {current_count} rays"
            )
    if len(walls) < ambient_rank:
        raise ValueError("a full-dimensional pointed wall cone needs at least d+1 walls")
    coefficients = seed.coefficients(dimension)
    functionals = tuple(
        tuple(coefficients[index] * wall[index] for index in range(ambient_rank))
        for wall in walls
    )
    initial = tuple(range(ambient_rank))
    if _row_rank(tuple(functionals[index] for index in initial), ambient_rank) != ambient_rank:
        raise ValueError(
            "the B_(d-1)xA1 walls plus first positive wall are not independent"
        )

    processed = list(initial)
    rays: dict[tuple[Cubic7, ...], tuple[Vector, frozenset[int]]] = {}
    for omitted in initial:
        support = tuple(index for index in initial if index != omitted)
        ray = _kernel(tuple(functionals[index] for index in support))
        if ray is None:
            raise AssertionError("an independent seed lost rank after one row was omitted")
        values = tuple(_functional_value(functionals[index], ray) for index in initial)
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-entry for entry in ray)
        else:  # pragma: no cover - impossible for a simplicial dual basis
            raise AssertionError("independent seed cone has inconsistent ray orientation")
        active = frozenset(
            index
            for index in initial
            if not _functional_value(functionals[index], oriented)
        )
        rays[_projective_key(oriented)] = (oriented, active)

    rank_cache: dict[frozenset[int], int] = {}
    for new_facet in range(ambient_rank, len(walls)):
        functional = functionals[new_facet]
        negative: list[tuple[tuple[Cubic7, ...], Vector, frozenset[int], Cubic7]] = []
        zero: list[tuple[tuple[Cubic7, ...], Vector, frozenset[int], Cubic7]] = []
        positive: list[tuple[tuple[Cubic7, ...], Vector, frozenset[int], Cubic7]] = []
        for key, (ray, active) in rays.items():
            value = _functional_value(functional, ray)
            destination = (
                negative
                if value.sign(0) < 0
                else positive
                if value.sign(0) > 0
                else zero
            )
            destination.append((key, ray, active, value))

        updated: dict[tuple[Cubic7, ...], tuple[Vector, frozenset[int]]] = {}
        # Retained rays already have canonical projective keys.  Reusing them
        # avoids repeatedly normalizing the increasingly large exact vectors.
        for key, ray, active, _value in negative:
            updated[key] = (ray, active)
            check_limits(len(updated), new_facet)
        for key, ray, active, _value in zero:
            updated[key] = (ray, active | frozenset((new_facet,)))
            check_limits(len(updated), new_facet)

        pair_count = 0
        for _negative_key, negative_ray, negative_active, negative_value in negative:
            for _positive_key, positive_ray, positive_active, positive_value in positive:
                pair_count += 1
                # Exact field row reduction and projective normalization can
                # become expensive long before a complete negative-ray pass
                # finishes.  Check bounded-search contracts inside that pass,
                # not only after it, so a caller's wall-clock cap is effective.
                if pair_count % 64 == 0:
                    check_limits(len(updated), new_facet)
                common = negative_active & positive_active
                if len(common) < ambient_rank - 2:
                    continue
                common_rank = rank_cache.get(common)
                if common_rank is None:
                    common_rank = _row_rank(
                        tuple(functionals[index] for index in common),
                        ambient_rank,
                    )
                    rank_cache[common] = common_rank
                if common_rank != ambient_rank - 2:
                    continue
                ray = tuple(
                    positive_value * negative_ray[index]
                    - negative_value * positive_ray[index]
                    for index in range(ambient_rank)
                )
                actual = frozenset(
                    index
                    for index in (*processed, new_facet)
                    if not _functional_value(functionals[index], ray)
                )
                updated[_projective_key(ray)] = (ray, actual)
            check_limits(len(updated), new_facet)
        rays = updated
        processed.append(new_facet)
        check_limits(len(rays), new_facet)

    if any(
        _row_rank(tuple(functionals[index] for index in active), ambient_rank)
        < dimension
        for _ray, active in rays.values()
    ):
        raise AssertionError("double description returned a non-extreme projective ray")
    return _classify_bna1_ray_pairs(
        seed,
        dimension,
        functionals,
        tuple((ray, tuple(active)) for ray, active in rays.values()),
    )


@dataclass(frozen=True, slots=True)
class Cubic7BnA1VinbergSearch:
    pool: Cubic7BnA1RootPool
    accepted_roots: tuple[Cubic7BnA1Root, ...]
    rejected_roots: tuple[Cubic7BnA1Root, ...]
    walls: tuple[Vector, ...]
    wall_gram: Matrix
    geometry: Cubic7BnA1Geometry | None


def run_bna1_vinberg_search(
    seed: Cubic7BnA1Seed,
    dimension: int,
    height_bound: int,
    *,
    analyze_geometry: bool = True,
    combination_limit: int = 2_000_000,
    geometry_engine: str = "support",
    geometry_ray_limit: int | None = None,
    geometry_time_limit_seconds: float | None = None,
) -> Cubic7BnA1VinbergSearch:
    if geometry_engine not in {"support", "incremental"}:
        raise ValueError("geometry_engine must be 'support' or 'incremental'")
    pool = enumerate_bna1_root_pool(seed, dimension, height_bound)
    walls = list(pool.cone_roots)
    accepted: list[Cubic7BnA1Root] = []
    rejected: list[Cubic7BnA1Root] = []
    # Vinberg's ordering is by distance/height, not by an arbitrary tie-break
    # among equal-height mirrors.  Select a whole height level against only
    # earlier walls, then prove its survivors are mutually acute before
    # adjoining the level.  This makes the result invariant under any
    # permutation of same-height candidates.
    for _height, level_iterator in groupby(pool.roots, key=lambda root: root.height):
        level = tuple(level_iterator)
        survivors = tuple(
            candidate
            for candidate in level
            if all(
                bna1_diagonal_pair(seed, dimension, candidate.vector, wall).sign(0) <= 0
                for wall in walls
            )
        )
        if any(
            bna1_diagonal_pair(seed, dimension, left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError(
                "same-height Vinberg survivors are not mutually acute; "
                "the normalization or height-zero chamber is incomplete"
            )
        survivor_set = set(survivors)
        accepted.extend(survivors)
        rejected.extend(candidate for candidate in level if candidate not in survivor_set)
        walls.extend(candidate.vector for candidate in survivors)
    if any(
        bna1_diagonal_pair(seed, dimension, left, right).sign(0) > 0
        for index, left in enumerate(walls)
        for right in walls[index + 1 :]
    ):
        raise AssertionError("the accepted Vinberg walls are not pairwise acute")
    wall_gram: Matrix = tuple(
        tuple(bna1_diagonal_pair(seed, dimension, left, right) for right in walls)
        for left in walls
    )
    if not analyze_geometry:
        geometry = None
    elif geometry_engine == "incremental":
        geometry = analyze_bna1_cone_incremental(
            seed,
            dimension,
            walls,
            ray_limit=geometry_ray_limit,
            time_limit_seconds=geometry_time_limit_seconds,
        )
    else:
        if geometry_ray_limit is not None or geometry_time_limit_seconds is not None:
            raise ValueError("ray/time limits are available only for incremental geometry")
        geometry = analyze_bna1_cone(
            seed,
            dimension,
            walls,
            combination_limit=combination_limit,
        )
    return Cubic7BnA1VinbergSearch(
        pool=pool,
        accepted_roots=tuple(accepted),
        rejected_roots=tuple(rejected),
        walls=tuple(walls),
        wall_gram=wall_gram,
        geometry=geometry,
    )


__all__ = [
    "DELTA7",
    "L2_BNA1_SEED",
    "L7_BNA1_SEED",
    "Cubic7BnA1ExtremeRay",
    "Cubic7BnA1Geometry",
    "Cubic7BnA1Root",
    "Cubic7BnA1RootPool",
    "Cubic7BnA1Seed",
    "Cubic7BnA1VinbergSearch",
    "HeightZeroPattern",
    "analyze_bna1_cone",
    "analyze_bna1_cone_incremental",
    "bna1_diagonal_pair",
    "dominant_bna1_representative",
    "enumerate_bna1_root_pool",
    "enumerate_height_zero_patterns",
    "height_zero_is_bn_minus_1_times_a1",
    "run_bna1_vinberg_search",
]
