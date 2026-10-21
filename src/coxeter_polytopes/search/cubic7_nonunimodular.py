"""Exact Vinberg prefixes for small-prime non-unimodular cubic-seven seeds.

The lattices in this module are

``L_2 = diag(-theta, 1, 1, 1, 2)`` and
``L_7 = diag(-theta, 1, 1, 1, 2-theta)``, and three forms whose
last coefficient generates one of the primes over 13,

where ``theta=2*cos(2*pi/7)``.  All arithmetic and all archimedean bounds are
exact.  More importantly, the root search does not guess a norm cutoff: it
uses the complete list of primitive mirror norms, modulo unit-square scaling.

Here is the local argument.  For a primitive vector ``r`` put

``I(r)=(r_0,r_1,r_2,r_3,t*r_4)``.

The reflection in ``r`` preserves the lattice exactly when
``(r,r) | 2*I(r)``.  Away from primes dividing ``t``, primitivity makes
``I(r)`` the unit ideal.  At a prime occurring once in ``t``, either
``I(r)=O`` or all first four coordinates are divisible by that prime and
``r_4`` is a unit.  In the latter case the diagonal norm has valuation
exactly one there.  Since 2 is inert, this gives

* ``L_2``: norm ideals ``(1)`` or ``(2)``;
* if ``t`` generates a prime distinct from 2: norm ideals
  ``(1)``, ``(2)``, ``(t)``, or ``(2*t)``.

The narrow class number is one and every totally positive unit is a square,
so the displayed totally positive generators are the complete canonical
element norms.  They are all squarefree ideals, hence a vector having one of
these norms is automatically primitive.  Coordinatewise reflection
integrality is nevertheless checked for every enumerated vector.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import Sequence

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    THETA7,
    Cubic7,
)
from coxeter_polytopes.search.cubic7_diagonal import bounded_cubic7_integers


Vector = tuple[Cubic7, ...]
DELTA7 = 2 - THETA7


@dataclass(frozen=True, slots=True)
class Cubic7NonunimodularSeed:
    seed_id: str
    last_coefficient: Cubic7
    canonical_mirror_norms: tuple[Cubic7, ...]
    local_norm_proof: str

    @property
    def coefficients(self) -> tuple[Cubic7, ...]:
        return (-THETA7, CUBIC7_ONE, CUBIC7_ONE, CUBIC7_ONE, self.last_coefficient)


L2_SEED = Cubic7NonunimodularSeed(
    seed_id="L2",
    last_coefficient=Cubic7(2),
    canonical_mirror_norms=(Cubic7(1), Cubic7(2)),
    local_norm_proof=(
        "I(r) is O or (2).  In the second case r0,...,r3 are in (2), "
        "r4 is a 2-adic unit, and q(r)=2*r4^2 mod (4), so v_2(q)=1. "
        "Thus (q) is (1) or (2)."
    ),
)

L7_SEED = Cubic7NonunimodularSeed(
    seed_id="L7",
    last_coefficient=DELTA7,
    canonical_mirror_norms=(Cubic7(1), Cubic7(2), DELTA7, 2 * DELTA7),
    local_norm_proof=(
        "For p=(2-theta), I(r) is O or p.  In the latter case "
        "q(r)=(2-theta)*r4^2 mod p^2, hence v_p(q)=1; the inert prime "
        "(2) can occur to exponent at most one."
    ),
)


def _single_prime_seed(seed_id: str, coefficient: Cubic7, prime_label: str) -> Cubic7NonunimodularSeed:
    if coefficient.norm() not in (7, 13):
        raise ValueError("a single-prime seed coefficient must have prime field norm")
    if any(coefficient.sign(embedding) <= 0 for embedding in range(3)):
        raise ValueError("a diagonal seed coefficient must be totally positive")
    return Cubic7NonunimodularSeed(
        seed_id=seed_id,
        last_coefficient=coefficient,
        canonical_mirror_norms=(Cubic7(1), Cubic7(2), coefficient, 2 * coefficient),
        local_norm_proof=(
            f"For p=({prime_label})=(t), N(p)={int(coefficient.norm())}, "
            "I(r) is O or p. If I(r)=p, then r0,...,r3 are in p and r4 "
            "is a p-adic unit, so q(r)=t*r4^2 mod p^2 and v_p(q)=1. "
            "If I(r)=O, reflection integrality forces v_p(q)=0. The inert "
            "prime (2) occurs to exponent at most one. Hence (q) is exactly "
            "one of (1),(2),p,2p; narrow class number one and squarehood of "
            "every totally positive unit give element representatives "
            "1,2,t,2t up to unit scaling."
        ),
    )


P13_R7_SEED = _single_prime_seed("P13-r7", Cubic7(0, -1, 2), "-theta+2theta^2")
P13_R10_SEED = _single_prime_seed("P13-r10", Cubic7(1, -1, 1), "1-theta+theta^2")
P13_R8_SEED = _single_prime_seed("P13-r8", Cubic7(1, 0, 1), "1+theta^2")
P13_SEEDS = (P13_R7_SEED, P13_R10_SEED, P13_R8_SEED)
SUPPORTED_NONUNIMODULAR_SEEDS = (L2_SEED, L7_SEED, *P13_SEEDS)


def nonunimodular_diagonal_pair(
    seed: Cubic7NonunimodularSeed,
    left: Sequence[Cubic7],
    right: Sequence[Cubic7],
) -> Cubic7:
    if len(left) != 5 or len(right) != 5:
        raise ValueError("the initial non-unimodular seeds have rank five")
    return sum(
        (coefficient * x * y for coefficient, x, y in zip(seed.coefficients, left, right)),
        start=CUBIC7_ZERO,
    )


def _embedding_integer_upper(value: Cubic7, embedding: int) -> int:
    """Return a certified integer upper bound at one real embedding."""

    if value.sign(embedding) < 0:
        raise ValueError("an archimedean square bound must be nonnegative")
    if not value:
        return 0
    upper = 1
    while (Cubic7(upper) - value).sign(embedding) < 0:
        upper *= 2
    return upper


@lru_cache(maxsize=None)
def _bounded_for_embedding_values_cached(
    values: tuple[Cubic7, Cubic7, Cubic7],
) -> tuple[Cubic7, ...]:
    return bounded_cubic7_integers(
        tuple(_embedding_integer_upper(value, embedding) for embedding, value in enumerate(values))
    )


def _bounded_for_embedding_values(values: Sequence[Cubic7]) -> tuple[Cubic7, ...]:
    if len(values) != 3:
        raise ValueError("one bound per real embedding is required")
    return _bounded_for_embedding_values_cached(tuple(values))  # type: ignore[arg-type]


def _reflection_is_integral(
    seed: Cubic7NonunimodularSeed, vector: Vector, norm: Cubic7
) -> bool:
    return all(
        (2 * coefficient * coordinate / norm).is_integral
        for coefficient, coordinate in zip(seed.coefficients, vector)
    )


def _three_square_solutions(
    alphabet: Sequence[Cubic7], target: Cubic7
) -> tuple[tuple[Cubic7, Cubic7, Cubic7], ...]:
    """Solve three squares exactly, retaining only the B3 dominant cone."""

    by_square: dict[Cubic7, tuple[Cubic7, ...]] = {}
    mutable: dict[Cubic7, list[Cubic7]] = {}
    for value in alphabet:
        mutable.setdefault(value * value, []).append(value)
    by_square = {square: tuple(values) for square, values in mutable.items()}

    solutions: list[tuple[Cubic7, Cubic7, Cubic7]] = []
    for first in alphabet:
        first_remainder = target - first * first
        if any(first_remainder.sign(embedding) < 0 for embedding in range(3)):
            continue
        for second in alphabet:
            if (first - second).sign(0) < 0:
                continue
            remainder = first_remainder - second * second
            if any(remainder.sign(embedding) < 0 for embedding in range(3)):
                continue
            for third in by_square.get(remainder, ()):
                if (second - third).sign(0) >= 0:
                    solutions.append((first, second, third))
    return tuple(solutions)


def _spatial_solutions(
    seed: Cubic7NonunimodularSeed,
    norm: Cubic7,
    height_bound: int,
    target: Cubic7,
) -> tuple[tuple[Cubic7, Cubic7, Cubic7, Cubic7], ...]:
    # At the distinguished place target <= norm*(H+1); at the two other
    # places theta is negative, so target <= norm.  These are field elements
    # used only at their indicated embeddings.
    target_caps = (norm * (height_bound + 1), norm, norm)
    ordinary_pool, last_pool = _spatial_alphabets(seed, norm, height_bound)
    ordinary = tuple(
        value
        for value in ordinary_pool
        if all((target - value * value).sign(i) >= 0 for i in range(3))
    )
    last_values = tuple(
        value
        for value in last_pool
        if all(
            (target - seed.last_coefficient * value * value).sign(i) >= 0
            for i in range(3)
        )
    )

    solutions: list[tuple[Cubic7, Cubic7, Cubic7, Cubic7]] = []
    for last in last_values:
        remainder = target - seed.last_coefficient * last * last
        for first, second, third in _three_square_solutions(ordinary, remainder):
            solutions.append((first, second, third, last))
    return tuple(solutions)


@lru_cache(maxsize=None)
def _spatial_alphabets(
    seed: Cubic7NonunimodularSeed, norm: Cubic7, height_bound: int
) -> tuple[tuple[Cubic7, ...], tuple[Cubic7, ...]]:
    target_caps = (norm * (height_bound + 1), norm, norm)
    ordinary = tuple(
        value for value in _bounded_for_embedding_values(target_caps) if value.sign(0) >= 0
    )
    weighted_caps = tuple(cap / seed.last_coefficient for cap in target_caps)
    last = tuple(
        value for value in _bounded_for_embedding_values(weighted_caps) if value.sign(0) >= 0
    )
    return ordinary, last


def _cone_roots(seed: Cubic7NonunimodularSeed) -> tuple[Vector, ...]:
    # Once height-zero exhaustion verifies B3 x A1, this is its dominant
    # chamber y1 >= y2 >= y3 >= 0 and y4 >= 0.
    return (
        (CUBIC7_ZERO, -CUBIC7_ONE, CUBIC7_ONE, CUBIC7_ZERO, CUBIC7_ZERO),
        (CUBIC7_ZERO, CUBIC7_ZERO, -CUBIC7_ONE, CUBIC7_ONE, CUBIC7_ZERO),
        (CUBIC7_ZERO, CUBIC7_ZERO, CUBIC7_ZERO, -CUBIC7_ONE, CUBIC7_ZERO),
        (CUBIC7_ZERO, CUBIC7_ZERO, CUBIC7_ZERO, CUBIC7_ZERO, -CUBIC7_ONE),
    )


@lru_cache(maxsize=None)
def enumerate_height_zero_mirrors(seed: Cubic7NonunimodularSeed) -> tuple[Vector, ...]:
    """Exhaust every normalized root perpendicular to the control vector."""

    roots: set[Vector] = set()
    for norm in seed.canonical_mirror_norms:
        ordinary = _bounded_for_embedding_values((norm, norm, norm))
        weighted = _bounded_for_embedding_values(
            tuple(norm / seed.last_coefficient for _ in range(3))
        )
        for last in weighted:
            remainder = norm - seed.last_coefficient * last * last
            if any(remainder.sign(i) < 0 for i in range(3)):
                continue
            # Height-zero pools are tiny, and retaining all signs here makes
            # the B3 x A1 verification independent of the chamber assumption.
            for first in ordinary:
                remainder_1 = remainder - first * first
                if any(remainder_1.sign(i) < 0 for i in range(3)):
                    continue
                for second in ordinary:
                    remainder_2 = remainder_1 - second * second
                    if any(remainder_2.sign(i) < 0 for i in range(3)):
                        continue
                    for third in ordinary:
                        if third * third != remainder_2:
                            continue
                        vector = (CUBIC7_ZERO, first, second, third, last)
                        if _reflection_is_integral(seed, vector, norm):
                            roots.add(vector)
    return tuple(sorted(roots))


def height_zero_is_b3_times_a1(seed: Cubic7NonunimodularSeed) -> bool:
    expected: set[Vector] = set()
    basis = [
        tuple(CUBIC7_ONE if i == j else CUBIC7_ZERO for i in range(5))
        for j in range(5)
    ]
    # B3: short coordinate roots and long signed pairs in slots 1,2,3.
    for i in (1, 2, 3):
        expected.add(basis[i])
        expected.add(tuple(-entry for entry in basis[i]))
    for i in (1, 2, 3):
        for j in range(i + 1, 4):
            for sign_i in (-1, 1):
                for sign_j in (-1, 1):
                    vector = [CUBIC7_ZERO] * 5
                    vector[i] = Cubic7(sign_i)
                    vector[j] = Cubic7(sign_j)
                    expected.add(tuple(vector))
    # The final coordinate gives A1, of norm 2 or 2-theta respectively.
    expected.add(basis[4])
    expected.add(tuple(-entry for entry in basis[4]))
    return set(enumerate_height_zero_mirrors(seed)) == expected


@dataclass(frozen=True, slots=True)
class Cubic7NonunimodularRoot:
    norm: Cubic7
    height: Cubic7
    vector: Vector

    @property
    def control(self) -> Cubic7:
        return self.vector[0]


@dataclass(frozen=True, slots=True)
class Cubic7NonunimodularRootPool:
    seed: Cubic7NonunimodularSeed
    height_bound: int
    cone_roots: tuple[Vector, ...]
    roots: tuple[Cubic7NonunimodularRoot, ...]
    control_counts_by_norm: tuple[tuple[Cubic7, int], ...]


def enumerate_nonunimodular_root_pool(
    seed: Cubic7NonunimodularSeed, height_bound: int
) -> Cubic7NonunimodularRootPool:
    """Exhaust all normalized mirrors with positive height at most ``H``."""

    if seed not in SUPPORTED_NONUNIMODULAR_SEEDS:
        raise ValueError("the seed has no registered complete local norm proof")
    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    if not height_zero_is_b3_times_a1(seed):
        raise AssertionError("the asserted height-zero Weyl chamber is incomplete")

    roots: set[Cubic7NonunimodularRoot] = set()
    control_counts: list[tuple[Cubic7, int]] = []
    for norm in seed.canonical_mirror_norms:
        control_caps = (
            norm * height_bound / THETA7,
            norm / (-THETA7),
            norm / (-THETA7),
        )
        controls = tuple(
            control
            for control in _bounded_for_embedding_values(control_caps)
            if control.sign(0) > 0
            and (norm * height_bound - THETA7 * control * control).sign(0) >= 0
            and all((norm + THETA7 * control * control).sign(i) >= 0 for i in (1, 2))
        )
        control_counts.append((norm, len(controls)))
        for control in controls:
            target = norm + THETA7 * control * control
            for spatial in _spatial_solutions(seed, norm, height_bound, target):
                vector = (control, *spatial)
                if nonunimodular_diagonal_pair(seed, vector, vector) != norm:
                    raise AssertionError("weighted enumerator returned the wrong norm")
                if not _reflection_is_integral(seed, vector, norm):
                    continue
                roots.add(
                    Cubic7NonunimodularRoot(
                        norm=norm,
                        height=THETA7 * control * control / norm,
                        vector=vector,
                    )
                )

    return Cubic7NonunimodularRootPool(
        seed=seed,
        height_bound=height_bound,
        cone_roots=_cone_roots(seed),
        roots=tuple(sorted(roots, key=lambda root: (root.height, root.norm, root.vector))),
        control_counts_by_norm=tuple(control_counts),
    )


@dataclass(frozen=True, slots=True)
class Cubic7NonunimodularVinbergSearch:
    pool: Cubic7NonunimodularRootPool
    accepted_roots: tuple[Cubic7NonunimodularRoot, ...]
    rejected_roots: tuple[Cubic7NonunimodularRoot, ...]
    walls: tuple[Vector, ...]
    wall_gram: tuple[tuple[Cubic7, ...], ...]


@dataclass(frozen=True, slots=True)
class Cubic7ConeRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Cubic7
    kind: str


@dataclass(frozen=True, slots=True)
class Cubic7ConeRayCensus:
    rays: tuple[Cubic7ConeRay, ...]

    @property
    def finite_count(self) -> int:
        return sum(ray.kind == "finite" for ray in self.rays)

    @property
    def cusp_count(self) -> int:
        return sum(ray.kind == "cusp" for ray in self.rays)

    @property
    def hyperideal_count(self) -> int:
        return sum(ray.kind == "hyperideal" for ray in self.rays)


def _dot(left: Sequence[Cubic7], right: Sequence[Cubic7]) -> Cubic7:
    return sum((x * y for x, y in zip(left, right)), start=CUBIC7_ZERO)


def _row_rank(rows: Sequence[Sequence[Cubic7]], width: int) -> int:
    work = [list(row) for row in rows]
    rank = 0
    for column in range(width):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for row in range(rank + 1, len(work)):
            if not work[row][column]:
                continue
            factor = work[row][column] / pivot_value
            for index in range(column, width):
                work[row][index] -= factor * work[rank][index]
        rank += 1
        if rank == len(work):
            break
    return rank


def _kernel_ray(rows: Sequence[Sequence[Cubic7]], width: int) -> Vector | None:
    """One exact kernel vector for a corank-one row matrix."""

    work = [list(row) for row in rows]
    pivot_columns: list[int] = []
    rank = 0
    for column in range(width):
        pivot = next((i for i in range(rank, len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [entry / pivot_value for entry in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[rank][index] for index in range(width)
            ]
        pivot_columns.append(column)
        rank += 1
        if rank == len(work):
            break
    if rank != width - 1:
        return None
    free = next(column for column in range(width) if column not in pivot_columns)
    result = [CUBIC7_ZERO for _ in range(width)]
    result[free] = CUBIC7_ONE
    for row, column in enumerate(pivot_columns):
        result[column] = -work[row][free]
    return tuple(result)


def _projective_key(vector: Sequence[Cubic7]) -> Vector:
    pivot = next((entry for entry in vector if entry), None)
    if pivot is None:
        raise ValueError("zero is not a projective ray")
    return tuple(entry / pivot for entry in vector)


def enumerate_cone_rays_incremental(
    search: Cubic7NonunimodularVinbergSearch,
) -> Cubic7ConeRayCensus:
    """Classify the current wall cone by exact double description.

    This is equivalent to intersecting every four-facet subset in rank five,
    but its running time is governed by the actual ray set instead of by
    ``binomial(number_of_walls, 4)``.  It is cross-checked against the brute
    cofactor engine in the test suite.
    """

    seed = search.pool.seed
    roots = search.walls
    width = 5
    functionals = tuple(
        tuple(coefficient * coordinate for coefficient, coordinate in zip(seed.coefficients, root))
        for root in roots
    )
    independent = next(
        (
            subset
            for subset in combinations(range(len(roots)), width)
            if _row_rank(tuple(functionals[index] for index in subset), width) == width
        ),
        None,
    )
    if independent is None:
        raise ValueError("wall functionals do not span the Lorentz space")
    independent_set = set(independent)
    order = (*independent, *(i for i in range(len(roots)) if i not in independent_set))
    processed = list(independent)
    rays: dict[Vector, tuple[Vector, frozenset[int]]] = {}
    for omitted in independent:
        active = tuple(index for index in independent if index != omitted)
        ray = _kernel_ray(tuple(functionals[index] for index in active), width)
        if ray is None:
            raise AssertionError("a proper subset of the independent seed was dependent")
        values = tuple(_dot(ray, functionals[index]) for index in independent)
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-entry for entry in ray)
        else:
            raise AssertionError("independent seed cone has inconsistent orientation")
        actual = frozenset(index for index in independent if not _dot(oriented, functionals[index]))
        rays[_projective_key(oriented)] = (oriented, actual)

    rank_cache: dict[frozenset[int], int] = {}
    for new_facet in order[width:]:
        functional = functionals[new_facet]
        negative = []
        zero = []
        positive = []
        for ray, active in rays.values():
            value = _dot(ray, functional)
            (negative if value.sign(0) < 0 else positive if value.sign(0) > 0 else zero).append(
                (ray, active, value)
            )
        updated: dict[Vector, tuple[Vector, frozenset[int]]] = {}
        for ray, active, _ in (*negative, *zero):
            updated[_projective_key(ray)] = (
                ray,
                active | (frozenset((new_facet,)) if not _dot(ray, functional) else frozenset()),
            )
        for negative_ray, negative_active, negative_value in negative:
            for positive_ray, positive_active, positive_value in positive:
                common = negative_active & positive_active
                common_rank = rank_cache.get(common)
                if common_rank is None:
                    common_rank = _row_rank(tuple(functionals[index] for index in common), width)
                    rank_cache[common] = common_rank
                if common_rank != width - 2:
                    continue
                ray = tuple(
                    positive_value * negative_ray[index] - negative_value * positive_ray[index]
                    for index in range(width)
                )
                actual = frozenset(
                    index
                    for index in (*processed, new_facet)
                    if not _dot(ray, functionals[index])
                )
                updated[_projective_key(ray)] = (ray, actual)
        rays = updated
        processed.append(new_facet)

    classified = []
    for ray, active in rays.values():
        norm = nonunimodular_diagonal_pair(seed, ray, ray)
        kind = "finite" if norm.sign(0) < 0 else "cusp" if not norm else "hyperideal"
        classified.append(Cubic7ConeRay(ray, tuple(sorted(active)), norm, kind))
    return Cubic7ConeRayCensus(
        tuple(sorted(classified, key=lambda ray: (_projective_key(ray.vector), ray.active_facets)))
    )


def run_nonunimodular_vinberg_search(
    seed: Cubic7NonunimodularSeed, height_bound: int
) -> Cubic7NonunimodularVinbergSearch:
    pool = enumerate_nonunimodular_root_pool(seed, height_bound)
    walls = list(pool.cone_roots)
    accepted: list[Cubic7NonunimodularRoot] = []
    rejected: list[Cubic7NonunimodularRoot] = []
    # Vinberg roots at one exact height must be selected simultaneously.
    # Using a secondary vector ordering would be mathematically invalid if
    # two survivors had positive mutual pairing.  Audit every height packet
    # against the lower-height walls first, then require pairwise acuteness.
    cursor = 0
    while cursor < len(pool.roots):
        height = pool.roots[cursor].height
        end = cursor + 1
        while end < len(pool.roots) and pool.roots[end].height == height:
            end += 1
        packet = pool.roots[cursor:end]
        survivors = tuple(
            candidate
            for candidate in packet
            if all(
                nonunimodular_diagonal_pair(seed, candidate.vector, wall).sign(0) <= 0
                for wall in walls
            )
        )
        if any(
            nonunimodular_diagonal_pair(seed, left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise ArithmeticError(
                "same-height Vinberg survivors are not pairwise acute; "
                "a Weyl-orbit resolution is required"
            )
        rejected.extend(candidate for candidate in packet if candidate not in survivors)
        accepted.extend(survivors)
        walls.extend(candidate.vector for candidate in survivors)
        cursor = end
    gram = tuple(
        tuple(nonunimodular_diagonal_pair(seed, left, right) for right in walls)
        for left in walls
    )
    return Cubic7NonunimodularVinbergSearch(
        pool=pool,
        accepted_roots=tuple(accepted),
        rejected_roots=tuple(rejected),
        walls=tuple(walls),
        wall_gram=gram,
    )


__all__ = [
    "DELTA7",
    "L2_SEED",
    "L7_SEED",
    "P13_R7_SEED",
    "P13_R10_SEED",
    "P13_R8_SEED",
    "P13_SEEDS",
    "SUPPORTED_NONUNIMODULAR_SEEDS",
    "Cubic7NonunimodularRoot",
    "Cubic7NonunimodularRootPool",
    "Cubic7NonunimodularSeed",
    "Cubic7NonunimodularVinbergSearch",
    "Cubic7ConeRay",
    "Cubic7ConeRayCensus",
    "enumerate_cone_rays_incremental",
    "enumerate_height_zero_mirrors",
    "enumerate_nonunimodular_root_pool",
    "height_zero_is_b3_times_a1",
    "nonunimodular_diagonal_pair",
    "run_nonunimodular_vinberg_search",
]
