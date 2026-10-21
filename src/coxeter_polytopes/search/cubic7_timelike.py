"""Exact Vinberg search for two timelike-scaled conductor-seven lattices.

Put ``theta = 2*cos(2*pi/7)`` and ``delta = 2-theta``.  This module studies

``T2 = diag(-2*theta, 1, 1, 1, 1)`` and
``T7 = diag(-delta*theta, 1, 1, 1, 1)``.

The mirror norms below are a complete quotient, not a search cutoff.  If
``r=(x,y1,...,y4)`` is primitive and ``q=(r,r)``, reflection integrality is

``q | 2 I(r)``,  where  ``I(r)=(-s*theta*x,y1,...,y4)``.

At a prime not dividing ``2s``, primitivity makes ``I(r)`` the unit ideal.
At the squarefree prime dividing ``s``, either a spatial coordinate is a
unit, or ``x`` is a unit and all spatial coordinates lie in that prime.  In
the second case the diagonal equation gives ``v(q)=1``.  The same argument
works at the inert prime ``(2)`` when ``s=2``.  Consequently the possible
norm ideals are ``(1),(2)`` for T2 and ``(1),(2),(delta),(2*delta)`` for T7.

For completeness, the number-theory normalization can be proved without a
database.  The Minkowski bound is

``(3!/3^3)*sqrt(49) = 14/9 < 2``,

so the class number is one.  The units ``-1``, ``theta`` and
``theta^2-2`` realize three independent sign vectors; hence the signature
map on ``O_K^*/O_K^{*2}`` is bijective.  Thus the narrow class number is one
and every totally positive unit is a square.  Scaling a root by the inverse
of that unit gives precisely the canonical element norms listed here.

All field signs, lattice points, root equations, wall pairings, and cone rays
are exact.  Equal-height Vinberg packets are selected simultaneously.  The
incremental double-description cone engine is accompanied by an independent
brute-force support enumerator for regression tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, groupby
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


def _signature_bits(value: Cubic7) -> tuple[int, int, int]:
    """Encode the three exact signs of a nonzero element as F2 bits."""

    signs = tuple(value.sign(embedding) for embedding in range(3))
    if 0 in signs:
        raise ValueError("a unit cannot vanish at a real embedding")
    return tuple(int(sign < 0) for sign in signs)  # type: ignore[return-value]


UNIT_SIGNATURE_GENERATORS = (-CUBIC7_ONE, THETA7, THETA7**2 - 2)
UNIT_SIGNATURE_BITS = tuple(_signature_bits(unit) for unit in UNIT_SIGNATURE_GENERATORS)
MINKOWSKI_BOUND_NUMERATOR = 14
MINKOWSKI_BOUND_DENOMINATOR = 9


def unit_signature_rank() -> int:
    """Rank over F2 of the displayed unit-signature matrix."""

    rows = [sum(bit << column for column, bit in enumerate(bits)) for bits in UNIT_SIGNATURE_BITS]
    rank = 0
    for column in range(3):
        pivot = next((index for index in range(rank, len(rows)) if rows[index] & (1 << column)), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and rows[index] & (1 << column):
                rows[index] ^= rows[rank]
        rank += 1
    return rank


@dataclass(frozen=True, slots=True)
class Cubic7TimelikeSeed:
    seed_id: str
    timelike_scale: Cubic7
    canonical_mirror_norms: tuple[Cubic7, ...]
    complete_norm_class_proof: str

    @property
    def coefficients(self) -> tuple[Cubic7, ...]:
        return (-self.timelike_scale * THETA7, CUBIC7_ONE, CUBIC7_ONE, CUBIC7_ONE, CUBIC7_ONE)


T2_TIMELIKE_SEED = Cubic7TimelikeSeed(
    seed_id="T2",
    timelike_scale=Cubic7(2),
    canonical_mirror_norms=(Cubic7(1), Cubic7(2)),
    complete_norm_class_proof=(
        "For primitive r, q(r) divides 2I(r), with "
        "I(r)=(-2*theta*x,y1,...,y4). Away from the inert prime (2), "
        "I(r) is locally O. At (2), either a spatial coordinate is a "
        "unit, giving v_2(q)<=1, or x is a unit and every yi lies in "
        "(2), in which case q=-2*theta*x^2 mod (4), so v_2(q)=1. "
        "Therefore (q) is (1) or (2)."
    ),
)

T7_TIMELIKE_SEED = Cubic7TimelikeSeed(
    seed_id="T7",
    timelike_scale=DELTA7,
    canonical_mirror_norms=(Cubic7(1), Cubic7(2), DELTA7, 2 * DELTA7),
    complete_norm_class_proof=(
        "For primitive r, q(r) divides 2I(r), with "
        "I(r)=(-delta*theta*x,y1,...,y4). At p=(delta), either a "
        "spatial coordinate is a unit, giving v_p(q)=0, or x is a unit "
        "and all yi lie in p, giving q=-delta*theta*x^2 mod p^2 and "
        "v_p(q)=1. Away from p and the inert prime (2), I(r) is locally "
        "O; at (2), v_2(q)<=1. Therefore (q) is one of "
        "(1),(2),(delta),(2*delta)."
    ),
)

SUPPORTED_TIMELIKE_SEEDS = (T2_TIMELIKE_SEED, T7_TIMELIKE_SEED)


def timelike_diagonal_pair(
    seed: Cubic7TimelikeSeed,
    left: Sequence[Cubic7],
    right: Sequence[Cubic7],
) -> Cubic7:
    if len(left) != 5 or len(right) != 5:
        raise ValueError("a rank-five timelike-scaled vector has five coordinates")
    return sum(
        (coefficient * x * y for coefficient, x, y in zip(seed.coefficients, left, right)),
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
    return bounded_cubic7_integers(
        tuple(_embedding_integer_upper(cap, embedding) for embedding, cap in enumerate(caps))
    )


def _reflection_is_integral(
    seed: Cubic7TimelikeSeed,
    vector: Vector,
    norm: Cubic7,
) -> bool:
    return all(
        (2 * coefficient * coordinate / norm).is_integral
        for coefficient, coordinate in zip(seed.coefficients, vector)
    )


def _dominant_square_partitions(
    values: Sequence[Cubic7], coordinate_count: int, target: Cubic7
) -> tuple[tuple[Cubic7, ...], ...]:
    """Solve a sum of squares modulo the signed-permutation group B_n."""

    if coordinate_count < 0:
        return ()
    positive = tuple(sorted((value for value in values if value.sign(0) > 0), reverse=True))
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
            if any(running.sign(embedding) < 0 for embedding in range(3)):
                break
            for suffix in solve(index + 1, slots - multiplicity, running):
                output.append((multiplicity, *suffix))
            running -= square
        return tuple(output)

    solutions: list[tuple[Cubic7, ...]] = []
    for multiplicities in solve(0, coordinate_count, target):
        used = sum(multiplicities)
        solutions.append(
            tuple(
                value
                for value, multiplicity in zip(positive, multiplicities)
                for _ in range(multiplicity)
            )
            + (CUBIC7_ZERO,) * (coordinate_count - used)
        )
    return tuple(solutions)


def _cone_roots() -> tuple[Vector, ...]:
    # Dominant B4 chamber: y1 >= y2 >= y3 >= y4 >= 0.
    return (
        (CUBIC7_ZERO, -CUBIC7_ONE, CUBIC7_ONE, CUBIC7_ZERO, CUBIC7_ZERO),
        (CUBIC7_ZERO, CUBIC7_ZERO, -CUBIC7_ONE, CUBIC7_ONE, CUBIC7_ZERO),
        (CUBIC7_ZERO, CUBIC7_ZERO, CUBIC7_ZERO, -CUBIC7_ONE, CUBIC7_ONE),
        (CUBIC7_ZERO, CUBIC7_ZERO, CUBIC7_ZERO, CUBIC7_ZERO, -CUBIC7_ONE),
    )


@lru_cache(maxsize=None)
def enumerate_height_zero_mirrors(seed: Cubic7TimelikeSeed) -> tuple[Vector, ...]:
    """Exhaust every normalized reflecting vector perpendicular to e0."""

    if seed not in SUPPORTED_TIMELIKE_SEEDS:
        raise ValueError("the seed has no registered complete norm-class proof")
    roots: set[Vector] = set()
    for norm in seed.canonical_mirror_norms:
        alphabet = _bounded_for_caps((norm, norm, norm))
        pairs: dict[Cubic7, list[tuple[Cubic7, Cubic7]]] = {}
        for first in alphabet:
            for second in alphabet:
                pairs.setdefault(first * first + second * second, []).append((first, second))
        for first_sum, first_pairs in pairs.items():
            second_pairs = pairs.get(norm - first_sum, ())
            for first, second in first_pairs:
                for third, fourth in second_pairs:
                    vector = (CUBIC7_ZERO, first, second, third, fourth)
                    if _reflection_is_integral(seed, vector, norm):
                        roots.add(vector)
    return tuple(sorted(roots))


def height_zero_is_b4(seed: Cubic7TimelikeSeed) -> bool:
    basis = [
        tuple(CUBIC7_ONE if index == coordinate else CUBIC7_ZERO for index in range(5))
        for coordinate in range(5)
    ]
    expected: set[Vector] = set()
    for coordinate in range(1, 5):
        expected.add(basis[coordinate])
        expected.add(tuple(-entry for entry in basis[coordinate]))
    for left in range(1, 5):
        for right in range(left + 1, 5):
            for left_sign in (-1, 1):
                for right_sign in (-1, 1):
                    vector = [CUBIC7_ZERO] * 5
                    vector[left] = Cubic7(left_sign)
                    vector[right] = Cubic7(right_sign)
                    expected.add(tuple(vector))
    return set(enumerate_height_zero_mirrors(seed)) == expected


@dataclass(frozen=True, slots=True)
class Cubic7TimelikeRoot:
    norm: Cubic7
    height: Cubic7
    vector: Vector

    @property
    def control(self) -> Cubic7:
        return self.vector[0]


@dataclass(frozen=True, slots=True)
class Cubic7TimelikeRootPool:
    seed: Cubic7TimelikeSeed
    height_bound: int
    cone_roots: tuple[Vector, ...]
    roots: tuple[Cubic7TimelikeRoot, ...]
    control_counts_by_norm: tuple[tuple[Cubic7, int], ...]
    spatial_alphabet_counts_by_norm: tuple[tuple[Cubic7, int], ...]


@lru_cache(maxsize=None)
def _controls(
    seed: Cubic7TimelikeSeed, norm: Cubic7, height_bound: int
) -> tuple[Cubic7, ...]:
    scale_theta = seed.timelike_scale * THETA7
    caps = (
        norm * height_bound / scale_theta,
        norm / (-scale_theta),
        norm / (-scale_theta),
    )
    return tuple(
        control
        for control in _bounded_for_caps(caps)
        if control.sign(0) > 0
        and (norm * height_bound - scale_theta * control * control).sign(0) >= 0
        and all(
            (norm + scale_theta * control * control).sign(embedding) >= 0
            for embedding in (1, 2)
        )
    )


def enumerate_timelike_root_pool(
    seed: Cubic7TimelikeSeed, height_bound: int
) -> Cubic7TimelikeRootPool:
    """Exhaust normalized positive-height roots through exact height H."""

    if seed not in SUPPORTED_TIMELIKE_SEEDS:
        raise ValueError("the seed has no registered complete norm-class proof")
    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    if not height_zero_is_b4(seed):
        raise AssertionError("the asserted B4 height-zero chamber is incomplete")

    roots: set[Cubic7TimelikeRoot] = set()
    control_counts: list[tuple[Cubic7, int]] = []
    spatial_counts: list[tuple[Cubic7, int]] = []
    scale_theta = seed.timelike_scale * THETA7
    for norm in seed.canonical_mirror_norms:
        controls = _controls(seed, norm, height_bound)
        control_counts.append((norm, len(controls)))
        target_caps = (norm * (height_bound + 1), norm, norm)
        spatial_pool = tuple(
            value for value in _bounded_for_caps(target_caps) if value.sign(0) >= 0
        )
        spatial_counts.append((norm, len(spatial_pool)))
        for control in controls:
            target = norm + scale_theta * control * control
            alphabet = tuple(
                value
                for value in spatial_pool
                if all(
                    (target - value * value).sign(embedding) >= 0
                    for embedding in range(3)
                )
            )
            for spatial in _dominant_square_partitions(alphabet, 4, target):
                vector = (control, *spatial)
                if timelike_diagonal_pair(seed, vector, vector) != norm:
                    raise AssertionError("the square-partition enumerator returned a wrong norm")
                if not _reflection_is_integral(seed, vector, norm):
                    continue
                roots.add(
                    Cubic7TimelikeRoot(
                        norm=norm,
                        height=scale_theta * control * control / norm,
                        vector=vector,
                    )
                )
    return Cubic7TimelikeRootPool(
        seed=seed,
        height_bound=height_bound,
        cone_roots=_cone_roots(),
        roots=tuple(sorted(roots, key=lambda root: (root.height, root.norm, root.vector))),
        control_counts_by_norm=tuple(control_counts),
        spatial_alphabet_counts_by_norm=tuple(spatial_counts),
    )


@dataclass(frozen=True, slots=True)
class Cubic7TimelikeVinbergSearch:
    pool: Cubic7TimelikeRootPool
    accepted_roots: tuple[Cubic7TimelikeRoot, ...]
    rejected_roots: tuple[Cubic7TimelikeRoot, ...]
    walls: tuple[Vector, ...]
    wall_gram: Matrix


def run_timelike_vinberg_search(
    seed: Cubic7TimelikeSeed, height_bound: int
) -> Cubic7TimelikeVinbergSearch:
    pool = enumerate_timelike_root_pool(seed, height_bound)
    walls = list(pool.cone_roots)
    accepted: list[Cubic7TimelikeRoot] = []
    rejected: list[Cubic7TimelikeRoot] = []
    for _height, level_iterator in groupby(pool.roots, key=lambda root: root.height):
        level = tuple(level_iterator)
        survivors = tuple(
            candidate
            for candidate in level
            if all(
                timelike_diagonal_pair(seed, candidate.vector, wall).sign(0) <= 0
                for wall in walls
            )
        )
        if any(
            timelike_diagonal_pair(seed, left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise ArithmeticError(
                "same-height Vinberg survivors are not pairwise acute; "
                "a Weyl-orbit resolution is required"
            )
        survivor_set = set(survivors)
        accepted.extend(survivors)
        rejected.extend(candidate for candidate in level if candidate not in survivor_set)
        walls.extend(root.vector for root in survivors)

    if any(
        timelike_diagonal_pair(seed, left, right).sign(0) > 0
        for index, left in enumerate(walls)
        for right in walls[index + 1 :]
    ):
        raise AssertionError("accepted walls are not pairwise acute")
    gram = tuple(
        tuple(timelike_diagonal_pair(seed, left, right) for right in walls)
        for left in walls
    )
    return Cubic7TimelikeVinbergSearch(
        pool=pool,
        accepted_roots=tuple(accepted),
        rejected_roots=tuple(rejected),
        walls=tuple(walls),
        wall_gram=gram,
    )


def _dot(left: Sequence[Cubic7], right: Sequence[Cubic7]) -> Cubic7:
    return sum((x * y for x, y in zip(left, right)), start=CUBIC7_ZERO)


def _row_rank(rows: Sequence[Sequence[Cubic7]], width: int) -> int:
    work = [list(row) for row in rows]
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
            factor = work[row][column] / pivot_value
            for index in range(column, width):
                work[row][index] -= factor * work[rank][index]
        rank += 1
        if rank == len(work):
            break
    return rank


def _kernel_ray(rows: Sequence[Sequence[Cubic7]], width: int) -> Vector | None:
    work = [list(row) for row in rows]
    pivots: list[int] = []
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
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
        pivots.append(column)
        rank += 1
        if rank == len(work):
            break
    if rank != width - 1:
        return None
    free = next(column for column in range(width) if column not in pivots)
    vector = [CUBIC7_ZERO for _ in range(width)]
    vector[free] = CUBIC7_ONE
    for row, column in enumerate(pivots):
        vector[column] = -work[row][free]
    return tuple(vector)


def _projective_key(vector: Sequence[Cubic7]) -> Vector:
    pivot = next((entry for entry in vector if entry), None)
    if pivot is None:
        raise ValueError("zero is not a projective ray")
    return tuple(entry / pivot for entry in vector)


@dataclass(frozen=True, slots=True)
class Cubic7TimelikeConeRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Cubic7
    kind: str


@dataclass(frozen=True, slots=True)
class Cubic7TimelikeConeGeometry:
    rays: tuple[Cubic7TimelikeConeRay, ...]
    finite_count: int
    cusp_count: int
    hyperideal_count: int
    irredundant: bool
    consistent_time_sheet: bool
    finite_volume: bool
    compact: bool


def _classify_rays(
    seed: Cubic7TimelikeSeed,
    functionals: Sequence[Vector],
    rays: Sequence[tuple[Vector, frozenset[int]]],
) -> Cubic7TimelikeConeGeometry:
    classified = []
    for vector, active in rays:
        norm = timelike_diagonal_pair(seed, vector, vector)
        sign = norm.sign(0)
        kind = "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal"
        classified.append(Cubic7TimelikeConeRay(vector, tuple(sorted(active)), norm, kind))
    ordered = tuple(
        sorted(classified, key=lambda ray: (_projective_key(ray.vector), ray.active_facets))
    )
    finite = sum(ray.kind == "finite" for ray in ordered)
    cusps = sum(ray.kind == "cusp" for ray in ordered)
    hyperideal = sum(ray.kind == "hyperideal" for ray in ordered)
    irredundant = bool(ordered) and all(
        _row_rank(
            tuple(ray.vector for ray in ordered if facet in ray.active_facets),
            5,
        )
        == 4
        for facet in range(len(functionals))
    )
    nonspacelike = tuple(ray for ray in ordered if ray.kind != "hyperideal")
    sheet_ok = all(
        timelike_diagonal_pair(seed, left.vector, right.vector).sign(0) <= 0
        for index, left in enumerate(nonspacelike)
        for right in nonspacelike[index + 1 :]
    )
    finite_volume = bool(ordered) and sheet_ok and hyperideal == 0
    return Cubic7TimelikeConeGeometry(
        rays=ordered,
        finite_count=finite,
        cusp_count=cusps,
        hyperideal_count=hyperideal,
        irredundant=irredundant,
        consistent_time_sheet=sheet_ok,
        finite_volume=finite_volume,
        compact=finite_volume and cusps == 0,
    )


def enumerate_timelike_cone_rays_incremental(
    search: Cubic7TimelikeVinbergSearch,
) -> Cubic7TimelikeConeGeometry:
    """Enumerate current extreme rays by exact double description."""

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
    order = (*independent, *(index for index in range(len(roots)) if index not in independent_set))
    processed = list(independent)
    rays: dict[Vector, tuple[Vector, frozenset[int]]] = {}
    for omitted in independent:
        active = tuple(index for index in independent if index != omitted)
        ray = _kernel_ray(tuple(functionals[index] for index in active), width)
        if ray is None:
            raise AssertionError("an independent seed lost rank after omitting one row")
        values = tuple(_dot(ray, functionals[index]) for index in independent)
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-entry for entry in ray)
        else:
            raise AssertionError("the independent seed cone has inconsistent orientation")
        actual = frozenset(
            index for index in independent if not _dot(oriented, functionals[index])
        )
        rays[_projective_key(oriented)] = (oriented, actual)

    rank_cache: dict[frozenset[int], int] = {}
    for new_facet in order[width:]:
        functional = functionals[new_facet]
        negative: list[tuple[Vector, frozenset[int], Cubic7]] = []
        zero: list[tuple[Vector, frozenset[int], Cubic7]] = []
        positive: list[tuple[Vector, frozenset[int], Cubic7]] = []
        for ray, active in rays.values():
            value = _dot(ray, functional)
            destination = negative if value.sign(0) < 0 else positive if value.sign(0) > 0 else zero
            destination.append((ray, active, value))
        updated: dict[Vector, tuple[Vector, frozenset[int]]] = {}
        for ray, active, _value in (*negative, *zero):
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
                    positive_value * negative_ray[index]
                    - negative_value * positive_ray[index]
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

    return _classify_rays(seed, functionals, tuple(rays.values()))


def enumerate_timelike_cone_rays_bruteforce(
    search: Cubic7TimelikeVinbergSearch,
    *,
    combination_limit: int = 1_000_000,
) -> Cubic7TimelikeConeGeometry:
    """Independent exact support enumeration, intended as a check oracle."""

    seed = search.pool.seed
    roots = search.walls
    width = 5
    support_count = len(tuple(combinations(range(len(roots)), width - 1)))
    if support_count > combination_limit:
        raise RuntimeError(
            f"brute cone check needs {support_count} supports, above {combination_limit}"
        )
    functionals = tuple(
        tuple(coefficient * coordinate for coefficient, coordinate in zip(seed.coefficients, root))
        for root in roots
    )
    rays: dict[Vector, tuple[Vector, frozenset[int]]] = {}
    for support in combinations(range(len(roots)), width - 1):
        ray = _kernel_ray(tuple(functionals[index] for index in support), width)
        if ray is None:
            continue
        values = tuple(_dot(ray, functional) for functional in functionals)
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
            oriented_values = values
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-entry for entry in ray)
            oriented_values = tuple(-value for value in values)
        else:
            continue
        active = frozenset(index for index, value in enumerate(oriented_values) if not value)
        rays[_projective_key(oriented)] = (oriented, active)
    return _classify_rays(seed, functionals, tuple(rays.values()))


__all__ = [
    "DELTA7",
    "MINKOWSKI_BOUND_DENOMINATOR",
    "MINKOWSKI_BOUND_NUMERATOR",
    "SUPPORTED_TIMELIKE_SEEDS",
    "T2_TIMELIKE_SEED",
    "T7_TIMELIKE_SEED",
    "UNIT_SIGNATURE_BITS",
    "UNIT_SIGNATURE_GENERATORS",
    "Cubic7TimelikeConeGeometry",
    "Cubic7TimelikeConeRay",
    "Cubic7TimelikeRoot",
    "Cubic7TimelikeRootPool",
    "Cubic7TimelikeSeed",
    "Cubic7TimelikeVinbergSearch",
    "enumerate_height_zero_mirrors",
    "enumerate_timelike_cone_rays_bruteforce",
    "enumerate_timelike_cone_rays_incremental",
    "enumerate_timelike_root_pool",
    "height_zero_is_b4",
    "run_timelike_vinberg_search",
    "timelike_diagonal_pair",
    "unit_signature_rank",
]
