"""Exact early Vinberg search for ``<-2phi> + D_n`` and ``D24+``.

Spatial vectors use the standard Euclidean model
``D_n={z in Z^n : sum(z_i) even}``, with dominant chamber
``y_1>=...>=y_(n-1)>=|y_n|``.  For ``D24+`` the spinor glue coset
``(1/2,...,1/2)+D24`` is included explicitly.  The first levels are bounded
simultaneously at the identity and conjugate embeddings, so the enumeration
is exact rather than a floating-point short-vector scan.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.golden import ONE, PHI, ZERO, Golden
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse


Vector: TypeAlias = tuple[Golden, ...]


def golden_dn_pair(left: Sequence[Golden], right: Sequence[Golden]) -> Golden:
    if len(left) != len(right) or len(left) < 5:
        raise ValueError("golden Dn vectors must have equal rank at least five")
    return -2 * PHI * left[0] * right[0] + sum(
        (left[index] * right[index] for index in range(1, len(left))),
        start=ZERO,
    )


@lru_cache(maxsize=None)
def golden_dn_height_zero_walls(dimension: int) -> tuple[Vector, ...]:
    if dimension < 4:
        raise ValueError("D_n requires n at least four")
    roots = []
    for index in range(dimension - 1):
        root = [ZERO] * (dimension + 1)
        root[index + 1] = -ONE
        root[index + 2] = ONE
        roots.append(tuple(root))
    final = [ZERO] * (dimension + 1)
    final[-2] = final[-1] = -ONE
    roots.append(tuple(final))
    return tuple(roots)


def _coefficient_box(target: Golden) -> tuple[Golden, ...]:
    """All coordinate values whose squares fit ``target`` at both places."""

    identity_upper = 0
    while (Golden(identity_upper * identity_upper) - target).sign(0) < 0:
        identity_upper += 1
    conjugate_upper = 0
    while (Golden(conjugate_upper * conjugate_upper) - target).sign(1) < 0:
        conjugate_upper += 1
    b_bound = (identity_upper + conjugate_upper) // 2 + 1
    a_bound = conjugate_upper + b_bound
    return tuple(
        sorted(
            Golden(a, b)
            for a in range(-a_bound, a_bound + 1)
            for b in range(-b_bound, b_bound + 1)
            if all(
                (target - Golden(a, b) ** 2).sign(place) >= 0
                for place in (0, 1)
            )
        )
    )


@dataclass(frozen=True, slots=True)
class GoldenDnRoot:
    control: Golden
    height: Golden
    vector: Vector
    spatial_coset: str


def _dominant(vector: Sequence[Golden]) -> bool:
    return all(
        (vector[index] - vector[index + 1]).sign(0) >= 0
        for index in range(len(vector) - 2)
    ) and (vector[-2] - vector[-1]).sign(0) >= 0 and (vector[-2] + vector[-1]).sign(0) >= 0


def _in_dn_root_coset(vector: Sequence[Golden]) -> bool:
    if not all(value.is_integral for value in vector):
        return False
    rational_sum = sum(value.a for value in vector)
    radical_sum = sum(value.b for value in vector)
    return rational_sum.denominator == radical_sum.denominator == 1 and rational_sum.numerator % 2 == 0 and radical_sum.numerator % 2 == 0


def _in_d24_plus_spinor_coset(vector: Sequence[Golden]) -> bool:
    """Membership in one of the three nonzero ``O_K/2`` spinor cosets.

    The same residue ``lambda in {1,phi,1+phi}`` multiplies the half-spinor
    in every Euclidean coordinate.  For each of the two Z-coefficients the
    doubled entries therefore have one common parity.  At least one parity
    is odd, and after subtracting that common half-spinor residue the D24
    translate has even coordinate sum coefficientwise.
    """

    parities = []
    translated_sums = []
    for coefficient in ("a", "b"):
        doubled = [2 * getattr(value, coefficient) for value in vector]
        if any(value.denominator != 1 for value in doubled):
            return False
        numerators = [value.numerator for value in doubled]
        parity = numerators[0] % 2
        if any(value % 2 != parity for value in numerators):
            return False
        parities.append(parity)
        translated_sums.append(sum((value - parity) // 2 for value in numerators))
    return any(parities) and all(value % 2 == 0 for value in translated_sums)


def _sum_squares_from_alphabet(
    dimension: int,
    target: Golden,
    *,
    alphabet: Sequence[Golden],
    membership,
    coset: str,
) -> tuple[tuple[Vector, str], ...]:
    # D_n dominance sorts absolute values and leaves only the sign of the
    # smallest coordinate.  Enumerating multiplicities of distinct squares
    # avoids the enormous zero-padding recursion in ranks 22--29.
    positive = tuple(
        sorted(
            {
                value if value.sign(0) >= 0 else -value
                for value in alphabet
                if value
            },
            reverse=True,
        )
    )
    squares = tuple(value * value for value in positive)

    @lru_cache(maxsize=None)
    def solve(index: int, slots: int, remainder: Golden) -> tuple[tuple[int, ...], ...]:
        if any(remainder.sign(place) < 0 for place in (0, 1)):
            return ()
        if index == len(positive):
            return ((),) if not remainder else ()
        output = []
        running = remainder
        for multiplicity in range(slots + 1):
            if all(running.sign(place) >= 0 for place in (0, 1)):
                for suffix in solve(index + 1, slots - multiplicity, running):
                    output.append((multiplicity, *suffix))
            running -= squares[index]
        return tuple(output)

    output = []
    for multiplicities in solve(0, dimension, target):
        used = sum(multiplicities)
        nonnegative = tuple(
            value
            for value, multiplicity in zip(positive, multiplicities)
            for _ in range(multiplicity)
        ) + (ZERO,) * (dimension - used)
        representatives = (nonnegative,)
        if nonnegative[-1]:
            representatives = (
                nonnegative,
                (*nonnegative[:-1], -nonnegative[-1]),
            )
        for representative in representatives:
            if _dominant(representative) and membership(representative):
                output.append((representative, coset))
    return tuple(output)


def _half_spinor_alphabet(target: Golden, residue: Golden) -> tuple[Golden, ...]:
    """All entries in ``residue/2 + O_K`` inside the square rectangle."""

    integral = _coefficient_box(4 * target)
    return tuple(
        sorted(
            value / 2
            for value in integral
            if (value - residue).is_integral
            and all(
                int(coefficient) % 2 == 0
                for coefficient in (
                    (value - residue).a,
                    (value - residue).b,
                )
            )
            and all((target - (value / 2) ** 2).sign(place) >= 0 for place in (0, 1))
        )
    )


def _sum_squares_dominant(
    dimension: int,
    target: Golden,
    *,
    include_d24_plus_glue: bool,
) -> tuple[tuple[Vector, str], ...]:
    output = list(
        _sum_squares_from_alphabet(
            dimension,
            target,
            alphabet=_coefficient_box(target),
            membership=_in_dn_root_coset,
            coset="D_n",
        )
    )
    if include_d24_plus_glue:
        if dimension != 24:
            raise ValueError("spinor glue is supported only for D24+")
        for residue in (ONE, PHI, ONE + PHI):
            output.extend(
                _sum_squares_from_alphabet(
                    dimension,
                    target,
                    alphabet=_half_spinor_alphabet(target, residue),
                    membership=_in_d24_plus_spinor_coset,
                    coset=f"spinor-{residue}",
                )
            )
    return tuple(sorted(set(output)))


@dataclass(frozen=True, slots=True)
class GoldenDnLevel:
    control: Golden
    height: Golden
    target: Golden
    candidate_count: int
    roots: tuple[GoldenDnRoot, ...]


@dataclass(frozen=True, slots=True)
class GoldenDnSearch:
    dimension: int
    d24_plus: bool
    controls: tuple[Golden, ...]
    levels: tuple[GoldenDnLevel, ...]
    accepted_roots: tuple[GoldenDnRoot, ...]
    walls: tuple[Vector, ...]


@lru_cache(maxsize=None)
def run_golden_dn_search(
    dimension: int,
    controls: tuple[Golden, ...] = (ONE, PHI),
    *,
    d24_plus: bool = False,
) -> GoldenDnSearch:
    if d24_plus and dimension != 24:
        raise ValueError("D24+ glue exists only in dimension 24")
    walls = list(golden_dn_height_zero_walls(dimension))
    accepted = []
    levels = []
    for control in controls:
        target = 2 + 2 * PHI * control * control
        candidates = tuple(
            GoldenDnRoot(control, PHI * control * control, (control, *spatial), coset)
            for spatial, coset in _sum_squares_dominant(
                dimension, target, include_d24_plus_glue=d24_plus
            )
        )
        survivors = tuple(
            root
            for root in candidates
            if all(golden_dn_pair(root.vector, wall).sign(0) <= 0 for wall in walls)
        )
        if any(
            golden_dn_pair(left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-control golden Dn survivors are not acute")
        accepted.extend(survivors)
        walls.extend(root.vector for root in survivors)
        levels.append(
            GoldenDnLevel(
                control,
                PHI * control * control,
                target,
                len(candidates),
                survivors,
            )
        )
    return GoldenDnSearch(
        dimension, d24_plus, controls, tuple(levels), tuple(accepted), tuple(walls)
    )


@dataclass(frozen=True, slots=True)
class GoldenDnFirstCapGeometry:
    extreme_ray_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int


@dataclass(frozen=True, slots=True)
class GoldenDnExtremeRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Golden
    kind: str


@dataclass(frozen=True, slots=True)
class GoldenDnGeometry:
    extreme_rays: tuple[GoldenDnExtremeRay, ...]
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    irredundant: bool
    finite_volume: bool
    compact: bool


def _projective_key(vector: Sequence[Golden]) -> Vector:
    pivot = next(value for value in vector if value)
    return tuple(value / pivot for value in vector)


def _functional_value(functional: Sequence[Golden], vector: Sequence[Golden]) -> Golden:
    return sum(
        (
            coefficient * coordinate
            for coefficient, coordinate in zip(functional, vector)
            if coefficient and coordinate
        ),
        start=ZERO,
    )


def _row_rank(rows: Sequence[Sequence[Golden]], width: int) -> int:
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


def _kernel(rows: Sequence[Sequence[Golden]]) -> Vector | None:
    if not rows:
        return None
    width = len(rows[0])
    work = [list(row) for row in rows]
    pivots = []
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        value = work[rank][column]
        work[rank] = [entry / value for entry in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[rank][index]
                for index in range(width)
            ]
        pivots.append(column)
        rank += 1
    if width - rank != 1:
        return None
    free = next(index for index in range(width) if index not in pivots)
    vector = [ZERO] * width
    vector[free] = ONE
    for row, column in enumerate(pivots):
        vector[column] = -work[row][free]
    return tuple(vector)


def analyze_golden_dn_cone(search: GoldenDnSearch) -> GoldenDnGeometry:
    """Enumerate all rays of the current wall cone by exact DD."""

    dimension = search.dimension
    ambient_rank = dimension + 1
    walls = search.walls
    if len(walls) < ambient_rank:
        raise ValueError("a pointed cone needs at least dimension+1 walls")
    functionals = tuple(
        (-2 * PHI * wall[0], *wall[1:]) for wall in walls
    )
    initial = tuple(range(ambient_rank))
    if _row_rank(tuple(functionals[index] for index in initial), ambient_rank) != ambient_rank:
        raise ValueError("height-zero Dn walls and first cap are not independent")
    processed = list(initial)
    rays: dict[Vector, tuple[Vector, frozenset[int]]] = {}
    for omitted in initial:
        support = tuple(index for index in initial if index != omitted)
        ray = _kernel(tuple(functionals[index] for index in support))
        if ray is None:
            raise AssertionError("independent seed lost a ray")
        values = tuple(_functional_value(functionals[index], ray) for index in initial)
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-value for value in ray)
        else:
            raise AssertionError("seed ray has inconsistent orientation")
        active = frozenset(
            index
            for index in initial
            if not _functional_value(functionals[index], oriented)
        )
        rays[_projective_key(oriented)] = (oriented, active)

    rank_cache: dict[frozenset[int], int] = {}
    for new_facet in range(ambient_rank, len(walls)):
        functional = functionals[new_facet]
        negative = []
        zero = []
        positive = []
        for key, (ray, active) in rays.items():
            value = _functional_value(functional, ray)
            destination = negative if value.sign(0) < 0 else positive if value.sign(0) > 0 else zero
            destination.append((key, ray, active, value))
        updated = {
            key: (ray, active)
            for key, ray, active, _value in negative
        }
        updated.update(
            {
                key: (ray, active | frozenset((new_facet,)))
                for key, ray, active, _value in zero
            }
        )
        for _negative_key, negative_ray, negative_active, negative_value in negative:
            for _positive_key, positive_ray, positive_active, positive_value in positive:
                common = negative_active & positive_active
                if len(common) < ambient_rank - 2:
                    continue
                rank = rank_cache.get(common)
                if rank is None:
                    rank = _row_rank(
                        tuple(functionals[index] for index in common), ambient_rank
                    )
                    rank_cache[common] = rank
                if rank != ambient_rank - 2:
                    continue
                ray = tuple(
                    positive_value * negative_ray[index]
                    - negative_value * positive_ray[index]
                    for index in range(ambient_rank)
                )
                active = frozenset(
                    index
                    for index in (*processed, new_facet)
                    if not _functional_value(functionals[index], ray)
                )
                updated[_projective_key(ray)] = (ray, active)
        rays = updated
        processed.append(new_facet)

    output = []
    for ray, active in rays.values():
        if _row_rank(tuple(functionals[index] for index in active), ambient_rank) != dimension:
            raise AssertionError("double description returned a non-extreme ray")
        norm = golden_dn_pair(ray, ray)
        sign = norm.sign(0)
        output.append(
            GoldenDnExtremeRay(
                ray,
                tuple(sorted(active)),
                norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
            )
        )
    ordered = tuple(sorted(output, key=lambda ray: _projective_key(ray.vector)))
    finite = sum(ray.kind == "finite" for ray in ordered)
    cusps = sum(ray.kind == "cusp" for ray in ordered)
    hyperideal = sum(ray.kind == "hyperideal" for ray in ordered)
    irredundant = bool(ordered) and all(
        any(_functional_value(functionals[facet], ray.vector).sign(0) < 0 for ray in ordered)
        for facet in range(len(functionals))
    )
    finite_volume = irredundant and hyperideal == 0
    return GoldenDnGeometry(
        ordered,
        finite,
        cusps,
        hyperideal,
        irredundant,
        finite_volume,
        finite_volume and cusps == 0,
    )


def golden_dn_first_cap_geometry(dimension: int) -> GoldenDnFirstCapGeometry:
    """Classify the vertices of the affine D_n highest-root alcove exactly."""

    squares = [Golden(0), Golden(4)]
    squares.extend(Golden(index) for index in range(2, dimension - 1))
    squares.extend((Golden(dimension), Golden(dimension)))
    if len(squares) != dimension + 1:
        raise AssertionError("D_n alcove vertex census has wrong size")
    norms = tuple(square - 2 * PHI for square in squares)
    return GoldenDnFirstCapGeometry(
        len(norms),
        sum(norm.sign(0) < 0 for norm in norms),
        sum(norm.sign(0) == 0 for norm in norms),
        sum(norm.sign(0) > 0 for norm in norms),
    )


__all__ = [
    "GoldenDnFirstCapGeometry",
    "GoldenDnExtremeRay",
    "GoldenDnGeometry",
    "GoldenDnLevel",
    "GoldenDnRoot",
    "GoldenDnSearch",
    "analyze_golden_dn_cone",
    "golden_dn_first_cap_geometry",
    "golden_dn_height_zero_walls",
    "golden_dn_pair",
    "run_golden_dn_search",
]
