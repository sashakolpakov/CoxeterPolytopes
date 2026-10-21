"""First exact Vinberg packet for ``<-2*theta> + D_n`` and ``D24+``.

The coordinate model is

``D_n = {y in O_K^n : sum(y_i) in 2 O_K}``.

For ``D24+`` we adjoin the spinor coset

``D24 + (1/2,...,1/2)``.

Since ``D24+`` is even unimodular, its primitive mirror-norm ideals for the
timelike coefficient ``2*theta`` are exactly ``(1)`` or ``(2)``: away from
the inert prime two, unimodularity and primitivity make the pairing ideal a
unit; if the spatial vector is locally divisible by two, primitive control
coordinate gives ``q=-2*theta*x^2 mod 4`` and valuation one.  The usual
cubic-seven class-number and unit-signature argument supplies canonical
element representatives 1 and 2.

At the first positive control level through height ten, exact parity-class
square partitions give a single dominant root for every D_n, 22<=n<=29,
and also for D24+: ``(x;s,s,0,...)`` with
``x=-1+theta+theta^2`` and ``s=1+theta``.  A spinor-coset coordinate must be
half-integral in *all* 24 orthonormal coordinates with one common nonzero
residue; zero padding is forbidden.  The exhaustive fixed-alphabet coset
census finds no spinor root on this level.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Iterable

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    THETA7,
    Cubic7,
)
from coxeter_polytopes.search.cubic7_bna1 import (
    _functional_value,
    _kernel,
    _row_rank,
    _square_partition_solutions,
)
from coxeter_polytopes.search.cubic7_root_blocks import _bounded_for_caps
from coxeter_polytopes.search.silver_d24_plus import (
    _integer_determinant,
    d24_plus_basis_euclidean,
    d24_plus_gram,
)


Vector = tuple[Cubic7, ...]
TIMELIKE = 2 * THETA7
FIRST_CONTROL = -1 + THETA7 + THETA7**2
FIRST_SCALE = 1 + THETA7
FIRST_HEIGHT = 2 * THETA7 + THETA7**2


COMPLETE_D24_PLUS_NORM_PROOF = (
    "D24+ is even unimodular. For primitive r, q(r)|2I(r), with "
    "I(r)=(-2*theta*x,y). Away from inert (2), I(r)=O. At (2), either "
    "I(r)=O or y lies in 2L and primitive x is a unit, giving "
    "q=-2*theta*x^2 mod 4 and v_2(q)=1. Hence (q)=(1) or (2). "
    "Class number one and full unit-signature rank make every totally "
    "positive unit a square, yielding canonical elements 1 and 2."
)


def dn_pair(left: Vector, right: Vector) -> Cubic7:
    if len(left) != len(right) or len(left) < 3:
        raise ValueError("D_n vectors have equal rank n+1")
    return -TIMELIKE * left[0] * right[0] + sum(
        (left[i] * right[i] for i in range(1, len(left))),
        start=CUBIC7_ZERO,
    )


def dn_height_zero_simple_roots(dimension: int) -> tuple[Vector, ...]:
    if dimension < 4:
        raise ValueError("D_n requires n>=4")
    roots: list[Vector] = []
    for index in range(dimension - 1):
        root = [CUBIC7_ZERO] * (dimension + 1)
        root[1 + index] = -CUBIC7_ONE
        root[2 + index] = CUBIC7_ONE
        roots.append(tuple(root))
    fork = [CUBIC7_ZERO] * (dimension + 1)
    fork[dimension - 1] = -CUBIC7_ONE
    fork[dimension] = -CUBIC7_ONE
    roots.append(tuple(fork))
    return tuple(roots)


def _coefficient_parity(value: Cubic7) -> tuple[int, int, int] | None:
    if not value.is_integral:
        return None
    return tuple(int(coefficient) % 2 for coefficient in (value.a, value.b, value.c))


def _sum_in_two_o(values: Iterable[Cubic7]) -> bool:
    total = sum(values, start=CUBIC7_ZERO)
    return total.is_integral and all(
        int(coefficient) % 2 == 0 for coefficient in (total.a, total.b, total.c)
    )


def d24_plus_membership(coordinates: Vector) -> str | None:
    """Return ``D``, ``spinor`` or ``None`` in orthonormal coordinates."""

    if len(coordinates) != 24:
        raise ValueError("D24+ has 24 orthonormal coordinates")
    if all(value.is_integral for value in coordinates):
        return "D" if _sum_in_two_o(coordinates) else None
    doubled = tuple(2 * value for value in coordinates)
    parities = tuple(_coefficient_parity(value) for value in doubled)
    first = parities[0]
    if (
        first is None
        or first == (0, 0, 0)
        or any(parity != first for parity in parities)
    ):
        return None
    return "spinor" if _sum_in_two_o(coordinates) else None


def d24_plus_integral_basis_coordinates(coordinates: Vector) -> Vector:
    """Coordinates in ``alpha_2,...,alpha_24,spinor`` over ``O_K``."""

    if len(coordinates) != 24:
        raise ValueError("D24+ has 24 orthonormal coordinates")
    basis = d24_plus_basis_euclidean()
    matrix = tuple(
        tuple(basis[column][row] for column in range(24))
        for row in range(24)
    )
    from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse

    inverse = _fraction_inverse(matrix)
    return tuple(
        sum(
            (Cubic7(inverse[i][j]) * coordinates[j] for j in range(24)),
            start=CUBIC7_ZERO,
        )
        for i in range(24)
    )


def d24_plus_basis_audit() -> bool:
    gram = d24_plus_gram()
    return (
        _integer_determinant(gram) == 1
        and all(entry.denominator == 1 for row in gram for entry in row)
        and all(gram[index][index].numerator % 2 == 0 for index in range(24))
    )


def _fixed_alphabet_square_partitions(
    values: tuple[Cubic7, ...], coordinate_count: int, target: Cubic7
) -> tuple[Vector, ...]:
    """Use exactly ``coordinate_count`` entries from ``values``.

    Unlike the B_n/D_n root-lattice helper, this never pads unused positions
    by zero.  That distinction is essential in a nonzero D24+ spinor coset.
    """

    ordered = tuple(sorted(values, reverse=True))
    squares = tuple(value * value for value in ordered)

    @lru_cache(maxsize=None)
    def solve(
        index: int, slots: int, remainder: Cubic7
    ) -> tuple[tuple[int, ...], ...]:
        if any(remainder.sign(place) < 0 for place in range(3)):
            return ()
        if index == len(ordered):
            return ((),) if slots == 0 and not remainder else ()
        square = squares[index]
        output: list[tuple[int, ...]] = []
        running = remainder
        for multiplicity in range(slots + 1):
            if any(running.sign(place) < 0 for place in range(3)):
                break
            for suffix in solve(index + 1, slots - multiplicity, running):
                output.append((multiplicity, *suffix))
            running -= square
        return tuple(output)

    return tuple(
        tuple(
            value
            for value, multiplicity in zip(ordered, multiplicities)
            for _ in range(multiplicity)
        )
        for multiplicities in solve(0, coordinate_count, target)
    )


@dataclass(frozen=True, slots=True)
class DominantDnPattern:
    coset: str
    coordinates: Vector


@lru_cache(maxsize=None)
def first_level_dominant_patterns(
    dimension: int, *, include_spinor_coset: bool
) -> tuple[DominantDnPattern, ...]:
    """Exhaust the first norm-two shell modulo signed permutations."""

    target = 2 + TIMELIKE * FIRST_CONTROL**2
    patterns: list[DominantDnPattern] = []

    integral_values = tuple(
        value
        for value in _bounded_for_caps((target, target, target))
        if value.sign(0) >= 0
        and all((target - value * value).sign(place) >= 0 for place in range(3))
    )
    for coordinates in _square_partition_solutions(integral_values, dimension, target):
        if _sum_in_two_o(coordinates):
            patterns.append(DominantDnPattern("D", coordinates))

    if include_spinor_coset:
        # Write a half-integral coordinate as z/2.  Every z has the same
        # nonzero coefficient-parity vector, and the coordinate sum lies in
        # 2O.  Iterating all seven classes is finite and exact.
        doubled = _bounded_for_caps((4 * target, 4 * target, 4 * target))
        for parity_code in range(1, 8):
            parity = tuple((parity_code >> bit) & 1 for bit in range(3))
            values = tuple(
                value / 2
                for value in doubled
                if _coefficient_parity(value) == parity
                and (value / 2).sign(0) >= 0
                and all(
                    (target - (value / 2) ** 2).sign(place) >= 0
                    for place in range(3)
                )
            )
            for coordinates in _fixed_alphabet_square_partitions(
                values, dimension, target
            ):
                if d24_plus_membership(coordinates) == "spinor":
                    patterns.append(DominantDnPattern(f"spinor-parity-{parity_code}", coordinates))

    return tuple(sorted(set(patterns), key=lambda pattern: (pattern.coset, pattern.coordinates)))


@dataclass(frozen=True, slots=True)
class DnExtremeRay:
    omitted_facet: int
    vector: Vector
    active_facets: tuple[int, ...]
    active_rank_verified: bool
    norm: Cubic7
    kind: str


@dataclass(frozen=True, slots=True)
class DnFirstPacket:
    dimension: int
    lattice: str
    canonical_mirror_norms: tuple[Cubic7, ...]
    complete_norm_class_proof: str | None
    patterns: tuple[DominantDnPattern, ...]
    walls: tuple[Vector, ...]
    wall_rank: int
    rays: tuple[DnExtremeRay, ...]
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    irredundant: bool
    finite_volume: bool
    compact: bool


@lru_cache(maxsize=None)
def run_dn_first_packet(
    dimension: int, *, include_spinor_coset: bool = False
) -> DnFirstPacket:
    if include_spinor_coset and dimension != 24:
        raise ValueError("the supported even unimodular spinor extension is D24+")
    if include_spinor_coset and not d24_plus_basis_audit():
        raise AssertionError("glue-aware D24+ basis is not even unimodular")
    patterns = first_level_dominant_patterns(
        dimension, include_spinor_coset=include_spinor_coset
    )
    if len(patterns) != 1:
        raise AssertionError(f"first D_n packet has {len(patterns)} dominant patterns")
    positive_roots = tuple((FIRST_CONTROL, *pattern.coordinates) for pattern in patterns)
    if any(dn_pair(root, root) != 2 for root in positive_roots):
        raise AssertionError("first D_n root has wrong norm")
    if include_spinor_coset:
        gram = d24_plus_gram()
        for root in positive_roots:
            basis_coordinates = d24_plus_integral_basis_coordinates(root[1:])
            if not all(value.is_integral for value in basis_coordinates):
                raise AssertionError("D24+ candidate is nonintegral in the glue basis")
            pairings = tuple(
                sum(
                    (
                        Cubic7(gram[i][j]) * basis_coordinates[j]
                        for j in range(24)
                        if gram[i][j]
                    ),
                    start=CUBIC7_ZERO,
                )
                for i in range(24)
            )
            if not all(value.is_integral for value in pairings):
                raise AssertionError("D24+ candidate reflection is nonintegral")
    if any(
        dn_pair(left, right).sign(0) > 0
        for index, left in enumerate(positive_roots)
        for right in positive_roots[index + 1 :]
    ):
        raise AssertionError("same-height D24+ glue packet is not mutually acute")
    walls = (*dn_height_zero_simple_roots(dimension), *positive_roots)
    functionals = tuple(
        (-TIMELIKE * wall[0], *wall[1:]) for wall in walls
    )
    wall_rank = _row_rank(functionals, dimension + 1)
    if wall_rank != dimension + 1:
        raise AssertionError("affine D_n wall system does not span")

    # The D_n simple roots plus the root-lattice affine root are an
    # independent simplicial seed.  Compute its n+1 rays once, then insert
    # the single D24+ glue half-space by exact double description.  Every
    # pair of rays of a simplicial cone is adjacent, so no support census is
    # needed for the one-step refinement.
    base_count = dimension + 1
    base_rays: list[tuple[Vector, frozenset[int]]] = []
    for omitted in range(base_count):
        support = tuple(index for index in range(base_count) if index != omitted)
        vector = _kernel(tuple(functionals[index] for index in support))
        if vector is None:
            raise AssertionError("affine D_n simplex lost a seed ray")
        values = tuple(
            _functional_value(functionals[index], vector)
            for index in range(base_count)
        )
        if all(value.sign(0) <= 0 for value in values):
            oriented = vector
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-entry for entry in vector)
        else:
            raise AssertionError("affine D_n simplex has inconsistent orientation")
        active = frozenset(
            index
            for index in range(base_count)
            if not _functional_value(functionals[index], oriented)
        )
        base_rays.append((oriented, active))

    refined = base_rays
    for new_facet in range(base_count, len(walls)):
        functional = functionals[new_facet]
        negative = []
        zero = []
        positive = []
        for ray, active in refined:
            value = _functional_value(functional, ray)
            (negative if value.sign(0) < 0 else positive if value.sign(0) > 0 else zero).append(
                (ray, active, value)
            )
        updated: dict[Vector, tuple[Vector, frozenset[int]]] = {}

        def projective_key(vector: Vector) -> Vector:
            pivot = next(entry for entry in vector if entry)
            return tuple(entry / pivot for entry in vector)

        for ray, active, _value in negative:
            updated[projective_key(ray)] = (ray, active)
        for ray, active, _value in zero:
            updated[projective_key(ray)] = (ray, active | frozenset((new_facet,)))
        for negative_ray, negative_active, negative_value in negative:
            for positive_ray, positive_active, positive_value in positive:
                intersection = tuple(
                    positive_value * negative_ray[index]
                    - negative_value * positive_ray[index]
                    for index in range(dimension + 1)
                )
                active = frozenset(
                    index
                    for index in range(new_facet + 1)
                    if not _functional_value(functionals[index], intersection)
                )
                updated[projective_key(intersection)] = (intersection, active)
        refined = list(updated.values())

    ray_map: dict[Vector, DnExtremeRay] = {}
    for oriented, _stored_active in refined:
        actual_active = tuple(
            index
            for index, functional in enumerate(functionals)
            if not _functional_value(functional, oriented)
        )
        active_rank_verified = (
            _row_rank(
                tuple(functionals[index] for index in actual_active), dimension + 1
            )
            == dimension
        )
        norm = dn_pair(oriented, oriented)
        sign = norm.sign(0)
        omitted = next(
            (index for index in range(len(walls)) if index not in actual_active), -1
        )
        candidate = DnExtremeRay(
                omitted,
                oriented,
                actual_active,
                active_rank_verified,
                norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
            )
        key_pivot = next(entry for entry in oriented if entry)
        key = tuple(entry / key_pivot for entry in oriented)
        ray_map.setdefault(key, candidate)
    ordered = tuple(sorted(ray_map.values(), key=lambda ray: ray.vector))
    if not ordered or not all(ray.active_rank_verified for ray in ordered):
        raise AssertionError("D_n support enumeration failed extremality")
    finite = sum(ray.kind == "finite" for ray in ordered)
    cusps = sum(ray.kind == "cusp" for ray in ordered)
    hyperideal = sum(ray.kind == "hyperideal" for ray in ordered)
    return DnFirstPacket(
        dimension,
        "D24+" if include_spinor_coset else f"D{dimension}",
        (Cubic7(1), Cubic7(2)),
        COMPLETE_D24_PLUS_NORM_PROOF if include_spinor_coset else None,
        patterns,
        walls,
        wall_rank,
        ordered,
        finite,
        cusps,
        hyperideal,
        True,
        hyperideal == 0,
        hyperideal == 0 and cusps == 0,
    )


def _cubic(value: Cubic7) -> dict[str, str]:
    return value.to_json()


def _vector(vector: Iterable[Cubic7]) -> list[dict[str, str]]:
    return [_cubic(value) for value in vector]


def build_dn_first_packet_transcript(
    dimension: int, *, include_spinor_coset: bool = False
) -> dict[str, Any]:
    result = run_dn_first_packet(
        dimension, include_spinor_coset=include_spinor_coset
    )
    return {
        "schema_version": "0.1.0",
        "transcript_kind": "cubic7_dn_first_vinberg_packet",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "lattice": f"<-2theta> + {result.lattice}",
        "hyperbolic_dimension": dimension,
        "height_bound": 10,
        "canonical_mirror_norms": [_cubic(norm) for norm in result.canonical_mirror_norms],
        "complete_norm_class_proof": result.complete_norm_class_proof,
        "scope_warning": None
        if result.complete_norm_class_proof
        else (
            "The D_n root-lattice entry is a first-shell screen only. Its "
            "nonunimodular discriminant-two local mirror-norm audit is not "
            "yet complete and no global completeness claim is made."
        ),
        "first_control": _cubic(FIRST_CONTROL),
        "first_height": _cubic(FIRST_HEIGHT),
        "dominant_patterns": [
            {"coset": pattern.coset, "coordinates": _vector(pattern.coordinates)}
            for pattern in result.patterns
        ],
        "wall_count": len(result.walls),
        "wall_rank": result.wall_rank,
        "geometry": {
            "engine": "exact_simplicial_omit_one_facet",
            "extreme_ray_count": len(result.rays),
            "finite_vertex_count": result.finite_vertex_count,
            "cusp_count": result.cusp_count,
            "hyperideal_ray_count": result.hyperideal_ray_count,
            "irredundant": result.irredundant,
            "finite_volume": result.finite_volume,
            "compact": result.compact,
            "rays": [
                {
                    "omitted_facet": ray.omitted_facet,
                    "vector": _vector(ray.vector),
                    "active_facets": list(ray.active_facets),
                    "active_rank_verified": ray.active_rank_verified,
                    "norm": _cubic(ray.norm),
                    "kind": ray.kind,
                }
                for ray in result.rays
            ],
        },
    }


__all__ = [
    "COMPLETE_D24_PLUS_NORM_PROOF",
    "FIRST_CONTROL",
    "FIRST_HEIGHT",
    "FIRST_SCALE",
    "build_dn_first_packet_transcript",
    "d24_plus_basis_audit",
    "d24_plus_integral_basis_coordinates",
    "d24_plus_membership",
    "dn_height_zero_simple_roots",
    "dn_pair",
    "first_level_dominant_patterns",
    "run_dn_first_packet",
]
