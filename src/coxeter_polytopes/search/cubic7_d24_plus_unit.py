"""Glue-complete H=100 prefix for ``<-(theta-1)> + D24+``.

The positive lattice is the even unimodular Niemeier lattice with root
system D24, in the independently audited glue basis from ``cubic7_dn``.
The first nonempty positive shell through H=100 is the norm-one affine
highest-root wall.  Its cone is the D24 alcove, so all 25 projective rays are
known exactly from the fundamental coweights and affine marks.  Later shells
are enumerated in orthonormal signed-permutation coordinates, including all
seven nonzero ``O_K/2O_K`` spinor residues, and inserted by rank-certified
exact double description.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import groupby
from typing import Iterable, Sequence, TypeAlias

from coxeter_polytopes.exact.cubic7 import CUBIC7_ONE, CUBIC7_ZERO, THETA7, Cubic7
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.cubic7_a1_golay import (
    A_UNIT,
    COMPLETE_MIRROR_NORM_PROOF,
    UNIT_A1_GOLAY_SEED,
    enumerate_control_levels,
)
from coxeter_polytopes.search.cubic7_dn import (
    _coefficient_parity,
    _fixed_alphabet_square_partitions,
    d24_plus_basis_audit,
    d24_plus_integral_basis_coordinates,
    d24_plus_membership,
)
from coxeter_polytopes.search.cubic7_bna1 import (
    _functional_value,
    _kernel,
    _projective_key,
    _row_rank,
    _square_partition_solutions,
)
from coxeter_polytopes.search.cubic7_root_blocks import _bounded_for_caps
from coxeter_polytopes.search.silver_d24_plus import d24_cartan_gram, d24_plus_gram


Vector: TypeAlias = tuple[Cubic7, ...]
D24_MARKS = (1, *(2 for _ in range(21)), 1, 1)
FIRST_CONTROL = 2 + 5 * THETA7 + 2 * THETA7**2
FIRST_SPATIAL = 2 * THETA7 + THETA7**2
FIRST_HEIGHT = 5 + 14 * THETA7 + 6 * THETA7**2


@dataclass(frozen=True, slots=True)
class D24UnitAlcoveRay:
    node: int
    mark: int
    inverse_cartan_diagonal: int
    norm: Cubic7
    kind: str


@dataclass(frozen=True, slots=True)
class D24PlusUnitFirstPacket:
    complete_control_height_bound: int
    control_level_count: int
    first_control: Cubic7
    first_height: Cubic7
    root_norm: Cubic7
    spatial_orthonormal: Vector
    spatial_coset: str | None
    glue_basis_coordinates: Vector
    glue_basis_integral: bool
    reflection_integral: bool
    affine_marks_verified: bool
    cap: Cubic7
    rays: tuple[D24UnitAlcoveRay, ...]
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    finite_volume: bool
    compact: bool


@lru_cache(maxsize=1)
def unit_d24_plus_first_packet() -> D24PlusUnitFirstPacket:
    if not d24_plus_basis_audit():
        raise AssertionError("the D24+ glue basis lost even unimodularity")
    levels = enumerate_control_levels(UNIT_A1_GOLAY_SEED, 100)
    level = next(
        item
        for item in levels
        if item.norm == 1 and item.control == FIRST_CONTROL
    )
    if level.height != FIRST_HEIGHT:
        raise AssertionError("the unit first height changed")

    spatial = (FIRST_SPATIAL, FIRST_SPATIAL, *(CUBIC7_ZERO for _ in range(22)))
    norm = -A_UNIT * FIRST_CONTROL**2 + sum(
        (value * value for value in spatial), start=CUBIC7_ZERO
    )
    if norm != 1 or d24_plus_membership(spatial) != "D":
        raise AssertionError("the affine candidate is not a norm-one D24 root")
    glue_coordinates = d24_plus_integral_basis_coordinates(spatial)
    reflection_integral = (
        (2 * A_UNIT * FIRST_CONTROL).is_integral
        and all((2 * value).is_integral for value in glue_coordinates)
    )

    roots = []
    for index in range(23):
        vector = [0] * 24
        vector[index] = 1
        vector[index + 1] = -1
        roots.append(tuple(vector))
    final = [0] * 24
    final[-2] = final[-1] = 1
    roots.append(tuple(final))
    reconstructed = tuple(
        sum(mark * roots[index][coordinate] for index, mark in enumerate(D24_MARKS))
        for coordinate in range(24)
    )
    marks_verified = reconstructed == (1, 1, *(0 for _ in range(22)))
    if not marks_verified:
        raise AssertionError("the D24 affine marks changed")

    cap = A_UNIT * FIRST_CONTROL / FIRST_SPATIAL
    inverse = _fraction_inverse(d24_cartan_gram())
    rays = [D24UnitAlcoveRay(0, 1, 0, -A_UNIT, "finite")]
    for index, mark in enumerate(D24_MARKS):
        diagonal = inverse[index][index]
        if diagonal.denominator != 1:
            raise AssertionError("the D24 inverse diagonal ceased to be integral")
        ray_norm = -A_UNIT + cap * cap * Cubic7(diagonal) / (mark * mark)
        sign = ray_norm.sign(0)
        rays.append(
            D24UnitAlcoveRay(
                index + 1,
                mark,
                diagonal.numerator,
                ray_norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
            )
        )
    ordered = tuple(rays)
    finite = sum(ray.kind == "finite" for ray in ordered)
    cusps = sum(ray.kind == "cusp" for ray in ordered)
    hyperideal = sum(ray.kind == "hyperideal" for ray in ordered)
    return D24PlusUnitFirstPacket(
        100,
        len(levels),
        FIRST_CONTROL,
        FIRST_HEIGHT,
        norm,
        spatial,
        d24_plus_membership(spatial),
        glue_coordinates,
        all(value.is_integral for value in glue_coordinates),
        reflection_integral,
        marks_verified,
        cap,
        ordered,
        finite,
        cusps,
        hyperideal,
        hyperideal == 0,
        hyperideal == 0 and cusps == 0,
    )


def unit_d24_pair(left: Sequence[Cubic7], right: Sequence[Cubic7]) -> Cubic7:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("a D24+ host vector has rank 25")
    return -A_UNIT * left[0] * right[0] + sum(
        (left[index] * right[index] for index in range(1, 25)),
        start=CUBIC7_ZERO,
    )


@lru_cache(maxsize=1)
def d24_height_zero_simple_roots_orthonormal() -> tuple[Vector, ...]:
    roots = []
    for index in range(23):
        root = [CUBIC7_ZERO] * 25
        root[1 + index] = -CUBIC7_ONE
        root[2 + index] = CUBIC7_ONE
        roots.append(tuple(root))
    fork = [CUBIC7_ZERO] * 25
    fork[23] = fork[24] = -CUBIC7_ONE
    roots.append(tuple(fork))
    return tuple(roots)


def _dominant_spatial_representatives(coordinates: Vector) -> tuple[Vector, ...]:
    if coordinates[-1]:
        return coordinates, (*coordinates[:-1], -coordinates[-1])
    return (coordinates,)


@dataclass(frozen=True, slots=True)
class Cubic7D24DominantPattern:
    coordinates: Vector
    coset: str


@lru_cache(maxsize=None)
def enumerate_dominant_d24_plus_spatial(
    target: Cubic7,
) -> tuple[Cubic7D24DominantPattern, ...]:
    """Exhaust D24 and all seven spinor residue classes at one square."""

    output: set[Cubic7D24DominantPattern] = set()
    integral_values = tuple(
        value
        for value in _bounded_for_caps((target, target, target))
        if value.sign(0) >= 0
        and all((target - value * value).sign(place) >= 0 for place in range(3))
    )
    for nonnegative in _square_partition_solutions(integral_values, 24, target):
        for representative in _dominant_spatial_representatives(nonnegative):
            if d24_plus_membership(representative) == "D":
                output.add(Cubic7D24DominantPattern(representative, "D"))

    doubled = _bounded_for_caps((4 * target, 4 * target, 4 * target))
    for parity_code in range(1, 8):
        parity = tuple((parity_code >> bit) & 1 for bit in range(3))
        alphabet = tuple(
            value / 2
            for value in doubled
            if _coefficient_parity(value) == parity
            and (value / 2).sign(0) >= 0
            and all(
                (target - (value / 2) ** 2).sign(place) >= 0
                for place in range(3)
            )
        )
        for nonnegative in _fixed_alphabet_square_partitions(alphabet, 24, target):
            for representative in _dominant_spatial_representatives(nonnegative):
                if d24_plus_membership(representative) == "spinor":
                    output.add(
                        Cubic7D24DominantPattern(
                            representative, f"spinor-parity-{parity_code}"
                        )
                    )
    return tuple(sorted(output, key=lambda item: (item.coordinates, item.coset)))


def _reflection_integral(
    norm: Cubic7, control: Cubic7, spatial: Vector
) -> bool:
    if not (2 * A_UNIT * control / norm).is_integral:
        return False
    coordinates = d24_plus_integral_basis_coordinates(spatial)
    if not all(value.is_integral for value in coordinates):
        return False
    gram = d24_plus_gram()
    pairings = tuple(
        sum(
            (
                Cubic7(gram[index][other]) * coordinates[other]
                for other in range(24)
                if gram[index][other]
            ),
            start=CUBIC7_ZERO,
        )
        for index in range(24)
    )
    return all((2 * value / norm).is_integral for value in pairings)


@dataclass(frozen=True, slots=True)
class Cubic7D24Candidate:
    norm: Cubic7
    height: Cubic7
    control: Cubic7
    target: Cubic7
    coset: str
    vector: Vector


@lru_cache(maxsize=None)
def enumerate_unit_d24_control_level(level) -> tuple[Cubic7D24Candidate, ...]:
    target = 2 * level.coordinate_square_target
    roots = []
    for pattern in enumerate_dominant_d24_plus_spatial(target):
        vector = (level.control, *pattern.coordinates)
        if unit_d24_pair(vector, vector) != level.norm:
            raise AssertionError("the D24 square partition has the wrong norm")
        if not _reflection_integral(level.norm, level.control, pattern.coordinates):
            raise AssertionError("a canonical D24+ mirror failed reflection integrality")
        roots.append(
            Cubic7D24Candidate(
                level.norm,
                level.height,
                level.control,
                target,
                pattern.coset,
                vector,
            )
        )
    return tuple(sorted(roots, key=lambda root: (root.norm, root.vector, root.coset)))


@dataclass(frozen=True, slots=True)
class Cubic7D24Packet:
    height: Cubic7
    candidate_count: int
    accepted_roots: tuple[Cubic7D24Candidate, ...]


@dataclass(frozen=True, slots=True)
class Cubic7D24Search:
    height_bound: int
    control_level_count: int
    packets: tuple[Cubic7D24Packet, ...]
    accepted_roots: tuple[Cubic7D24Candidate, ...]
    walls: tuple[Vector, ...]


@lru_cache(maxsize=None)
def run_unit_d24_plus_prefix(height_bound: int = 100) -> Cubic7D24Search:
    """Run simultaneous-height Vinberg selection through ``height_bound``."""

    if not d24_plus_basis_audit():
        raise AssertionError("the D24+ glue basis lost even unimodularity")
    levels = enumerate_control_levels(UNIT_A1_GOLAY_SEED, height_bound)
    walls = list(d24_height_zero_simple_roots_orthonormal())
    accepted: list[Cubic7D24Candidate] = []
    packets = []
    for height, iterator in groupby(levels, key=lambda level: level.height):
        packet_levels = tuple(iterator)
        candidates = tuple(
            root
            for level in packet_levels
            for root in enumerate_unit_d24_control_level(level)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(unit_d24_pair(root.vector, wall).sign(0) <= 0 for wall in walls)
        )
        if any(
            unit_d24_pair(left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height D24+ survivors are not mutually acute")
        accepted.extend(survivors)
        walls.extend(root.vector for root in survivors)
        packets.append(Cubic7D24Packet(height, len(candidates), survivors))
    return Cubic7D24Search(
        height_bound,
        len(levels),
        tuple(packets),
        tuple(accepted),
        tuple(walls),
    )


@dataclass(frozen=True, slots=True)
class Cubic7D24ExtremeRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Cubic7
    kind: str


@dataclass(frozen=True, slots=True)
class Cubic7D24InsertionProfile:
    facet_index: int
    negative_ray_count: int
    zero_ray_count: int
    positive_ray_count: int
    output_ray_count: int


@dataclass(frozen=True, slots=True)
class Cubic7D24PacketGeometry:
    height: Cubic7
    wall_count: int
    extreme_ray_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int


@dataclass(frozen=True, slots=True)
class Cubic7D24Geometry:
    packet_geometries: tuple[Cubic7D24PacketGeometry, ...]
    insertion_profiles: tuple[Cubic7D24InsertionProfile, ...]
    extreme_rays: tuple[Cubic7D24ExtremeRay, ...]
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    finite_volume: bool
    compact: bool


def _functional(wall: Vector) -> Vector:
    return (-A_UNIT * wall[0], *wall[1:])


def _classify_rays(
    rays: Iterable[tuple[Vector, frozenset[int]]],
    functionals: Sequence[Vector],
) -> tuple[Cubic7D24ExtremeRay, ...]:
    output = []
    for ray, active in rays:
        if _row_rank(tuple(functionals[index] for index in active), 25) != 24:
            raise AssertionError("double description returned a non-extreme ray")
        norm = unit_d24_pair(ray, ray)
        sign = norm.sign(0)
        output.append(
            Cubic7D24ExtremeRay(
                ray,
                tuple(sorted(active)),
                norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
            )
        )
    return tuple(sorted(output, key=lambda item: _projective_key(item.vector)))


@lru_cache(maxsize=None)
def analyze_unit_d24_plus_prefix(search: Cubic7D24Search) -> Cubic7D24Geometry:
    """Insert every nonempty equal-height packet by exact DD."""

    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    if not nonempty or len(nonempty[0].accepted_roots) != 1:
        raise ValueError("the unit D24+ prefix lost its affine seed packet")
    functionals = tuple(_functional(wall) for wall in search.walls)
    if _row_rank(functionals[:25], 25) != 25:
        raise AssertionError("the D24 affine seed is not independent")

    rays: dict[Vector, tuple[Vector, frozenset[int]]] = {}
    for omitted in range(25):
        ray = _kernel(tuple(functionals[index] for index in range(25) if index != omitted))
        if ray is None:
            raise AssertionError("the D24 affine cone lost a seed ray")
        values = tuple(_functional_value(functional, ray) for functional in functionals[:25])
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-value for value in ray)
        else:
            raise AssertionError("the D24 affine ray has inconsistent orientation")
        active = frozenset(
            index
            for index in range(25)
            if not _functional_value(functionals[index], oriented)
        )
        rays[_projective_key(oriented)] = (oriented, active)

    classified = _classify_rays(rays.values(), functionals)
    packet_geometries = [
        Cubic7D24PacketGeometry(
            nonempty[0].height,
            25,
            len(classified),
            sum(ray.kind == "finite" for ray in classified),
            sum(ray.kind == "cusp" for ray in classified),
            sum(ray.kind == "hyperideal" for ray in classified),
        )
    ]
    profiles = []
    processed = list(range(25))
    next_facet = 25
    rank_cache: dict[frozenset[int], int] = {}
    for packet in nonempty[1:]:
        for _root in packet.accepted_roots:
            functional = functionals[next_facet]
            negative = []
            zero = []
            positive = []
            for key, (ray, active) in rays.items():
                value = _functional_value(functional, ray)
                destination = (
                    negative if value.sign(0) < 0 else positive if value.sign(0) > 0 else zero
                )
                destination.append((key, ray, active, value))
            updated = {key: (ray, active) for key, ray, active, _value in negative}
            updated.update(
                {
                    key: (ray, active | frozenset((next_facet,)))
                    for key, ray, active, _value in zero
                }
            )
            for _nk, negative_ray, negative_active, negative_value in negative:
                for _pk, positive_ray, positive_active, positive_value in positive:
                    common = negative_active & positive_active
                    if len(common) < 23:
                        continue
                    rank = rank_cache.get(common)
                    if rank is None:
                        rank = _row_rank(
                            tuple(functionals[index] for index in common), 25
                        )
                        rank_cache[common] = rank
                    if rank != 23:
                        continue
                    ray = tuple(
                        positive_value * negative_ray[index]
                        - negative_value * positive_ray[index]
                        for index in range(25)
                    )
                    active = frozenset(
                        index
                        for index in (*processed, next_facet)
                        if not _functional_value(functionals[index], ray)
                    )
                    updated[_projective_key(ray)] = (ray, active)
            rays = updated
            profiles.append(
                Cubic7D24InsertionProfile(
                    next_facet,
                    len(negative),
                    len(zero),
                    len(positive),
                    len(rays),
                )
            )
            processed.append(next_facet)
            next_facet += 1
        classified = _classify_rays(rays.values(), functionals)
        packet_geometries.append(
            Cubic7D24PacketGeometry(
                packet.height,
                next_facet,
                len(classified),
                sum(ray.kind == "finite" for ray in classified),
                sum(ray.kind == "cusp" for ray in classified),
                sum(ray.kind == "hyperideal" for ray in classified),
            )
        )

    classified = _classify_rays(rays.values(), functionals)
    finite = sum(ray.kind == "finite" for ray in classified)
    cusps = sum(ray.kind == "cusp" for ray in classified)
    hyperideal = sum(ray.kind == "hyperideal" for ray in classified)
    return Cubic7D24Geometry(
        tuple(packet_geometries),
        tuple(profiles),
        classified,
        finite,
        cusps,
        hyperideal,
        hyperideal == 0,
        hyperideal == 0 and cusps == 0,
    )


__all__ = [
    "Cubic7D24Candidate",
    "Cubic7D24Geometry",
    "Cubic7D24Packet",
    "Cubic7D24Search",
    "D24_MARKS",
    "FIRST_CONTROL",
    "FIRST_HEIGHT",
    "FIRST_SPATIAL",
    "analyze_unit_d24_plus_prefix",
    "d24_height_zero_simple_roots_orthonormal",
    "enumerate_dominant_d24_plus_spatial",
    "enumerate_unit_d24_control_level",
    "run_unit_d24_plus_prefix",
    "unit_d24_pair",
    "unit_d24_plus_first_packet",
]
