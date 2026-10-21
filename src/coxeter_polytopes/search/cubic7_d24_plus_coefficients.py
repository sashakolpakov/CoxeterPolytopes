"""Exact coefficient screen for cubic-seven ``D24+`` hosts.

The lattices are

``L(a) = < -a > + D24+``

over ``O_K = Z[theta]``, ``theta=2*cos(2*pi/7)``.  This module compares the
three genuinely different small squareclasses supported at the inert prime
two and the ramified prime ``delta=2-theta``.  Unit associates are not
rescanned: explicitly

``(theta^2-2)^2*theta = theta-1``.

The spatial enumerator is the glue-complete orthonormal ``D24+`` enumerator:
it treats the root coset and all seven nonzero classes in ``O_K/2O_K`` for
the spinor coset.  Equal-height Vinberg candidates are selected
simultaneously.  Cone snapshots use exact double description and certify
adjacency by rank 23 of common active functionals.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import groupby
from typing import Iterable, Sequence, TypeAlias

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ZERO,
    THETA7,
    Cubic7,
)
from coxeter_polytopes.search.cubic7_a1_golay import A_UNIT
from coxeter_polytopes.search.cubic7_bna1 import (
    _functional_value,
    _kernel,
    _projective_key,
    _row_rank,
)
from coxeter_polytopes.search.cubic7_d24_plus_unit import (
    d24_height_zero_simple_roots_orthonormal,
    enumerate_dominant_d24_plus_spatial,
)
from coxeter_polytopes.search.cubic7_dn import (
    d24_plus_basis_audit,
    d24_plus_gram,
    d24_plus_integral_basis_coordinates,
)
from coxeter_polytopes.search.cubic7_root_blocks import _bounded_for_caps


Vector: TypeAlias = tuple[Cubic7, ...]
DELTA7 = 2 - THETA7
UNIT_ISOMETRY_MULTIPLIER = THETA7**2 - 2


_LOCAL_NORM_ARGUMENT = (
    "For primitive r=(x,y), let I(r)=(a*x,(y,D24+)). Since D24+ is even "
    "unimodular, I(r) is a unit away from primes dividing a. Reflection "
    "integrality gives (q(r))|2I(r), so no other prime can divide q. At the "
    "odd squarefree prime p=(2-theta), either y is primitive and v_p(q)=0, "
    "or y lies in pD24+, x is a unit, and q=-a*x^2 mod p^2 has v_p(q)=1. "
    "At inert (2), evenness gives y^2 in (2). If a is odd, x a 2-adic "
    "unit gives v_2(q)=0, while x in (2) and primitive y gives v_2(q)=1. "
    "If a has v_2(a)=1, the whole form is even and the same congruence plus "
    "q|2I forces v_2(q)=1. Class number one and U^+=U^2 turn these ideal "
    "classes into the displayed totally-positive element representatives. "
)


@dataclass(frozen=True, slots=True)
class Cubic7D24CoefficientSeed:
    seed_id: str
    timelike_coefficient: Cubic7
    canonical_mirror_norms: tuple[Cubic7, ...]
    support: str
    complete_norm_class_proof: str

    @property
    def admissible(self) -> bool:
        return self.timelike_coefficient.sign(0) > 0 and all(
            self.timelike_coefficient.sign(place) < 0 for place in (1, 2)
        )


TWO_THETA_D24_SEED = Cubic7D24CoefficientSeed(
    "cubic7-H24-D24+-a-2theta",
    2 * THETA7,
    (Cubic7(2),),
    "(2)",
    _LOCAL_NORM_ARGUMENT
    + "Here v_2(a)=1 and there is no odd support, so every mirror norm ideal "
    "is exactly (2); canonical norm 2 is complete.",
)

P7_D24_SEED = Cubic7D24CoefficientSeed(
    "cubic7-H24-D24+-a-theta-delta",
    THETA7 * DELTA7,
    (Cubic7(1), Cubic7(2), DELTA7, 2 * DELTA7),
    "(2-theta)",
    _LOCAL_NORM_ARGUMENT
    + "Here a is odd with squarefree support p=(2-theta), so the complete "
    "canonical norms are 1, 2, delta, and 2*delta.",
)

TWO_P7_D24_SEED = Cubic7D24CoefficientSeed(
    "cubic7-H24-D24+-a-2theta-delta",
    2 * THETA7 * DELTA7,
    (Cubic7(2), 2 * DELTA7),
    "(2)*(2-theta)",
    _LOCAL_NORM_ARGUMENT
    + "Here v_2(a)=v_p(a)=1. The 2-adic valuation of every mirror norm is "
    "one and the p-adic valuation is optional, so the complete canonical "
    "norms are 2 and 2*delta.",
)

SUPPORTED_D24_COEFFICIENT_SEEDS = (
    TWO_THETA_D24_SEED,
    P7_D24_SEED,
    TWO_P7_D24_SEED,
)


def unit_coefficient_isometry_audit() -> bool:
    """Verify that ``theta`` and ``theta-1`` are time-line isometric."""

    multiplier = UNIT_ISOMETRY_MULTIPLIER
    return (
        multiplier.is_integral
        and abs(multiplier.norm()) == 1
        and multiplier * multiplier * THETA7 == A_UNIT
    )


@dataclass(frozen=True, slots=True)
class Cubic7D24CoefficientLevel:
    norm: Cubic7
    height: Cubic7
    control: Cubic7
    target: Cubic7


@lru_cache(maxsize=None)
def enumerate_coefficient_control_levels(
    seed: Cubic7D24CoefficientSeed, height_bound: int
) -> tuple[Cubic7D24CoefficientLevel, ...]:
    """Exhaust every reflective control compatible with the archimedean box."""

    if seed not in SUPPORTED_D24_COEFFICIENT_SEEDS or not seed.admissible:
        raise ValueError("unsupported or inadmissible cubic-seven D24+ seed")
    if height_bound < 0:
        raise ValueError("height bound must be nonnegative")
    coefficient = seed.timelike_coefficient
    levels: set[Cubic7D24CoefficientLevel] = set()
    for norm in seed.canonical_mirror_norms:
        caps = (
            norm * height_bound / coefficient,
            norm / (-coefficient),
            norm / (-coefficient),
        )
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0:
                continue
            if not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if any(target.sign(place) < 0 for place in (1, 2)):
                continue
            height = coefficient * control * control / norm
            if (Cubic7(height_bound) - height).sign(0) < 0:
                continue
            levels.add(Cubic7D24CoefficientLevel(norm, height, control, target))
    return tuple(
        sorted(levels, key=lambda level: (level.height, level.norm, level.control))
    )


def d24_coefficient_pair(
    seed: Cubic7D24CoefficientSeed,
    left: Sequence[Cubic7],
    right: Sequence[Cubic7],
) -> Cubic7:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("a D24+ host vector has rank 25")
    return -seed.timelike_coefficient * left[0] * right[0] + sum(
        (left[index] * right[index] for index in range(1, 25)),
        start=CUBIC7_ZERO,
    )


def _reflection_integral(
    seed: Cubic7D24CoefficientSeed,
    norm: Cubic7,
    control: Cubic7,
    spatial: Vector,
) -> bool:
    if not (2 * seed.timelike_coefficient * control / norm).is_integral:
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
class Cubic7D24CoefficientCandidate:
    norm: Cubic7
    height: Cubic7
    control: Cubic7
    target: Cubic7
    coset: str
    vector: Vector


@lru_cache(maxsize=None)
def enumerate_coefficient_level(
    seed: Cubic7D24CoefficientSeed,
    level: Cubic7D24CoefficientLevel,
) -> tuple[Cubic7D24CoefficientCandidate, ...]:
    output = []
    # If delta divides q, unimodularity makes reflection integrality
    # equivalent to y lying in delta*D24+: for q=delta the factor 2 is a
    # p-adic unit, and for q=2*delta it cancels the displayed 2.  Enumerating
    # z=y/delta at target/delta^2 is both exact and dramatically smaller than
    # generating nonreflective vectors and filtering them afterwards.
    delta_divisible = level.norm in (DELTA7, 2 * DELTA7)
    spatial_scale = DELTA7 if delta_divisible else Cubic7(1)
    reduced_target = (
        level.target / (DELTA7 * DELTA7)
        if delta_divisible
        else level.target
    )
    if delta_divisible and not reduced_target.is_integral:
        return ()
    for pattern in enumerate_dominant_d24_plus_spatial(reduced_target):
        spatial = tuple(spatial_scale * value for value in pattern.coordinates)
        if not _reflection_integral(
            seed, level.norm, level.control, spatial
        ):
            continue
        vector = (level.control, *spatial)
        if d24_coefficient_pair(seed, vector, vector) != level.norm:
            raise AssertionError("a D24+ coefficient candidate has wrong norm")
        output.append(
            Cubic7D24CoefficientCandidate(
                level.norm,
                level.height,
                level.control,
                level.target,
                pattern.coset,
                vector,
            )
        )
    return tuple(sorted(output, key=lambda root: (root.norm, root.vector, root.coset)))


@dataclass(frozen=True, slots=True)
class Cubic7D24CoefficientPacket:
    height: Cubic7
    candidate_count: int
    accepted_roots: tuple[Cubic7D24CoefficientCandidate, ...]


@dataclass(frozen=True, slots=True)
class Cubic7D24CoefficientSearch:
    seed: Cubic7D24CoefficientSeed
    height_bound: int
    control_level_count: int
    packets: tuple[Cubic7D24CoefficientPacket, ...]
    accepted_roots: tuple[Cubic7D24CoefficientCandidate, ...]


@lru_cache(maxsize=None)
def run_coefficient_prefix(
    seed: Cubic7D24CoefficientSeed, height_bound: int = 100
) -> Cubic7D24CoefficientSearch:
    """Run glue-complete simultaneous-height Vinberg selection."""

    if not d24_plus_basis_audit():
        raise AssertionError("the D24+ glue basis lost even unimodularity")
    levels = enumerate_coefficient_control_levels(seed, height_bound)
    walls = list(d24_height_zero_simple_roots_orthonormal())
    accepted: list[Cubic7D24CoefficientCandidate] = []
    packets = []
    for height, iterator in groupby(levels, key=lambda level: level.height):
        candidates = tuple(
            root
            for level in iterator
            for root in enumerate_coefficient_level(seed, level)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(
                d24_coefficient_pair(seed, root.vector, wall).sign(0) <= 0
                for wall in walls
            )
        )
        if any(
            d24_coefficient_pair(seed, left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height D24+ coefficient roots are not acute")
        accepted.extend(survivors)
        walls.extend(root.vector for root in survivors)
        packets.append(
            Cubic7D24CoefficientPacket(height, len(candidates), survivors)
        )
    return Cubic7D24CoefficientSearch(
        seed, height_bound, len(levels), tuple(packets), tuple(accepted)
    )


@dataclass(frozen=True, slots=True)
class Cubic7D24CoefficientRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Cubic7
    kind: str


@dataclass(frozen=True, slots=True)
class Cubic7D24CoefficientInsertion:
    facet_index: int
    negative_ray_count: int
    zero_ray_count: int
    positive_ray_count: int
    output_ray_count: int


@dataclass(frozen=True, slots=True)
class Cubic7D24CoefficientSnapshot:
    height: Cubic7
    wall_count: int
    extreme_ray_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int


@dataclass(frozen=True, slots=True)
class Cubic7D24CoefficientGeometry:
    packet_limit: int | None
    snapshots: tuple[Cubic7D24CoefficientSnapshot, ...]
    insertions: tuple[Cubic7D24CoefficientInsertion, ...]
    rays: tuple[Cubic7D24CoefficientRay, ...]
    ray_stream_sha256: str

    @property
    def finite_vertex_count(self) -> int:
        return sum(ray.kind == "finite" for ray in self.rays)

    @property
    def cusp_count(self) -> int:
        return sum(ray.kind == "cusp" for ray in self.rays)

    @property
    def hyperideal_ray_count(self) -> int:
        return sum(ray.kind == "hyperideal" for ray in self.rays)

    @property
    def finite_volume(self) -> bool:
        return self.hyperideal_ray_count == 0

    @property
    def compact(self) -> bool:
        return self.finite_volume and self.cusp_count == 0


def _functional(seed: Cubic7D24CoefficientSeed, wall: Vector) -> Vector:
    return (-seed.timelike_coefficient * wall[0], *wall[1:])


def _classify(
    seed: Cubic7D24CoefficientSeed,
    rays: Iterable[tuple[Vector, frozenset[int]]],
    functionals: Sequence[Vector],
) -> tuple[Cubic7D24CoefficientRay, ...]:
    output = []
    for ray, active in rays:
        if _row_rank(tuple(functionals[index] for index in active), 25) != 24:
            raise AssertionError("double description returned a non-extreme ray")
        norm = d24_coefficient_pair(seed, ray, ray)
        sign = norm.sign(0)
        output.append(
            Cubic7D24CoefficientRay(
                ray,
                tuple(sorted(active)),
                norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
            )
        )
    return tuple(sorted(output, key=lambda item: _projective_key(item.vector)))


def _ray_digest(rays: Sequence[Cubic7D24CoefficientRay]) -> str:
    digest = sha256()
    for ray in rays:
        for value in ray.vector:
            digest.update(
                f"{value.a}/{value.b}/{value.c};".encode("ascii")
            )
        digest.update(
            (",".join(map(str, ray.active_facets)) + ";" + ray.kind + "\n").encode(
                "ascii"
            )
        )
    return digest.hexdigest()


@lru_cache(maxsize=None)
def analyze_coefficient_prefix(
    search: Cubic7D24CoefficientSearch,
    packet_limit: int | None = None,
) -> Cubic7D24CoefficientGeometry:
    """Classify the cone after each selected nonempty packet exactly."""

    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    if packet_limit is not None:
        if packet_limit < 1:
            raise ValueError("packet_limit must be positive")
        nonempty = nonempty[:packet_limit]
    if not nonempty or len(nonempty[0].accepted_roots) != 1:
        raise ValueError("the coefficient screen lost its affine seed packet")
    seed = search.seed
    positive_roots = tuple(
        root.vector for packet in nonempty for root in packet.accepted_roots
    )
    walls = (*d24_height_zero_simple_roots_orthonormal(), *positive_roots)
    functionals = tuple(_functional(seed, wall) for wall in walls)
    if _row_rank(functionals[:25], 25) != 25:
        raise AssertionError("the D24 affine seed is not independent")

    rays: dict[Vector, tuple[Vector, frozenset[int]]] = {}
    for omitted in range(25):
        ray = _kernel(
            tuple(functionals[index] for index in range(25) if index != omitted)
        )
        if ray is None:
            raise AssertionError("the D24 affine cone lost a seed ray")
        values = tuple(
            _functional_value(functional, ray) for functional in functionals[:25]
        )
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

    def snapshot(height: Cubic7, wall_count: int) -> Cubic7D24CoefficientSnapshot:
        classified = _classify(seed, rays.values(), functionals)
        return Cubic7D24CoefficientSnapshot(
            height,
            wall_count,
            len(classified),
            sum(ray.kind == "finite" for ray in classified),
            sum(ray.kind == "cusp" for ray in classified),
            sum(ray.kind == "hyperideal" for ray in classified),
        )

    snapshots = [snapshot(nonempty[0].height, 25)]
    insertions = []
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
                    negative
                    if value.sign(0) < 0
                    else positive
                    if value.sign(0) > 0
                    else zero
                )
                destination.append((key, ray, active, value))
            updated = {key: (ray, active) for key, ray, active, _ in negative}
            updated.update(
                {
                    key: (ray, active | frozenset((next_facet,)))
                    for key, ray, active, _ in zero
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
            insertions.append(
                Cubic7D24CoefficientInsertion(
                    next_facet,
                    len(negative),
                    len(zero),
                    len(positive),
                    len(rays),
                )
            )
            processed.append(next_facet)
            next_facet += 1
        snapshots.append(snapshot(packet.height, next_facet))

    classified = _classify(seed, rays.values(), functionals)
    return Cubic7D24CoefficientGeometry(
        packet_limit,
        tuple(snapshots),
        tuple(insertions),
        classified,
        _ray_digest(classified),
    )


__all__ = [
    "DELTA7",
    "P7_D24_SEED",
    "SUPPORTED_D24_COEFFICIENT_SEEDS",
    "TWO_P7_D24_SEED",
    "TWO_THETA_D24_SEED",
    "UNIT_ISOMETRY_MULTIPLIER",
    "Cubic7D24CoefficientCandidate",
    "Cubic7D24CoefficientGeometry",
    "Cubic7D24CoefficientLevel",
    "Cubic7D24CoefficientSearch",
    "Cubic7D24CoefficientSeed",
    "analyze_coefficient_prefix",
    "d24_coefficient_pair",
    "enumerate_coefficient_control_levels",
    "enumerate_coefficient_level",
    "run_coefficient_prefix",
    "unit_coefficient_isometry_audit",
]
