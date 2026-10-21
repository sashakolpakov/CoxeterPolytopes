"""Exact H24 Vinberg search for ``<-(2*sqrt(2)-1)> + D24+``.

Put ``pi=sqrt(2)``, ``epsilon=1+pi`` and ``delta=3+pi``.  The
timelike coefficient ``a=2*pi-1`` has signs ``(+,-)`` and

``epsilon*a = delta,    N(delta)=7``.

Thus its only odd ideal support is the split prime ``p=(delta)``.  The
complete primitive mirror norms, modulo totally-positive unit squares, are
``1, 2, delta, 2*delta``.  A mirror whose norm contains ``delta`` has
spatial component in ``delta*D24+``; we impose that local condition before
enumerating the glue lattice.  This turns the apparently large prime-seven
targets into small, exact square-partition problems.

The spatial enumeration is the orthonormal, glue-complete ``D24+``
enumerator from :mod:`silver_d24_plus`: it covers the root coset and all
three nonzero residue classes in ``O_K/2O_K`` for the spinor coset.
Primitivity is checked in the integral glue basis.  In particular, a raw
norm-two vector divisible by ``pi`` is removed because it represents the
same mirror as its primitive norm-one quotient.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import groupby
from typing import Iterable, Sequence, TypeAlias

from coxeter_polytopes.exact.silver import (
    SILVER_ONE,
    SILVER_ZERO,
    SQRT2,
    Silver,
)
from coxeter_polytopes.search.silver_d24_plus import (
    _integer_determinant,
    _silver_functional_value,
    _silver_kernel,
    _silver_projective_key,
    _silver_row_rank,
    d24_plus_gram,
    d24_plus_integral_basis_coordinates,
    enumerate_dominant_d24_plus_euclidean,
)
from coxeter_polytopes.search.silver_e8_cubes import _bounded_for_caps


Vector: TypeAlias = tuple[Silver, ...]
EPSILON = SILVER_ONE + SQRT2
P7_DELTA = Silver(3, 1)
P7_TIMELIKE = Silver(-1, 2)
H100_FROZEN_WITNESS_VECTOR: Vector = (
    Silver(20, 12) / 7,
    Silver(2, 2),
    Silver(2, 2),
    Silver(2),
    Silver(2),
    *(SILVER_ZERO for _ in range(20)),
)
H300_REPLACEMENT_WITNESS_VECTOR: Vector = (
    SILVER_ONE,
    *(Silver(172, -120) / 7 for _ in range(19)),
    *(SILVER_ZERO for _ in range(5)),
)
H500_REPLACEMENT_WITNESS_VECTOR: Vector = (
    SILVER_ONE,
    *(Silver(4, Fraction(-12, 5)) for _ in range(5)),
    *(SILVER_ZERO for _ in range(19)),
)


COMPLETE_P7_D24_NORM_PROOF = (
    "Put pi=sqrt(2), epsilon=1+pi, delta=3+pi, and a=2*pi-1, so "
    "epsilon*a=delta and (delta) is a prime of norm seven. For primitive "
    "r=(x,y), reflection integrality says (q(r)) divides 2I(r), where "
    "I(r)=(a*x,(y,D24+)). Even unimodularity excludes every prime away "
    "from (pi) and (delta). At (delta), either y is primitive and the "
    "valuation is zero, or reflection integrality gives y in delta*D24+; "
    "then x is a unit and q=-a*x^2 mod delta^2 has valuation one. At "
    "(pi), if x is a unit then q=-a*x^2 mod pi^2 has valuation zero. "
    "Otherwise primitivity makes the spatial pairing ideal a unit; "
    "evenness gives valuation at least two and q|2I gives valuation "
    "exactly two, namely the ideal (2). Since Z[sqrt(2)] has class number "
    "one and every totally-positive unit is a square, the complete "
    "canonical norms are 1, 2, delta, and 2*delta."
)


@dataclass(frozen=True, slots=True)
class SilverP7D24Seed:
    seed_id: str
    timelike_coefficient: Silver
    canonical_mirror_norms: tuple[Silver, ...]
    complete_norm_class_proof: str

    @property
    def admissible(self) -> bool:
        return (
            self.timelike_coefficient.sign(0) > 0
            and self.timelike_coefficient.sign(1) < 0
        )


SILVER_P7_D24_SEED = SilverP7D24Seed(
    "silver-H24-D24+-a-2sqrt2-minus1",
    P7_TIMELIKE,
    (SILVER_ONE, Silver(2), P7_DELTA, 2 * P7_DELTA),
    COMPLETE_P7_D24_NORM_PROOF,
)


def p7_factorization_audit() -> bool:
    """Check the sign and prime-norm identities used by the local proof."""

    return (
        SILVER_P7_D24_SEED.admissible
        and EPSILON * P7_TIMELIKE == P7_DELTA
        and P7_TIMELIKE.norm() == -7
        and P7_DELTA.norm() == 7
        and P7_DELTA.sign(0) > 0
        and P7_DELTA.sign(1) > 0
    )


@dataclass(frozen=True, slots=True)
class SilverP7D24ControlLevel:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_p7_d24_control_levels(
    height_bound: int = 100,
) -> tuple[SilverP7D24ControlLevel, ...]:
    """Exhaust all reflective control coordinates through a height bound."""

    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    if not p7_factorization_audit():
        raise AssertionError("the silver prime-seven seed identities changed")
    coefficient = P7_TIMELIKE
    levels: set[SilverP7D24ControlLevel] = set()
    for norm in SILVER_P7_D24_SEED.canonical_mirror_norms:
        caps = (norm * height_bound / coefficient, norm / (-coefficient))
        for control in _bounded_for_caps(caps):
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
            levels.add(SilverP7D24ControlLevel(norm, height, control, target))
    return tuple(
        sorted(levels, key=lambda level: (level.height, level.norm, level.control))
    )


def silver_p7_d24_pair(left: Sequence[Silver], right: Sequence[Silver]) -> Silver:
    """Lorentzian pairing in orthonormal spatial coordinates."""

    if len(left) != 25 or len(right) != 25:
        raise ValueError("a silver prime-seven D24+ vector has rank 25")
    return -P7_TIMELIKE * left[0] * right[0] + sum(
        (left[index] * right[index] for index in range(1, 25)),
        start=SILVER_ZERO,
    )


@lru_cache(maxsize=1)
def silver_d24_height_zero_simple_roots() -> tuple[Vector, ...]:
    """Negative D24 simple roots compatible with dominant coordinates."""

    roots = []
    for index in range(23):
        root = [SILVER_ZERO] * 25
        root[1 + index] = -SILVER_ONE
        root[2 + index] = SILVER_ONE
        roots.append(tuple(root))
    fork = [SILVER_ZERO] * 25
    fork[23] = fork[24] = -SILVER_ONE
    roots.append(tuple(fork))
    expected = d24_plus_gram()
    # The basis-independent check is determinant/even-unimodularity below;
    # here only the actual root norms are needed for the orthonormal model.
    if any(silver_p7_d24_pair(root, root) != 2 for root in roots):
        raise AssertionError("a height-zero D24 wall lost norm two")
    if _integer_determinant(expected) != 1:
        raise AssertionError("the D24+ glue basis lost unimodularity")
    return tuple(roots)


def _glue_basis_coordinates(spatial: Sequence[Silver]) -> Vector:
    coordinates = d24_plus_integral_basis_coordinates(spatial)
    if not all(value.is_integral for value in coordinates):
        raise AssertionError("a D24+ pattern left the integral glue basis")
    return coordinates


def _spatial_pairings(spatial: Sequence[Silver]) -> Vector:
    coordinates = _glue_basis_coordinates(spatial)
    gram = d24_plus_gram()
    return tuple(
        sum(
            (
                Silver(gram[index][other]) * coordinates[other]
                for other in range(24)
                if gram[index][other]
            ),
            start=SILVER_ZERO,
        )
        for index in range(24)
    )


def _reflection_integral(
    norm: Silver, control: Silver, spatial: Sequence[Silver]
) -> bool:
    if not (2 * P7_TIMELIKE * control / norm).is_integral:
        return False
    return all((2 * pairing / norm).is_integral for pairing in _spatial_pairings(spatial))


def _divisible_by_sqrt2(control: Silver, spatial: Sequence[Silver]) -> bool:
    """Test divisibility in ``O*x + D24+`` rather than coordinatewise."""

    if not (control / SQRT2).is_integral:
        return False
    return all(
        (coordinate / SQRT2).is_integral
        for coordinate in _glue_basis_coordinates(spatial)
    )


def silver_p7_d24_vector_divisible_by_sqrt2(vector: Sequence[Silver]) -> bool:
    if len(vector) != 25:
        raise ValueError("a silver prime-seven D24+ vector has rank 25")
    return _divisible_by_sqrt2(vector[0], vector[1:])


@dataclass(frozen=True, slots=True)
class SilverP7D24Candidate:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver
    coset: str
    vector: Vector


@dataclass(frozen=True, slots=True)
class SilverP7D24LevelEnumeration:
    level: SilverP7D24ControlLevel
    reduced_target: Silver | None
    glue_pattern_count: int
    reflection_integral_count: int
    imprimitive_sqrt2_count: int
    imprimitive_vectors: tuple[Vector, ...]
    candidates: tuple[SilverP7D24Candidate, ...]


def _has_p7_norm(norm: Silver) -> bool:
    return norm in (P7_DELTA, 2 * P7_DELTA)


@lru_cache(maxsize=None)
def enumerate_p7_d24_level(
    level: SilverP7D24ControlLevel,
) -> SilverP7D24LevelEnumeration:
    """Enumerate one glue-complete level, imposing local divisibility first."""

    if (
        level.norm not in SILVER_P7_D24_SEED.canonical_mirror_norms
        or level.control.sign(0) <= 0
        or level.target != level.norm + P7_TIMELIKE * level.control**2
        or level.height != P7_TIMELIKE * level.control**2 / level.norm
        or not (2 * P7_TIMELIKE * level.control / level.norm).is_integral
        or level.target.sign(1) < 0
    ):
        raise ValueError("unsupported prime-seven D24+ control level")
    scale = P7_DELTA if _has_p7_norm(level.norm) else SILVER_ONE
    reduced_target = level.target / (scale * scale)
    if not reduced_target.is_integral:
        return SilverP7D24LevelEnumeration(level, None, 0, 0, 0, (), ())
    if any(reduced_target.sign(place) < 0 for place in (0, 1)):
        return SilverP7D24LevelEnumeration(
            level, reduced_target, 0, 0, 0, (), ()
        )

    patterns = enumerate_dominant_d24_plus_euclidean(reduced_target)
    reflection_integral = 0
    imprimitive = []
    output = []
    for reduced_spatial, coset in patterns:
        spatial = tuple(scale * value for value in reduced_spatial)
        if not _reflection_integral(level.norm, level.control, spatial):
            continue
        reflection_integral += 1
        vector = (level.control, *spatial)
        if silver_p7_d24_pair(vector, vector) != level.norm:
            raise AssertionError("a prime-seven D24+ candidate has wrong norm")
        if _divisible_by_sqrt2(level.control, spatial):
            imprimitive.append(vector)
            continue
        output.append(
            SilverP7D24Candidate(
                level.norm,
                level.height,
                level.control,
                level.target,
                coset,
                vector,
            )
        )
    return SilverP7D24LevelEnumeration(
        level,
        reduced_target,
        len(patterns),
        reflection_integral,
        len(imprimitive),
        tuple(imprimitive),
        tuple(sorted(output, key=lambda root: (root.norm, root.vector, root.coset))),
    )


@dataclass(frozen=True, slots=True)
class SilverP7D24Packet:
    height: Silver
    raw_reflective_count: int
    imprimitive_sqrt2_count: int
    candidate_count: int
    accepted_roots: tuple[SilverP7D24Candidate, ...]


@dataclass(frozen=True, slots=True)
class SilverP7D24Search:
    height_bound: int
    control_levels: tuple[SilverP7D24ControlLevel, ...]
    level_enumerations: tuple[SilverP7D24LevelEnumeration, ...]
    packets: tuple[SilverP7D24Packet, ...]
    accepted_roots: tuple[SilverP7D24Candidate, ...]
    walls: tuple[Vector, ...]


@lru_cache(maxsize=None)
def run_silver_p7_d24_prefix(height_bound: int = 100) -> SilverP7D24Search:
    """Run simultaneous-height Vinberg selection through ``height_bound``."""

    levels = enumerate_p7_d24_control_levels(height_bound)
    enumerations = tuple(enumerate_p7_d24_level(level) for level in levels)
    walls = list(silver_d24_height_zero_simple_roots())
    accepted: list[SilverP7D24Candidate] = []
    packets = []
    for height, iterator in groupby(enumerations, key=lambda item: item.level.height):
        packet_levels = tuple(iterator)
        candidates = tuple(
            root for enumeration in packet_levels for root in enumeration.candidates
        )
        survivors = tuple(
            root
            for root in candidates
            if all(
                silver_p7_d24_pair(root.vector, wall).sign(0) <= 0
                for wall in walls
            )
        )
        if any(
            silver_p7_d24_pair(left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height prime-seven roots are not acute")
        accepted.extend(survivors)
        walls.extend(root.vector for root in survivors)
        packets.append(
            SilverP7D24Packet(
                height,
                sum(item.reflection_integral_count for item in packet_levels),
                sum(item.imprimitive_sqrt2_count for item in packet_levels),
                len(candidates),
                survivors,
            )
        )
    return SilverP7D24Search(
        height_bound,
        levels,
        enumerations,
        tuple(packets),
        tuple(accepted),
        tuple(walls),
    )


@dataclass(frozen=True, slots=True)
class SilverP7D24Ray:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Silver
    kind: str


@dataclass(frozen=True, slots=True)
class SilverP7D24Insertion:
    facet_index: int
    negative_ray_count: int
    zero_ray_count: int
    positive_ray_count: int
    output_ray_count: int


@dataclass(frozen=True, slots=True)
class SilverP7D24Snapshot:
    height: Silver
    wall_count: int
    extreme_ray_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int


@dataclass(frozen=True, slots=True)
class SilverP7D24Geometry:
    packet_limit: int | None
    snapshots: tuple[SilverP7D24Snapshot, ...]
    insertions: tuple[SilverP7D24Insertion, ...]
    rays: tuple[SilverP7D24Ray, ...]
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


@dataclass(frozen=True, slots=True)
class SilverP7D24HyperidealWitness:
    """An old extreme ray that survives every wall of a longer prefix."""

    ray: SilverP7D24Ray
    later_wall_count: int
    later_pairings: tuple[Silver, ...]
    new_zero_facets: tuple[int, ...]
    final_active_rank: int
    surviving_hyperideal_ray_count: int


@dataclass(frozen=True, slots=True)
class SilverP7D24WitnessPacketAudit:
    height: Silver
    wall_count: int
    negative_pairing_count: int
    zero_pairing_count: int
    positive_pairing_count: int
    minimum_pairing: Silver
    maximum_pairing: Silver


@dataclass(frozen=True, slots=True)
class SilverP7D24WitnessCutAudit:
    vector: Vector
    norm: Silver
    packets: tuple[SilverP7D24WitnessPacketAudit, ...]
    first_cut_height: Silver | None
    first_cut_packet_positive_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class SilverP7D24InitialRaySurvival:
    height: Silver
    cumulative_later_wall_count: int
    surviving_hyperideal_ray_count: int


@dataclass(frozen=True, slots=True)
class SilverP7D24DirectPrefixWitness:
    """A directly specified exact extreme ray of a complete wall prefix."""

    vector: Vector
    norm: Silver
    wall_count: int
    pairings: tuple[Silver, ...]
    active_facets: tuple[int, ...]
    active_rank: int


def _functional(wall: Vector) -> Vector:
    return (-P7_TIMELIKE * wall[0], *wall[1:])


def _fraction_mod_prime(value: Fraction, prime: int) -> int | None:
    denominator = value.denominator % prime
    if not denominator:
        return None
    return value.numerator % prime * pow(denominator, -1, prime) % prime


def _silver_row_rank_specialization(
    rows: Sequence[Sequence[Silver]],
    width: int,
    *,
    prime: int,
    radical: int,
) -> int | None:
    """Rank after the exact specialization ``sqrt(2) -> radical mod prime``.

    A nonzero minor after specialization is a certificate that the same
    minor is nonzero over ``Q(sqrt(2))``.  A rank drop is inconclusive.
    """

    if radical * radical % prime != 2 % prime:
        raise ValueError("radical is not a square root of two modulo prime")
    work = []
    for row in rows:
        converted = []
        for value in row:
            rational = _fraction_mod_prime(value.a, prime)
            irrational = _fraction_mod_prime(value.b, prime)
            if rational is None or irrational is None:
                return None
            converted.append((rational + radical * irrational) % prime)
        work.append(converted)
    rank = 0
    for column in range(width):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        for row in range(rank + 1, len(work)):
            if not work[row][column]:
                continue
            factor = work[row][column] * inverse % prime
            for index in range(column, width):
                work[row][index] = (
                    work[row][index] - factor * work[rank][index]
                ) % prime
        rank += 1
        if rank == len(work):
            break
    return rank


_RANK_SPECIALIZATIONS = ((17, 6), (17, 11), (23, 5), (23, 18), (31, 8), (31, 23))


def _silver_rank_equals_certified(
    rows: Sequence[Sequence[Silver]], width: int, expected: int
) -> bool:
    """Decide a rank equality, using modular nonvanishing when conclusive."""

    maximum = min(len(rows), width)
    for prime, radical in _RANK_SPECIALIZATIONS:
        rank = _silver_row_rank_specialization(
            rows, width, prime=prime, radical=radical
        )
        if rank is None:
            continue
        if rank > expected:
            return False
        if rank == expected == maximum:
            return True
    return _silver_row_rank(rows, width) == expected


def _classify_rays(
    rays: Iterable[tuple[Vector, frozenset[int]]],
    functionals: Sequence[Vector],
) -> tuple[SilverP7D24Ray, ...]:
    output = []
    for ray, active in rays:
        if not _silver_rank_equals_certified(
            tuple(functionals[index] for index in active), 25, 24
        ):
            raise AssertionError("double description returned a non-extreme ray")
        norm = silver_p7_d24_pair(ray, ray)
        sign = norm.sign(0)
        output.append(
            SilverP7D24Ray(
                ray,
                tuple(sorted(active)),
                norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
            )
        )
    return tuple(sorted(output, key=lambda item: _silver_projective_key(item.vector)))


def _ray_digest(rays: Sequence[SilverP7D24Ray]) -> str:
    digest = sha256()
    for ray in rays:
        for value in ray.vector:
            digest.update(
                (
                    f"{value.a.numerator}/{value.a.denominator}/"
                    f"{value.b.numerator}/{value.b.denominator};"
                ).encode("ascii")
            )
        digest.update(
            (",".join(map(str, ray.active_facets)) + ";" + ray.kind + "\n").encode(
                "ascii"
            )
        )
    return digest.hexdigest()


@lru_cache(maxsize=None)
def analyze_silver_p7_d24_prefix(
    search: SilverP7D24Search,
    packet_limit: int | None = None,
) -> SilverP7D24Geometry:
    """Classify exact cone rays through the requested nonempty packets."""

    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    if packet_limit is not None:
        if packet_limit < 1:
            raise ValueError("packet_limit must be positive")
        nonempty = nonempty[:packet_limit]
    if not nonempty:
        raise ValueError("the prefix has no positive-height walls")
    first_roots = nonempty[0].accepted_roots
    affine = tuple(root for root in first_roots if root.norm == SILVER_ONE)
    if len(affine) != 1:
        raise ValueError("the first packet lost its unique affine norm-one wall")

    positive_roots = [affine[0]]
    positive_roots.extend(root for root in first_roots if root is not affine[0])
    for packet in nonempty[1:]:
        positive_roots.extend(packet.accepted_roots)
    walls = (*silver_d24_height_zero_simple_roots(), *(root.vector for root in positive_roots))
    functionals = tuple(_functional(wall) for wall in walls)
    if not _silver_rank_equals_certified(functionals[:25], 25, 25):
        raise AssertionError("the prime-seven D24 affine seed is not independent")

    rays: dict[Vector, tuple[Vector, frozenset[int]]] = {}
    for omitted in range(25):
        support = tuple(index for index in range(25) if index != omitted)
        ray = _silver_kernel(tuple(functionals[index] for index in support))
        if ray is None:
            raise AssertionError("the affine seed cone lost a ray")
        values = tuple(
            _silver_functional_value(functionals[index], ray) for index in range(25)
        )
        if all(value.sign(0) <= 0 for value in values):
            oriented = ray
        elif all(value.sign(0) >= 0 for value in values):
            oriented = tuple(-value for value in ray)
        else:
            raise AssertionError("the affine seed ray has inconsistent orientation")
        active = frozenset(
            index
            for index in range(25)
            if not _silver_functional_value(functionals[index], oriented)
        )
        rays[_silver_projective_key(oriented)] = (oriented, active)

    insertions = []
    snapshots = []
    processed = list(range(25))
    rank_cache: dict[frozenset[int], int] = {}
    next_facet = 25

    def insert_one() -> None:
        nonlocal rays, next_facet
        functional = functionals[next_facet]
        negative = []
        zero = []
        positive = []
        for key, (ray, active) in rays.items():
            value = _silver_functional_value(functional, ray)
            destination = (
                negative if value.sign(0) < 0 else positive if value.sign(0) > 0 else zero
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
                    rank = (
                        23
                        if _silver_rank_equals_certified(
                            tuple(functionals[index] for index in common), 25, 23
                        )
                        else -1
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
                    if not _silver_functional_value(functionals[index], ray)
                )
                updated[_silver_projective_key(ray)] = (ray, active)
        rays = updated
        insertions.append(
            SilverP7D24Insertion(
                next_facet, len(negative), len(zero), len(positive), len(rays)
            )
        )
        processed.append(next_facet)
        next_facet += 1

    def snapshot(height: Silver) -> None:
        classified = _classify_rays(rays.values(), functionals)
        snapshots.append(
            SilverP7D24Snapshot(
                height,
                next_facet,
                len(classified),
                sum(ray.kind == "finite" for ray in classified),
                sum(ray.kind == "cusp" for ray in classified),
                sum(ray.kind == "hyperideal" for ray in classified),
            )
        )

    # Insert the rest of the simultaneous first packet before recording it.
    for _root in first_roots:
        if _root is affine[0]:
            continue
        insert_one()
    snapshot(nonempty[0].height)
    for packet in nonempty[1:]:
        for _root in packet.accepted_roots:
            insert_one()
        snapshot(packet.height)

    classified = _classify_rays(rays.values(), functionals)
    return SilverP7D24Geometry(
        packet_limit,
        tuple(snapshots),
        tuple(insertions),
        classified,
        _ray_digest(classified),
    )


@lru_cache(maxsize=None)
def certify_full_prefix_hyperideal_ray(
    search: SilverP7D24Search,
) -> SilverP7D24HyperidealWitness:
    """Certify non-finite volume without expanding a large later packet.

    A rank-24 extreme ray of the first complete packet remains extreme after
    refinement if it satisfies every later half-space.  A positive-norm such
    ray is therefore an exact hyperideal ray of the full enumerated prefix.
    """

    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    if not nonempty:
        raise ValueError("the prefix has no nonempty packet")
    first = nonempty[0]
    geometry = analyze_silver_p7_d24_prefix(search, 1)
    later_roots = tuple(
        root for packet in nonempty[1:] for root in packet.accepted_roots
    )
    survivors = []
    for ray in geometry.rays:
        if ray.kind != "hyperideal":
            continue
        pairings = tuple(
            silver_p7_d24_pair(ray.vector, root.vector) for root in later_roots
        )
        if all(value.sign(0) <= 0 for value in pairings):
            survivors.append((ray, pairings))
    if not survivors:
        raise AssertionError("no first-packet hyperideal ray survives the full prefix")

    ray, pairings = survivors[0]
    first_wall_count = 24 + len(first.accepted_roots)
    new_zero_facets = tuple(
        first_wall_count + index for index, value in enumerate(pairings) if not value
    )
    functionals = tuple(_functional(wall) for wall in search.walls)
    active = (*ray.active_facets, *new_zero_facets)
    final_rank = (
        24
        if _silver_rank_equals_certified(
            tuple(functionals[index] for index in active), 25, 24
        )
        else -1
    )
    if final_rank != 24 or ray.norm.sign(0) <= 0:
        raise AssertionError("the surviving ray is not rank-24 hyperideal")
    return SilverP7D24HyperidealWitness(
        ray,
        len(later_roots),
        pairings,
        new_zero_facets,
        final_rank,
        len(survivors),
    )


@lru_cache(maxsize=None)
def audit_frozen_h100_witness(
    search: SilverP7D24Search,
) -> SilverP7D24WitnessCutAudit:
    """Pair the frozen H=100 witness with every later nonempty packet."""

    vector = H100_FROZEN_WITNESS_VECTOR
    norm = silver_p7_d24_pair(vector, vector)
    if norm != Silver(48, -16) / 7 or norm.sign(0) <= 0:
        raise AssertionError("the frozen H=100 hyperideal witness changed")
    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    if not nonempty or nonempty[0].height != Silver(5, 4):
        raise AssertionError("the first prime-seven packet changed")
    output = []
    first_cut_height = None
    first_cut_indices: tuple[int, ...] = ()
    for packet in nonempty[1:]:
        pairings = tuple(
            silver_p7_d24_pair(vector, root.vector)
            for root in packet.accepted_roots
        )
        positive_indices = tuple(
            index for index, value in enumerate(pairings) if value.sign(0) > 0
        )
        if first_cut_height is None and positive_indices:
            first_cut_height = packet.height
            first_cut_indices = positive_indices
        output.append(
            SilverP7D24WitnessPacketAudit(
                packet.height,
                len(pairings),
                sum(value.sign(0) < 0 for value in pairings),
                sum(not value for value in pairings),
                len(positive_indices),
                min(pairings),
                max(pairings),
            )
        )
    return SilverP7D24WitnessCutAudit(
        vector,
        norm,
        tuple(output),
        first_cut_height,
        first_cut_indices,
    )


@lru_cache(maxsize=None)
def audit_first_packet_hyperideal_survival(
    search: SilverP7D24Search,
) -> tuple[SilverP7D24InitialRaySurvival, ...]:
    """Count first-packet hyperideal rays surviving cumulative later walls."""

    geometry = analyze_silver_p7_d24_prefix(search, 1)
    rays = tuple(ray for ray in geometry.rays if ray.kind == "hyperideal")
    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    later_roots: list[SilverP7D24Candidate] = []
    output = []
    for packet in nonempty[1:]:
        later_roots.extend(packet.accepted_roots)
        surviving = sum(
            all(
                silver_p7_d24_pair(ray.vector, root.vector).sign(0) <= 0
                for root in later_roots
            )
            for ray in rays
        )
        output.append(
            SilverP7D24InitialRaySurvival(
                packet.height, len(later_roots), surviving
            )
        )
    return tuple(output)


@lru_cache(maxsize=None)
def certify_h300_replacement_hyperideal_ray(
    search: SilverP7D24Search,
) -> SilverP7D24DirectPrefixWitness:
    """Certify a replacement hyperideal extreme ray of the H=300 prefix.

    This certificate does not rely on completing double description through
    the dense later packets.  The displayed vector is paired with every wall
    directly; its zero-wall functionals have exact rank 24.
    """

    if len(search.walls) != 348 or len(search.accepted_roots) != 324:
        raise ValueError("the direct witness is pinned to the 348-wall H=300 prefix")
    vector = H300_REPLACEMENT_WITNESS_VECTOR
    pairings = tuple(silver_p7_d24_pair(vector, wall) for wall in search.walls)
    if any(value.sign(0) > 0 for value in pairings):
        raise AssertionError("the H=300 replacement ray leaves the chamber")
    active = tuple(index for index, value in enumerate(pairings) if not value)
    functionals = tuple(_functional(wall) for wall in search.walls)
    active_rank = _silver_row_rank(
        tuple(functionals[index] for index in active), 25
    )
    norm = silver_p7_d24_pair(vector, vector)
    if active_rank != 24 or norm.sign(0) <= 0:
        raise AssertionError("the replacement ray is not rank-24 hyperideal")
    return SilverP7D24DirectPrefixWitness(
        vector,
        norm,
        len(search.walls),
        pairings,
        active,
        active_rank,
    )


def certify_h500_replacement_hyperideal_ray(
    walls: Sequence[Vector],
) -> SilverP7D24DirectPrefixWitness:
    """Certify the frozen rank-24 hyperideal ray of the H=500 prefix.

    The dense H=500 extension is generated from the independently frozen
    H=300 wall list by ``audit_silver_p7_d24_plus_h500.py``.  Accepting the
    wall sequence here keeps this small certificate independent of the
    accelerator used for the 4,003-root simultaneous-height packet.
    """

    if len(walls) != 4374:
        raise ValueError("the direct witness is pinned to the 4,374-wall H=500 prefix")
    vector = H500_REPLACEMENT_WITNESS_VECTOR
    pairings = tuple(silver_p7_d24_pair(vector, wall) for wall in walls)
    if any(value.sign(0) > 0 for value in pairings):
        raise AssertionError("the H=500 replacement ray leaves the chamber")
    active = tuple(index for index, value in enumerate(pairings) if not value)
    functionals = tuple(_functional(wall) for wall in walls)
    active_rank = _silver_row_rank(
        tuple(functionals[index] for index in active), 25
    )
    norm = silver_p7_d24_pair(vector, vector)
    if active_rank != 24 or norm.sign(0) <= 0:
        raise AssertionError("the H=500 replacement ray is not rank-24 hyperideal")
    return SilverP7D24DirectPrefixWitness(
        vector,
        norm,
        len(walls),
        pairings,
        active,
        active_rank,
    )


__all__ = [
    "COMPLETE_P7_D24_NORM_PROOF",
    "EPSILON",
    "H100_FROZEN_WITNESS_VECTOR",
    "H300_REPLACEMENT_WITNESS_VECTOR",
    "H500_REPLACEMENT_WITNESS_VECTOR",
    "P7_DELTA",
    "P7_TIMELIKE",
    "SILVER_P7_D24_SEED",
    "SilverP7D24Candidate",
    "SilverP7D24ControlLevel",
    "SilverP7D24DirectPrefixWitness",
    "SilverP7D24Geometry",
    "SilverP7D24HyperidealWitness",
    "SilverP7D24InitialRaySurvival",
    "SilverP7D24LevelEnumeration",
    "SilverP7D24Packet",
    "SilverP7D24Search",
    "SilverP7D24WitnessCutAudit",
    "SilverP7D24WitnessPacketAudit",
    "analyze_silver_p7_d24_prefix",
    "audit_first_packet_hyperideal_survival",
    "audit_frozen_h100_witness",
    "certify_full_prefix_hyperideal_ray",
    "certify_h300_replacement_hyperideal_ray",
    "certify_h500_replacement_hyperideal_ray",
    "enumerate_p7_d24_control_levels",
    "enumerate_p7_d24_level",
    "p7_factorization_audit",
    "run_silver_p7_d24_prefix",
    "silver_d24_height_zero_simple_roots",
    "silver_p7_d24_pair",
    "silver_p7_d24_vector_divisible_by_sqrt2",
]
