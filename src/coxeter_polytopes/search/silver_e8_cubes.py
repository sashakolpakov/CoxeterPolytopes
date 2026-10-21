"""Exact H24 Vinberg prefixes for ``<-a> + E8^3`` over ``Z[sqrt(2)]``.

Only coefficients ``a`` with signs ``(+,-)`` are admissible: ``-a`` is then
timelike at the distinguished embedding and positive at the conjugate one.
Thus ``sqrt(2)`` and ``2sqrt(2)`` are supported, whereas ``2+sqrt(2)`` is
explicitly rejected because it is positive at both embeddings.

The complete primitive mirror-norm classes follow from the pairing ideal
``I(r)=(a*x,G*y)`` and the ramified uniformizer ``pi=sqrt(2)``, with
``(2)=(pi)^2``.  The even unimodular ``E8^3`` block makes every spatial
square divisible by ``pi^2``.  For ``a=pi`` the only possible norm
valuations are one and two; for ``a=pi^3`` they are two and three.  The unit
``epsilon=1+pi`` has signs ``(+,-)``, so the totally-positive canonical
generators are respectively ``epsilon*pi=2+pi, pi^2=2`` and
``pi^2=2, epsilon*pi^3=4+2pi``.  Since the totally-positive units are exactly
``epsilon^(2Z)``, these lists are complete modulo unit-square rescaling.

Positive-height roots are enumerated modulo the full height-zero ``E8^3``
Weyl group.  Dynkin pairings in each block are dominant, exact Schur bounds
are imposed at both embeddings, and equal-height Vinberg packets are selected
simultaneously.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby, product
from math import isqrt
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import (
    SILVER_ONE,
    SILVER_ZERO,
    SQRT2,
    Silver,
)
from coxeter_polytopes.lattices import e8_cartan_gram
from coxeter_polytopes.search.cubic7_root_blocks import _block_data


Vector: TypeAlias = tuple[Silver, ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]
EPSILON = SILVER_ONE + SQRT2


@dataclass(frozen=True, slots=True)
class SilverE8CubeSeed:
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


SQRT2_E8_CUBE_SEED = SilverE8CubeSeed(
    seed_id="silver-H24-E8^3-a-sqrt2",
    timelike_coefficient=SQRT2,
    canonical_mirror_norms=(Silver(2), Silver(2, 1)),
    complete_norm_class_proof=(
        "Put pi=sqrt(2), so (2)=(pi)^2. For primitive r=(x,y), "
        "reflection integrality is (q(r))|2I(r), I(r)=(pi*x,G*y). "
        "The E8^3 block is even unimodular. Away from (pi), local "
        "unimodularity and primitivity make I(r)=O. If x is a unit, "
        "q=-pi*x^2 mod pi^2 and v_pi(q)=1. If x is not a unit, some "
        "spatial coordinate is a unit, I(r)=O, and evenness gives "
        "v_pi(q)>=2; divisibility by 2 forces v_pi(q)=2. Thus the norm "
        "ideals are (pi),(pi^2). The unit epsilon=1+pi has signs (+,-), "
        "and O*=+/-epsilon^Z, so every totally-positive unit is a square. "
        "The totally-positive representatives are epsilon*pi=2+pi and "
        "pi^2=2."
    ),
)


TWO_SQRT2_E8_CUBE_SEED = SilverE8CubeSeed(
    seed_id="silver-H24-E8^3-a-2sqrt2",
    timelike_coefficient=2 * SQRT2,
    canonical_mirror_norms=(Silver(2), Silver(4, 2)),
    complete_norm_class_proof=(
        "Put pi=sqrt(2), so a=2sqrt(2)=pi^3 and (2)=(pi)^2. For "
        "primitive r=(x,y), reflection integrality is (q(r))|2I(r), "
        "I(r)=(pi^3*x,G*y). The E8^3 block is even unimodular. Away "
        "from (pi), local unimodularity and primitivity make I(r)=O. "
        "If some spatial coordinate is a unit, I(r)=O and q is divisible "
        "by pi^2; integrality forces v_pi(q)=2. Otherwise x is a unit "
        "and y is divisible by pi, so q=-pi^3*x^2 mod pi^4 and "
        "v_pi(q)=3. Thus the norm ideals are (pi^2),(pi^3). The unit "
        "epsilon=1+pi has signs (+,-), and O*=+/-epsilon^Z, so every "
        "totally-positive unit is a square. The totally-positive "
        "representatives are pi^2=2 and epsilon*pi^3=4+2pi."
    ),
)


INADMISSIBLE_TWO_PLUS_SQRT2 = Silver(2, 1)
SUPPORTED_SILVER_E8_CUBE_SEEDS = (
    SQRT2_E8_CUBE_SEED,
    TWO_SQRT2_E8_CUBE_SEED,
)


def coefficient_is_arithmetically_admissible(value: Silver) -> bool:
    """Whether ``<-value> + E8^3`` is hyperbolic at exactly one place."""

    coefficient = Silver.coerce(value)
    return coefficient.sign(0) > 0 and coefficient.sign(1) < 0


def _ceiling_square_root(value: Fraction) -> int:
    if value < 0:
        raise ValueError("a square cap cannot be negative")
    floor = isqrt(value.numerator // value.denominator)
    while Fraction(floor * floor) < value:
        floor += 1
    return floor


@lru_cache(maxsize=None)
def bounded_silver_integers(square_caps: tuple[Fraction, Fraction]) -> tuple[Silver, ...]:
    """Enumerate all ``m+n*sqrt(2)`` inside an exact Minkowski rectangle."""

    caps = tuple(Fraction(cap) for cap in square_caps)
    if len(caps) != 2 or any(cap < 0 for cap in caps):
        return ()
    embedding_bounds = tuple(_ceiling_square_root(cap) for cap in caps)
    # If s0=m+n*pi and s1=m-n*pi, then 2m=s0+s1 and
    # 2n*pi=s0-s1.  The same ceil((B0+B1)/2) is a certified (slightly
    # slack) integral bound for both m and n because pi>1.
    coefficient_bound = (sum(embedding_bounds) + 1) // 2
    values = []
    for rational in range(-coefficient_bound, coefficient_bound + 1):
        for radical in range(-coefficient_bound, coefficient_bound + 1):
            value = Silver(rational, radical)
            square = value * value
            if all(
                (Silver(caps[place]) - square).sign(place) >= 0
                for place in (0, 1)
            ):
                values.append(value)
    return tuple(sorted(values))


def _embedding_integer_upper(value: Silver, embedding: int) -> int:
    if value.sign(embedding) < 0:
        raise ValueError("the selected embedding cap is negative")
    upper = 0 if not value else 1
    while (Silver(upper) - value).sign(embedding) < 0:
        upper *= 2
    return upper


def _bounded_for_caps(caps: tuple[Silver, Silver]) -> tuple[Silver, ...]:
    return bounded_silver_integers(
        tuple(
            Fraction(_embedding_integer_upper(cap, place))
            for place, cap in enumerate(caps)
        )
    )


def _quadratic(vector: Sequence[Silver], matrix: RationalMatrix) -> Silver:
    return sum(
        (
            vector[row] * Silver(matrix[row][column]) * vector[column]
            for row in range(len(vector))
            for column in range(len(vector))
            if matrix[row][column]
        ),
        start=SILVER_ZERO,
    )


def _mat_vec(matrix: RationalMatrix, vector: Sequence[Silver]) -> Vector:
    return tuple(
        sum(
            (
                Silver(matrix[row][column]) * vector[column]
                for column in range(len(vector))
                if matrix[row][column]
            ),
            start=SILVER_ZERO,
        )
        for row in range(len(matrix))
    )


@dataclass(frozen=True, slots=True)
class SilverDominantE8Weight:
    square: Silver
    weights: Vector
    coordinates: Vector


@dataclass(frozen=True, slots=True)
class SilverDominantE8Enumeration:
    target: Silver
    root_norm: Silver
    alphabet: tuple[Silver, ...]
    visited_by_depth: tuple[int, ...]
    schur_pruned_by_depth: tuple[int, ...]
    identity_pruned_by_depth: tuple[int, ...]
    vectors: tuple[SilverDominantE8Weight, ...]


@lru_cache(maxsize=None)
def enumerate_dominant_e8_weights(
    target: Silver, root_norm: Silver
) -> SilverDominantE8Enumeration:
    """Exhaust integral dominant E8 vectors with square at most ``target``."""

    _gram, inverse, schur = _block_data("E8")
    twice_target = 2 * target
    alphabet = tuple(
        value
        for value in _bounded_for_caps((twice_target, twice_target))
        if value.sign(0) >= 0
        and all((twice_target - value * value).sign(place) >= 0 for place in (0, 1))
        and (2 * value / root_norm).is_integral
    )
    visited = [0] * 9
    schur_pruned = [0] * 9
    identity_pruned = [0] * 9
    output: list[SilverDominantE8Weight] = []

    def extend(prefix: Vector) -> None:
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            minimum = _quadratic(prefix, schur[depth - 1])
            if any((target - minimum).sign(place) < 0 for place in (0, 1)):
                schur_pruned[depth] += 1
                return
            prefix_energy = _quadratic(
                prefix,
                tuple(tuple(inverse[i][j] for j in range(depth)) for i in range(depth)),
            )
            if (target - prefix_energy).sign(0) < 0:
                identity_pruned[depth] += 1
                return
        if depth == 8:
            square = _quadratic(prefix, inverse)
            coordinates = _mat_vec(inverse, prefix)
            if not all(entry.is_integral for entry in coordinates):
                raise AssertionError("the unimodular E8 inverse lost integrality")
            if any((target - square).sign(place) < 0 for place in (0, 1)):
                raise AssertionError("Schur pruning admitted a block above target")
            output.append(SilverDominantE8Weight(square, prefix, coordinates))
            return
        for entry in alphabet:
            extend((*prefix, entry))

    extend(())
    vectors = tuple(sorted(output, key=lambda item: (item.square, item.weights)))
    return SilverDominantE8Enumeration(
        target,
        root_norm,
        alphabet,
        tuple(visited),
        tuple(schur_pruned),
        tuple(identity_pruned),
        vectors,
    )


def silver_e8_cube_pair(
    seed: SilverE8CubeSeed, left: Sequence[Silver], right: Sequence[Silver]
) -> Silver:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("a silver H24 E8-cube vector has rank 25")
    gram = e8_cartan_gram()
    result = -seed.timelike_coefficient * left[0] * right[0]
    for block in range(3):
        offset = 1 + 8 * block
        result += sum(
            (
                left[offset + row]
                * Silver(gram[row][column])
                * right[offset + column]
                for row in range(8)
                for column in range(8)
                if gram[row][column]
            ),
            start=SILVER_ZERO,
        )
    return result


@dataclass(frozen=True, slots=True)
class SilverH24ControlLevel:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_control_levels(
    seed: SilverE8CubeSeed, height_bound: int
) -> tuple[SilverH24ControlLevel, ...]:
    if seed not in SUPPORTED_SILVER_E8_CUBE_SEEDS or not seed.admissible:
        raise ValueError("unsupported or arithmetically inadmissible silver seed")
    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    levels: set[SilverH24ControlLevel] = set()
    coefficient = seed.timelike_coefficient
    for norm in seed.canonical_mirror_norms:
        caps = (norm * height_bound / coefficient, norm / (-coefficient))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0 or not (2 * coefficient * control / norm).is_integral:
                continue
            target = norm + coefficient * control * control
            if target.sign(1) < 0:
                continue
            height = coefficient * control * control / norm
            if (Silver(height_bound) - height).sign(0) < 0:
                continue
            levels.add(SilverH24ControlLevel(norm, height, control, target))
    return tuple(sorted(levels, key=lambda item: (item.height, item.norm, item.control)))


@dataclass(frozen=True, slots=True)
class SilverH24Candidate:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver
    vector: Vector


@lru_cache(maxsize=None)
def enumerate_control_level(
    seed: SilverE8CubeSeed, level: SilverH24ControlLevel
) -> tuple[SilverH24Candidate, ...]:
    e8 = enumerate_dominant_e8_weights(level.target, level.norm)
    buckets: dict[Silver, list[SilverDominantE8Weight]] = {}
    for vector in e8.vectors:
        buckets.setdefault(vector.square, []).append(vector)
    candidates: list[SilverH24Candidate] = []
    for first in e8.vectors:
        for second in e8.vectors:
            remainder = level.target - first.square - second.square
            for third in buckets.get(remainder, ()):
                vector = (
                    level.control,
                    *first.coordinates,
                    *second.coordinates,
                    *third.coordinates,
                )
                if silver_e8_cube_pair(seed, vector, vector) != level.norm:
                    raise AssertionError("E8-cube reconstruction returned the wrong norm")
                candidates.append(
                    SilverH24Candidate(
                        level.norm,
                        level.height,
                        level.control,
                        level.target,
                        vector,
                    )
                )
    return tuple(sorted(candidates, key=lambda root: root.vector))


@lru_cache(maxsize=None)
def height_zero_dominant_roots(
    seed: SilverE8CubeSeed,
) -> tuple[SilverH24Candidate, ...]:
    roots = tuple(
        root
        for norm in seed.canonical_mirror_norms
        for root in enumerate_control_level(
            seed, SilverH24ControlLevel(norm, SILVER_ZERO, SILVER_ZERO, norm)
        )
    )
    return tuple(sorted(roots, key=lambda root: (root.norm, root.vector)))


def height_zero_is_e8_cubed(seed: SilverE8CubeSeed) -> bool:
    roots = height_zero_dominant_roots(seed)
    return len(roots) == 3 and all(root.norm == 2 for root in roots)


@lru_cache(maxsize=1)
def height_zero_simple_roots() -> tuple[Vector, ...]:
    roots = []
    for index in range(24):
        root = [SILVER_ZERO] * 25
        root[index + 1] = -SILVER_ONE
        roots.append(tuple(root))
    return tuple(roots)


@dataclass(frozen=True, slots=True)
class SilverH24VinbergPacket:
    height: Silver
    controls: tuple[SilverH24ControlLevel, ...]
    candidate_count: int
    lower_wall_rejected_count: int
    roots: tuple[SilverH24Candidate, ...]


@dataclass(frozen=True, slots=True)
class SilverH24E8CubeSearch:
    seed: SilverE8CubeSeed
    height_bound: int
    control_levels: tuple[SilverH24ControlLevel, ...]
    packets: tuple[SilverH24VinbergPacket, ...]
    cone_roots: tuple[Vector, ...]
    accepted_roots: tuple[SilverH24Candidate, ...]
    rejected_roots: tuple[SilverH24Candidate, ...]
    walls: tuple[Vector, ...]


@dataclass(frozen=True, slots=True)
class SilverH24FirstPacketGeometry:
    cap: Silver
    highest_root: tuple[Silver, ...]
    extreme_ray_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    finite_volume: bool
    compact: bool
    exhaustive_argument: str


@lru_cache(maxsize=None)
def run_silver_h24_e8_cube_vinberg(
    seed: SilverE8CubeSeed, height_bound: int
) -> SilverH24E8CubeSearch:
    """Exhaust and select all silver H24 walls through exact height ``H``."""

    if seed not in SUPPORTED_SILVER_E8_CUBE_SEEDS or not seed.admissible:
        raise ValueError("unsupported or arithmetically inadmissible silver seed")
    if not height_zero_is_e8_cubed(seed):
        raise AssertionError("height-zero dominant census is larger than E8^3")
    levels = enumerate_control_levels(seed, height_bound)
    walls = list(height_zero_simple_roots())
    accepted: list[SilverH24Candidate] = []
    rejected: list[SilverH24Candidate] = []
    packets: list[SilverH24VinbergPacket] = []
    for height, iterator in groupby(levels, key=lambda item: item.height):
        controls = tuple(iterator)
        candidates = tuple(
            root for level in controls for root in enumerate_control_level(seed, level)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(
                silver_e8_cube_pair(seed, root.vector, wall).sign(0) <= 0
                for wall in walls
            )
        )
        if any(
            silver_e8_cube_pair(seed, left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height silver H24 survivors are not mutually acute")
        survivor_set = set(survivors)
        accepted.extend(survivors)
        rejected.extend(root for root in candidates if root not in survivor_set)
        walls.extend(root.vector for root in survivors)
        packets.append(
            SilverH24VinbergPacket(
                height,
                controls,
                len(candidates),
                len(candidates) - len(survivors),
                survivors,
            )
        )
    return SilverH24E8CubeSearch(
        seed,
        height_bound,
        levels,
        tuple(packets),
        height_zero_simple_roots(),
        tuple(accepted),
        tuple(rejected),
        tuple(walls),
    )


def first_packet_product_geometry(
    search: SilverH24E8CubeSearch,
) -> SilverH24FirstPacketGeometry:
    """Certify the first three caps as a product of three E8 simplices."""

    accepted = search.accepted_roots
    if len(accepted) != 3 or len(search.walls) != 27:
        raise ValueError("the current wall set is not exactly the first three-wall packet")
    gram = e8_cartan_gram()
    _block_gram, inverse, _schur = _block_data("E8")
    controls = {root.control for root in accepted}
    if len(controls) != 1:
        raise AssertionError("first-packet controls differ")
    control = next(iter(controls))
    data: list[tuple[int, Silver, Vector]] = []
    for root in accepted:
        nonzero_blocks = []
        for block in range(3):
            coordinates = root.vector[1 + 8 * block : 1 + 8 * (block + 1)]
            weights = tuple(
                sum(
                    (
                        Silver(gram[row][column]) * coordinates[column]
                        for column in range(8)
                        if gram[row][column]
                    ),
                    start=SILVER_ZERO,
                )
                for row in range(8)
            )
            if any(weights):
                nonzero_blocks.append((block, weights, coordinates))
        if len(nonzero_blocks) != 1:
            raise AssertionError("a first-packet wall meets more than one E8 block")
        block, weights, coordinates = nonzero_blocks[0]
        nonzero = tuple((index, value) for index, value in enumerate(weights) if value)
        if len(nonzero) != 1:
            raise AssertionError("a first-packet wall is not a highest-root cap")
        _index, scale = nonzero[0]
        data.append((block, scale, tuple(entry / scale for entry in coordinates)))
    if {block for block, _scale, _highest in data} != {0, 1, 2}:
        raise AssertionError("the first packet does not cap all three E8 factors")
    scales = {scale for _block, scale, _highest in data}
    highest_roots = {highest for _block, _scale, highest in data}
    if len(scales) != 1 or len(highest_roots) != 1:
        raise AssertionError("the three E8 caps are not Weyl-identical")
    scale = next(iter(scales))
    highest = next(iter(highest_roots))
    cap = search.seed.timelike_coefficient * control / scale
    axis_energies = [SILVER_ZERO]
    for index in range(8):
        weight = cap / highest[index]
        axis_energies.append(weight * Silver(inverse[index][index]) * weight)
    finite = cusps = hyperideal = 0
    for choices in product(range(9), repeat=3):
        norm = -search.seed.timelike_coefficient + sum(
            (axis_energies[choice] for choice in choices), start=SILVER_ZERO
        )
        sign = norm.sign(0)
        if sign < 0:
            finite += 1
        elif sign == 0:
            cusps += 1
        else:
            hyperideal += 1
    return SilverH24FirstPacketGeometry(
        cap,
        highest,
        9**3,
        finite,
        cusps,
        hyperideal,
        hyperideal == 0,
        hyperideal == 0 and cusps == 0,
        (
            "In dominant Dynkin weights each E8 factor has z_i>=0 and one "
            "highest-root inequality highest.z<=cap*t. Its t=1 section is "
            "an 8-simplex. The three inequalities use disjoint blocks, so "
            "the whole section is their Cartesian product and its 9^3 "
            "vertices exhaust all projective cone rays."
        ),
    )


__all__ = [
    "EPSILON",
    "INADMISSIBLE_TWO_PLUS_SQRT2",
    "SQRT2_E8_CUBE_SEED",
    "SUPPORTED_SILVER_E8_CUBE_SEEDS",
    "TWO_SQRT2_E8_CUBE_SEED",
    "SilverDominantE8Enumeration",
    "SilverDominantE8Weight",
    "SilverE8CubeSeed",
    "SilverH24Candidate",
    "SilverH24ControlLevel",
    "SilverH24E8CubeSearch",
    "SilverH24FirstPacketGeometry",
    "SilverH24VinbergPacket",
    "bounded_silver_integers",
    "coefficient_is_arithmetically_admissible",
    "enumerate_control_level",
    "enumerate_control_levels",
    "enumerate_dominant_e8_weights",
    "first_packet_product_geometry",
    "height_zero_dominant_roots",
    "height_zero_is_e8_cubed",
    "height_zero_simple_roots",
    "run_silver_h24_e8_cube_vinberg",
    "silver_e8_cube_pair",
]
