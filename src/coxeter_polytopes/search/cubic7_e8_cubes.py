"""Exact H24 Vinberg prefixes for ``<-a> + E8^3`` over cubic seven.

The positive block is unimodular.  Hence, for ``a=theta``, every primitive
pairing ideal is the unit ideal and reflection integrality leaves exactly the
norm representatives ``1,2``.  For ``a=2*theta`` (respectively
``theta*(2-theta)``), the same local argument permits one squarefree factor
at ``(2)`` (respectively ``(2-theta)``), giving the complete lists below.
Class number one and full unit-signature rank identify totally positive unit
associates modulo root rescaling by unit squares.

The search uses the exact cubic dominant-E8 enumerator implemented in
``cubic7_root_blocks``.  Three block squares are joined exactly, and Vinberg
packets of equal height are selected simultaneously.  For the first
``a=2*theta`` packet, the cone admits an especially small exact certificate:
in Dynkin-weight coordinates it is a product of three simplices.  The specialized
ray enumerator records every projective ray without a generic rank-25
support census.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import groupby, product
from typing import Any, Iterable, Sequence

from coxeter_polytopes.exact.cubic7 import (
    CUBIC7_ONE,
    CUBIC7_ZERO,
    THETA7,
    Cubic7,
)
from coxeter_polytopes.lattices import e8_cartan_gram
from coxeter_polytopes.search.cubic7_root_blocks import (
    DominantBlockWeight,
    _bounded_for_caps,
    _block_data,
    _mat_vec,
    enumerate_dominant_block_weights,
)


Vector = tuple[Cubic7, ...]
DELTA7 = 2 - THETA7


@dataclass(frozen=True, slots=True)
class Cubic7E8CubeSeed:
    seed_id: str
    timelike_coefficient: Cubic7
    canonical_mirror_norms: tuple[Cubic7, ...]
    complete_norm_class_proof: str


_BASE_PROOF = (
    "For primitive r=(x,y), reflection integrality is q(r)|2I(r), with "
    "I(r)=(a*x,G*y). The E8^3 block is unimodular. Away from primes "
    "dividing 2a, primitivity therefore makes I(r)=O. The inert prime (2) "
    "can occur with exponent at most one. Class number one and the full "
    "unit-signature map imply every totally positive unit is a square."
)


THETA_E8_CUBE_SEED = Cubic7E8CubeSeed(
    "H24-E8^3-a-theta",
    THETA7,
    (Cubic7(1), Cubic7(2)),
    _BASE_PROOF + " Since theta is a unit, (q) is (1) or (2).",
)

TWO_THETA_E8_CUBE_SEED = Cubic7E8CubeSeed(
    "H24-E8^3-a-2theta",
    2 * THETA7,
    (Cubic7(1), Cubic7(2)),
    _BASE_PROOF
    + " If y is divisible by (2), primitive x is a unit and "
    "q=-2*theta*x^2 mod (4), so v_2(q)=1. Thus (q) is (1) or (2).",
)

DELTA_THETA_E8_CUBE_SEED = Cubic7E8CubeSeed(
    "H24-E8^3-a-theta-delta",
    THETA7 * DELTA7,
    (Cubic7(1), Cubic7(2), DELTA7, 2 * DELTA7),
    _BASE_PROOF
    + " At p=(2-theta), either I(r)=O or y is in p and primitive x is "
    "a unit, giving q=-theta*(2-theta)*x^2 mod p^2 and v_p(q)=1. "
    "Thus (q) is (1),(2),p or 2p.",
)

SUPPORTED_E8_CUBE_SEEDS = (
    THETA_E8_CUBE_SEED,
    TWO_THETA_E8_CUBE_SEED,
    DELTA_THETA_E8_CUBE_SEED,
)


def e8_cube_pair(
    seed: Cubic7E8CubeSeed, left: Vector, right: Vector
) -> Cubic7:
    if len(left) != 25 or len(right) != 25:
        raise ValueError("an H24 E8-cube vector has rank 25")
    e8 = e8_cartan_gram()
    result = -seed.timelike_coefficient * left[0] * right[0]
    for block in range(3):
        offset = 1 + 8 * block
        result += sum(
            (
                left[offset + i]
                * Cubic7(e8[i][j])
                * right[offset + j]
                for i in range(8)
                for j in range(8)
                if e8[i][j]
            ),
            start=CUBIC7_ZERO,
        )
    return result


@dataclass(frozen=True, slots=True)
class H24ControlLevel:
    norm: Cubic7
    height: Cubic7
    control: Cubic7
    target: Cubic7


def _control_levels(
    seed: Cubic7E8CubeSeed, height_bound: int
) -> tuple[H24ControlLevel, ...]:
    if height_bound < 0:
        raise ValueError("height_bound must be nonnegative")
    a = seed.timelike_coefficient
    levels: set[H24ControlLevel] = set()
    for norm in seed.canonical_mirror_norms:
        caps = (norm * height_bound / a, norm / (-a), norm / (-a))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0 or not (2 * a * control / norm).is_integral:
                continue
            target = norm + a * control * control
            if any(target.sign(place) < 0 for place in (1, 2)):
                continue
            height = a * control * control / norm
            if (Cubic7(height_bound) - height).sign(0) < 0:
                continue
            levels.add(H24ControlLevel(norm, height, control, target))
    return tuple(sorted(levels, key=lambda item: (item.height, item.norm, item.control)))


@dataclass(frozen=True, slots=True)
class H24Candidate:
    norm: Cubic7
    height: Cubic7
    control: Cubic7
    target: Cubic7
    vector: Vector


@lru_cache(maxsize=None)
def enumerate_h24_control_level(
    seed: Cubic7E8CubeSeed, level: H24ControlLevel
) -> tuple[H24Candidate, ...]:
    e8 = enumerate_dominant_block_weights("E8", level.target, level.norm)
    buckets: dict[Cubic7, list[DominantBlockWeight]] = {}
    for vector in e8.vectors:
        buckets.setdefault(vector.square, []).append(vector)
    output: list[H24Candidate] = []
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
                if e8_cube_pair(seed, vector, vector) != level.norm:
                    raise AssertionError("E8-cube reconstruction has the wrong norm")
                output.append(
                    H24Candidate(
                        level.norm, level.height, level.control, level.target, vector
                    )
                )
    return tuple(sorted(output, key=lambda candidate: candidate.vector))


@lru_cache(maxsize=None)
def height_zero_dominant_roots(seed: Cubic7E8CubeSeed) -> tuple[H24Candidate, ...]:
    roots: list[H24Candidate] = []
    for norm in seed.canonical_mirror_norms:
        roots.extend(
            enumerate_h24_control_level(
                seed, H24ControlLevel(norm, CUBIC7_ZERO, CUBIC7_ZERO, norm)
            )
        )
    return tuple(sorted(roots, key=lambda root: (root.norm, root.vector)))


def height_zero_is_e8_cubed(seed: Cubic7E8CubeSeed) -> bool:
    roots = height_zero_dominant_roots(seed)
    return len(roots) == 3 and all(root.norm == 2 for root in roots)


@lru_cache(maxsize=1)
def h24_height_zero_simple_roots() -> tuple[Vector, ...]:
    roots: list[Vector] = []
    for index in range(24):
        root = [CUBIC7_ZERO] * 25
        root[index + 1] = -CUBIC7_ONE
        roots.append(tuple(root))
    return tuple(roots)


@dataclass(frozen=True, slots=True)
class H24VinbergPacket:
    height: Cubic7
    controls: tuple[H24ControlLevel, ...]
    candidate_count: int
    lower_wall_rejected_count: int
    roots: tuple[H24Candidate, ...]


@dataclass(frozen=True, slots=True)
class H24ExtremeRay:
    vector: Vector
    active_facets: tuple[int, ...]
    norm: Cubic7
    kind: str
    active_rank_verified: bool


@dataclass(frozen=True, slots=True)
class H24Geometry:
    engine: str
    exhaustive_argument: str
    extreme_rays: tuple[H24ExtremeRay, ...]
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    irredundant: bool
    finite_volume: bool
    compact: bool


def _e8_weights(coordinates: Sequence[Cubic7]) -> Vector:
    gram = e8_cartan_gram()
    return tuple(
        sum(
            (Cubic7(gram[i][j]) * coordinates[j] for j in range(8) if gram[i][j]),
            start=CUBIC7_ZERO,
        )
        for i in range(8)
    )


def _from_block_weights(control: Cubic7, block_weights: Sequence[Vector]) -> Vector:
    _gram, inverse, _schur = _block_data("E8")
    return (control, *tuple(entry for weights in block_weights for entry in _mat_vec(inverse, weights)))


def _structured_first_packet_geometry(
    seed: Cubic7E8CubeSeed,
    accepted: Sequence[H24Candidate],
    walls: Sequence[Vector],
) -> H24Geometry | None:
    """Recognize and exhaust the three affine-block inequalities exactly."""

    if len(accepted) != 3 or len(walls) != 27:
        return None
    controls = {root.control for root in accepted}
    if len(controls) != 1:
        return None
    control = next(iter(controls))
    special_data: list[tuple[int, Cubic7, Vector]] = []
    seen_blocks: set[int] = set()
    new_facet_for_block: dict[int, int] = {}
    for root_offset, root in enumerate(accepted):
        nonzero_blocks = []
        for block in range(3):
            coordinates = root.vector[1 + 8 * block : 1 + 8 * (block + 1)]
            weights = _e8_weights(coordinates)
            if any(weights):
                nonzero_blocks.append((block, weights, coordinates))
        if len(nonzero_blocks) != 1:
            return None
        block, weights, coordinates = nonzero_blocks[0]
        nonzero = [(index, value) for index, value in enumerate(weights) if value]
        if len(nonzero) != 1 or block in seen_blocks:
            return None
        seen_blocks.add(block)
        new_facet_for_block[block] = 24 + root_offset
        special_index, special_value = nonzero[0]
        highest = tuple(entry / special_value for entry in coordinates)
        if not all(entry.is_rational and entry.a > 0 for entry in highest):
            return None
        special_data.append((special_index, special_value, highest))
    if seen_blocks != {0, 1, 2}:
        return None
    special_indices = {index for index, _value, _highest in special_data}
    special_values = {value for _index, value, _highest in special_data}
    highest_vectors = {highest for _index, _value, highest in special_data}
    if len(special_indices) != 1 or len(special_values) != 1 or len(highest_vectors) != 1:
        return None
    special = next(iter(special_indices))
    scale = next(iter(special_values))
    highest = next(iter(highest_vectors))
    cap = seed.timelike_coefficient * control / scale
    if cap != 2:
        return None

    raw_rays: list[Vector] = []
    # In each block the new inequality is highest dot z <= 2t, so its
    # t=1 section is an 8-simplex with vertices 0 and
    # (2/highest_i)e_i.  The three blocks are independent at fixed t;
    # therefore the complete section is their Cartesian product (9^3
    # vertices), and its homogenization has exactly those projective rays.
    for choices in product(range(9), repeat=3):
        weights = [[CUBIC7_ZERO] * 8 for _ in range(3)]
        for block, choice in enumerate(choices):
            if choice:
                coordinate = choice - 1
                weights[block][coordinate] = cap / highest[coordinate]
        raw_rays.append(
            _from_block_weights(CUBIC7_ONE, tuple(tuple(row) for row in weights))
        )

    _gram, inverse, _schur = _block_data("E8")
    rays: list[H24ExtremeRay] = []
    all_choices = tuple(product(range(9), repeat=3))
    for choices, vector in zip(all_choices, raw_rays):
        active_list: list[int] = []
        spatial_square = CUBIC7_ZERO
        for block, choice in enumerate(choices):
            for coordinate in range(8):
                if choice != coordinate + 1:
                    active_list.append(8 * block + coordinate)
            if choice:
                coordinate = choice - 1
                active_list.append(new_facet_for_block[block])
                weight = cap / highest[coordinate]
                spatial_square += weight * Cubic7(inverse[coordinate][coordinate]) * weight
        active = tuple(sorted(active_list))
        expected_active: set[int] = set()
        for block, choice in enumerate(choices):
            expected_active.update(
                8 * block + coordinate
                for coordinate in range(8)
                if choice != coordinate + 1
            )
            if choice:
                expected_active.add(new_facet_for_block[block])
                if not highest[choice - 1]:
                    raise AssertionError("affine replacement has zero missing coefficient")
        active_rank_verified = len(active) == 24 and set(active) == expected_active
        if not active_rank_verified:
            raise AssertionError("product-simplex active basis has the wrong rank pattern")
        norm = -seed.timelike_coefficient + spatial_square
        sign = norm.sign(0)
        rays.append(
            H24ExtremeRay(
                vector,
                active,
                norm,
                "finite" if sign < 0 else "cusp" if sign == 0 else "hyperideal",
                True,
            )
        )
    ordered = tuple(sorted(rays, key=lambda ray: ray.vector))
    if len(ordered) != 9**3 or len({ray.vector for ray in ordered}) != 9**3:
        raise AssertionError("the product of three E8 simplices lost a projective ray")
    finite = sum(ray.kind == "finite" for ray in ordered)
    cusps = sum(ray.kind == "cusp" for ray in ordered)
    hyperideal = sum(ray.kind == "hyperideal" for ray in ordered)
    # Each simple facet is strict at the vertex using its coordinate axis;
    # each new facet is strict at a vertex using zero in that block.
    irredundant = True
    finite_volume = irredundant and hyperideal == 0
    return H24Geometry(
        engine="exact_dynkin_weight_product_of_simplices",
        exhaustive_argument=(
            "The simple E8 walls give z_(b,i)>=0. Each new root is scale "
            "times the highest E8 root and a*x/scale=2, hence its block "
            "inequality is highest dot z_b <= 2t. At t=1 each block is an "
            "8-simplex with vertices 0 and (2/highest_i)e_i. The whole "
            "section is the Cartesian product of three such simplices, so "
            "its 9^3 vertices are exactly all projective cone rays."
        ),
        extreme_rays=ordered,
        finite_vertex_count=finite,
        cusp_count=cusps,
        hyperideal_ray_count=hyperideal,
        irredundant=irredundant,
        finite_volume=finite_volume,
        compact=finite_volume and cusps == 0,
    )


@dataclass(frozen=True, slots=True)
class H24E8CubeSearch:
    seed: Cubic7E8CubeSeed
    height_bound: int
    control_levels: tuple[H24ControlLevel, ...]
    packets: tuple[H24VinbergPacket, ...]
    cone_roots: tuple[Vector, ...]
    accepted_roots: tuple[H24Candidate, ...]
    rejected_roots: tuple[H24Candidate, ...]
    walls: tuple[Vector, ...]
    geometry: H24Geometry | None


@lru_cache(maxsize=None)
def run_h24_e8_cube_vinberg(
    seed: Cubic7E8CubeSeed, height_bound: int, *, analyze_geometry: bool = True
) -> H24E8CubeSearch:
    if seed not in SUPPORTED_E8_CUBE_SEEDS:
        raise ValueError("unsupported E8-cube seed")
    if not height_zero_is_e8_cubed(seed):
        raise AssertionError("height-zero dominant census found an extra mirror orbit")
    controls = _control_levels(seed, height_bound)
    walls = list(h24_height_zero_simple_roots())
    accepted: list[H24Candidate] = []
    rejected: list[H24Candidate] = []
    packets: list[H24VinbergPacket] = []
    for height, iterator in groupby(controls, key=lambda item: item.height):
        packet_controls = tuple(iterator)
        candidates = tuple(
            root
            for level in packet_controls
            for root in enumerate_h24_control_level(seed, level)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(e8_cube_pair(seed, root.vector, wall).sign(0) <= 0 for wall in walls)
        )
        if any(
            e8_cube_pair(seed, left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height H24 survivors are not mutually acute")
        survivor_set = set(survivors)
        accepted.extend(survivors)
        rejected.extend(root for root in candidates if root not in survivor_set)
        walls.extend(root.vector for root in survivors)
        packets.append(
            H24VinbergPacket(
                height,
                packet_controls,
                len(candidates),
                len(candidates) - len(survivors),
                survivors,
            )
        )
    geometry = (
        _structured_first_packet_geometry(seed, accepted, walls)
        if analyze_geometry
        else None
    )
    return H24E8CubeSearch(
        seed,
        height_bound,
        controls,
        tuple(packets),
        h24_height_zero_simple_roots(),
        tuple(accepted),
        tuple(rejected),
        tuple(walls),
        geometry,
    )


def _cubic(value: Cubic7) -> dict[str, str]:
    return value.to_json()


def _vector(vector: Iterable[Cubic7]) -> list[dict[str, str]]:
    return [_cubic(entry) for entry in vector]


def build_h24_e8_cube_transcript(
    seed: Cubic7E8CubeSeed, height_bound: int, *, analyze_geometry: bool = True
) -> dict[str, Any]:
    search = run_h24_e8_cube_vinberg(seed, height_bound, analyze_geometry=analyze_geometry)
    geometry = search.geometry
    return {
        "schema_version": "0.1.0",
        "transcript_kind": "cubic7_h24_e8_cubed_vinberg_prefix",
        "seed_id": seed.seed_id,
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "lattice": "<-a> + E8^3",
        "hyperbolic_dimension": 24,
        "timelike_coefficient": _cubic(seed.timelike_coefficient),
        "height_bound": height_bound,
        "canonical_mirror_norms": [_cubic(norm) for norm in seed.canonical_mirror_norms],
        "complete_norm_class_proof": seed.complete_norm_class_proof,
        "height_zero": {
            "type": "E8^3",
            "all_root_count": 720,
            "dominant_orbit_count": len(height_zero_dominant_roots(seed)),
            "dominant_orbits": [
                {"norm": _cubic(root.norm), "vector": _vector(root.vector)}
                for root in height_zero_dominant_roots(seed)
            ],
            "simple_wall_count": len(search.cone_roots),
        },
        "control_level_count": len(search.control_levels),
        "packets": [
            {
                "height": _cubic(packet.height),
                "controls": [
                    {
                        "norm": _cubic(level.norm),
                        "control": _cubic(level.control),
                        "target": _cubic(level.target),
                    }
                    for level in packet.controls
                ],
                "candidate_count": packet.candidate_count,
                "lower_wall_rejected_count": packet.lower_wall_rejected_count,
                "accepted_root_count": len(packet.roots),
                "accepted_roots": [
                    {"norm": _cubic(root.norm), "vector": _vector(root.vector)}
                    for root in packet.roots
                ],
            }
            for packet in search.packets
        ],
        "output": {
            "accepted_positive_root_count": len(search.accepted_roots),
            "rejected_positive_root_count": len(search.rejected_roots),
            "wall_count": len(search.walls),
            "walls": [_vector(wall) for wall in search.walls],
        },
        "geometry": None
        if geometry is None
        else {
            "engine": geometry.engine,
            "exhaustive_argument": geometry.exhaustive_argument,
            "extreme_ray_count": len(geometry.extreme_rays),
            "finite_vertex_count": geometry.finite_vertex_count,
            "cusp_count": geometry.cusp_count,
            "hyperideal_ray_count": geometry.hyperideal_ray_count,
            "irredundant": geometry.irredundant,
            "finite_volume": geometry.finite_volume,
            "compact": geometry.compact,
            "rays": [
                {
                    "vector": _vector(ray.vector),
                    "active_facets": list(ray.active_facets),
                    "norm": _cubic(ray.norm),
                    "kind": ray.kind,
                    "active_rank_verified": ray.active_rank_verified,
                }
                for ray in geometry.extreme_rays
            ],
        },
        "scope": (
            "The norm classes and every Vinberg level through the stated exact "
            "height are exhaustive. A hyperideal prefix is not a nonreflectivity "
            "proof because later walls may remove its rays."
        ),
    }


__all__ = [
    "DELTA_THETA_E8_CUBE_SEED",
    "SUPPORTED_E8_CUBE_SEEDS",
    "THETA_E8_CUBE_SEED",
    "TWO_THETA_E8_CUBE_SEED",
    "build_h24_e8_cube_transcript",
    "e8_cube_pair",
    "enumerate_h24_control_level",
    "height_zero_dominant_roots",
    "height_zero_is_e8_cubed",
    "run_h24_e8_cube_vinberg",
]
