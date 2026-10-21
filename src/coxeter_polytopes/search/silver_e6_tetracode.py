"""Exact silver Vinberg prefix for the Niemeier lattice ``N(E6^4)``.

The spatial lattice is constructed inside ``(E6^*)^4`` by the ternary
tetracode.  If ``z=C_E6*y`` is the Dynkin-weight vector of an E6 weight,
its discriminant residue is

``rho(z) = (1,2,0,1,2,0).z mod 3``.

For ``y=u+sqrt(2)*v``, membership in ``N(E6^4) tensor Z[sqrt(2)]`` means
that the four residues of both ``u`` and ``v`` are tetracode words.  This
module imposes that condition exactly; it never replaces the Niemeier
lattice by the root lattice ``E6^4``.

The supported Lorentzian seed is

``<-(1+sqrt(2))> + N(E6^4)``.

The coefficient is a unit, the spatial lattice is even unimodular, and the
complete primitive mirror-norm representatives are therefore ``1`` and
``2``.  Exact two-embedding Schur enumeration gives four factorwise affine
walls in the first nonempty packet and 44 walls in the next nonempty packet.
Twenty-four of the latter have nonzero tetracode residue and hence genuinely
couple three E6 factors.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby, product
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import (
    SILVER_ONE,
    SILVER_ZERO,
    SQRT2,
    Silver,
)
from coxeter_polytopes.search.bugaenko_vinberg import (
    _fraction_inverse,
    _prefix_schur_complements,
)
from coxeter_polytopes.search.silver_e8_cubes import (
    _bounded_for_caps,
    _mat_vec,
    _quadratic,
)


Vector: TypeAlias = tuple[Silver, ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]

EPSILON = SILVER_ONE + SQRT2
UNIT_E6_TETRACODE_NORMS = (SILVER_ONE, Silver(2))
E6_DISCRIMINANT_FUNCTIONAL = (1, 2, 0, 1, 2, 0)
E6_HIGHEST_ROOT = (1, 2, 3, 2, 1, 2)
SECOND_PACKET_HEIGHT = Silver(7, 5)
HYPERIDEAL_WITNESS_WEIGHT = Fraction(3, 5)

COMPLETE_UNIT_NORM_PROOF = (
    "Put pi=sqrt(2) and epsilon=1+pi.  The coefficient epsilon is a unit "
    "and N(E6^4) is even unimodular.  For a primitive r=(x,y), reflection "
    "integrality says q(r) divides 2I(r), where I(r)=(epsilon*x,(y,N)). "
    "Away from (pi), local unimodularity and primitivity make I(r) the unit "
    "ideal.  At (pi), a unit control gives v_pi(q)=0; otherwise a primitive "
    "spatial coordinate makes I(r) a unit, while evenness and q|2I force "
    "v_pi(q)=2.  Since every totally-positive unit of Z[sqrt(2)] is a "
    "square, the complete canonical positive norms are 1 and 2."
)


@lru_cache(maxsize=1)
def e6_cartan_gram() -> RationalMatrix:
    """Return E6 with arms of lengths two, two, and one at vertex two."""

    size = 6
    gram = [[Fraction(2 if i == j else 0) for j in range(size)] for i in range(size)]
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 4), (2, 5)):
        gram[i][j] = gram[j][i] = Fraction(-1)
    return tuple(tuple(row) for row in gram)


@lru_cache(maxsize=1)
def e6_inverse_cartan() -> RationalMatrix:
    return _fraction_inverse(e6_cartan_gram())


@lru_cache(maxsize=1)
def e6_prefix_schur_complements() -> tuple[RationalMatrix, ...]:
    return _prefix_schur_complements(e6_inverse_cartan())


@lru_cache(maxsize=1)
def ternary_tetracode() -> frozenset[tuple[int, int, int, int]]:
    """The self-dual ternary ``[4,2,3]`` tetracode.

    The chosen generators are ``(1,0,1,1)`` and ``(0,1,1,-1)``.  Changing
    signs or permuting the four E6 factors gives an equivalent Niemeier glue.
    """

    return frozenset(
        (left, right, (left + right) % 3, (left - right) % 3)
        for left in range(3)
        for right in range(3)
    )


def e6_discriminant_residue(weights: Sequence[Silver]) -> tuple[int, int]:
    """Return the rational and radical E6 discriminant residues."""

    if len(weights) != 6 or not all(value.is_integral for value in weights):
        raise ValueError("E6 Dynkin weights must be six silver integers")
    rational = sum(
        coefficient * int(value.a)
        for coefficient, value in zip(E6_DISCRIMINANT_FUNCTIONAL, weights)
    )
    radical = sum(
        coefficient * int(value.b)
        for coefficient, value in zip(E6_DISCRIMINANT_FUNCTIONAL, weights)
    )
    return rational % 3, radical % 3


def tetracode_membership(block_weights: Sequence[Sequence[Silver]]) -> bool:
    """Test membership in the full silver extension of ``N(E6^4)``."""

    if len(block_weights) != 4:
        return False
    try:
        residues = tuple(e6_discriminant_residue(block) for block in block_weights)
    except ValueError:
        return False
    code = ternary_tetracode()
    return (
        tuple(residue[0] for residue in residues) in code
        and tuple(residue[1] for residue in residues) in code
    )


@dataclass(frozen=True, slots=True)
class SilverE6Weight:
    square: Silver
    weights: Vector
    coordinates: Vector
    residue: tuple[int, int]


@dataclass(frozen=True, slots=True)
class SilverE6Enumeration:
    target: Silver
    root_norm: Silver
    alphabet: tuple[Silver, ...]
    visited_by_depth: tuple[int, ...]
    pruned_by_depth: tuple[int, ...]
    vectors: tuple[SilverE6Weight, ...]


@lru_cache(maxsize=None)
def enumerate_dominant_e6_weights(
    target: Silver, root_norm: Silver
) -> SilverE6Enumeration:
    """Exhaust dominant E6 weights below ``target`` at both embeddings."""

    inverse = e6_inverse_cartan()
    schur = e6_prefix_schur_complements()
    twice_target = 2 * target
    alphabet = tuple(
        value
        for value in _bounded_for_caps((twice_target, twice_target))
        if value.sign(0) >= 0
        and all((twice_target - value * value).sign(place) >= 0 for place in (0, 1))
        and (2 * value / root_norm).is_integral
    )
    visited = [0] * 7
    pruned = [0] * 7
    output: list[SilverE6Weight] = []

    def extend(prefix: Vector) -> None:
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            minimum = _quadratic(prefix, schur[depth - 1])
            if any((target - minimum).sign(place) < 0 for place in (0, 1)):
                pruned[depth] += 1
                return
        if depth == 6:
            square = _quadratic(prefix, inverse)
            coordinates = _mat_vec(inverse, prefix)
            output.append(
                SilverE6Weight(
                    square,
                    prefix,
                    coordinates,
                    e6_discriminant_residue(prefix),
                )
            )
            return
        for value in alphabet:
            extend((*prefix, value))

    extend(())
    return SilverE6Enumeration(
        target,
        root_norm,
        alphabet,
        tuple(visited),
        tuple(pruned),
        tuple(sorted(output, key=lambda item: (item.square, item.weights))),
    )


@dataclass(frozen=True, slots=True)
class SilverE6ControlLevel:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_control_levels(height_bound: int = 15) -> tuple[SilverE6ControlLevel, ...]:
    """Exhaust reflective control coordinates through a rational bound."""

    if height_bound < 0:
        raise ValueError("height bound must be nonnegative")
    levels: set[SilverE6ControlLevel] = set()
    for norm in UNIT_E6_TETRACODE_NORMS:
        caps = (norm * height_bound / EPSILON, norm / (-EPSILON))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0 or not (2 * EPSILON * control / norm).is_integral:
                continue
            target = norm + EPSILON * control * control
            if target.sign(1) < 0:
                continue
            height = EPSILON * control * control / norm
            if (Silver(height_bound) - height).sign(0) < 0:
                continue
            levels.add(SilverE6ControlLevel(norm, height, control, target))
    return tuple(sorted(levels, key=lambda level: (level.height, level.norm, level.control)))


@dataclass(frozen=True, slots=True)
class SilverE6Candidate:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver
    blocks: tuple[SilverE6Weight, ...]
    vector: Vector

    @property
    def glue_coupled(self) -> bool:
        return any(block.residue != (0, 0) for block in self.blocks)


def _candidate_is_primitive(candidate: SilverE6Candidate) -> bool:
    if candidate.norm != 2 or not (candidate.control / SQRT2).is_integral:
        return True
    divided_weights = tuple(
        tuple(value / SQRT2 for value in block.weights) for block in candidate.blocks
    )
    if not all(value.is_integral for block in divided_weights for value in block):
        return True
    return not tetracode_membership(divided_weights)


@lru_cache(maxsize=None)
def enumerate_control_level(level: SilverE6ControlLevel) -> tuple[SilverE6Candidate, ...]:
    """Enumerate the complete glue-aware dominant shell at one control."""

    blocks = enumerate_dominant_e6_weights(level.target, level.norm).vectors
    pairs: dict[Silver, list[tuple[SilverE6Weight, SilverE6Weight]]] = defaultdict(list)
    for left in blocks:
        for right in blocks:
            pairs[left.square + right.square].append((left, right))
    roots: list[SilverE6Candidate] = []
    for square, left_pairs in pairs.items():
        for left in left_pairs:
            for right in pairs.get(level.target - square, ()):
                four = (*left, *right)
                if not tetracode_membership(tuple(block.weights for block in four)):
                    continue
                vector = (
                    level.control,
                    *(entry for block in four for entry in block.coordinates),
                )
                candidate = SilverE6Candidate(
                    level.norm,
                    level.height,
                    level.control,
                    level.target,
                    four,
                    vector,
                )
                if _candidate_is_primitive(candidate):
                    roots.append(candidate)
    return tuple(sorted(roots, key=lambda root: root.vector))


def silver_e6_pair(left: Sequence[Silver], right: Sequence[Silver]) -> Silver:
    """The exact Lorentzian pairing in four E6 simple-root bases."""

    if len(left) != 25 or len(right) != 25:
        raise ValueError("an H24 E6-tetracode vector has rank 25")
    cartan = e6_cartan_gram()
    result = -EPSILON * left[0] * right[0]
    for block in range(4):
        offset = 1 + 6 * block
        result += sum(
            (
                left[offset + i] * Silver(cartan[i][j]) * right[offset + j]
                for i in range(6)
                for j in range(6)
                if cartan[i][j]
            ),
            start=SILVER_ZERO,
        )
    return result


@lru_cache(maxsize=1)
def height_zero_simple_roots() -> tuple[Vector, ...]:
    roots = []
    for index in range(24):
        root = [SILVER_ZERO] * 25
        root[index + 1] = -SILVER_ONE
        roots.append(tuple(root))
    return tuple(roots)


@lru_cache(maxsize=1)
def height_zero_dominant_roots() -> tuple[SilverE6Candidate, ...]:
    roots = []
    for norm in UNIT_E6_TETRACODE_NORMS:
        roots.extend(
            enumerate_control_level(
                SilverE6ControlLevel(norm, SILVER_ZERO, SILVER_ZERO, norm)
            )
        )
    return tuple(sorted(roots, key=lambda root: (root.norm, root.vector)))


@dataclass(frozen=True, slots=True)
class SilverE6Packet:
    height: Silver
    controls: tuple[SilverE6ControlLevel, ...]
    candidate_count: int
    rejected_count: int
    roots: tuple[SilverE6Candidate, ...]


@dataclass(frozen=True, slots=True)
class SilverE6Prefix:
    height_bound: int
    controls: tuple[SilverE6ControlLevel, ...]
    packets: tuple[SilverE6Packet, ...]
    accepted_roots: tuple[SilverE6Candidate, ...]
    walls: tuple[Vector, ...]


@lru_cache(maxsize=None)
def run_e6_tetracode_prefix(height_bound: int = 15) -> SilverE6Prefix:
    """Run the exact chronological Vinberg prefix through the second packet."""

    if len(height_zero_dominant_roots()) != 4 or any(
        root.norm != 2 for root in height_zero_dominant_roots()
    ):
        raise AssertionError("height zero is not exactly E6^4")
    controls = enumerate_control_levels(height_bound)
    walls = list(height_zero_simple_roots())
    accepted: list[SilverE6Candidate] = []
    packets: list[SilverE6Packet] = []
    for height, iterator in groupby(controls, key=lambda level: level.height):
        packet_controls = tuple(iterator)
        candidates = tuple(
            root
            for level in packet_controls
            for root in enumerate_control_level(level)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(silver_e6_pair(root.vector, wall).sign(0) <= 0 for wall in walls)
        )
        if any(
            silver_e6_pair(left.vector, right.vector).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height E6-tetracode roots are not mutually acute")
        accepted.extend(survivors)
        walls.extend(root.vector for root in survivors)
        packets.append(
            SilverE6Packet(
                height,
                packet_controls,
                len(candidates),
                len(candidates) - len(survivors),
                survivors,
            )
        )
    return SilverE6Prefix(
        height_bound,
        controls,
        tuple(packets),
        tuple(accepted),
        tuple(walls),
    )


@dataclass(frozen=True, slots=True)
class SilverE6FirstPacketGeometry:
    cap: Silver
    extreme_ray_count: int
    finite_vertex_count: int
    cusp_count: int
    hyperideal_ray_count: int
    compact: bool


@lru_cache(maxsize=1)
def first_packet_geometry() -> SilverE6FirstPacketGeometry:
    """Classify every ray of the product of four affine E6 alcoves."""

    inverse = e6_inverse_cartan()
    energies = [Fraction(0)]
    for index, mark in enumerate(E6_HIGHEST_ROOT):
        energies.append(Fraction(4) * inverse[index][index] / (mark * mark))
    finite = cusps = hyperideal = 0
    for choices in product(range(7), repeat=4):
        norm = Silver(sum((energies[index] for index in choices), start=Fraction(0))) - EPSILON
        sign = norm.sign(0)
        finite += sign < 0
        cusps += sign == 0
        hyperideal += sign > 0
    return SilverE6FirstPacketGeometry(
        Silver(2),
        7**4,
        finite,
        cusps,
        hyperideal,
        hyperideal == 0 and cusps == 0,
    )


@dataclass(frozen=True, slots=True)
class SilverE6WitnessAudit:
    dynkin_weights: tuple[Vector, ...]
    norm: Silver
    height_zero_nonpositive: bool
    positive_wall_pairings: tuple[Silver, ...]
    all_positive_walls_strict: bool
    hyperideal: bool


@lru_cache(maxsize=1)
def second_packet_hyperideal_witness() -> SilverE6WitnessAudit:
    """Audit the symmetric hyperideal ray through both nonempty packets."""

    block = (
        SILVER_ZERO,
        SILVER_ZERO,
        SILVER_ZERO,
        SILVER_ZERO,
        SILVER_ZERO,
        Silver(HYPERIDEAL_WITNESS_WEIGHT),
    )
    weights = (block, block, block, block)
    inverse = e6_inverse_cartan()
    spatial_norm = sum(
        (
            _quadratic(block_weights, inverse)
            for block_weights in weights
        ),
        start=SILVER_ZERO,
    )
    norm = spatial_norm - EPSILON
    prefix = run_e6_tetracode_prefix(15)
    pairings = tuple(
        -EPSILON * root.control
        + sum(
            (
                block_weights[index] * root.blocks[block_index].coordinates[index]
                for block_index, block_weights in enumerate(weights)
                for index in range(6)
            ),
            start=SILVER_ZERO,
        )
        for root in prefix.accepted_roots
    )
    return SilverE6WitnessAudit(
        weights,
        norm,
        all(value.sign(0) >= 0 for block_weights in weights for value in block_weights),
        pairings,
        all(pairing.sign(0) < 0 for pairing in pairings),
        norm.sign(0) > 0,
    )


__all__ = [
    "COMPLETE_UNIT_NORM_PROOF",
    "E6_DISCRIMINANT_FUNCTIONAL",
    "E6_HIGHEST_ROOT",
    "EPSILON",
    "SECOND_PACKET_HEIGHT",
    "SilverE6Candidate",
    "SilverE6ControlLevel",
    "SilverE6Enumeration",
    "SilverE6FirstPacketGeometry",
    "SilverE6Packet",
    "SilverE6Prefix",
    "SilverE6Weight",
    "SilverE6WitnessAudit",
    "e6_cartan_gram",
    "e6_discriminant_residue",
    "e6_inverse_cartan",
    "enumerate_control_level",
    "enumerate_control_levels",
    "enumerate_dominant_e6_weights",
    "first_packet_geometry",
    "height_zero_dominant_roots",
    "height_zero_simple_roots",
    "run_e6_tetracode_prefix",
    "second_packet_hyperideal_witness",
    "silver_e6_pair",
    "ternary_tetracode",
    "tetracode_membership",
]
