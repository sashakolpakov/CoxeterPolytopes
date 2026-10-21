"""Exact silver Vinberg packets for the Niemeier lattice ``N(D12^2)``.

The root lattice ``D12 + D12`` has discriminant group
``{0,v,s,c}^2``.  The Niemeier glue is the graph which exchanges the
vector and one spinor class,

``{(0,0), (v,s), (s,v), (c,c)}``.

This distinction is essential: diagonal glue would contain ``(v,v)`` of
square two and therefore would not have root system ``D12^2``.  Membership
in the full silver extension is imposed independently on the rational and
``sqrt(2)`` coefficients.

The primary Lorentzian seed is

``<-(1+sqrt(2))> + N(D12^2)``.

Since its spatial block is even unimodular and the coefficient is a unit,
the complete canonical primitive mirror norms are ``1`` and ``2``.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby
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
from coxeter_polytopes.search.niemeier_host_triage import (
    affine_vertex_energies,
    cartan_gram,
    highest_root_marks,
)
from coxeter_polytopes.search.silver_e8_cubes import (
    _bounded_for_caps,
    _mat_vec,
    _quadratic,
)


Vector: TypeAlias = tuple[Silver, ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]

EPSILON = SILVER_ONE + SQRT2
UNIT_NORM_CLASSES = (SILVER_ONE, Silver(2))
D12_COMPONENT = ("D", 12)
FIRST_PACKET_HEIGHT = Silver(2, 2)
SECOND_PACKET_HEIGHT = Silver(7, 5)

COMPLETE_UNIT_NORM_PROOF = (
    "Put pi=sqrt(2) and epsilon=1+pi. The coefficient epsilon is a unit "
    "and N(D12^2) is even unimodular. For primitive r=(x,y), reflection "
    "integrality says q(r) divides 2I(r), I(r)=(epsilon*x,(y,N)). Away "
    "from (pi), local unimodularity and primitivity make I(r) the unit "
    "ideal. At (pi), a unit control gives v_pi(q)=0; otherwise a primitive "
    "spatial coordinate makes I(r) a unit, while evenness and q|2I force "
    "v_pi(q)=2. Every totally-positive unit of Z[sqrt(2)] is a square, so "
    "the complete canonical positive norms are 1 and 2."
)

COMPLETE_SQRT2_NORM_PROOF = (
    "Put pi=sqrt(2), so (2)=(pi)^2. For primitive r=(x,y), reflection "
    "integrality is q(r)|2I(r), I(r)=(pi*x,(y,N)). The spatial lattice "
    "N(D12^2) is even unimodular. Away from (pi), local unimodularity and "
    "primitivity make I(r) the unit ideal. If x is a unit then v_pi(q)=1; "
    "otherwise a primitive spatial coordinate is a unit and evenness plus "
    "q|2I gives v_pi(q)=2. Totally-positive units are squares, so the "
    "complete canonical norms are 2 and 2+pi."
)

COMPLETE_TWO_SQRT2_NORM_PROOF = (
    "Put pi=sqrt(2), so the coefficient 2sqrt(2)=pi^3. For primitive "
    "r=(x,y), q(r)|2I(r), I(r)=(pi^3*x,(y,N)), and N(D12^2) is even "
    "unimodular. If a spatial coordinate is a unit, evenness and q|2I "
    "give v_pi(q)=2. Otherwise x is a unit and y is pi-divisible, giving "
    "v_pi(q)=3. Totally-positive units are squares, so the complete "
    "canonical norms are 2 and 4+2pi."
)


@dataclass(frozen=True, slots=True)
class SilverD12SquaredSeed:
    seed_id: str
    timelike_coefficient: Silver
    canonical_mirror_norms: tuple[Silver, ...]
    affine_weight_cap: Silver
    complete_norm_class_proof: str


UNIT_D12_SQUARED_SEED = SilverD12SquaredSeed(
    "silver-H24-D12^2-a-1+sqrt2",
    EPSILON,
    UNIT_NORM_CLASSES,
    Silver(2),
    COMPLETE_UNIT_NORM_PROOF,
)

SQRT2_D12_SQUARED_SEED = SilverD12SquaredSeed(
    "silver-H24-D12^2-a-sqrt2",
    SQRT2,
    (Silver(2), Silver(2, 1)),
    SQRT2,
    COMPLETE_SQRT2_NORM_PROOF,
)

TWO_SQRT2_D12_SQUARED_SEED = SilverD12SquaredSeed(
    "silver-H24-D12^2-a-2sqrt2",
    2 * SQRT2,
    (Silver(2), Silver(4, 2)),
    Silver(2),
    COMPLETE_TWO_SQRT2_NORM_PROOF,
)

SUPPORTED_D12_SQUARED_SEEDS = (
    UNIT_D12_SQUARED_SEED,
    SQRT2_D12_SQUARED_SEED,
    TWO_SQRT2_D12_SQUARED_SEED,
)


def _determinant(matrix: RationalMatrix) -> Fraction:
    work = [list(row) for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next(row for row in range(column, len(work)) if work[row][column])
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result *= -1
        value = work[column][column]
        result *= value
        for row in range(column + 1, len(work)):
            factor = work[row][column] / value
            for index in range(column + 1, len(work)):
                work[row][index] -= factor * work[column][index]
    return result


@lru_cache(maxsize=1)
def d12_squared_niemeier_basis() -> RationalMatrix:
    """Return a Euclidean basis for the vector--spinor glued lattice.

    The vectors are stored directly in the orthonormal coordinates of
    ``R^12 + R^12``.  Twice this basis is an integer HNF basis of the module
    generated by ``D12^2``, ``(v,s)``, and ``(s,v)``.
    """

    scaled: list[list[int]] = []
    row = [0] * 24
    row[0] = 4
    scaled.append(row)
    for index in range(1, 12):
        row = [0] * 24
        row[0] = row[index] = 2
        scaled.append(row)
    for index in range(12, 23):
        row = [1] * 12 + [0] * 12
        row[index] = 2
        scaled.append(row)
    scaled.append([2] + [0] * 11 + [1] * 12)
    return tuple(tuple(Fraction(value, 2) for value in row) for row in scaled)


@lru_cache(maxsize=1)
def d12_squared_niemeier_gram() -> RationalMatrix:
    basis = d12_squared_niemeier_basis()
    return tuple(
        tuple(sum(left * right for left, right in zip(row, column)) for column in basis)
        for row in basis
    )


def _euclidean_d12_class(coordinates: Sequence[Fraction]) -> tuple[int, int]:
    """Classify a vector of ``D12*`` as 0, v, s, or c in bit notation."""

    if len(coordinates) != 12:
        raise ValueError("a D12 Euclidean block has twelve coordinates")
    values = tuple(Fraction(value) for value in coordinates)
    if all(value.denominator == 1 for value in values):
        return (sum(value.numerator for value in values) % 2, 0)
    shifted = tuple(value - Fraction(1, 2) for value in values)
    if all(value.denominator == 1 for value in shifted):
        return (sum(value.numerator for value in shifted) % 2, 1)
    raise ValueError("the vector does not belong to D12 dual")


def d12_squared_basis_glue_audit() -> bool:
    """Tie the displayed HNF basis to the vector--spinor graph code."""

    allowed = {
        ((0, 0), (0, 0)),
        ((1, 0), (0, 1)),
        ((0, 1), (1, 0)),
        ((1, 1), (1, 1)),
    }
    try:
        classes = tuple(
            (
                _euclidean_d12_class(row[:12]),
                _euclidean_d12_class(row[12:]),
            )
            for row in d12_squared_niemeier_basis()
        )
    except ValueError:
        return False
    # The four displayed classes form an index-four subgroup of the
    # discriminant group.  Since det(D12^2)=16, its overlattice has
    # determinant 16/4^2=1; containment plus the basis Gram determinant one
    # therefore proves equality with that graph-glued overlattice.
    generated = {((0, 0), (0, 0))}
    for left, right in classes:
        generated |= {
            (
                (base_left[0] ^ left[0], base_left[1] ^ left[1]),
                (base_right[0] ^ right[0], base_right[1] ^ right[1]),
            )
            for base_left, base_right in tuple(generated)
        }
    return all(value in allowed for value in classes) and generated == allowed


def d12_squared_basis_audit() -> bool:
    gram = d12_squared_niemeier_gram()
    return (
        d12_squared_basis_glue_audit()
        and len(d12_squared_niemeier_basis()) == 24
        and _determinant(gram) == 1
        and all(value.denominator == 1 for row in gram for value in row)
        and all(gram[index][index].numerator % 2 == 0 for index in range(24))
    )


def d12_residue(weights: Sequence[Silver]) -> tuple[tuple[int, int], tuple[int, int]]:
    """Return the vector/spinor bits in both silver coefficients.

    We label ``v=(1,0)``, ``s=(0,1)``, and ``c=(1,1)``.  For the standard
    D12 numbering, odd fundamental weights along the length-ten chain have
    class ``v``; the last two weights have classes ``s`` and ``c``.
    """

    if len(weights) != 12 or not all(value.is_integral for value in weights):
        raise ValueError("D12 Dynkin weights must be twelve silver integers")
    output = []
    for coefficient in ("a", "b"):
        values = [int(getattr(value, coefficient)) for value in weights]
        vector_bit = (sum(values[index] for index in (0, 2, 4, 6, 8)) + values[11]) % 2
        spinor_bit = (values[10] + values[11]) % 2
        output.append((vector_bit, spinor_bit))
    return output[0], output[1]


def d12_squared_glue_membership(
    left: Sequence[Silver], right: Sequence[Silver]
) -> bool:
    """Test ``right = swap(v,s)(left)`` in both coefficient residues."""

    try:
        left_residue = d12_residue(left)
        right_residue = d12_residue(right)
    except ValueError:
        return False
    return all(
        right_residue[place] == tuple(reversed(left_residue[place]))
        for place in (0, 1)
    )


@dataclass(frozen=True, slots=True)
class DominantD12Weight:
    square: Silver
    weights: Vector
    coordinates: Vector
    residue: tuple[tuple[int, int], tuple[int, int]]


@dataclass(frozen=True, slots=True)
class D12Enumeration:
    target: Silver
    root_norm: Silver
    control_cap: Silver | None
    alphabet: tuple[Silver, ...]
    visited_by_depth: tuple[int, ...]
    cap_pruned_by_depth: tuple[int, ...]
    schur_pruned_by_depth: tuple[int, ...]
    vectors: tuple[DominantD12Weight, ...]


def _dominant_alphabet(target: Silver, root_norm: Silver) -> tuple[Silver, ...]:
    twice_target = 2 * target
    return tuple(
        value
        for value in _bounded_for_caps((twice_target, twice_target))
        if value.sign(0) >= 0
        and all((twice_target - value * value).sign(place) >= 0 for place in (0, 1))
        and (2 * value / root_norm).is_integral
    )


@lru_cache(maxsize=None)
def enumerate_dominant_d12_weights(
    target: Silver,
    root_norm: Silver,
    control_cap: Silver | None,
) -> D12Enumeration:
    """Exhaust dominant D12 weights below ``target`` at both embeddings."""

    gram = cartan_gram(D12_COMPONENT)
    inverse = _fraction_inverse(gram)
    schur = _prefix_schur_complements(inverse)
    marks = highest_root_marks(D12_COMPONENT)
    alphabet = _dominant_alphabet(target, root_norm)
    visited = [0] * 13
    cap_pruned = [0] * 13
    schur_pruned = [0] * 13
    output: list[DominantD12Weight] = []

    def extend(prefix: Vector) -> None:
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            if control_cap is not None and (
                2 * control_cap
                - sum(
                    (marks[index] * prefix[index] for index in range(depth)),
                    start=SILVER_ZERO,
                )
            ).sign(0) < 0:
                cap_pruned[depth] += 1
                return
            minimum = _quadratic(prefix, schur[depth - 1])
            if any((target - minimum).sign(place) < 0 for place in (0, 1)):
                schur_pruned[depth] += 1
                return
        if depth == 12:
            square = _quadratic(prefix, inverse)
            output.append(
                DominantD12Weight(
                    square,
                    prefix,
                    _mat_vec(inverse, prefix),
                    d12_residue(prefix),
                )
            )
            return
        for value in alphabet:
            extend((*prefix, value))

    extend(())
    return D12Enumeration(
        target,
        root_norm,
        control_cap,
        alphabet,
        tuple(visited),
        tuple(cap_pruned),
        tuple(schur_pruned),
        tuple(sorted(output, key=lambda item: (item.square, item.weights))),
    )


@dataclass(frozen=True, slots=True)
class SilverD12SquaredControlLevel:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_seed_control_levels(
    seed: SilverD12SquaredSeed,
    height_bound: int,
) -> tuple[SilverD12SquaredControlLevel, ...]:
    """Exhaust reflective controls for any supported D12-squared seed."""

    if seed not in SUPPORTED_D12_SQUARED_SEEDS:
        raise ValueError("unsupported D12^2 seed")
    if height_bound < 0:
        raise ValueError("height bound must be nonnegative")
    levels: set[SilverD12SquaredControlLevel] = set()
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
            levels.add(SilverD12SquaredControlLevel(norm, height, control, target))
    return tuple(sorted(levels, key=lambda level: (level.height, level.norm, level.control)))


@lru_cache(maxsize=None)
def enumerate_control_levels(height_bound: int = 15) -> tuple[SilverD12SquaredControlLevel, ...]:
    return enumerate_seed_control_levels(UNIT_D12_SQUARED_SEED, height_bound)


@dataclass(frozen=True, slots=True)
class SilverD12SquaredCandidate:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver
    left: DominantD12Weight
    right: DominantD12Weight

    @property
    def glue_nontrivial(self) -> bool:
        zero = ((0, 0), (0, 0))
        return self.left.residue != zero or self.right.residue != zero

    @property
    def cross_component_glue(self) -> bool:
        return self.glue_nontrivial


def _candidate_is_primitive(candidate: SilverD12SquaredCandidate) -> bool:
    if not (candidate.control / SQRT2).is_integral:
        return True
    divided_left = tuple(value / SQRT2 for value in candidate.left.weights)
    divided_right = tuple(value / SQRT2 for value in candidate.right.weights)
    if not all(value.is_integral for value in (*divided_left, *divided_right)):
        return True
    return not d12_squared_glue_membership(divided_left, divided_right)


@lru_cache(maxsize=None)
def enumerate_control_level(
    level: SilverD12SquaredControlLevel,
    affine_caps_active: bool,
) -> tuple[SilverD12SquaredCandidate, ...]:
    """Enumerate one complete dominant shell with exact D12^2 glue."""

    cap = level.control if affine_caps_active else None
    vectors = enumerate_dominant_d12_weights(level.target, level.norm, cap).vectors
    by_square: dict[Silver, list[DominantD12Weight]] = defaultdict(list)
    for vector in vectors:
        by_square[vector.square].append(vector)
    candidates = []
    for left in vectors:
        for right in by_square.get(level.target - left.square, ()):
            if not d12_squared_glue_membership(left.weights, right.weights):
                continue
            candidate = SilverD12SquaredCandidate(
                level.norm,
                level.height,
                level.control,
                level.target,
                left,
                right,
            )
            if _candidate_is_primitive(candidate):
                if d12_squared_pair(candidate, candidate) != level.norm:
                    raise AssertionError("D12^2 reconstruction has the wrong norm")
                candidates.append(candidate)
    return tuple(
        sorted(candidates, key=lambda root: (root.control, root.left.weights, root.right.weights))
    )


@lru_cache(maxsize=None)
def enumerate_seed_control_level(
    seed: SilverD12SquaredSeed,
    level: SilverD12SquaredControlLevel,
    affine_caps_active: bool,
) -> tuple[SilverD12SquaredCandidate, ...]:
    """Enumerate a complete shell for a supported coefficient.

    ``enumerate_dominant_d12_weights`` expresses its optional cap as
    ``2*control_cap``.  The actual affine inequality is
    ``sum(mark_i*z_i) <= seed.affine_weight_cap*level.control``.
    """

    if seed not in SUPPORTED_D12_SQUARED_SEEDS:
        raise ValueError("unsupported D12^2 seed")
    cap = (
        seed.affine_weight_cap * level.control / 2
        if affine_caps_active
        else None
    )
    vectors = enumerate_dominant_d12_weights(level.target, level.norm, cap).vectors
    by_square: dict[Silver, list[DominantD12Weight]] = defaultdict(list)
    for vector in vectors:
        by_square[vector.square].append(vector)
    candidates = []
    for left in vectors:
        for right in by_square.get(level.target - left.square, ()):
            if not d12_squared_glue_membership(left.weights, right.weights):
                continue
            candidate = SilverD12SquaredCandidate(
                level.norm,
                level.height,
                level.control,
                level.target,
                left,
                right,
            )
            if _candidate_is_primitive(candidate):
                if d12_squared_seed_pair(seed, candidate, candidate) != level.norm:
                    raise AssertionError("D12^2 reconstruction has the wrong norm")
                candidates.append(candidate)
    return tuple(
        sorted(candidates, key=lambda root: (root.control, root.left.weights, root.right.weights))
    )


def d12_squared_pair(
    left: SilverD12SquaredCandidate,
    right: SilverD12SquaredCandidate,
) -> Silver:
    return (
        -EPSILON * left.control * right.control
        + sum(
            (
                left.left.weights[index] * right.left.coordinates[index]
                for index in range(12)
            ),
            start=SILVER_ZERO,
        )
        + sum(
            (
                left.right.weights[index] * right.right.coordinates[index]
                for index in range(12)
            ),
            start=SILVER_ZERO,
        )
    )


def d12_squared_seed_pair(
    seed: SilverD12SquaredSeed,
    left: SilverD12SquaredCandidate,
    right: SilverD12SquaredCandidate,
) -> Silver:
    return (
        -seed.timelike_coefficient * left.control * right.control
        + sum(
            (
                left.left.weights[index] * right.left.coordinates[index]
                for index in range(12)
            ),
            start=SILVER_ZERO,
        )
        + sum(
            (
                left.right.weights[index] * right.right.coordinates[index]
                for index in range(12)
            ),
            start=SILVER_ZERO,
        )
    )


@lru_cache(maxsize=None)
def height_zero_seed_dominant_roots(
    seed: SilverD12SquaredSeed,
) -> tuple[SilverD12SquaredCandidate, ...]:
    roots = []
    for norm in seed.canonical_mirror_norms:
        roots.extend(
            enumerate_seed_control_level(
                seed,
                SilverD12SquaredControlLevel(norm, SILVER_ZERO, SILVER_ZERO, norm),
                False,
            )
        )
    return tuple(sorted(roots, key=lambda root: (root.norm, root.left.weights, root.right.weights)))


@lru_cache(maxsize=1)
def height_zero_dominant_roots() -> tuple[SilverD12SquaredCandidate, ...]:
    roots = []
    for norm in UNIT_NORM_CLASSES:
        roots.extend(
            enumerate_control_level(
                SilverD12SquaredControlLevel(norm, SILVER_ZERO, SILVER_ZERO, norm),
                False,
            )
        )
    return tuple(sorted(roots, key=lambda root: (root.norm, root.left.weights, root.right.weights)))


@dataclass(frozen=True, slots=True)
class SilverD12SquaredPacket:
    height: Silver
    controls: tuple[SilverD12SquaredControlLevel, ...]
    candidate_count: int
    rejected_count: int
    roots: tuple[SilverD12SquaredCandidate, ...]


@dataclass(frozen=True, slots=True)
class SilverD12SquaredPrefix:
    height_bound: int
    controls: tuple[SilverD12SquaredControlLevel, ...]
    packets: tuple[SilverD12SquaredPacket, ...]
    accepted_roots: tuple[SilverD12SquaredCandidate, ...]


def _affine_caps_present(
    seed: SilverD12SquaredSeed,
    roots: Sequence[SilverD12SquaredCandidate],
) -> bool:
    """Check the two factorwise caps which justify later weight pruning."""

    marks = highest_root_marks(D12_COMPONENT)
    zero_weights = (SILVER_ZERO,) * 12
    for selected_side in ("left", "right"):
        found = False
        for root in roots:
            selected = getattr(root, selected_side)
            other = root.right if selected_side == "left" else root.left
            if other.weights != zero_weights:
                continue
            scale = selected.coordinates[0] / marks[0]
            if scale.sign(0) <= 0:
                continue
            if selected.coordinates != tuple(scale * mark for mark in marks):
                continue
            if seed.timelike_coefficient * root.control != seed.affine_weight_cap * scale:
                continue
            found = True
            break
        if not found:
            return False
    return True


@lru_cache(maxsize=None)
def run_d12_squared_seed_prefix(
    seed: SilverD12SquaredSeed,
    height_bound: int,
    nonempty_packet_limit: int | None = None,
) -> SilverD12SquaredPrefix:
    """Run a complete chronological glue-aware prefix for ``seed``."""

    if seed not in SUPPORTED_D12_SQUARED_SEEDS:
        raise ValueError("unsupported D12^2 seed")
    if not d12_squared_basis_audit():
        raise AssertionError("the displayed D12^2 glue basis is not even unimodular")
    height_zero = height_zero_seed_dominant_roots(seed)
    if len(height_zero) != 2 or any(root.norm != 2 for root in height_zero):
        raise AssertionError("height zero is not exactly D12^2")
    controls = enumerate_seed_control_levels(seed, height_bound)
    accepted: list[SilverD12SquaredCandidate] = []
    packets: list[SilverD12SquaredPacket] = []
    nonempty_packets = 0
    caps_certified = False
    for height, iterator in groupby(controls, key=lambda level: level.height):
        packet_controls = tuple(iterator)
        candidates = tuple(
            root
            for level in packet_controls
            for root in enumerate_seed_control_level(seed, level, caps_certified)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(d12_squared_seed_pair(seed, root, wall).sign(0) <= 0 for wall in accepted)
        )
        if any(
            d12_squared_seed_pair(seed, left, right).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height D12^2 survivors are not mutually acute")
        accepted.extend(survivors)
        packets.append(
            SilverD12SquaredPacket(
                height,
                packet_controls,
                len(candidates),
                len(candidates) - len(survivors),
                survivors,
            )
        )
        if survivors:
            nonempty_packets += 1
            if not caps_certified:
                if not _affine_caps_present(seed, survivors):
                    raise AssertionError("the first nonempty packet lacks its declared affine caps")
                caps_certified = True
        if nonempty_packet_limit is not None and nonempty_packets == nonempty_packet_limit:
            break
    return SilverD12SquaredPrefix(height_bound, controls, tuple(packets), tuple(accepted))


@lru_cache(maxsize=None)
def run_d12_squared_prefix(
    height_bound: int = 15,
    nonempty_packet_limit: int | None = 2,
) -> SilverD12SquaredPrefix:
    """Run the exact chronological unit-coefficient Vinberg prefix."""

    if not d12_squared_basis_audit():
        raise AssertionError("the displayed D12^2 glue basis is not even unimodular")
    height_zero = height_zero_dominant_roots()
    if len(height_zero) != 2 or any(root.norm != 2 for root in height_zero):
        raise AssertionError("height zero is not exactly D12^2")
    controls = enumerate_control_levels(height_bound)
    accepted: list[SilverD12SquaredCandidate] = []
    packets: list[SilverD12SquaredPacket] = []
    nonempty_packets = 0
    for height, iterator in groupby(controls, key=lambda level: level.height):
        packet_controls = tuple(iterator)
        candidates = tuple(
            root
            for level in packet_controls
            for root in enumerate_control_level(level, nonempty_packets > 0)
        )
        survivors = tuple(
            root
            for root in candidates
            if all(d12_squared_pair(root, wall).sign(0) <= 0 for wall in accepted)
        )
        if any(
            d12_squared_pair(left, right).sign(0) > 0
            for index, left in enumerate(survivors)
            for right in survivors[index + 1 :]
        ):
            raise AssertionError("same-height D12^2 survivors are not mutually acute")
        accepted.extend(survivors)
        packets.append(
            SilverD12SquaredPacket(
                height,
                packet_controls,
                len(candidates),
                len(candidates) - len(survivors),
                survivors,
            )
        )
        if survivors:
            nonempty_packets += 1
        if nonempty_packet_limit is not None and nonempty_packets == nonempty_packet_limit:
            break
    return SilverD12SquaredPrefix(height_bound, controls, tuple(packets), tuple(accepted))


@dataclass(frozen=True, slots=True)
class SilverD12SquaredWitness:
    seed_id: str
    height_bound: int
    block: int
    dynkin_index: int
    weight: Fraction
    norm: Silver
    positive_wall_pairings: tuple[Silver, ...]
    height_zero_nonpositive: bool
    all_positive_walls_strict: bool
    hyperideal: bool


_WITNESS_DATA = {
    UNIT_D12_SQUARED_SEED: (30, 0, 0, Fraction(8, 5)),
    SQRT2_D12_SQUARED_SEED: (40, 0, 2, Fraction(7, 10)),
    TWO_SQRT2_D12_SQUARED_SEED: (27, 0, 4, Fraction(7, 9)),
}


@lru_cache(maxsize=None)
def strict_hyperideal_witness(
    seed: SilverD12SquaredSeed,
) -> SilverD12SquaredWitness:
    """Return and audit a sparse strict witness for the frozen prefix."""

    if seed not in _WITNESS_DATA:
        raise ValueError("unsupported D12^2 seed")
    height_bound, block, dynkin_index, rational_weight = _WITNESS_DATA[seed]
    prefix = run_d12_squared_seed_prefix(seed, height_bound)
    inverse = _fraction_inverse(cartan_gram(D12_COMPONENT))
    weight = Silver(rational_weight)
    norm = weight * weight * Silver(inverse[dynkin_index][dynkin_index]) - seed.timelike_coefficient
    pairings = tuple(
        -seed.timelike_coefficient * root.control
        + weight
        * (root.left if block == 0 else root.right).coordinates[dynkin_index]
        for root in prefix.accepted_roots
    )
    return SilverD12SquaredWitness(
        seed.seed_id,
        height_bound,
        block,
        dynkin_index,
        rational_weight,
        norm,
        pairings,
        rational_weight >= 0,
        all(pairing.sign(0) < 0 for pairing in pairings),
        norm.sign(0) > 0,
    )


@dataclass(frozen=True, slots=True)
class SilverD12SquaredFirstGeometry:
    extreme_ray_count: int
    finite_count: int
    cusp_count: int
    hyperideal_count: int


@lru_cache(maxsize=1)
def first_packet_geometry() -> SilverD12SquaredFirstGeometry:
    distribution = {Fraction(0): 1}
    for _ in range(2):
        updated: dict[Fraction, int] = {}
        for energy, multiplicity in distribution.items():
            for vertex_energy in affine_vertex_energies(D12_COMPONENT):
                total = energy + vertex_energy
                updated[total] = updated.get(total, 0) + multiplicity
        distribution = updated
    finite = sum(
        multiplicity
        for energy, multiplicity in distribution.items()
        if (Silver(energy) - EPSILON).sign(0) < 0
    )
    cusps = sum(
        multiplicity
        for energy, multiplicity in distribution.items()
        if (Silver(energy) - EPSILON).sign(0) == 0
    )
    total = sum(distribution.values())
    return SilverD12SquaredFirstGeometry(total, finite, cusps, total - finite - cusps)


__all__ = [
    "COMPLETE_SQRT2_NORM_PROOF",
    "COMPLETE_TWO_SQRT2_NORM_PROOF",
    "COMPLETE_UNIT_NORM_PROOF",
    "D12_COMPONENT",
    "EPSILON",
    "FIRST_PACKET_HEIGHT",
    "SECOND_PACKET_HEIGHT",
    "SQRT2_D12_SQUARED_SEED",
    "SUPPORTED_D12_SQUARED_SEEDS",
    "TWO_SQRT2_D12_SQUARED_SEED",
    "UNIT_D12_SQUARED_SEED",
    "D12Enumeration",
    "DominantD12Weight",
    "SilverD12SquaredCandidate",
    "SilverD12SquaredControlLevel",
    "SilverD12SquaredFirstGeometry",
    "SilverD12SquaredPacket",
    "SilverD12SquaredPrefix",
    "SilverD12SquaredSeed",
    "SilverD12SquaredWitness",
    "d12_residue",
    "d12_squared_basis_audit",
    "d12_squared_basis_glue_audit",
    "d12_squared_glue_membership",
    "d12_squared_niemeier_basis",
    "d12_squared_niemeier_gram",
    "d12_squared_pair",
    "d12_squared_seed_pair",
    "enumerate_control_level",
    "enumerate_control_levels",
    "enumerate_seed_control_level",
    "enumerate_seed_control_levels",
    "enumerate_dominant_d12_weights",
    "first_packet_geometry",
    "height_zero_dominant_roots",
    "height_zero_seed_dominant_roots",
    "run_d12_squared_prefix",
    "run_d12_squared_seed_prefix",
    "strict_hyperideal_witness",
]
