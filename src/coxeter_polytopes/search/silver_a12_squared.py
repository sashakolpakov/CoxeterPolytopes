"""Exact silver-unit Vinberg packets for the Niemeier lattice ``N(A12^2)``.

The discriminant group of each A12 factor is Z/13.  In the standard
numbering, the Niemeier glue is the graph ``{(k,5k):k in Z/13}``, since
``5*8/13 + 1*12/13 = 4``.  Membership is imposed separately on the two
coefficients of Z[sqrt(2)].
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import groupby
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.silver import SILVER_ONE, SILVER_ZERO, SQRT2, Silver
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.niemeier_host_triage import (
    affine_vertex_energies,
    cartan_gram,
)
from coxeter_polytopes.search.silver_a17_e7 import (
    A17Enumeration,
    DominantWeight,
    _determinant,
    _pair_add,
    _pair_conjugate,
    _pair_mul,
    _pair_scale,
    _pair_sign,
    _pair_sub,
)
from coxeter_polytopes.search.silver_e8_cubes import (
    _bounded_for_caps,
    _mat_vec,
    _quadratic,
)


Vector: TypeAlias = tuple[Silver, ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]
IntegerPair: TypeAlias = tuple[int, int]

EPSILON = SILVER_ONE + SQRT2
UNIT_NORM_CLASSES = (SILVER_ONE, Silver(2))
A12_COMPONENT = ("A", 12)
GLUE_MULTIPLIER = 8

COMPLETE_UNIT_NORM_PROOF = (
    "Put pi=sqrt(2) and epsilon=1+pi. The coefficient is a unit and "
    "N(A12^2) is even unimodular. For primitive r=(x,y), reflection "
    "integrality says q(r)|2I(r), I(r)=(epsilon*x,(y,N)). Away from "
    "(pi), local primitivity makes I(r) a unit. At (pi), the same even "
    "unimodular valuation argument leaves valuations zero or two. Every "
    "totally-positive unit is a square, so the complete canonical norms "
    "are 1 and 2."
)


@lru_cache(maxsize=1)
def a12_squared_niemeier_gram() -> RationalMatrix:
    """Gram of alpha_1..alpha_12, beta_1..beta_11, g=(omega_5,omega_1)."""

    a12 = cartan_gram(A12_COMPONENT)
    gram = [[Fraction(0) for _ in range(24)] for _ in range(24)]
    for i in range(12):
        for j in range(12):
            gram[i][j] = a12[i][j]
    for i in range(11):
        for j in range(11):
            gram[12 + i][12 + j] = a12[i][j]
    # The final basis vector is g=omega_5 + omega_1.  Its pairings with
    # simple roots are one at nodes 5 and 1, and g^2=40/13+12/13=4.
    gram[4][23] = gram[23][4] = Fraction(1)
    gram[12][23] = gram[23][12] = Fraction(1)
    gram[23][23] = Fraction(4)
    return tuple(tuple(row) for row in gram)


def a12_squared_basis_audit() -> bool:
    gram = a12_squared_niemeier_gram()
    return (
        _determinant(gram) == 1
        and all(value.denominator == 1 for row in gram for value in row)
        and all(gram[index][index].numerator % 2 == 0 for index in range(24))
    )


def a12_residue(weights: Sequence[Silver]) -> tuple[int, int]:
    if len(weights) != 12 or not all(value.is_integral for value in weights):
        raise ValueError("A12 Dynkin weights must be twelve silver integers")
    return tuple(
        sum(
            (index + 1) * int(getattr(value, coefficient))
            for index, value in enumerate(weights)
        )
        % 13
        for coefficient in ("a", "b")
    )


def a12_squared_glue_membership(
    left: Sequence[Silver], right: Sequence[Silver]
) -> bool:
    try:
        left_residue = a12_residue(left)
        right_residue = a12_residue(right)
    except ValueError:
        return False
    return all(
        right_residue[place] == GLUE_MULTIPLIER * left_residue[place] % 13
        for place in (0, 1)
    )


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
def enumerate_dominant_a12_weights(
    target_energies: tuple[Silver, ...],
    root_norm: Silver,
    control_cap: Silver | None,
) -> A17Enumeration:
    """Exhaust A12 weights at requested exact energies by chain DP."""

    targets = tuple(sorted(set(target_energies)))
    zero_counts = (0,) * 13
    if not targets or any(target.sign(place) < 0 for target in targets for place in (0, 1)):
        return A17Enumeration(targets, root_norm, control_cap, (), zero_counts, zero_counts, zero_counts, 0, ())
    largest = max(targets)
    alphabet = _dominant_alphabet(largest, root_norm)
    if any(not value.is_integral for value in alphabet):
        raise AssertionError("the A12 Dynkin alphabet left Z[sqrt(2)]")
    integer_alphabet = tuple((int(value.a), int(value.b)) for value in alphabet)
    target_pairs = []
    for target in targets:
        scaled = 13 * target
        if not scaled.is_integral:
            raise AssertionError("an A12 target does not lie in (1/13)O")
        target_pairs.append((int(scaled.a), int(scaled.b)))
    target_pair_set = set(target_pairs)
    cap_pair = None
    if control_cap is not None:
        cap = 2 * control_cap
        if not cap.is_integral:
            raise AssertionError("the A12 affine cap is not integral")
        cap_pair = int(cap.a), int(cap.b)
    visited = [0] * 13
    cap_pruned = [0] * 13
    energy_pruned = [0] * 13
    leaves = 0
    raw: list[tuple[IntegerPair, tuple[IntegerPair, ...], IntegerPair]] = []
    zero = (0, 0)

    def extend(
        depth: int,
        weight_sum: IntegerPair,
        weighted_sum: IntegerPair,
        suffix_square_sum: IntegerPair,
        prefix: tuple[IntegerPair, ...],
    ) -> None:
        nonlocal leaves
        visited[depth] += 1
        if depth:
            minimum_numerator = _pair_sub(
                _pair_scale(depth + 1, suffix_square_sum),
                _pair_mul(weighted_sum, weighted_sum),
            )
            possible = False
            for target_pair in target_pairs:
                difference = _pair_sub(
                    _pair_scale(depth + 1, target_pair),
                    _pair_scale(13, minimum_numerator),
                )
                if _pair_sign(difference) >= 0 and _pair_sign(_pair_conjugate(difference)) >= 0:
                    possible = True
                    break
            if not possible:
                energy_pruned[depth] += 1
                return
        if depth == 12:
            leaves += 1
            energy_times_thirteen = _pair_sub(
                _pair_scale(13, suffix_square_sum),
                _pair_mul(weighted_sum, weighted_sum),
            )
            if energy_times_thirteen in target_pair_set:
                raw.append((energy_times_thirteen, prefix, weighted_sum))
            return
        position = depth + 1
        for value in integer_alphabet:
            new_sum = _pair_add(weight_sum, value)
            if cap_pair is not None and _pair_sign(_pair_sub(cap_pair, new_sum)) < 0:
                cap_pruned[depth + 1] += 1
                continue
            new_weighted = _pair_add(weighted_sum, _pair_scale(position, value))
            new_suffix_squares = _pair_add(
                suffix_square_sum,
                _pair_add(
                    _pair_scale(2, _pair_mul(value, weighted_sum)),
                    _pair_scale(position, _pair_mul(value, value)),
                ),
            )
            extend(depth + 1, new_sum, new_weighted, new_suffix_squares, (*prefix, value))

    extend(0, zero, zero, zero, ())
    inverse = _fraction_inverse(cartan_gram(A12_COMPONENT))
    output = []
    for energy_pair, integer_weights, weighted_sum in raw:
        weights = tuple(Silver(*value) for value in integer_weights)
        square = Silver(*energy_pair) / 13
        coordinates = _mat_vec(inverse, weights)
        if _quadratic(weights, inverse) != square:
            raise AssertionError("the A12 recurrence disagrees with C^-1")
        residue = weighted_sum[0] % 13, weighted_sum[1] % 13
        output.append(DominantWeight(square, weights, coordinates, residue))
    return A17Enumeration(
        targets,
        root_norm,
        control_cap,
        alphabet,
        tuple(visited),
        tuple(cap_pruned),
        tuple(energy_pruned),
        leaves,
        tuple(sorted(output, key=lambda item: (item.square, item.weights))),
    )


@dataclass(frozen=True, slots=True)
class ControlLevel:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver


@lru_cache(maxsize=None)
def enumerate_control_levels(height_bound: int = 15) -> tuple[ControlLevel, ...]:
    levels = set()
    for norm in UNIT_NORM_CLASSES:
        caps = (norm * height_bound / EPSILON, norm / (-EPSILON))
        for control in _bounded_for_caps(caps):
            if control.sign(0) <= 0 or not (2 * EPSILON * control / norm).is_integral:
                continue
            target = norm + EPSILON * control * control
            height = EPSILON * control * control / norm
            if target.sign(1) >= 0 and (Silver(height_bound) - height).sign(0) >= 0:
                levels.add(ControlLevel(norm, height, control, target))
    return tuple(sorted(levels, key=lambda item: (item.height, item.norm, item.control)))


@dataclass(frozen=True, slots=True)
class Candidate:
    norm: Silver
    height: Silver
    control: Silver
    target: Silver
    left: DominantWeight
    right: DominantWeight

    @property
    def glue_nontrivial(self) -> bool:
        return self.left.residue != (0, 0) or self.right.residue != (0, 0)


def candidate_is_primitive(candidate: Candidate) -> bool:
    if not (candidate.control / SQRT2).is_integral:
        return True
    left = tuple(value / SQRT2 for value in candidate.left.weights)
    right = tuple(value / SQRT2 for value in candidate.right.weights)
    if not all(value.is_integral for value in (*left, *right)):
        return True
    return not a12_squared_glue_membership(left, right)


@lru_cache(maxsize=None)
def enumerate_control_level(level: ControlLevel, caps_active: bool) -> tuple[Candidate, ...]:
    cap = level.control if caps_active else None
    provisional = enumerate_dominant_a12_weights((level.target,), level.norm, cap).vectors
    targets = tuple(sorted({level.target - item.square for item in provisional if all((level.target-item.square).sign(p)>=0 for p in (0,1))}))
    complete = enumerate_dominant_a12_weights(
        tuple(sorted(set(targets) | {level.target})), level.norm, cap
    ).vectors
    by_square: dict[Silver, list[DominantWeight]] = defaultdict(list)
    for item in complete:
        by_square[item.square].append(item)
    output = []
    for left in complete:
        for right in by_square.get(level.target - left.square, ()):
            if not a12_squared_glue_membership(left.weights, right.weights):
                continue
            candidate = Candidate(level.norm, level.height, level.control, level.target, left, right)
            if candidate_is_primitive(candidate):
                output.append(candidate)
    return tuple(sorted(output, key=lambda item: (item.left.weights, item.right.weights)))


def pair(left: Candidate, right: Candidate) -> Silver:
    return (
        -EPSILON * left.control * right.control
        + sum((left.left.weights[i] * right.left.coordinates[i] for i in range(12)), start=SILVER_ZERO)
        + sum((left.right.weights[i] * right.right.coordinates[i] for i in range(12)), start=SILVER_ZERO)
    )


@dataclass(frozen=True, slots=True)
class Packet:
    height: Silver
    controls: tuple[ControlLevel, ...]
    candidate_count: int
    rejected_count: int
    roots: tuple[Candidate, ...]


@dataclass(frozen=True, slots=True)
class Prefix:
    height_bound: int
    controls: tuple[ControlLevel, ...]
    packets: tuple[Packet, ...]
    accepted_roots: tuple[Candidate, ...]


@lru_cache(maxsize=1)
def height_zero_dominant_roots() -> tuple[Candidate, ...]:
    roots = []
    for norm in UNIT_NORM_CLASSES:
        roots.extend(enumerate_control_level(ControlLevel(norm, SILVER_ZERO, SILVER_ZERO, norm), False))
    return tuple(roots)


@lru_cache(maxsize=None)
def run_prefix(height_bound: int = 15, nonempty_limit: int | None = 2) -> Prefix:
    if not a12_squared_basis_audit():
        raise AssertionError("A12^2 basis failed")
    h0 = height_zero_dominant_roots()
    if len(h0) != 2 or any(root.norm != 2 for root in h0):
        raise AssertionError("height zero is not A12^2")
    controls = enumerate_control_levels(height_bound)
    accepted = []
    packets = []
    nonempty = 0
    for height, iterator in groupby(controls, key=lambda item: item.height):
        packet_controls = tuple(iterator)
        candidates = tuple(root for level in packet_controls for root in enumerate_control_level(level, nonempty > 0))
        survivors = tuple(root for root in candidates if all(pair(root, wall).sign(0) <= 0 for wall in accepted))
        if any(pair(left, right).sign(0) > 0 for i,left in enumerate(survivors) for right in survivors[i+1:]):
            raise AssertionError("same-height A12^2 roots are not acute")
        accepted.extend(survivors)
        packets.append(Packet(height, packet_controls, len(candidates), len(candidates)-len(survivors), survivors))
        if survivors:
            nonempty += 1
        if nonempty_limit is not None and nonempty == nonempty_limit:
            break
    return Prefix(height_bound, controls, tuple(packets), tuple(accepted))


@dataclass(frozen=True, slots=True)
class Geometry:
    extreme_ray_count: int
    finite_count: int
    cusp_count: int
    hyperideal_count: int


def first_packet_geometry() -> Geometry:
    distribution = {Fraction(0): 1}
    for _ in range(2):
        updated = {}
        for energy, multiplicity in distribution.items():
            for value in affine_vertex_energies(A12_COMPONENT):
                updated[energy+value] = updated.get(energy+value,0)+multiplicity
        distribution=updated
    finite=sum(m for e,m in distribution.items() if (Silver(e)-EPSILON).sign(0)<0)
    cusp=sum(m for e,m in distribution.items() if (Silver(e)-EPSILON).sign(0)==0)
    total=sum(distribution.values())
    return Geometry(total,finite,cusp,total-finite-cusp)


@dataclass(frozen=True, slots=True)
class Witness:
    height_bound: int
    dynkin_index: int
    weight: Fraction
    norm: Silver
    pairings: tuple[Silver, ...]
    all_positive_walls_strict: bool
    hyperideal: bool


@lru_cache(maxsize=1)
def strict_hyperideal_witness() -> Witness:
    prefix = run_prefix(150, None)
    inverse = _fraction_inverse(cartan_gram(A12_COMPONENT))
    weight = Fraction(5, 3)
    norm = Silver(weight * weight * inverse[0][0]) - EPSILON
    pairings = tuple(
        -EPSILON * root.control + Silver(weight) * root.left.coordinates[0]
        for root in prefix.accepted_roots
    )
    return Witness(
        150,
        0,
        weight,
        norm,
        pairings,
        all(value.sign(0) < 0 for value in pairings),
        norm.sign(0) > 0,
    )


__all__ = [name for name in globals() if not name.startswith("_")]
