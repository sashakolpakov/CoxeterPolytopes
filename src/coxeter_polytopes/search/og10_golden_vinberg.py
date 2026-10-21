"""Exact norm-two Vinberg search for the golden OG10 carrier.

The field lattice is

    H_OG = [[-2*phi,1],[1,2]] + E8 + <2>

over ``Z[phi]``.  It is the discriminant-11 H11 form with one positive
``<2>`` line.  For the same control vector ``p=(2,-1,0,...,0)`` put

    D = 1 + 4*phi,  t = x + 2*y,  w = E8*z.

For ``r=(x,y,z,u)`` the norm-two equation becomes

    t^2 + 2*w^T*E8^{-1}*w + 4*u^2 = 4 + D*x^2,
    t == x (mod 2),  t,w_i,u >= 0.

The height-zero root system is ``A1 + E8 + A1``.  This module certifies a
bounded prefix only; termination and finite covolume require separate tests.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO

from .bugaenko_vinberg import (
    RationalMatrix,
    Vector,
    _fraction_inverse,
    _mat_vec,
    _pair,
    _prefix_schur_complements,
)
from .h11_vinberg import (
    CONJUGATE_X_SQUARE_BOUND,
    D,
    CoordinateAlphabet,
    RingInteger,
    _coordinate_alphabet,
    _from_ring_integer,
    _integral_form,
    _ri_add,
    _ri_multiply,
    _ri_scale,
    _ri_sign,
    _ri_subtract,
    _ring_integer,
    _standard_e8_gram,
    _unit_vector,
    admissible_control_coordinates,
    standard_e8_roots_from_simple_reflections,
)


@lru_cache(maxsize=1)
def og10_field_gram() -> tuple[tuple[Golden, ...], ...]:
    e8 = _standard_e8_gram()
    gram = [[ZERO for _ in range(11)] for _ in range(11)]
    gram[0][0] = -2 * PHI
    gram[0][1] = gram[1][0] = ONE
    gram[1][1] = Golden(2)
    for i in range(8):
        for j in range(8):
            gram[i + 2][j + 2] = Golden(e8[i][j])
    gram[10][10] = Golden(2)
    return tuple(tuple(row) for row in gram)


def _cone_roots() -> tuple[Vector, ...]:
    standard_e8_roots_from_simple_reflections()
    first_a1 = (ZERO, Golden(-1), *(ZERO for _ in range(9)))
    e8 = tuple(
        (
            ZERO,
            ZERO,
            *tuple(Golden(entry) for entry in _unit_vector(i, sign=-1)),
            ZERO,
        )
        for i in range(8)
    )
    second_a1 = (ZERO, ZERO, *(ZERO for _ in range(8)), Golden(-1))
    return (first_a1, *e8, second_a1)


@lru_cache(maxsize=1)
def _positive_definite_enumeration_data() -> tuple[
    RationalMatrix, tuple[RationalMatrix, ...]
]:
    e8 = _standard_e8_gram()
    inverse_e8 = _fraction_inverse(e8)
    if any(entry.denominator != 1 for row in inverse_e8 for entry in row):
        raise AssertionError("E8 inverse must be integral")
    positive_form = (
        (Fraction(1), *(Fraction(0) for _ in range(9))),
        *tuple(
            (
                Fraction(0),
                *(2 * inverse_e8[i][j] for j in range(8)),
                Fraction(0),
            )
            for i in range(8)
        ),
        (*(Fraction(0) for _ in range(9)), Fraction(4)),
    )
    return inverse_e8, _prefix_schur_complements(positive_form)


def _reconstruct_root(
    x: Golden, prefix: Sequence[Golden], inverse_e8: RationalMatrix
) -> Vector:
    t = prefix[0]
    weights = prefix[1:9]
    u = prefix[9]
    y = (t - x) / 2
    z = _mat_vec(inverse_e8, weights)
    if not y.is_integral or any(not entry.is_integral for entry in z) or not u.is_integral:
        raise AssertionError("positive-definite coordinates did not reconstruct an integral root")
    return (x, y, *z, u)


@dataclass(frozen=True, slots=True)
class OG10GoldenHeightLevel:
    control_coordinate: Golden
    height: Golden
    target: Golden
    alphabet: CoordinateAlphabet
    visited_by_depth: tuple[int, ...]
    parity_pruned: int
    schur_pruned_by_depth: tuple[int, ...]
    wall_pruned_by_depth: tuple[int, ...]
    leaf_norm_mismatches: int
    roots: tuple[Vector, ...]


def _enumerate_height_level(
    control_coordinate: Golden,
    previous_nonzero_roots: Sequence[Vector],
    inverse_e8: RationalMatrix,
    schur: Sequence[RationalMatrix],
) -> OG10GoldenHeightLevel:
    x = control_coordinate
    target = 4 + D * x * x
    alphabet = _coordinate_alphabet(target)
    visited = [0] * 11
    schur_pruned = [0] * 11
    wall_pruned = [0] * 11
    parity_pruned = 0
    leaf_norm_mismatches = 0
    solutions: list[Vector] = []

    x_integer = _ring_integer(x)
    target_integer = _ring_integer(target)
    alphabet_integers = tuple(_ring_integer(entry) for entry in alphabet.values)
    wall_coefficients: list[tuple[RingInteger, ...]] = []
    wall_bounds: list[RingInteger] = []
    for wall in previous_nonzero_roots:
        wall_t = wall[0] + 2 * wall[1]
        if (
            wall_t.sign(0) < 0
            or any(entry.sign(0) < 0 for entry in wall[2:10])
            or wall[10].sign(0) < 0
        ):
            raise AssertionError("dominant-wall pruning requires nonnegative coefficients")
        wall_coefficients.append(
            tuple(
                _ring_integer(entry)
                for entry in (
                    wall_t,
                    *(2 * entry for entry in wall[2:10]),
                    4 * wall[10],
                )
            )
        )
        wall_bounds.append(_ring_integer(D * x * wall[0]))

    integral_schur = tuple(_integral_form(matrix) for matrix in schur)

    def extend(
        prefix: tuple[RingInteger, ...],
        minimum: tuple[RingInteger, int] | None,
        parent_partial_products: tuple[RingInteger, ...],
    ) -> None:
        nonlocal parity_pruned, leaf_norm_mismatches
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            assert minimum is not None
            minimum_numerator, minimum_denominator = minimum
            remainder = _ri_subtract(
                _ri_scale(target_integer, minimum_denominator), minimum_numerator
            )
            if _ri_sign(remainder, 0) < 0 or _ri_sign(remainder, 1) < 0:
                schur_pruned[depth] += 1
                return
            entry = prefix[-1]
            partial_products = tuple(
                _ri_add(
                    parent_partial_products[i],
                    _ri_multiply(entry, wall_coefficients[i][depth - 1]),
                )
                for i in range(len(wall_coefficients))
            )
            for partial, bound in zip(partial_products, wall_bounds, strict=True):
                if _ri_sign(_ri_subtract(bound, partial), 0) < 0:
                    wall_pruned[depth] += 1
                    return
        else:
            partial_products = parent_partial_products

        if depth == 10:
            assert minimum is not None
            minimum_numerator, minimum_denominator = minimum
            if minimum_numerator != _ri_scale(target_integer, minimum_denominator):
                leaf_norm_mismatches += 1
                return
            golden_prefix = tuple(_from_ring_integer(entry) for entry in prefix)
            solutions.append(_reconstruct_root(x, golden_prefix, inverse_e8))
            return

        denominator, next_form = integral_schur[depth]
        base: RingInteger = (0, 0)
        for i in range(depth):
            for j in range(depth):
                base = _ri_add(
                    base,
                    _ri_scale(
                        _ri_multiply(prefix[i], prefix[j]), next_form[i][j]
                    ),
                )
        linear: RingInteger = (0, 0)
        for i in range(depth):
            linear = _ri_add(
                linear,
                _ri_scale(prefix[i], 2 * next_form[i][depth]),
            )
        diagonal = next_form[depth][depth]
        for entry in alphabet_integers:
            if depth == 0 and (
                (entry[0] - x_integer[0]) % 2
                or (entry[1] - x_integer[1]) % 2
            ):
                parity_pruned += 1
                continue
            child_minimum = _ri_add(
                _ri_add(base, _ri_multiply(linear, entry)),
                _ri_scale(_ri_multiply(entry, entry), diagonal),
            )
            extend((*prefix, entry), (child_minimum, denominator), partial_products)

    extend((), None, tuple((0, 0) for _ in wall_coefficients))
    gram = og10_field_gram()
    if any(_pair(root, gram, root) != 2 for root in solutions):
        raise AssertionError("height enumerator produced a vector of the wrong norm")
    if any(
        _pair(left, gram, right).sign(0) > 0
        for i, left in enumerate(solutions)
        for right in solutions[i + 1 :]
    ):
        raise AssertionError("same-height Vinberg roots do not meet acutely")
    return OG10GoldenHeightLevel(
        control_coordinate=x,
        height=D * x * x / 4,
        target=target,
        alphabet=alphabet,
        visited_by_depth=tuple(visited),
        parity_pruned=parity_pruned,
        schur_pruned_by_depth=tuple(schur_pruned),
        wall_pruned_by_depth=tuple(wall_pruned),
        leaf_norm_mismatches=leaf_norm_mismatches,
        roots=tuple(solutions),
    )


@dataclass(frozen=True, slots=True)
class OG10GoldenVinbergSearch:
    control_vector: Vector
    identity_bound: int
    cone_roots: tuple[Vector, ...]
    candidate_coordinates: tuple[Golden, ...]
    levels: tuple[OG10GoldenHeightLevel, ...]
    roots: tuple[Vector, ...]


@dataclass(frozen=True, slots=True)
class OG10GoldenVinbergExtension:
    start_after: Golden
    identity_bound: int
    previous_nonzero_root_count: int
    candidate_coordinates: tuple[Golden, ...]
    levels: tuple[OG10GoldenHeightLevel, ...]
    new_roots: tuple[Vector, ...]


def run_og10_golden_vinberg_extension(
    previous_nonzero_roots: Sequence[Vector],
    *,
    start_after: Golden,
    identity_bound: int,
    max_new_levels: int | None = None,
) -> OG10GoldenVinbergExtension:
    """Resume a certified positive-height prefix strictly after one level."""

    if start_after.sign(0) < 0:
        raise ValueError("start_after must be nonnegative")
    if identity_bound < 1:
        raise ValueError("identity_bound must be positive")
    if max_new_levels is not None and max_new_levels < 1:
        raise ValueError("max_new_levels must be positive when supplied")

    gram = og10_field_gram()
    previous = tuple(tuple(entry for entry in root) for root in previous_nonzero_roots)
    if any(len(root) != 11 for root in previous):
        raise ValueError("every previous OG10 root must have eleven coordinates")
    if len(previous) != len(set(previous)):
        raise ValueError("the previous OG10 prefix contains duplicate roots")
    if any(_pair(root, gram, root) != 2 for root in previous):
        raise ValueError("the previous OG10 prefix contains a non-root")
    if any(
        root[0].sign(0) <= 0 or (root[0] - start_after).sign(0) > 0
        for root in previous
    ):
        raise ValueError("previous roots do not lie in the declared height prefix")
    if any(
        (root[0] + 2 * root[1]).sign(0) < 0
        or any(entry.sign(0) < 0 for entry in root[2:10])
        or root[10].sign(0) < 0
        for root in previous
    ):
        raise ValueError("previous roots are not height-zero dominant")
    if any(
        _pair(left, gram, right).sign(0) > 0
        for index, left in enumerate(previous)
        for right in previous[index + 1 :]
    ):
        raise ValueError("previous roots are not pairwise acute")

    candidates = tuple(
        value
        for value in admissible_control_coordinates(identity_bound)
        if (value - start_after).sign(0) > 0
    )
    if max_new_levels is not None:
        candidates = candidates[:max_new_levels]

    inverse_e8, schur = _positive_definite_enumeration_data()
    roots = list(previous)
    levels: list[OG10GoldenHeightLevel] = []
    for value in candidates:
        level = _enumerate_height_level(value, roots, inverse_e8, schur)
        levels.append(level)
        roots.extend(level.roots)
    return OG10GoldenVinbergExtension(
        start_after=start_after,
        identity_bound=identity_bound,
        previous_nonzero_root_count=len(previous),
        candidate_coordinates=candidates,
        levels=tuple(levels),
        new_roots=tuple(root for level in levels for root in level.roots),
    )


@lru_cache(maxsize=None)
def run_og10_golden_vinberg_search(identity_bound: int = 1) -> OG10GoldenVinbergSearch:
    """Exhaust every norm-two Vinberg height with ``0 < x <= bound``."""

    gram = og10_field_gram()
    control = (Golden(2), Golden(-1), *(ZERO for _ in range(9)))
    if _pair(control, gram, control) != -2 * D:
        raise AssertionError("the control vector has the wrong norm")
    control_covector = tuple(
        sum((gram[i][j] * control[j] for j in range(11)), start=ZERO)
        for i in range(11)
    )
    if control_covector != (-D, *(ZERO for _ in range(10))):
        raise AssertionError("the control vector does not diagonalize the binary block")

    inverse_e8, schur = _positive_definite_enumeration_data()
    cone = _cone_roots()
    candidates = admissible_control_coordinates(identity_bound)
    levels: list[OG10GoldenHeightLevel] = []
    nonzero_roots: list[Vector] = []
    for x in candidates:
        level = _enumerate_height_level(x, nonzero_roots, inverse_e8, schur)
        levels.append(level)
        nonzero_roots.extend(level.roots)
    roots = (*cone, *nonzero_roots)
    if any(
        _pair(left, gram, right).sign(0) > 0
        for i, left in enumerate(roots)
        for right in roots[i + 1 :]
    ):
        raise AssertionError("enumerated simple roots are not pairwise acute")
    return OG10GoldenVinbergSearch(
        control_vector=control,
        identity_bound=identity_bound,
        cone_roots=cone,
        candidate_coordinates=candidates,
        levels=tuple(levels),
        roots=tuple(roots),
    )


__all__ = [
    "CONJUGATE_X_SQUARE_BOUND",
    "OG10GoldenHeightLevel",
    "OG10GoldenVinbergSearch",
    "og10_field_gram",
    "run_og10_golden_vinberg_extension",
    "run_og10_golden_vinberg_search",
]
