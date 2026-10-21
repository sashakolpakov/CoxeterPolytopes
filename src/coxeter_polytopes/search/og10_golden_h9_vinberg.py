"""Direct exact Vinberg search on the golden OG10 fixed lattice.

The fixed field lattice is

    M = <-2*phi> + E8 + <2>

over ``Z[phi]``.  Write a vector as ``r=(x,z,u)`` and put ``w=E8*z``.
For the control vector ``p=(1,0,...,0)``, the norm-two equation is

    2*w^T*E8^{-1}*w + 4*u^2 = 4 + 4*phi*x^2,

and dominance for the height-zero ``E8+A1`` cone says ``w_i,u >= 0``.
If ``R=(X,Z,U)`` is an earlier simple root, the remaining Vinberg wall
condition is

    2*w.Z + 4*u*U <= 4*phi*x*X.

Thus a height level is a nine-variable positive-definite recursion.  This is
strictly smaller than filtering the eleven-coordinate ambient OG10 search.
Every comparison is exact at both real embeddings.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Sequence

from coxeter_polytopes.exact.golden import Golden, PHI, ZERO

from .bugaenko_vinberg import (
    RationalMatrix,
    Vector,
    _fraction_inverse,
    _mat_vec,
    _pair,
    _prefix_schur_complements,
)
from .h11_vinberg import (
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


def admissible_h9_control_coordinates(identity_bound: int) -> tuple[Golden, ...]:
    """Return every possible positive control coordinate through ``bound``.

    At the conjugate embedding the norm equation gives

        2 + 2*phi'*x'^2 >= 0,

    or ``x'^2 <= phi``.  The ambient candidate generator uses the weaker
    bound ``x'^2 <= 4*(1+4*phi)/11`` and an independently proved exhaustive
    coefficient box, so exact filtering by the sharper bound remains
    exhaustive.
    """

    return tuple(
        x
        for x in admissible_control_coordinates(identity_bound)
        if (PHI - x.conjugate() * x.conjugate()).sign(0) >= 0
    )


@lru_cache(maxsize=1)
def h9_field_gram() -> tuple[tuple[Golden, ...], ...]:
    e8 = _standard_e8_gram()
    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][0] = -2 * PHI
    for i in range(8):
        for j in range(8):
            gram[i + 1][j + 1] = Golden(e8[i][j])
    gram[9][9] = Golden(2)
    return tuple(tuple(row) for row in gram)


@lru_cache(maxsize=1)
def h9_cone_roots() -> tuple[Vector, ...]:
    standard_e8_roots_from_simple_reflections()
    e8 = tuple(
        (
            ZERO,
            *tuple(Golden(entry) for entry in _unit_vector(i, sign=-1)),
            ZERO,
        )
        for i in range(8)
    )
    a1 = (ZERO, *(ZERO for _ in range(8)), Golden(-1))
    return (*e8, a1)


@lru_cache(maxsize=2)
def _positive_definite_enumeration_data(a1_norm: int = 2) -> tuple[
    RationalMatrix, tuple[RationalMatrix, ...]
]:
    e8 = _standard_e8_gram()
    inverse_e8 = _fraction_inverse(e8)
    if any(entry.denominator != 1 for row in inverse_e8 for entry in row):
        raise AssertionError("E8 inverse must be integral")
    # The doubled norm equation is intentional.  The successive conditional
    # coefficients of ``2*E8^-1`` are at least one, so the shared coordinate
    # alphabet ``entry^2 <= target`` is an exhaustive per-coordinate box.
    positive_form = (
        *tuple(
            (*(2 * entry for entry in inverse_e8[i]), Fraction(0))
            for i in range(8)
        ),
        (*(Fraction(0) for _ in range(8)), Fraction(2 * a1_norm)),
    )
    return inverse_e8, _prefix_schur_complements(positive_form)


def _reconstruct_root(
    x: Golden, prefix: Sequence[Golden], inverse_e8: RationalMatrix
) -> Vector:
    weights = prefix[:8]
    u = prefix[8]
    z = _mat_vec(inverse_e8, weights)
    if any(not entry.is_integral for entry in z) or not u.is_integral:
        raise AssertionError(
            "positive-definite coordinates did not reconstruct an integral root"
        )
    return (x, *z, u)


@dataclass(frozen=True, slots=True)
class OG10GoldenH9HeightLevel:
    control_coordinate: Golden
    height: Golden
    target: Golden
    alphabet: CoordinateAlphabet
    visited_by_depth: tuple[int, ...]
    schur_pruned_by_depth: tuple[int, ...]
    wall_pruned_by_depth: tuple[int, ...]
    leaf_norm_mismatches: int
    roots: tuple[Vector, ...]


def _enumerate_height_level(
    control_coordinate: Golden,
    previous_nonzero_roots: Sequence[Vector],
    inverse_e8: RationalMatrix,
    schur: Sequence[RationalMatrix],
    norm: int = 2,
    a1_norm: int = 2,
) -> OG10GoldenH9HeightLevel:
    if a1_norm not in (1, 2) or norm not in ((2, 4) if a1_norm == 2 else (1, 2)):
        raise ValueError("unsupported golden H9 host/norm combination")
    x = control_coordinate
    target = 2 * norm + 4 * PHI * x * x
    if target.sign(0) < 0 or target.sign(1) < 0:
        raise ValueError("the control coordinate violates conjugate positivity")
    alphabet = _coordinate_alphabet(target)
    visited = [0] * 10
    schur_pruned = [0] * 10
    wall_pruned = [0] * 10
    leaf_norm_mismatches = 0
    solutions: list[Vector] = []

    target_integer = _ring_integer(target)
    alphabet_integers = tuple(_ring_integer(entry) for entry in alphabet.values)
    wall_coefficients: list[tuple[RingInteger, ...]] = []
    wall_bounds: list[RingInteger] = []
    e8 = _standard_e8_gram()
    for wall in previous_nonzero_roots:
        wall_weights = _mat_vec(e8, wall[1:9])
        if (
            any(entry.sign(0) < 0 for entry in wall_weights)
            or any(entry.sign(0) < 0 for entry in wall[1:9])
            or wall[9].sign(0) < 0
        ):
            raise AssertionError(
                "dominant-wall pruning requires nonnegative coefficients"
            )
        # In w-coordinates, (r,wall)<=0 is
        #     2*w.wall_z + 4*u*wall_u <= 4*phi*x*wall_x.
        wall_coefficients.append(
            tuple(
                _ring_integer(entry)
                for entry in (
                    *(2 * entry for entry in wall[1:9]),
                    2 * a1_norm * wall[9],
                )
            )
        )
        wall_bounds.append(_ring_integer(4 * PHI * x * wall[0]))

    integral_schur = tuple(_integral_form(matrix) for matrix in schur)

    def extend(
        prefix: tuple[RingInteger, ...],
        minimum: tuple[RingInteger, int] | None,
        parent_partial_products: tuple[RingInteger, ...],
    ) -> None:
        nonlocal leaf_norm_mismatches
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

        if depth == 9:
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
            # A norm-four reflection is integral exactly when its E8
            # pairings are even; E8 is unimodular.  The two outer
            # coordinates already have pairings divisible by two.
            if norm == 4 and depth < 8 and any(value % 2 for value in entry):
                continue
            child_minimum = _ri_add(
                _ri_add(base, _ri_multiply(linear, entry)),
                _ri_scale(_ri_multiply(entry, entry), diagonal),
            )
            extend((*prefix, entry), (child_minimum, denominator), partial_products)

    extend((), None, tuple((0, 0) for _ in wall_coefficients))
    gram = h9_field_gram()
    if a1_norm != 2:
        gram = tuple(tuple(Golden(a1_norm) if i == j == 9 else gram[i][j] for j in range(10)) for i in range(10))
    if any(_pair(root, gram, root) != norm for root in solutions):
        raise AssertionError("height enumerator produced a vector of the wrong norm")
    if any(
        _pair(left, gram, right).sign(0) > 0
        for i, left in enumerate(solutions)
        for right in solutions[i + 1 :]
    ):
        raise AssertionError("same-height Vinberg roots do not meet acutely")
    return OG10GoldenH9HeightLevel(
        control_coordinate=x,
        height=2 * PHI * x * x / norm,
        target=target,
        alphabet=alphabet,
        visited_by_depth=tuple(visited),
        schur_pruned_by_depth=tuple(schur_pruned),
        wall_pruned_by_depth=tuple(wall_pruned),
        leaf_norm_mismatches=leaf_norm_mismatches,
        roots=tuple(solutions),
    )


@dataclass(frozen=True, slots=True)
class OG10GoldenH9VinbergExtension:
    start_after: Golden
    identity_bound: int
    previous_nonzero_root_count: int
    candidate_coordinates: tuple[Golden, ...]
    levels: tuple[OG10GoldenH9HeightLevel, ...]
    new_roots: tuple[Vector, ...]


def run_og10_golden_h9_vinberg_extension(
    previous_nonzero_roots: Sequence[Vector],
    *,
    start_after: Golden,
    identity_bound: int,
    max_new_levels: int | None = None,
) -> OG10GoldenH9VinbergExtension:
    """Resume a certified direct H9 prefix strictly after one level."""

    if start_after.sign(0) < 0:
        raise ValueError("start_after must be nonnegative")
    if identity_bound < 1:
        raise ValueError("identity_bound must be positive")
    if max_new_levels is not None and max_new_levels < 1:
        raise ValueError("max_new_levels must be positive when supplied")

    gram = h9_field_gram()
    e8 = _standard_e8_gram()
    previous = tuple(tuple(entry for entry in root) for root in previous_nonzero_roots)
    if any(len(root) != 10 for root in previous):
        raise ValueError("every previous H9 root must have ten coordinates")
    if len(previous) != len(set(previous)):
        raise ValueError("the previous H9 prefix contains duplicate roots")
    if any(_pair(root, gram, root) != 2 for root in previous):
        raise ValueError("the previous H9 prefix contains a non-root")
    if any(
        root[0].sign(0) <= 0 or (root[0] - start_after).sign(0) > 0
        for root in previous
    ):
        raise ValueError("previous roots do not lie in the declared height prefix")
    if any(
        any(entry.sign(0) < 0 for entry in _mat_vec(e8, root[1:9]))
        or any(entry.sign(0) < 0 for entry in root[1:9])
        or root[9].sign(0) < 0
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
        for value in admissible_h9_control_coordinates(identity_bound)
        if (value - start_after).sign(0) > 0
    )
    if max_new_levels is not None:
        candidates = candidates[:max_new_levels]

    inverse_e8, schur = _positive_definite_enumeration_data()
    roots = list(previous)
    levels: list[OG10GoldenH9HeightLevel] = []
    for value in candidates:
        level = _enumerate_height_level(value, roots, inverse_e8, schur)
        levels.append(level)
        roots.extend(level.roots)
    return OG10GoldenH9VinbergExtension(
        start_after=start_after,
        identity_bound=identity_bound,
        previous_nonzero_root_count=len(previous),
        candidate_coordinates=candidates,
        levels=tuple(levels),
        new_roots=tuple(root for level in levels for root in level.roots),
    )


__all__ = [
    "OG10GoldenH9HeightLevel",
    "OG10GoldenH9VinbergExtension",
    "admissible_h9_control_coordinates",
    "h9_cone_roots",
    "h9_field_gram",
    "run_og10_golden_h9_vinberg_extension",
]
