"""Exact norm-two Vinberg enumeration for the discriminant-11 field lattice.

The lattice is

    H11 = [[-2*phi, 1], [1, 2]] + E8

over ``Z[phi]``.  The apparently coupled binary summand becomes especially
simple for the control vector ``p=(2,-1,0,...,0)``.  Put

    D = 1 + 4*phi,       t = x + 2*y,       w = E8*z

for a vector ``r=(x,y,z)``.  Then

    <p,p> = -2*D,        <p,r> = -D*x,
    2*<r,r> = t^2 - D*x^2 + 2*w^T*E8^{-1}*w.

Consequently norm-two roots at fixed positive ``x`` are exactly the integral
solutions of

    t^2 + 2*w^T*E8^{-1}*w = 4 + D*x^2,
    t == x (mod 2),       t,w_i >= 0.

The inequalities are the chamber conditions for the height-zero root system
``A1 + E8``.  This is the same positive-definite Schur-complement recursion as
the Bugaenko enumerator, with one extra ``A1`` coordinate and a parity gate.
Every comparison is exact at both real embeddings.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import lcm
from typing import Any, Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.lattices import e8_cartan_gram

from .bugaenko_vinberg import (
    RationalMatrix,
    Vector,
    _fraction_inverse,
    _integer_upper_square_root,
    _mat_vec,
    _pair,
    _prefix_schur_complements,
    _quadratic,
)


D = ONE + 4 * PHI
CONJUGATE_X_SQUARE_BOUND = 4 * D / 11
RingInteger = tuple[int, int]


def _ring_integer(value: Golden) -> RingInteger:
    if not value.is_integral:
        raise ValueError(f"{value} is not in Z[phi]")
    return value.a.numerator, value.b.numerator


def _from_ring_integer(value: RingInteger) -> Golden:
    return Golden(value[0], value[1])


def _ri_add(left: RingInteger, right: RingInteger) -> RingInteger:
    return left[0] + right[0], left[1] + right[1]


def _ri_subtract(left: RingInteger, right: RingInteger) -> RingInteger:
    return left[0] - right[0], left[1] - right[1]


def _ri_scale(value: RingInteger, scalar: int) -> RingInteger:
    return scalar * value[0], scalar * value[1]


def _ri_multiply(left: RingInteger, right: RingInteger) -> RingInteger:
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def _ri_sign(value: RingInteger, embedding: int = 0) -> int:
    if embedding == 1:
        value = value[0] + value[1], -value[1]
    elif embedding != 0:
        raise ValueError("embedding must be 0 or 1")
    rational_part = 2 * value[0] + value[1]
    radical_coefficient = value[1]
    rational_sign = (rational_part > 0) - (rational_part < 0)
    radical_sign = (radical_coefficient > 0) - (radical_coefficient < 0)
    if rational_sign == 0:
        return radical_sign
    if radical_sign == 0 or rational_sign == radical_sign:
        return rational_sign
    comparison = rational_part * rational_part - 5 * radical_coefficient * radical_coefficient
    return rational_sign if comparison > 0 else radical_sign


def _integral_form(matrix: RationalMatrix) -> tuple[int, tuple[tuple[int, ...], ...]]:
    denominator = lcm(*(entry.denominator for row in matrix for entry in row))
    integral = tuple(
        tuple(int(entry * denominator) for entry in row) for row in matrix
    )
    return denominator, integral


@lru_cache(maxsize=1)
def h11_basis_gram() -> tuple[tuple[Golden, ...], ...]:
    """Return ``[[-2*phi,1],[1,2]] + E8`` in the standard E8 basis."""

    e8 = e8_cartan_gram()
    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][0] = -2 * PHI
    gram[0][1] = gram[1][0] = ONE
    gram[1][1] = Golden(2)
    for i in range(8):
        for j in range(8):
            gram[i + 2][j + 2] = Golden(e8[i][j])
    return tuple(tuple(row) for row in gram)


@lru_cache(maxsize=1)
def _standard_e8_gram() -> RationalMatrix:
    return tuple(
        tuple(Fraction(entry) for entry in row) for row in e8_cartan_gram()
    )


def _unit_vector(index: int, size: int = 8, sign: int = 1) -> tuple[int, ...]:
    return tuple(sign if i == index else 0 for i in range(size))


@lru_cache(maxsize=1)
def standard_e8_roots_from_simple_reflections() -> tuple[tuple[int, ...], ...]:
    """Generate all 240 E8 roots in the standard simple-root basis."""

    e8 = _standard_e8_gram()
    roots = {_unit_vector(i, sign=sign) for i in range(8) for sign in (-1, 1)}
    frontier = list(roots)
    while frontier:
        root = frontier.pop()
        for i in range(8):
            product = sum(
                (Fraction(root[j]) * e8[j][i] for j in range(8)),
                start=Fraction(0),
            )
            reflected = list(root)
            reflected[i] -= int(product)
            candidate = tuple(reflected)
            if candidate not in roots:
                roots.add(candidate)
                frontier.append(candidate)
    ordered = tuple(sorted(roots))
    if len(ordered) != 240:
        raise AssertionError(f"E8 reflection closure has {len(ordered)} roots, not 240")
    if any(
        _quadratic(tuple(Golden(entry) for entry in root), e8) != 2
        for root in ordered
    ):
        raise AssertionError("E8 reflection closure contains a non-root")
    if any(
        not (all(entry >= 0 for entry in root) or all(entry <= 0 for entry in root))
        for root in ordered
    ):
        raise AssertionError("the standard E8 basis is not a simple-root basis")
    return ordered


def _cone_roots() -> tuple[Vector, ...]:
    # Execute the finite reflection closure rather than assuming the E8 claim.
    standard_e8_roots_from_simple_reflections()
    a1 = (ZERO, Golden(-1), *(ZERO for _ in range(8)))
    e8 = tuple(
        (ZERO, ZERO, *tuple(Golden(entry) for entry in _unit_vector(i, sign=-1)))
        for i in range(8)
    )
    return (a1, *e8)


def admissible_control_coordinates(identity_bound: int) -> tuple[Golden, ...]:
    """Enumerate every ``x`` with ``0 < x <= B`` allowed by conjugate norm.

    The conjugate norm equation gives

        x'^2 <= 4/(4*phi-5) = 4*(1+4*phi)/11 < 4.

    Thus ``|x'|<2``.  If ``x=a+b*phi``, then
    ``|b|=|x-x'|/sqrt(5)<(B+2)/2`` and
    ``|a|<=|x'|+|b|*|1-phi|<2+|b|``.  The coefficient box below is
    deliberately one integer wider, followed by exact filters.
    """

    if identity_bound < 1:
        return ()
    b_bound = (identity_bound + 2) // 2 + 1
    a_bound = b_bound + 2
    values = {
        x
        for a in range(-a_bound, a_bound + 1)
        for b in range(-b_bound, b_bound + 1)
        if (x := Golden(a, b)).sign(0) > 0
        and (Golden(identity_bound) - x).sign(0) >= 0
        and (CONJUGATE_X_SQUARE_BOUND - x.conjugate() ** 2).sign(0) >= 0
    }
    return tuple(sorted(values))


@dataclass(frozen=True, slots=True)
class CoordinateAlphabet:
    values: tuple[Golden, ...]
    identity_upper: int
    conjugate_upper: int
    a_bound: int
    b_bound: int


def _coordinate_alphabet(radius: Golden) -> CoordinateAlphabet:
    """All dominant coordinates whose square is bounded by ``radius``."""

    identity_upper = _integer_upper_square_root(radius, 0)
    conjugate_upper = _integer_upper_square_root(radius, 1)
    b_bound = (identity_upper + conjugate_upper) // 2 + 1
    a_bound = conjugate_upper + b_bound
    values = {
        coordinate
        for a in range(-a_bound, a_bound + 1)
        for b in range(-b_bound, b_bound + 1)
        if (coordinate := Golden(a, b)).sign(0) >= 0
        and (radius - coordinate * coordinate).sign(0) >= 0
        and (radius - coordinate * coordinate).sign(1) >= 0
    }
    return CoordinateAlphabet(
        tuple(sorted(values)), identity_upper, conjugate_upper, a_bound, b_bound
    )


@dataclass(frozen=True, slots=True)
class H11HeightLevel:
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


def _reconstruct_root(x: Golden, prefix: Sequence[Golden], inverse_e8: RationalMatrix) -> Vector:
    t, *weights = prefix
    y = (t - x) / 2
    z = _mat_vec(inverse_e8, weights)
    if not y.is_integral or any(not entry.is_integral for entry in z):
        raise AssertionError("positive-definite coordinates did not reconstruct an integral root")
    return (x, y, *z)


def _enumerate_height_level(
    control_coordinate: Golden,
    previous_nonzero_roots: Sequence[Vector],
    inverse_e8: RationalMatrix,
    schur: Sequence[RationalMatrix],
) -> H11HeightLevel:
    x = control_coordinate
    target = 4 + D * x * x
    alphabet = _coordinate_alphabet(target)
    visited = [0] * 10
    schur_pruned = [0] * 10
    wall_pruned = [0] * 10
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
        if wall_t.sign(0) < 0 or any(entry.sign(0) < 0 for entry in wall[2:]):
            raise AssertionError("dominant-wall pruning requires nonnegative coefficients")
        wall_coefficients.append(
            tuple(_ring_integer(entry) for entry in (wall_t, *(2 * entry for entry in wall[2:])))
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
    gram = h11_basis_gram()
    if any(_pair(root, gram, root) != 2 for root in solutions):
        raise AssertionError("height enumerator produced a vector of the wrong norm")
    if any(
        _pair(left, gram, right).sign(0) > 0
        for i, left in enumerate(solutions)
        for right in solutions[i + 1 :]
    ):
        raise AssertionError("same-height Vinberg roots do not meet acutely")
    return H11HeightLevel(
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
class H11VinbergSearch:
    control_vector: Vector
    identity_bound: int
    cone_root_count: int
    cone_roots: tuple[Vector, ...]
    candidate_coordinates: tuple[Golden, ...]
    levels: tuple[H11HeightLevel, ...]
    roots: tuple[Vector, ...]


@dataclass(frozen=True, slots=True)
class H11VinbergExtension:
    """An exact continuation of a previously certified nonzero-height prefix.

    ``previous_nonzero_roots`` is deliberately explicit: an archived prefix
    can be resumed without rerunning its (potentially very large) recursion
    tree.  The caller remains responsible for pinning or independently
    verifying that prefix; this object certifies every newly enumerated
    height and its pruning counters.
    """

    start_after: Golden
    identity_bound: int
    previous_nonzero_root_count: int
    candidate_coordinates: tuple[Golden, ...]
    levels: tuple[H11HeightLevel, ...]
    new_roots: tuple[Vector, ...]


@lru_cache(maxsize=1)
def _positive_definite_enumeration_data() -> tuple[
    RationalMatrix, tuple[RationalMatrix, ...]
]:
    e8 = _standard_e8_gram()
    inverse_e8 = _fraction_inverse(e8)
    if any(entry.denominator != 1 for row in inverse_e8 for entry in row):
        raise AssertionError("E8 inverse must be integral")
    positive_form = (
        (Fraction(1), *(Fraction(0) for _ in range(8))),
        *tuple(
            (
                Fraction(0),
                *(2 * inverse_e8[i][j] for j in range(8)),
            )
            for i in range(8)
        ),
    )
    return inverse_e8, _prefix_schur_complements(positive_form)


def run_h11_vinberg_extension(
    previous_nonzero_roots: Sequence[Vector],
    *,
    start_after: Golden,
    identity_bound: int,
    max_new_levels: int | None = None,
) -> H11VinbergExtension:
    """Resume the exact H11 search strictly after ``start_after``.

    The input roots must be the complete positive-height simple-root prefix
    through ``start_after``.  Basic lattice, dominance, height, and acuteness
    conditions are checked here; completeness of the archived prefix is a
    separate certificate obligation.  ``max_new_levels`` is useful for
    checkpointing a long continuation without weakening any completed level.
    """

    if start_after.sign(0) < 0:
        raise ValueError("start_after must be nonnegative")
    if identity_bound < 1:
        raise ValueError("identity_bound must be positive")
    if max_new_levels is not None and max_new_levels < 1:
        raise ValueError("max_new_levels must be positive when supplied")

    gram = h11_basis_gram()
    previous = tuple(tuple(entry for entry in root) for root in previous_nonzero_roots)
    if any(len(root) != 10 for root in previous):
        raise ValueError("every previous H11 root must have ten coordinates")
    if len(previous) != len(set(previous)):
        raise ValueError("previous H11 root prefix contains duplicates")
    if any(_pair(root, gram, root) != 2 for root in previous):
        raise ValueError("previous H11 prefix contains a non-root")
    if any(root[0].sign(0) <= 0 or (root[0] - start_after).sign(0) > 0 for root in previous):
        raise ValueError("previous H11 roots do not lie in the declared height prefix")
    if any(
        (root[0] + 2 * root[1]).sign(0) < 0
        or any(entry.sign(0) < 0 for entry in root[2:])
        for root in previous
    ):
        raise ValueError("previous H11 roots are not height-zero dominant")
    if any(
        _pair(left, gram, right).sign(0) > 0
        for i, left in enumerate(previous)
        for right in previous[i + 1 :]
    ):
        raise ValueError("previous H11 roots are not pairwise acute")

    candidates = tuple(
        x
        for x in admissible_control_coordinates(identity_bound)
        if (x - start_after).sign(0) > 0
    )
    if max_new_levels is not None:
        candidates = candidates[:max_new_levels]

    inverse_e8, schur = _positive_definite_enumeration_data()
    levels: list[H11HeightLevel] = []
    roots = list(previous)
    for x in candidates:
        level = _enumerate_height_level(x, roots, inverse_e8, schur)
        levels.append(level)
        roots.extend(level.roots)
    new_roots = tuple(root for level in levels for root in level.roots)
    return H11VinbergExtension(
        start_after=start_after,
        identity_bound=identity_bound,
        previous_nonzero_root_count=len(previous),
        candidate_coordinates=candidates,
        levels=tuple(levels),
        new_roots=new_roots,
    )


@lru_cache(maxsize=None)
def run_h11_vinberg_search(identity_bound: int = 2) -> H11VinbergSearch:
    """Exhaust every norm-two Vinberg height with ``0 < x <= identity_bound``."""

    gram = h11_basis_gram()
    control = (Golden(2), Golden(-1), *(ZERO for _ in range(8)))
    if _pair(control, gram, control) != -2 * D:
        raise AssertionError("the control vector must have norm -2*(1+4*phi)")
    control_covector = tuple(
        sum((gram[i][j] * control[j] for j in range(10)), start=ZERO)
        for i in range(10)
    )
    if control_covector != (-D, *(ZERO for _ in range(9))):
        raise AssertionError("the control vector does not diagonalize the binary block")

    inverse_e8, schur = _positive_definite_enumeration_data()
    cone = _cone_roots()
    candidates = admissible_control_coordinates(identity_bound)
    levels: list[H11HeightLevel] = []
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
    return H11VinbergSearch(
        control_vector=control,
        identity_bound=identity_bound,
        cone_root_count=242,
        cone_roots=cone,
        candidate_coordinates=candidates,
        levels=tuple(levels),
        roots=tuple(roots),
    )


def _golden(value: Golden) -> dict[str, str]:
    return value.to_json()


def _vector(vector: Iterable[Golden]) -> list[dict[str, str]]:
    return [_golden(entry) for entry in vector]


def build_h11_vinberg_transcript(identity_bound: int = 2) -> dict[str, Any]:
    """Build a deterministic JSON transcript for the bounded exact search."""

    search = run_h11_vinberg_search(identity_bound)
    return {
        "schema_version": "0.1.0",
        "transcript_kind": "exact_norm_two_vinberg_height_enumeration",
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "lattice": "[[-2*phi,1],[1,2]] + E8",
        "control_vector": _vector(search.control_vector),
        "control_norm": _golden(-2 * D),
        "coordinate_change": {
            "D": _golden(D),
            "t": "x+2*y",
            "w": "E8*z",
            "norm_equation": "t^2+2*w^T*E8^{-1}*w=4+D*x^2",
            "integrality_condition": "t=x modulo 2*Z[phi]",
        },
        "height_formula": "delta=D*x^2/4 for norm-two roots r=(x,y,z)",
        "height_zero": {
            "root_system": "A1+E8",
            "all_root_count": search.cone_root_count,
            "exhaustion": (
                "2*y^2+E8(z)=2 splits into totally nonnegative algebraic "
                "integers y^2 and E8(z)/2 summing to 1; the algebraic norm "
                "forces one summand to be 1 and the other 0, and coefficient "
                "comparison then reduces E8(z)=2 to the 240 ordinary E8 roots"
            ),
            "simple_roots": [_vector(root) for root in search.cone_roots],
        },
        "control_coordinate_scan": {
            "identity_bound": search.identity_bound,
            "conjugate_constraint": "x'^2 <= 4*(1+4*phi)/11",
            "candidates": [_golden(x) for x in search.candidate_coordinates],
        },
        "levels": [
            {
                "control_coordinate": _golden(level.control_coordinate),
                "height": _golden(level.height),
                "positive_definite_target": _golden(level.target),
                "coordinate_alphabet": {
                    "size": len(level.alphabet.values),
                    "identity_upper": level.alphabet.identity_upper,
                    "conjugate_upper": level.alphabet.conjugate_upper,
                    "a_bound": level.alphabet.a_bound,
                    "b_bound": level.alphabet.b_bound,
                },
                "visited_by_depth": list(level.visited_by_depth),
                "parity_pruned": level.parity_pruned,
                "schur_pruned_by_depth": list(level.schur_pruned_by_depth),
                "wall_pruned_by_depth": list(level.wall_pruned_by_depth),
                "leaf_norm_mismatches": level.leaf_norm_mismatches,
                "solution_count": len(level.roots),
                "roots": [_vector(root) for root in level.roots],
            }
            for level in search.levels
        ],
        "output": {
            "wall_count": len(search.roots),
            "root_order": (
                "height_zero.simple_roots followed by the roots arrays of "
                "levels in listed order"
            ),
        },
        "status": {
            "kind": "bounded_prefix",
            "last_processed_control_coordinate": (
                _golden(search.candidate_coordinates[-1])
                if search.candidate_coordinates
                else None
            ),
            "gap": (
                "a finite-volume chamber certificate is required before "
                "Vinberg's stopping criterion applies"
            ),
        },
    }
