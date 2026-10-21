"""Exact norm-two Vinberg prefix for the discriminant-19 field lattice.

The lattice is

    H19 = [[-2*phi, phi], [phi, 2]] + E8

over ``Z[phi]``.  The integral timelike control vector
``p=(2,-phi,0,...,0)`` diagonalizes the coupled binary block.  Put

    D = 1 + 5*phi,       t = phi*x + 2*y,       w = E8*z

for a vector ``r=(x,y,z)``.  Then

    <p,p> = -2*D,        <p,r> = -D*x,
    2*<r,r> = t^2 - D*x^2 + 2*w^T*E8^{-1}*w.

Thus norm-two roots of fixed positive control coordinate ``x`` are exactly
the integral solutions of

    t^2 + 2*w^T*E8^{-1}*w = 4 + D*x^2,
    t == phi*x (mod 2),       t,w_i >= 0.

The inequalities are the dominant-chamber conditions for the height-zero
root system ``A1 + E8``.  The enumeration below is exhaustive through the
declared identity-embedding bound and makes no termination claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Any, Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, PHI, ZERO

from .bugaenko_vinberg import (
    RationalMatrix,
    Vector,
    _fraction_inverse,
    _mat_vec,
    _pair,
    _prefix_schur_complements,
    _quadratic,
)
from .h11_vinberg import (
    CoordinateAlphabet,
    _coordinate_alphabet,
    _standard_e8_gram,
    _unit_vector,
    standard_e8_roots_from_simple_reflections,
)


D = Golden(1) + 5 * PHI
CONJUGATE_X_SQUARE_BOUND = 4 * D / 19


@lru_cache(maxsize=1)
def h19_basis_gram() -> tuple[tuple[Golden, ...], ...]:
    """Return ``[[-2*phi,phi],[phi,2]] + E8`` in the standard basis."""

    e8 = _standard_e8_gram()
    gram = [[ZERO for _ in range(10)] for _ in range(10)]
    gram[0][0] = -2 * PHI
    gram[0][1] = gram[1][0] = PHI
    gram[1][1] = Golden(2)
    for i in range(8):
        for j in range(8):
            gram[i + 2][j + 2] = Golden(e8[i][j])
    return tuple(tuple(row) for row in gram)


@lru_cache(maxsize=1)
def height_zero_roots() -> tuple[Vector, ...]:
    """Construct all 242 roots in the complete height-zero ``A1 + E8``."""

    e8_roots = standard_e8_roots_from_simple_reflections()
    roots = (
        (ZERO, Golden(-1), *(ZERO for _ in range(8))),
        (ZERO, Golden(1), *(ZERO for _ in range(8))),
        *tuple(
            (ZERO, ZERO, *tuple(Golden(entry) for entry in root))
            for root in e8_roots
        ),
    )
    gram = h19_basis_gram()
    if len(roots) != 242 or len(set(roots)) != 242:
        raise AssertionError("height-zero A1+E8 root construction is incomplete")
    if any(_pair(root, gram, root) != 2 for root in roots):
        raise AssertionError("height-zero construction contains a non-root")
    return roots


def _cone_roots() -> tuple[Vector, ...]:
    """Simple roots of the complete height-zero system ``A1 + E8``."""

    # Construct the entire system before selecting its dominant simple basis.
    height_zero_roots()
    a1 = (ZERO, Golden(-1), *(ZERO for _ in range(8)))
    e8 = tuple(
        (ZERO, ZERO, *tuple(Golden(entry) for entry in _unit_vector(i, sign=-1)))
        for i in range(8)
    )
    return (a1, *e8)


def admissible_control_coordinates(identity_bound: int) -> tuple[Golden, ...]:
    """Enumerate every admissible ``x`` with ``0 < x <= identity_bound``.

    At the conjugate embedding the positive-definite norm equation gives

        x'^2 <= 4/(5*phi-6) = 4*(1+5*phi)/19 < 2.

    The coefficient box is deliberately wider than the resulting bounds and
    is followed by exact two-embedding filters.
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
class H19HeightLevel:
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


def _reconstruct_root(
    x: Golden, prefix: Sequence[Golden], inverse_e8: RationalMatrix
) -> Vector:
    t, *weights = prefix
    y = (t - PHI * x) / 2
    z = _mat_vec(inverse_e8, weights)
    if not y.is_integral or any(not entry.is_integral for entry in z):
        raise AssertionError("positive-definite coordinates did not reconstruct an integral root")
    return (x, y, *z)


def _enumerate_height_level(
    control_coordinate: Golden,
    previous_nonzero_roots: Sequence[Vector],
    inverse_e8: RationalMatrix,
    schur: Sequence[RationalMatrix],
) -> H19HeightLevel:
    x = control_coordinate
    target = 4 + D * x * x
    alphabet = _coordinate_alphabet(target)
    visited = [0] * 10
    schur_pruned = [0] * 10
    wall_pruned = [0] * 10
    parity_pruned = 0
    leaf_norm_mismatches = 0
    solutions: list[Vector] = []

    for wall in previous_nonzero_roots:
        wall_t = PHI * wall[0] + 2 * wall[1]
        if wall_t.sign(0) < 0 or any(entry.sign(0) < 0 for entry in wall[2:]):
            raise AssertionError("dominant-wall pruning requires nonnegative coordinates")

    def extend(prefix: tuple[Golden, ...]) -> None:
        nonlocal parity_pruned, leaf_norm_mismatches
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            minimum = _quadratic(prefix, schur[depth - 1])
            remainder = target - minimum
            if remainder.sign(0) < 0 or remainder.sign(1) < 0:
                schur_pruned[depth] += 1
                return
            for wall in previous_nonzero_roots:
                wall_t = PHI * wall[0] + 2 * wall[1]
                partial = prefix[0] * wall_t + 2 * sum(
                    (prefix[i] * wall[i + 1] for i in range(1, depth)),
                    start=ZERO,
                )
                bound = D * x * wall[0]
                if (bound - partial).sign(0) < 0:
                    wall_pruned[depth] += 1
                    return
        if depth == 9:
            if _quadratic(prefix, schur[-1]) != target:
                leaf_norm_mismatches += 1
                return
            solutions.append(_reconstruct_root(x, prefix, inverse_e8))
            return
        for entry in alphabet.values:
            if depth == 0 and not ((entry - PHI * x) / 2).is_integral:
                parity_pruned += 1
                continue
            extend((*prefix, entry))

    extend(())
    gram = h19_basis_gram()
    if any(_pair(root, gram, root) != 2 for root in solutions):
        raise AssertionError("height enumerator produced a vector of the wrong norm")
    if any(
        _pair(left, gram, right).sign(0) > 0
        for i, left in enumerate(solutions)
        for right in solutions[i + 1 :]
    ):
        raise AssertionError("same-height Vinberg roots do not meet acutely")
    return H19HeightLevel(
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
class H19VinbergSearch:
    control_vector: Vector
    identity_bound: int
    cone_root_count: int
    cone_roots: tuple[Vector, ...]
    candidate_coordinates: tuple[Golden, ...]
    levels: tuple[H19HeightLevel, ...]
    roots: tuple[Vector, ...]


@lru_cache(maxsize=None)
def run_h19_vinberg_search(identity_bound: int = 2) -> H19VinbergSearch:
    """Exhaust every norm-two Vinberg height with ``0 < x <= identity_bound``."""

    gram = h19_basis_gram()
    control = (Golden(2), -PHI, *(ZERO for _ in range(8)))
    if _pair(control, gram, control) != -2 * D:
        raise AssertionError("the control vector must have norm -2*(1+5*phi)")
    control_covector = tuple(
        sum((gram[i][j] * control[j] for j in range(10)), start=ZERO)
        for i in range(10)
    )
    if control_covector != (-D, *(ZERO for _ in range(9))):
        raise AssertionError("the control vector does not diagonalize the binary block")

    e8 = _standard_e8_gram()
    inverse_e8 = _fraction_inverse(e8)
    if any(entry.denominator != 1 for row in inverse_e8 for entry in row):
        raise AssertionError("E8 inverse must be integral")
    positive_form = (
        (Fraction(1), *(Fraction(0) for _ in range(8))),
        *tuple(
            (Fraction(0), *(2 * inverse_e8[i][j] for j in range(8)))
            for i in range(8)
        ),
    )
    schur = _prefix_schur_complements(positive_form)
    cone = _cone_roots()
    candidates = admissible_control_coordinates(identity_bound)
    levels: list[H19HeightLevel] = []
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
    return H19VinbergSearch(
        control_vector=control,
        identity_bound=identity_bound,
        cone_root_count=len(height_zero_roots()),
        cone_roots=cone,
        candidate_coordinates=candidates,
        levels=tuple(levels),
        roots=tuple(roots),
    )


def _golden(value: Golden) -> dict[str, str]:
    return value.to_json()


def _vector(vector: Iterable[Golden]) -> list[dict[str, str]]:
    return [_golden(entry) for entry in vector]


def build_h19_vinberg_transcript(identity_bound: int = 2) -> dict[str, Any]:
    """Build a deterministic JSON transcript for the bounded exact search."""

    search = run_h19_vinberg_search(identity_bound)
    return {
        "schema_version": "0.1.0",
        "transcript_kind": "exact_norm_two_vinberg_height_enumeration",
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "lattice": "[[-2*phi,phi],[phi,2]] + E8",
        "control_vector": _vector(search.control_vector),
        "control_norm": _golden(-2 * D),
        "coordinate_change": {
            "D": _golden(D),
            "t": "phi*x+2*y",
            "w": "E8*z",
            "norm_equation": "t^2+2*w^T*E8^{-1}*w=4+D*x^2",
            "integrality_condition": "t=phi*x modulo 2*Z[phi]",
        },
        "height_formula": "delta=D*x^2/4 for norm-two roots r=(x,y,z)",
        "height_zero": {
            "root_system": "A1+E8",
            "all_root_count": search.cone_root_count,
            "exhaustion": (
                "2*y^2+E8(z)=2 splits into totally nonnegative algebraic "
                "integers y^2 and E8(z)/2 summing to 1; the algebraic norm "
                "forces one summand to be 1 and the other 0, and coefficient "
                "comparison reduces E8(z)=2 to the 240 ordinary E8 roots"
            ),
            "simple_roots": [_vector(root) for root in search.cone_roots],
        },
        "control_coordinate_scan": {
            "identity_bound": search.identity_bound,
            "conjugate_constraint": "x'^2 <= 4*(1+5*phi)/19",
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
            "roots_in_vinberg_order": [_vector(root) for root in search.roots],
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
