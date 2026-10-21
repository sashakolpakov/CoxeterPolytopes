"""Bounded exact Vinberg enumeration for ``U + 2 E8(-1) + A2(-1)``.

This is the integral rank-20 lattice of signature ``(1, 19)`` that occurs in
Vinberg's and Looijenga's cubic-fourfold chamber.  We use coordinates

    r = (a, b; v_1, v_2, z),

where ``(a,b)`` is the hyperbolic-plane coordinate, the ``v_i`` are in the
standard simple-root basis of ``E8``, and ``z`` is in the simple-root basis
of ``A2``.  Thus

    r^2 = 2*a*b - v_1^2 - v_2^2 - z^2.

There are two integral reflection types:

* primitive roots of square ``-2`` and divisibility one;
* primitive roots of square ``-6`` and divisibility three.

For the control point ``p=e+f``, put ``n=-(p,r)=-(a+b)>0`` and
``d=a-b``.  At a fixed root square ``-q`` (``q`` is 2 or 6), the positive
definite equation is

    d^2 + 2*(v_1^2 + v_2^2 + z^2) = n^2 + 2*q,
    d == n (mod 2).

The complete height-zero reflection system is ``A1 + E8 + E8 + G2``.
Choosing its displayed simple chamber gives ``d <= 0``, dominant ``E8``
weights, and ``A2`` weights ``0 <= w_1 <= w_2``.  Consequently every
positive level is a finite enumeration of dominant weights.  Vinberg levels
are ordered by the exact distance invariant ``n^2/q``; this is essential
because the two reflection types interleave.

The module deliberately separates bounded enumeration from stopping.  An
exact incremental extreme-ray computation is used to decide whether the
current outer chamber already has finite volume.  If it does, Vinberg's
stopping theorem makes the accepted prefix complete.  Failure at a bounded
prefix is only diagnostic: later roots may cut off its hyperideal rays.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd, isqrt
from typing import Any, Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.lattices import block_diagonal, e8_cartan_gram, hyperbolic_plane
from coxeter_polytopes.search.bugaenko_vinberg import _fraction_inverse
from coxeter_polytopes.search.h11_vinberg import (
    standard_e8_roots_from_simple_reflections,
)
from coxeter_polytopes.verify.chamber import ChamberVerification, diagnose_chamber_cone


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def _a2_cartan_gram() -> list[list[int]]:
    return [[2, -1], [-1, 2]]


@lru_cache(maxsize=1)
def s19_basis_gram() -> Matrix:
    """Return ``U + 2 E8(-1) + A2(-1)`` in the declared basis."""

    e8_negative = [[-entry for entry in row] for row in e8_cartan_gram()]
    a2_negative = [[-entry for entry in row] for row in _a2_cartan_gram()]
    return tuple(
        tuple(row)
        for row in block_diagonal(
            hyperbolic_plane(), e8_negative, e8_negative, a2_negative
        )
    )


def _mat_vec(matrix: Sequence[Sequence[int]], vector: Sequence[int]) -> tuple[int, ...]:
    return tuple(sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix)


def pair(left: Sequence[int], right: Sequence[int]) -> int:
    """Pair two vectors in the integral S19 lattice."""

    gram = s19_basis_gram()
    return sum(left[i] * gram[i][j] * right[j] for i in range(20) for j in range(20))


def root_divisibility(root: Sequence[int]) -> int:
    """Return the positive generator of ``(root,S19)``."""

    covector = _mat_vec(s19_basis_gram(), root)
    result = 0
    for entry in covector:
        result = gcd(result, abs(entry))
    return result


def _primitive(vector: Sequence[int]) -> bool:
    content = 0
    for entry in vector:
        content = gcd(content, abs(entry))
    return content == 1


def _unit(index: int, size: int) -> tuple[int, ...]:
    return tuple(1 if coordinate == index else 0 for coordinate in range(size))


@lru_cache(maxsize=1)
def height_zero_roots() -> tuple[Vector, ...]:
    """Construct all 494 roots in the finite ``A1+E8+E8+G2`` system."""

    zero18 = (0,) * 18
    roots: list[Vector] = [(1, -1, *zero18), (-1, 1, *zero18)]
    e8_roots = standard_e8_roots_from_simple_reflections()
    for block in (0, 1):
        for root in e8_roots:
            negative = [0] * 18
            negative[8 * block : 8 * block + 8] = root
            roots.append((0, 0, *negative))

    a2 = _a2_cartan_gram()
    for x in range(-3, 4):
        for y in range(-3, 4):
            z = (x, y)
            square = sum(z[i] * a2[i][j] * z[j] for i in range(2) for j in range(2))
            weights = _mat_vec(a2, z)
            if square == 2 or (square == 6 and gcd(abs(weights[0]), abs(weights[1])) == 3):
                roots.append((0, 0, *(0,) * 16, x, y))

    result = tuple(sorted(set(roots)))
    if len(result) != 494:
        raise AssertionError(f"height-zero closure has {len(result)} roots, not 494")
    type_counts = {
        (-pair(root, root), root_divisibility(root)) for root in result
    }
    if type_counts != {(2, 1), (6, 3)}:
        raise AssertionError(f"unexpected height-zero reflection types {type_counts}")
    return result


@lru_cache(maxsize=1)
def height_zero_simple_roots() -> tuple[Vector, ...]:
    """Return a simple system of type ``A1+E8+E8+G2``.

    In the final ``A2`` block the short simple root is ``-alpha_1`` and the
    long simple root is ``alpha_1-alpha_2``.  Their lattice Gram is
    ``[[-2,3],[3,-6]]``, hence Coxeter label six after sign reversal.
    """

    roots: list[Vector] = [(1, -1, *(0,) * 18)]
    for block in (0, 1):
        for index in range(8):
            finite = [0] * 18
            finite[8 * block + index] = -1
            roots.append((0, 0, *finite))
    roots.append((0, 0, *(0,) * 16, -1, 0))
    roots.append((0, 0, *(0,) * 16, 1, -1))
    result = tuple(roots)
    if len(result) != 19 or any(root not in height_zero_roots() for root in result):
        raise AssertionError("declared height-zero simple system is invalid")
    if any(pair(result[i], result[j]) < 0 for i in range(len(result)) for j in range(i)):
        raise AssertionError("height-zero simple roots are not pairwise acute")
    return result


@dataclass(frozen=True, slots=True)
class DominantVector:
    square: int
    coordinates: tuple[int, ...]
    weights: tuple[int, ...]


@lru_cache(maxsize=None)
def dominant_e8_vectors(max_square: int) -> tuple[DominantVector, ...]:
    """Enumerate all dominant E8 vectors through an exact norm bound.

    E8 is unimodular, so every nonnegative integral Dynkin-weight vector is
    represented by an integral lattice vector.  The inverse Cartan matrix has
    nonnegative integral entries; this makes the recursion monotone and the
    finite coefficient bounds exact.
    """

    if max_square < 0:
        return ()
    cartan = e8_cartan_gram()
    inverse_fraction = _fraction_inverse(
        tuple(tuple(Fraction(entry) for entry in row) for row in cartan)
    )
    if any(entry.denominator != 1 or entry < 0 for row in inverse_fraction for entry in row):
        raise AssertionError("the standard E8 inverse must be nonnegative integral")
    inverse = tuple(tuple(int(entry) for entry in row) for row in inverse_fraction)
    bounds = tuple(isqrt(max_square // inverse[i][i]) for i in range(8))
    output: list[DominantVector] = []

    def extend(weights: tuple[int, ...]) -> None:
        depth = len(weights)
        partial_square = sum(
            weights[i] * inverse[i][j] * weights[j]
            for i in range(depth)
            for j in range(depth)
        )
        if partial_square > max_square:
            return
        if depth == 8:
            coordinates = tuple(
                sum(inverse[i][j] * weights[j] for j in range(8)) for i in range(8)
            )
            output.append(DominantVector(partial_square, coordinates, weights))
            return
        for entry in range(bounds[depth] + 1):
            extend((*weights, entry))

    extend(())
    result = tuple(sorted(output, key=lambda vector: (vector.square, vector.weights)))
    if any(
        _mat_vec(cartan, vector.coordinates) != vector.weights
        or sum(
            vector.coordinates[i] * cartan[i][j] * vector.coordinates[j]
            for i in range(8)
            for j in range(8)
        )
        != vector.square
        for vector in result
    ):
        raise AssertionError("dominant E8 enumerator failed reconstruction")
    return result


@lru_cache(maxsize=None)
def dominant_a2_vectors(max_square: int) -> tuple[DominantVector, ...]:
    """Enumerate integral A2 vectors in the selected G2 chamber."""

    if max_square < 0:
        return ()
    # z^2 = (2/3)(w1^2+w1*w2+w2^2), with w1==w2 mod 3.
    weight_bound = isqrt(3 * max_square // 2) + 1
    result: list[DominantVector] = []
    for w1 in range(weight_bound + 1):
        for w2 in range(w1, weight_bound + 1):
            if (w1 - w2) % 3:
                continue
            numerator = 2 * (w1 * w1 + w1 * w2 + w2 * w2)
            if numerator % 3:
                raise AssertionError("A2 congruence did not make the norm integral")
            square = numerator // 3
            if square > max_square:
                continue
            coordinates = ((2 * w1 + w2) // 3, (w1 + 2 * w2) // 3)
            result.append(DominantVector(square, coordinates, (w1, w2)))
    return tuple(sorted(result, key=lambda vector: (vector.square, vector.weights)))


@lru_cache(maxsize=None)
def dominant_a2_divisible_three_vectors(
    max_square: int,
) -> tuple[DominantVector, ...]:
    """Enumerate dominant A2 vectors whose pairing ideal is divisible by 3."""

    if max_square < 0:
        return ()
    # Write the Dynkin weights as (3*s,3*t).  Then the simple-root
    # coordinates are (2*s+t,s+2*t), and the square is
    # 6*(s^2+s*t+t^2).  This avoids enumerating the overwhelmingly irrelevant
    # non-divisible A2 weights on late norm-six Vinberg levels.
    reduced_bound = isqrt(max_square // 6) + 1
    result: list[DominantVector] = []
    for s in range(reduced_bound + 1):
        for t in range(s, reduced_bound + 1):
            square = 6 * (s * s + s * t + t * t)
            if square > max_square:
                continue
            result.append(
                DominantVector(
                    square=square,
                    coordinates=(2 * s + t, s + 2 * t),
                    weights=(3 * s, 3 * t),
                )
            )
    return tuple(sorted(result, key=lambda vector: (vector.square, vector.weights)))


def _buckets(
    vectors: Iterable[DominantVector],
) -> dict[int, tuple[DominantVector, ...]]:
    mutable: dict[int, list[DominantVector]] = {}
    for vector in vectors:
        mutable.setdefault(vector.square, []).append(vector)
    return {square: tuple(values) for square, values in mutable.items()}


@dataclass(frozen=True, slots=True)
class S19Candidate:
    root: Vector
    norm_abs: int
    divisibility: int
    control: int
    difference: int
    finite_square: int

    @property
    def height(self) -> Fraction:
        return Fraction(self.control * self.control, self.norm_abs)


@lru_cache(maxsize=None)
def enumerate_level_candidates(control: int, norm_abs: int) -> tuple[S19Candidate, ...]:
    """Exhaust one ``(n,q)`` level before lower-wall rejection."""

    if control <= 0:
        raise ValueError("control must be positive")
    if norm_abs not in (2, 6):
        raise ValueError("the supported reflection norms are -2 and -6")
    if norm_abs == 6 and control % 3:
        return ()  # divisibility three forces a+b, hence n, to be divisible by 3

    radius = control * control + 2 * norm_abs
    maximum_finite_square = radius // 2
    if norm_abs == 6:
        # Divisibility three and E8 unimodularity force each E8 vector to be
        # three times an E8 vector.  Enforcing that before the product census
        # reduces a late shell by several orders of magnitude.
        e8 = _buckets(
            DominantVector(
                square=9 * vector.square,
                coordinates=tuple(3 * entry for entry in vector.coordinates),
                weights=tuple(3 * entry for entry in vector.weights),
            )
            for vector in dominant_e8_vectors(maximum_finite_square // 9)
        )
        a2 = _buckets(dominant_a2_divisible_three_vectors(maximum_finite_square))
    else:
        e8 = _buckets(dominant_e8_vectors(maximum_finite_square))
        a2 = _buckets(dominant_a2_vectors(maximum_finite_square))
    candidates: list[S19Candidate] = []

    for difference in range(-isqrt(radius), 1):
        if (difference - control) % 2:
            continue
        if norm_abs == 6 and difference % 3:
            continue
        numerator = radius - difference * difference
        if numerator < 0 or numerator % 2:
            continue
        finite_square = numerator // 2
        a = (-control + difference) // 2
        b = (-control - difference) // 2
        for square1, vectors1 in e8.items():
            for square2, vectors2 in e8.items():
                square_a2 = finite_square - square1 - square2
                if square_a2 not in a2:
                    continue
                for vector1 in vectors1:
                    for vector2 in vectors2:
                        for vector_a2 in a2[square_a2]:
                            root = (
                                a,
                                b,
                                *vector1.coordinates,
                                *vector2.coordinates,
                                *vector_a2.coordinates,
                            )
                            if not _primitive(root):
                                continue
                            divisibility = root_divisibility(root)
                            required_divisibility = 1 if norm_abs == 2 else 3
                            if divisibility != required_divisibility:
                                continue
                            if pair(root, root) != -norm_abs:
                                raise AssertionError("positive-definite equation reconstructed wrong norm")
                            candidates.append(
                                S19Candidate(
                                    root=root,
                                    norm_abs=norm_abs,
                                    divisibility=divisibility,
                                    control=control,
                                    difference=difference,
                                    finite_square=finite_square,
                                )
                            )

    result = tuple(sorted(candidates, key=lambda candidate: candidate.root))
    cone = height_zero_simple_roots()
    if any(pair(candidate.root, wall) < 0 for candidate in result for wall in cone):
        raise AssertionError("dominant enumeration violated a height-zero wall")
    return result


@dataclass(frozen=True, slots=True)
class S19HeightLevel:
    control: int
    norm_abs: int
    height: Fraction
    candidate_count: int
    lower_wall_rejected_count: int
    roots: tuple[Vector, ...]


@dataclass(frozen=True, slots=True)
class S19StoppingAudit:
    attempted: bool
    fired: bool
    finite_vertex_count: int | None
    cusp_count: int | None
    hyperideal_ray_count: int | None
    opposite_sheet_ray_count: int | None
    geometry: ChamberVerification | None


@dataclass(frozen=True, slots=True)
class S19VinbergSearch:
    control_bound: int
    height_zero_root_count: int
    cone_roots: tuple[Vector, ...]
    levels: tuple[S19HeightLevel, ...]
    roots: tuple[Vector, ...]
    stopping: S19StoppingAudit


def _ordered_level_keys(control_bound: int) -> tuple[tuple[int, int], ...]:
    """All mixed-type levels through height ``control_bound^2/2``.

    The public bound is the largest square-minus-two control coordinate.
    Square-minus-six roots of control as large as ``floor(sqrt(3)*B)`` can
    occur earlier in Vinberg order and therefore must also be exhausted.
    """

    keys = [(control, 2) for control in range(1, control_bound + 1)]
    long_control_bound = isqrt(3 * control_bound * control_bound)
    keys.extend((control, 6) for control in range(3, long_control_bound + 1, 3))
    return tuple(sorted(keys, key=lambda item: (Fraction(item[0] ** 2, item[1]), item)))


def _diagnostic_interior(roots: Sequence[Vector]) -> tuple[Golden, ...]:
    """Construct an exact strict interior point near ``p=e+f``."""

    inverse_e8 = _fraction_inverse(
        tuple(tuple(Fraction(entry) for entry in row) for row in e8_cartan_gram())
    )
    anti_weyl = tuple(-sum(row) for row in inverse_e8)
    if any(entry.denominator != 1 for entry in anti_weyl):
        raise AssertionError("E8 anti-Weyl vector should be integral")
    h = (1, -1, *(int(entry) for entry in anti_weyl), *(int(entry) for entry in anti_weyl), -2, -3)

    # The chamber verifier uses the sign-reversed lattice and roots.  Search
    # for a small exact perturbation that remains timelike and is strict on
    # every wall; no floating-point feasibility decision is involved.
    diagnostic_gram = tuple(
        tuple(Golden(-entry) for entry in row) for row in s19_basis_gram()
    )
    diagnostic_roots = tuple(tuple(Golden(-entry) for entry in root) for root in roots)
    for denominator in (10, 100, 1_000, 10_000, 100_000):
        point = (Golden(1) + Golden(h[0]) / denominator, Golden(1) + Golden(h[1]) / denominator, *tuple(Golden(entry) / denominator for entry in h[2:]))
        norm = sum(point[i] * diagnostic_gram[i][j] * point[j] for i in range(20) for j in range(20))
        wall_values = tuple(
            sum(point[i] * diagnostic_gram[i][j] * root[j] for i in range(20) for j in range(20))
            for root in diagnostic_roots
        )
        if norm.sign() < 0 and all(value.sign() < 0 for value in wall_values):
            return point
    raise AssertionError("failed to construct a strict rational chamber point")


def audit_stopping(roots: Sequence[Vector]) -> S19StoppingAudit:
    """Run the exact finite-volume outer-chamber test on an accepted prefix."""

    if len(roots) < 20:
        return S19StoppingAudit(False, False, None, None, None, None, None)
    gram = tuple(tuple(Golden(-entry) for entry in row) for row in s19_basis_gram())
    normals = tuple(tuple(Golden(-entry) for entry in root) for root in roots)
    interior = _diagnostic_interior(roots)
    geometry = diagnose_chamber_cone(normals, gram, interior, dimension=19)
    fired = not geometry.hyperideal_rays and not geometry.opposite_sheet_rays
    return S19StoppingAudit(
        attempted=True,
        fired=fired,
        finite_vertex_count=len(geometry.finite_vertices),
        cusp_count=len(geometry.cusps),
        hyperideal_ray_count=len(geometry.hyperideal_rays),
        opposite_sheet_ray_count=len(geometry.opposite_sheet_rays),
        geometry=geometry,
    )


@lru_cache(maxsize=None)
def run_s19_vinberg_search(
    control_bound: int = 3, *, attempt_stopping: bool = False
) -> S19VinbergSearch:
    """Enumerate the complete mixed prefix through height ``B^2/2``.

    Here ``B=control_bound`` is the largest processed square-minus-two
    control.  Earlier square-minus-six levels with larger control are included
    automatically.
    """

    if control_bound < 0:
        raise ValueError("control_bound must be nonnegative")
    roots = list(height_zero_simple_roots())
    levels: list[S19HeightLevel] = []
    for control, norm_abs in _ordered_level_keys(control_bound):
        candidates = enumerate_level_candidates(control, norm_abs)
        accepted = tuple(
            candidate.root
            for candidate in candidates
            if all(pair(candidate.root, wall) >= 0 for wall in roots)
        )
        if any(pair(left, right) < 0 for i, left in enumerate(accepted) for right in accepted[i + 1 :]):
            raise AssertionError("same-height accepted roots are not pairwise acute")
        roots.extend(accepted)
        levels.append(
            S19HeightLevel(
                control=control,
                norm_abs=norm_abs,
                height=Fraction(control * control, norm_abs),
                candidate_count=len(candidates),
                lower_wall_rejected_count=len(candidates) - len(accepted),
                roots=accepted,
            )
        )
    stopping = audit_stopping(roots) if attempt_stopping else S19StoppingAudit(
        False, False, None, None, None, None, None
    )
    return S19VinbergSearch(
        control_bound=control_bound,
        height_zero_root_count=len(height_zero_roots()),
        cone_roots=height_zero_simple_roots(),
        levels=tuple(levels),
        roots=tuple(roots),
        stopping=stopping,
    )


def _fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def build_s19_vinberg_transcript(
    control_bound: int = 3, *, attempt_stopping: bool = False
) -> dict[str, Any]:
    """Build a compact deterministic transcript of the bounded search."""

    search = run_s19_vinberg_search(control_bound, attempt_stopping=attempt_stopping)
    norm_counts = {
        str(norm_abs): sum(pair(root, root) == -norm_abs for root in search.roots)
        for norm_abs in (2, 6)
    }
    return {
        "schema_version": "0.1.0",
        "transcript_kind": "bounded_exact_integral_mixed_root_vinberg_enumeration",
        "lattice": "U + 2E8(-1) + A2(-1)",
        "dimension": 19,
        "control_vector": [1, 1, *([0] * 18)],
        "height_zero": {
            "root_system": "A1 + E8 + E8 + G2",
            "all_root_count": search.height_zero_root_count,
            "simple_root_count": len(search.cone_roots),
            "simple_roots": [list(root) for root in search.cone_roots],
        },
        "reflection_types": [
            {"square": -2, "divisibility": 1},
            {"square": -6, "divisibility": 3},
        ],
        "height_equation": "d^2+2*K(v)=n^2+2*q; height=n^2/q",
        "control_bound": search.control_bound,
        "height_cutoff": _fraction(Fraction(search.control_bound**2, 2)),
        "levels": [
            {
                "control": level.control,
                "root_square": -level.norm_abs,
                "height": _fraction(level.height),
                "candidate_count": level.candidate_count,
                "lower_wall_rejected_count": level.lower_wall_rejected_count,
                "accepted_count": len(level.roots),
                "roots": [list(root) for root in level.roots],
            }
            for level in search.levels
        ],
        "output": {
            "wall_count": len(search.roots),
            "wall_count_by_abs_norm": norm_counts,
            "roots_in_vinberg_order": [list(root) for root in search.roots],
        },
        "stopping": {
            "attempted": search.stopping.attempted,
            "fired": search.stopping.fired,
            "finite_vertex_count": search.stopping.finite_vertex_count,
            "cusp_count": search.stopping.cusp_count,
            "hyperideal_ray_count": search.stopping.hyperideal_ray_count,
            "opposite_sheet_ray_count": search.stopping.opposite_sheet_ray_count,
            "interpretation": (
                "finite-volume outer chamber; Vinberg stopping applies"
                if search.stopping.fired
                else "bounded prefix only; later roots may cut current escape rays"
            ),
        },
    }
