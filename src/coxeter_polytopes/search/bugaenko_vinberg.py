"""A small exact Vinberg enumeration for Bugaenko's rank-nine lattice.

The search uses only the lattice Gram matrix, not the archived simple roots.
Its control vector is the generator ``p`` of the ``<-2*phi>`` summand.  After
normalizing reflective vectors by squared units, every root hyperplane has a
norm-two representative ``r=(x,y)``.  The Vinberg distance parameter is

    delta(r) = -<p,r>^2 / (<p,p><r,r>) = phi*x^2.

At height zero the roots form the finite E8 root system.  Orienting its simple
roots as ``-e_i`` changes the remaining acute-angle conditions into
coordinatewise nonnegativity for ``z=E8*y``.  For fixed ``x`` the norm equation
is then

    z^T E8^{-1} z = 2 + 2*phi*x^2.

The right-hand side and its conjugate bound every coordinate.  The recursive
enumerator below is exhaustive: a rational Schur complement gives the exact
minimum remaining quadratic energy at both real embeddings, while already
found walls give exact linear pruning.  No floating-point comparison occurs.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from typing import Any, Iterable, Sequence

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.known_examples.bugaenko import bugaenko_basis_gram


Vector = tuple[Golden, ...]
RationalMatrix = tuple[tuple[Fraction, ...], ...]


def _pair(
    left: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    right: Sequence[Golden],
) -> Golden:
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        start=ZERO,
    )


def _fraction_inverse(matrix: RationalMatrix) -> RationalMatrix:
    size = len(matrix)
    work = [
        [*matrix[i], *(Fraction(int(i == j)) for j in range(size))]
        for i in range(size)
    ]
    for column in range(size):
        pivot = next((i for i in range(column, size) if work[i][column]), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            multiplier = work[row][column]
            work[row] = [
                work[row][j] - multiplier * work[column][j]
                for j in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


@lru_cache(maxsize=1)
def _e8_gram() -> RationalMatrix:
    ambient = bugaenko_basis_gram()
    return tuple(
        tuple(Fraction(ambient[i + 1][j + 1].a) for j in range(8))
        for i in range(8)
    )


def _matmul(
    left: RationalMatrix, right: RationalMatrix
) -> RationalMatrix:
    return tuple(
        tuple(
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                start=Fraction(0),
            )
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def _transpose(matrix: RationalMatrix) -> RationalMatrix:
    return tuple(tuple(matrix[i][j] for i in range(len(matrix))) for j in range(len(matrix[0])))


def _subtract(left: RationalMatrix, right: RationalMatrix) -> RationalMatrix:
    return tuple(
        tuple(left[i][j] - right[i][j] for j in range(len(left[0])))
        for i in range(len(left))
    )


def _prefix_schur_complements(matrix: RationalMatrix) -> tuple[RationalMatrix, ...]:
    """Minimum-energy forms after fixing each nonempty coordinate prefix."""

    size = len(matrix)
    output: list[RationalMatrix] = []
    for prefix in range(1, size + 1):
        a = tuple(tuple(matrix[i][j] for j in range(prefix)) for i in range(prefix))
        if prefix == size:
            output.append(a)
            continue
        b = tuple(tuple(matrix[i][j] for j in range(prefix, size)) for i in range(prefix))
        c = tuple(
            tuple(matrix[i][j] for j in range(prefix, size))
            for i in range(prefix, size)
        )
        output.append(_subtract(a, _matmul(_matmul(b, _fraction_inverse(c)), _transpose(b))))
    return tuple(output)


def _quadratic(vector: Sequence[Golden], matrix: RationalMatrix) -> Golden:
    return sum(
        (
            vector[i] * Golden(matrix[i][j]) * vector[j]
            for i in range(len(vector))
            for j in range(len(vector))
        ),
        start=ZERO,
    )


def _mat_vec(matrix: RationalMatrix, vector: Sequence[Golden]) -> Vector:
    return tuple(
        sum(
            (Golden(matrix[i][j]) * vector[j] for j in range(len(vector))),
            start=ZERO,
        )
        for i in range(len(matrix))
    )


def _unit_vector(index: int, size: int = 8, sign: int = 1) -> tuple[int, ...]:
    return tuple(sign if i == index else 0 for i in range(size))


@lru_cache(maxsize=1)
def e8_roots_from_simple_reflections() -> tuple[tuple[int, ...], ...]:
    """Generate all E8 roots by exact closure under the eight simple reflections."""

    e8 = _e8_gram()
    roots = {_unit_vector(i, sign=sign) for i in range(8) for sign in (-1, 1)}
    frontier = list(roots)
    while frontier:
        root = frontier.pop()
        for i in range(8):
            product = sum((Fraction(root[j]) * e8[j][i] for j in range(8)), Fraction(0))
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
        raise AssertionError("the chosen E8 basis is not a simple-root basis")
    return ordered


def _cone_roots() -> tuple[Vector, ...]:
    # The E8 calculation above proves that these eight walls form a chamber
    # for every height-zero root; it is deliberately executed, not assumed.
    e8_roots_from_simple_reflections()
    return tuple(
        (ZERO, *tuple(Golden(entry) for entry in _unit_vector(i, sign=-1)))
        for i in range(8)
    )


def admissible_control_coordinates(identity_bound: int) -> tuple[Golden, ...]:
    """Enumerate every positive ``x`` possible below an identity bound.

    The conjugate norm equation gives ``x'^2 <= phi``.  Since ``sqrt(5)>2``,
    ``phi<2``, and ``|1-phi|<1``, the coefficient box below follows from
    ``b=(x-x')/sqrt(5)`` and ``a=x'-b(1-phi)``.  Exact filters then decide the
    inequalities, so the deliberately loose integral box cannot lose a value.
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
        and (PHI - x.conjugate() ** 2).sign(0) >= 0
    }
    return tuple(sorted(values))


def _integer_upper_square_root(value: Golden, embedding: int) -> int:
    upper = 0
    while (Golden(upper * upper) - value).sign(embedding) < 0:
        upper += 1
    return upper


@dataclass(frozen=True, slots=True)
class CoordinateAlphabet:
    values: tuple[Golden, ...]
    identity_upper: int
    conjugate_upper: int
    a_bound: int
    b_bound: int


def _coordinate_alphabet(target: Golden) -> CoordinateAlphabet:
    """All possible entries of ``z=E8*y`` for a fixed norm target."""

    twice_target = 2 * target
    identity_upper = _integer_upper_square_root(twice_target, 0)
    conjugate_upper = _integer_upper_square_root(twice_target, 1)
    # |b|=|z-z'|/sqrt(5)<=(U+V)/sqrt(5)<(U+V)/2.
    b_bound = (identity_upper + conjugate_upper) // 2 + 1
    # a=z'-b(1-phi), with |1-phi|<1.
    a_bound = conjugate_upper + b_bound
    values = {
        z
        for a in range(-a_bound, a_bound + 1)
        for b in range(-b_bound, b_bound + 1)
        if (z := Golden(a, b)).sign(0) >= 0
        and (twice_target - z * z).sign(0) >= 0
        and (twice_target - z * z).sign(1) >= 0
    }
    return CoordinateAlphabet(
        tuple(sorted(values)), identity_upper, conjugate_upper, a_bound, b_bound
    )


@dataclass(frozen=True, slots=True)
class HeightLevel:
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
) -> HeightLevel:
    x = control_coordinate
    target = 2 + 2 * PHI * x * x
    alphabet = _coordinate_alphabet(target)
    visited = [0] * 9
    schur_pruned = [0] * 9
    wall_pruned = [0] * 9
    leaf_norm_mismatches = 0
    solutions: list[Vector] = []

    for wall in previous_nonzero_roots:
        if any(entry.sign(0) < 0 for entry in wall[1:]):
            raise AssertionError("dominant-wall pruning requires nonnegative E8 coordinates")

    def extend(prefix: tuple[Golden, ...]) -> None:
        nonlocal leaf_norm_mismatches
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            minimum = _quadratic(prefix, schur[depth - 1])
            remainder = target - minimum
            if remainder.sign(0) < 0 or remainder.sign(1) < 0:
                schur_pruned[depth] += 1
                return
            for wall in previous_nonzero_roots:
                partial_product = sum(
                    (prefix[i] * wall[i + 1] for i in range(depth)),
                    start=ZERO,
                )
                bound = 2 * PHI * x * wall[0]
                if (bound - partial_product).sign(0) < 0:
                    wall_pruned[depth] += 1
                    return
        if depth == 8:
            if _quadratic(prefix, inverse_e8) != target:
                leaf_norm_mismatches += 1
                return
            y = _mat_vec(inverse_e8, prefix)
            solutions.append((x, *y))
            return
        for entry in alphabet.values:
            extend((*prefix, entry))

    extend(())
    gram = bugaenko_basis_gram()
    if any(_pair(root, gram, root) != 2 for root in solutions):
        raise AssertionError("height enumerator produced a vector of the wrong norm")
    if any(
        _pair(left, gram, right).sign(0) > 0
        for i, left in enumerate(solutions)
        for right in solutions[i + 1 :]
    ):
        raise AssertionError("same-height Vinberg roots do not meet acutely")
    return HeightLevel(
        control_coordinate=x,
        height=PHI * x * x,
        target=target,
        alphabet=alphabet,
        visited_by_depth=tuple(visited),
        schur_pruned_by_depth=tuple(schur_pruned),
        wall_pruned_by_depth=tuple(wall_pruned),
        leaf_norm_mismatches=leaf_norm_mismatches,
        roots=tuple(solutions),
    )


@dataclass(frozen=True, slots=True)
class BugaenkoVinbergSearch:
    control_vector: Vector
    identity_bound: int
    cone_roots: tuple[Vector, ...]
    candidate_coordinates: tuple[Golden, ...]
    levels: tuple[HeightLevel, ...]
    roots: tuple[Vector, ...]


@lru_cache(maxsize=None)
def run_bugaenko_vinberg_search(identity_bound: int = 7) -> BugaenkoVinbergSearch:
    """Exhaust every Vinberg height with ``0 < x <= identity_bound``.

    Bound seven is the first integer shell containing the wall after which the
    independently verified chamber is compact.  It is a search budget, not a
    list of expected roots: every admissible ``x`` and every dominant E8 weight
    below it is generated from the norm and acute-angle equations.
    """

    gram = bugaenko_basis_gram()
    control = (ONE, *(ZERO for _ in range(8)))
    if _pair(control, gram, control) != -2 * PHI:
        raise AssertionError("the chosen control vector is not the negative summand generator")
    e8 = _e8_gram()
    inverse_e8 = _fraction_inverse(e8)
    if any(entry.denominator != 1 for row in inverse_e8 for entry in row):
        raise AssertionError("E8 inverse must be integral")
    schur = _prefix_schur_complements(inverse_e8)
    cone = _cone_roots()
    candidates = admissible_control_coordinates(identity_bound)
    levels: list[HeightLevel] = []
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
    return BugaenkoVinbergSearch(
        control_vector=control,
        identity_bound=identity_bound,
        cone_roots=cone,
        candidate_coordinates=candidates,
        levels=tuple(levels),
        roots=tuple(roots),
    )


def _golden(value: Golden) -> dict[str, str]:
    return value.to_json()


def _vector(vector: Iterable[Golden]) -> list[dict[str, str]]:
    return [_golden(entry) for entry in vector]


def build_bugaenko_vinberg_transcript(identity_bound: int = 7) -> dict[str, Any]:
    """Build the deterministic, JSON-serializable wall-completeness transcript."""

    search = run_bugaenko_vinberg_search(identity_bound)
    return {
        "schema_version": "0.1.0",
        "transcript_kind": "exact_vinberg_height_enumeration",
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "lattice": "<-2*phi> + E8",
        "control_vector": _vector(search.control_vector),
        "control_norm": _golden(-2 * PHI),
        "height_formula": "delta=phi*x^2 for norm-two roots r=(x,y)",
        "root_norm_normalization": {
            "smith_nonunit_invariant": "2",
            "candidates_modulo_unit_squares": ["1", "2", "4"],
            "norm_1_rejection": "the lattice is even",
            "norm_4_rejection": (
                "reflection divisibility forces y in 2*E8; the norm equation "
                "modulo 4 and inertness of (2) force x in 2*Z[phi], "
                "contradicting primitivity"
            ),
            "remaining_representative": "2",
        },
        "height_zero": {
            "all_root_count": 240,
            "exhaustion": (
                "writing y=a+b*phi gives E8(a)+E8(b)=2 and "
                "2<a,b>+E8(b)=0; even positive-definiteness forces b=0 "
                "and E8(a)=2"
            ),
            "simple_roots": [_vector(root) for root in search.cone_roots],
        },
        "control_coordinate_scan": {
            "identity_bound": search.identity_bound,
            "conjugate_constraint": "x'^2 <= phi",
            "candidates": [_golden(x) for x in search.candidate_coordinates],
        },
        "levels": [
            {
                "control_coordinate": _golden(level.control_coordinate),
                "height": _golden(level.height),
                "e8_target": _golden(level.target),
                "coordinate_alphabet": {
                    "size": len(level.alphabet.values),
                    "identity_upper": level.alphabet.identity_upper,
                    "conjugate_upper": level.alphabet.conjugate_upper,
                    "a_bound": level.alphabet.a_bound,
                    "b_bound": level.alphabet.b_bound,
                },
                "visited_by_depth": list(level.visited_by_depth),
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
        "termination": {
            "last_processed_control_coordinate": _golden(
                search.candidate_coordinates[-1]
            ),
            "criterion": "independent exact compactness verification of the current chamber",
            "geometry_certificate": "data/known_examples/bugaenko_compact_h8.json",
            "reason": (
                "Vinberg height ordering exhausts every earlier admissible root; "
                "once the current acute chamber is compact, no later root can be a "
                "new simple wall"
            ),
        },
    }
