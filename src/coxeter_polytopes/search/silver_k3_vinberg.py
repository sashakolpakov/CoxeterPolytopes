"""Exact norm-two Vinberg prefix for the silver rank-20 K3 carrier.

The Alexeev-aligned carrier has a free rank-ten presentation over
``O = Z[sqrt(2)]``.  Its first five generators are the saturated one-4
``C5`` block and the remaining generators are ``F4 + A1``.  In the basis

``e0, e1, e2, e3, (sqrt(2)/2)(e1+e2+e4)``

of the saturated block, the vector dual to ``e0`` is integral.  We use its
negative as Vinberg control vector.  If ``x`` is the first coordinate of a
norm-two root, then

``(p,r)=-x`` and ``delta=(1+sqrt(2))*x^2/2``.

This module freezes the original prefix ``0 < x <= 2`` and ten exact
continuations, first through ``x=1+sqrt(2)``, then through
``x=2+sqrt(2)``, through the empty shell ``x=1+2sqrt(2)``, and finally
through ``x=3+sqrt(2)``, ``x=2+2sqrt(2)``, ``x=3+2sqrt(2)``,
``x=4+2sqrt(2)``, ``x=3+3sqrt(2)``, and the empty shell
``x=5+2sqrt(2)``, ``x=4+3sqrt(2)``, ``x=5+3sqrt(2)``, and
``x=4+4sqrt(2)``, and ``x=6+3sqrt(2)``.  It deliberately does not claim
termination: the bounded wall lists do not satisfy a global stopping theorem,
so compactness and finite covolume remain undecidable from these prefixes.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
from math import isqrt
from typing import Sequence, TypeAlias

from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia
from coxeter_polytopes.exact.silver import (
    SILVER_ZERO,
    SQRT2,
    Silver,
    silver_determinant,
    silver_inertia,
    silver_trace_gram,
)
from coxeter_polytopes.exact.trace_lattice import integral_gram
from coxeter_polytopes.geometric.silver_k3_carrier import (
    audit_silver_c5_saturation,
    silver_finite_completion_audits,
)


Vector: TypeAlias = tuple[Silver, ...]
SilverMatrix: TypeAlias = tuple[tuple[Silver, ...], ...]
IntegerMatrix: TypeAlias = tuple[tuple[int, ...], ...]
RationalMatrix: TypeAlias = tuple[tuple[Fraction, ...], ...]


def _silver_inverse(matrix: SilverMatrix) -> SilverMatrix:
    size = len(matrix)
    work = [
        [
            *(Silver.coerce(entry) for entry in matrix[row]),
            *(Silver(1 if row == column else 0) for column in range(size)),
        ]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
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
                work[row][index] - multiplier * work[column][index]
                for index in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


def _fraction_inverse(matrix: RationalMatrix) -> RationalMatrix:
    size = len(matrix)
    work = [
        [
            *matrix[row],
            *(Fraction(row == column) for column in range(size)),
        ]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
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
                work[row][index] - multiplier * work[column][index]
                for index in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


def _fraction_prefix_schur(matrix: IntegerMatrix) -> tuple[RationalMatrix, ...]:
    rational = tuple(tuple(Fraction(entry) for entry in row) for row in matrix)
    size = len(rational)
    output: list[RationalMatrix] = []
    for prefix in range(1, size + 1):
        leading = tuple(
            tuple(rational[i][j] for j in range(prefix)) for i in range(prefix)
        )
        if prefix == size:
            output.append(leading)
            continue
        cross = tuple(
            tuple(rational[i][j] for j in range(prefix, size))
            for i in range(prefix)
        )
        tail = tuple(
            tuple(rational[i][j] for j in range(prefix, size))
            for i in range(prefix, size)
        )
        inverse_tail = _fraction_inverse(tail)
        schur = []
        for i in range(prefix):
            row = []
            for j in range(prefix):
                correction = sum(
                    cross[i][a] * inverse_tail[a][b] * cross[j][b]
                    for a in range(size - prefix)
                    for b in range(size - prefix)
                )
                row.append(leading[i][j] - correction)
            schur.append(tuple(row))
        output.append(tuple(schur))
    return tuple(output)


def _silver_mat_vec(matrix: SilverMatrix, vector: Sequence[Silver]) -> Vector:
    return tuple(
        sum(
            (
                matrix[row][column] * vector[column]
                for column in range(len(vector))
                if matrix[row][column]
            ),
            SILVER_ZERO,
        )
        for row in range(len(matrix))
    )


def _silver_dot(left: Sequence[Silver], right: Sequence[Silver]) -> Silver:
    if len(left) != len(right):
        raise ValueError("silver dot-product dimensions do not match")
    return sum(
        (left[index] * right[index] for index in range(len(left))),
        SILVER_ZERO,
    )


def _silver_pair(
    left: Sequence[Silver], matrix: SilverMatrix, right: Sequence[Silver]
) -> Silver:
    return sum(
        (
            left[row] * matrix[row][column] * right[column]
            for row in range(len(left))
            for column in range(len(right))
            if matrix[row][column]
        ),
        SILVER_ZERO,
    )


def _silver_quadratic(vector: Sequence[Silver], matrix: SilverMatrix) -> Silver:
    # Every form passed here is symmetric.  Evaluating only the diagonal and
    # strict upper triangle roughly halves the exact Fraction arithmetic in
    # the shell enumerators, which matters once the finite alphabet reaches
    # hundreds of thousands of component tuples.
    size = len(vector)
    return sum(
        (vector[index] * matrix[index][index] * vector[index] for index in range(size)),
        SILVER_ZERO,
    ) + 2 * sum(
        (
            vector[row] * matrix[row][column] * vector[column]
            for row in range(size)
            for column in range(row + 1, size)
            if matrix[row][column]
        ),
        SILVER_ZERO,
    )


def _fraction_quadratic(
    vector: Sequence[int], matrix: RationalMatrix
) -> Fraction:
    return sum(
        Fraction(vector[row]) * matrix[row][column] * Fraction(vector[column])
        for row in range(len(vector))
        for column in range(len(vector))
    )


def _transpose(matrix: SilverMatrix) -> SilverMatrix:
    return tuple(
        tuple(matrix[row][column] for row in range(len(matrix)))
        for column in range(len(matrix[0]))
    )


def _silver_congruence(change: SilverMatrix, gram: SilverMatrix) -> SilverMatrix:
    return tuple(
        tuple(
            sum(
                (
                    change[row][left]
                    * gram[row][column]
                    * change[column][right]
                    for row in range(len(gram))
                    for column in range(len(gram))
                ),
                SILVER_ZERO,
            )
            for right in range(len(change[0]))
        )
        for left in range(len(change[0]))
    )


def _block_diagonal(*blocks: SilverMatrix) -> SilverMatrix:
    size = sum(len(block) for block in blocks)
    matrix = [[SILVER_ZERO for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        for row in range(len(block)):
            for column in range(len(block)):
                matrix[offset + row][offset + column] = block[row][column]
        offset += len(block)
    return tuple(tuple(row) for row in matrix)


@lru_cache(maxsize=1)
def silver_k3_alexeev_field_gram() -> SilverMatrix:
    """Return the saturated ``C5 + (F4+A1)`` field Gram matrix."""

    seed = audit_silver_c5_saturation()
    old = seed.field_gram
    columns: list[Vector] = []
    for index in range(4):
        column = [SILVER_ZERO for _ in range(5)]
        column[index] = Silver(1)
        columns.append(tuple(column))
    columns.append(
        (
            SILVER_ZERO,
            SQRT2 / 2,
            SQRT2 / 2,
            SILVER_ZERO,
            SQRT2 / 2,
        )
    )
    change = tuple(
        tuple(columns[column][row] for column in range(5)) for row in range(5)
    )
    saturated_c5 = _silver_congruence(change, old)
    completion = next(
        item for item in silver_finite_completion_audits() if item.name == "F4 + A1"
    )
    return _block_diagonal(saturated_c5, completion.field_gram)


def silver_k3_control_vector() -> Vector:
    """The integral timelike vector ``-e0^vee``."""

    return (
        SQRT2 - 1,
        SQRT2 - 1,
        2 - SQRT2,
        Silver(1),
        Silver(1),
        *(SILVER_ZERO for _ in range(5)),
    )


def _height_zero_simple_matrix() -> SilverMatrix:
    """Columns are positive finite simple roots in coordinates 1,...,9."""

    matrix = [
        [Silver(1 if row == column else 0) for column in range(9)]
        for row in range(9)
    ]
    # The fourth B4 root is old e4 = sqrt(2)*f4-e1-e2.
    matrix[0][3] = Silver(-1)
    matrix[1][3] = Silver(-1)
    matrix[2][3] = SILVER_ZERO
    matrix[3][3] = SQRT2
    return tuple(tuple(row) for row in matrix)


def silver_k3_height_zero_walls() -> tuple[Vector, ...]:
    """Return the nine oriented walls of type ``B4 + F4 + A1``."""

    simple = _height_zero_simple_matrix()
    return tuple(
        (
            SILVER_ZERO,
            *tuple(-simple[row][column] for row in range(9)),
        )
        for column in range(9)
    )


@dataclass(frozen=True, slots=True)
class IntegerShortVectorAudit:
    rank: int
    determinant: int
    target: int
    coordinate_bounds: tuple[int, ...]
    visited_by_depth: tuple[int, ...]
    schur_pruned_by_depth: tuple[int, ...]
    vectors: tuple[tuple[int, ...], ...]


def _floor_sqrt_fraction(value: Fraction) -> int:
    if value < 0:
        raise ValueError("square-root bound must be nonnegative")
    candidate = isqrt(value.numerator // value.denominator)
    while Fraction((candidate + 1) ** 2) <= value:
        candidate += 1
    while Fraction(candidate**2) > value:
        candidate -= 1
    return candidate


def _integer_short_vectors(
    gram: IntegerMatrix, target: int = 2
) -> IntegerShortVectorAudit:
    """Exhaust integral vectors of a positive rational quadratic form."""

    if inertia(gram) != Inertia(len(gram), 0, 0):
        raise ValueError("short-vector enumeration requires a positive form")
    inverse = _fraction_inverse(
        tuple(tuple(Fraction(entry) for entry in row) for row in gram)
    )
    bounds = tuple(
        _floor_sqrt_fraction(Fraction(target) * inverse[index][index])
        for index in range(len(gram))
    )
    schur = _fraction_prefix_schur(gram)
    visited = [0 for _ in range(len(gram) + 1)]
    pruned = [0 for _ in range(len(gram) + 1)]
    vectors: list[tuple[int, ...]] = []

    def extend(prefix: tuple[int, ...]) -> None:
        depth = len(prefix)
        visited[depth] += 1
        if depth and _fraction_quadratic(prefix, schur[depth - 1]) > target:
            pruned[depth] += 1
            return
        if depth == len(gram):
            if _fraction_quadratic(prefix, schur[-1]) == target:
                vectors.append(prefix)
            return
        for entry in range(-bounds[depth], bounds[depth] + 1):
            extend((*prefix, entry))

    extend(())
    determinant_value = determinant(gram)
    if determinant_value.b or determinant_value.a.denominator != 1:
        raise AssertionError("an integral Gram matrix had nonintegral determinant")
    return IntegerShortVectorAudit(
        rank=len(gram),
        determinant=determinant_value.a.numerator,
        target=target,
        coordinate_bounds=bounds,
        visited_by_depth=tuple(visited),
        schur_pruned_by_depth=tuple(pruned),
        vectors=tuple(vectors),
    )


def _field_norm_two_vectors(
    block: SilverMatrix,
) -> tuple[IntegerShortVectorAudit, tuple[Vector, ...]]:
    trace = integral_gram(silver_trace_gram(block))
    positive_trace = tuple(tuple(-entry for entry in row) for row in trace)
    if any(positive_trace[index][index] % 2 for index in range(len(trace))):
        raise AssertionError("a height-zero trace block is not even")
    enumeration = _integer_short_vectors(positive_trace)
    field_vectors = []
    for vector in enumeration.vectors:
        field = tuple(
            Silver(vector[2 * index], vector[2 * index + 1])
            for index in range(len(block))
        )
        if _silver_quadratic(field, block) == 2:
            field_vectors.append(field)
    return enumeration, tuple(field_vectors)


def _reflect(vector: Vector, root: Vector, gram: SilverMatrix) -> Vector:
    # Every root used here has square two.
    coefficient = _silver_pair(vector, gram, root)
    return tuple(vector[index] - coefficient * root[index] for index in range(len(vector)))


def _root_closure(simple_roots: Sequence[Vector], gram: SilverMatrix) -> tuple[Vector, ...]:
    roots = set(simple_roots)
    roots.update(tuple(-entry for entry in root) for root in simple_roots)
    frontier = list(roots)
    while frontier:
        vector = frontier.pop()
        for root in simple_roots:
            reflected = _reflect(vector, root, gram)
            if reflected not in roots:
                roots.add(reflected)
                frontier.append(reflected)
    return tuple(sorted(roots, key=lambda vector: tuple((entry.a, entry.b) for entry in vector)))


@dataclass(frozen=True, slots=True)
class SilverHeightZeroAudit:
    finite_type: str
    trace_short_vector_counts: tuple[int, int, int]
    field_norm_two_counts: tuple[int, int, int]
    total_field_norm_two_roots: int
    reflection_closure_count: int
    reflection_closure_matches_exhaustion: bool
    trace_determinants: tuple[int, int, int]

    @property
    def passes(self) -> bool:
        return (
            self.finite_type == "B4 + F4 + A1"
            and self.trace_short_vector_counts == (240, 240, 4)
            and self.field_norm_two_counts == (32, 48, 2)
            and self.total_field_norm_two_roots == 82
            and self.reflection_closure_count == 82
            and self.reflection_closure_matches_exhaustion
            and self.trace_determinants == (1, 1, 4)
        )


@lru_cache(maxsize=1)
def audit_silver_k3_height_zero() -> SilverHeightZeroAudit:
    gram = silver_k3_alexeev_field_gram()
    blocks = (
        tuple(tuple(gram[i][j] for j in range(1, 5)) for i in range(1, 5)),
        tuple(tuple(gram[i][j] for j in range(5, 9)) for i in range(5, 9)),
        ((gram[9][9],),),
    )
    exhaustive: set[Vector] = set()
    trace_counts: list[int] = []
    field_counts: list[int] = []
    trace_determinants: list[int] = []
    offsets = (1, 5, 9)
    for block, offset in zip(blocks, offsets, strict=True):
        trace_audit, field_roots = _field_norm_two_vectors(block)
        trace_counts.append(len(trace_audit.vectors))
        field_counts.append(len(field_roots))
        trace_determinants.append(trace_audit.determinant)
        for local in field_roots:
            vector = [SILVER_ZERO for _ in range(10)]
            vector[offset : offset + len(local)] = local
            exhaustive.add(tuple(vector))

    positive_simple = tuple(
        tuple(-entry for entry in wall) for wall in silver_k3_height_zero_walls()
    )
    closure = set(_root_closure(positive_simple, gram))
    audit = SilverHeightZeroAudit(
        finite_type="B4 + F4 + A1",
        trace_short_vector_counts=tuple(trace_counts),
        field_norm_two_counts=tuple(field_counts),
        total_field_norm_two_roots=len(exhaustive),
        reflection_closure_count=len(closure),
        reflection_closure_matches_exhaustion=closure == exhaustive,
        trace_determinants=tuple(trace_determinants),
    )
    if not audit.passes:
        raise AssertionError("the silver height-zero audit failed")
    return audit


@dataclass(frozen=True, slots=True)
class SilverCoordinateAlphabet:
    values: tuple[Silver, ...]
    identity_upper: int
    conjugate_upper: int
    coefficient_bound: int


def _integer_upper_square_root(value: Silver, embedding: int) -> int:
    upper = 0
    while (Silver(upper * upper) - value).sign(embedding) < 0:
        upper += 1
    return upper


def _weight_alphabet(target: Silver) -> SilverCoordinateAlphabet:
    twice_target = 2 * target
    identity_upper = _integer_upper_square_root(twice_target, 0)
    conjugate_upper = _integer_upper_square_root(twice_target, 1)
    coefficient_bound = identity_upper + conjugate_upper
    values = {
        value
        for a in range(-coefficient_bound, coefficient_bound + 1)
        for b in range(-coefficient_bound, coefficient_bound + 1)
        if (value := Silver(a, b)).sign(0) >= 0
        and (twice_target - value * value).sign(0) >= 0
        and (twice_target - value * value).sign(1) >= 0
    }
    return SilverCoordinateAlphabet(
        tuple(sorted(values)),
        identity_upper,
        conjugate_upper,
        coefficient_bound,
    )


def admissible_silver_control_coordinates(identity_bound: int) -> tuple[Silver, ...]:
    """All ``x`` with ``0<x<=bound`` allowed by the conjugate norm equation."""

    if identity_bound < 1:
        return ()
    # The conjugate equation gives x'^2 <= 2(1+sqrt(2)) < 9.  Hence
    # |a|,|b| < (bound+3)/2; the deliberately looser box is exhaustive.
    coefficient_bound = identity_bound + 3
    values = {
        value
        for a in range(-coefficient_bound, coefficient_bound + 1)
        for b in range(-coefficient_bound, coefficient_bound + 1)
        if (value := Silver(a, b)).sign(0) > 0
        and (Silver(identity_bound) - value).sign(0) >= 0
        and (2 * (1 + SQRT2) - value.conjugate() ** 2).sign(0) >= 0
    }
    return tuple(sorted(values))


@dataclass(frozen=True, slots=True)
class SilverK3HeightLevel:
    control_coordinate: Silver
    height: Silver
    target: Silver
    alphabet: SilverCoordinateAlphabet
    visited_by_depth: tuple[int, ...]
    schur_pruned_by_depth: tuple[int, ...]
    exact_weight_solutions: int
    nonintegral_reconstructions: int
    previous_wall_rejections: int
    previous_wall_rejection_witnesses: tuple[
        tuple[Vector, tuple[int, ...]], ...
    ]
    roots: tuple[Vector, ...]


def _enumeration_data() -> tuple[
    SilverMatrix,
    SilverMatrix,
    Vector,
    SilverMatrix,
    SilverMatrix,
    tuple[SilverMatrix, ...],
]:
    gram = silver_k3_alexeev_field_gram()
    finite = tuple(tuple(gram[i][j] for j in range(1, 10)) for i in range(1, 10))
    cross = tuple(gram[0][j] for j in range(1, 10))
    simple = _height_zero_simple_matrix()
    finite_simple_gram = _silver_congruence(simple, finite)
    prefixes = tuple(
        _silver_inverse(
            tuple(tuple(finite_simple_gram[i][j] for j in range(size)) for i in range(size))
        )
        for size in range(1, 10)
    )
    return (
        gram,
        _silver_inverse(finite),
        cross,
        _transpose(_silver_inverse(simple)),
        _silver_inverse(finite_simple_gram),
        prefixes,
    )


def _enumerate_height_level(
    control_coordinate: Silver, previous_roots: Sequence[Vector]
) -> SilverK3HeightLevel:
    gram, finite_inverse, cross, inverse_simple_transpose, weight_form, prefixes = (
        _enumeration_data()
    )
    x = control_coordinate
    target = 2 + (1 + SQRT2) * x * x
    alphabet = _weight_alphabet(target)
    visited = [0 for _ in range(10)]
    pruned = [0 for _ in range(10)]
    exact_weight_solutions = 0
    nonintegral_reconstructions = 0
    previous_wall_rejections = 0
    previous_wall_rejection_witnesses: list[
        tuple[Vector, tuple[int, ...]]
    ] = []
    roots: list[Vector] = []

    def extend(prefix: tuple[Silver, ...]) -> None:
        nonlocal exact_weight_solutions
        nonlocal nonintegral_reconstructions
        nonlocal previous_wall_rejections
        depth = len(prefix)
        visited[depth] += 1
        if depth:
            remainder = target - _silver_quadratic(prefix, prefixes[depth - 1])
            if remainder.sign(0) < 0 or remainder.sign(1) < 0:
                pruned[depth] += 1
                return
        if depth == 9:
            if _silver_quadratic(prefix, weight_form) != target:
                return
            exact_weight_solutions += 1
            weight = _silver_mat_vec(inverse_simple_transpose, prefix)
            finite_coordinates = _silver_mat_vec(
                finite_inverse,
                tuple(weight[index] - cross[index] * x for index in range(9)),
            )
            if any(not entry.is_integral for entry in finite_coordinates):
                nonintegral_reconstructions += 1
                return
            root = (x, *finite_coordinates)
            if _silver_quadratic(root, gram) != 2:
                raise AssertionError("a reconstructed silver root has the wrong norm")
            # Symmetry gives B(root, wall)=wall^T(G root).  Computing G root
            # once makes the complete blocker census exact but substantially
            # cheaper than repeating the sparse matrix product for each wall.
            gram_root = _silver_mat_vec(gram, root)
            blockers = tuple(
                index
                for index, wall in enumerate(previous_roots)
                if _silver_dot(wall, gram_root).sign(0) > 0
            )
            if blockers:
                previous_wall_rejections += 1
                previous_wall_rejection_witnesses.append((root, blockers))
                return
            roots.append(root)
            return
        for entry in alphabet.values:
            extend((*prefix, entry))

    extend(())
    if any(
        _silver_pair(left, gram, right).sign(0) > 0
        for index, left in enumerate(roots)
        for right in roots[index + 1 :]
    ):
        raise AssertionError("same-height silver Vinberg candidates are not acute")
    return SilverK3HeightLevel(
        control_coordinate=x,
        height=(1 + SQRT2) * x * x / 2,
        target=target,
        alphabet=alphabet,
        visited_by_depth=tuple(visited),
        schur_pruned_by_depth=tuple(pruned),
        exact_weight_solutions=exact_weight_solutions,
        nonintegral_reconstructions=nonintegral_reconstructions,
        previous_wall_rejections=previous_wall_rejections,
        previous_wall_rejection_witnesses=tuple(
            previous_wall_rejection_witnesses
        ),
        roots=tuple(roots),
    )


@dataclass(frozen=True, slots=True)
class SilverPrefixDiagramAudit:
    elliptic_rank_eight_count: int
    elliptic_rank_nine_count: int
    finite_endpoint_histogram: tuple[tuple[int, int], ...]
    open_edge_count: int
    connected_affine_subdiagrams: tuple[tuple[int, ...], ...]
    compact_stopping_available: bool
    finite_volume_stopping_available: bool


def _connected(matrix: SilverMatrix, vertices: tuple[int, ...]) -> bool:
    unseen = set(vertices)
    start = min(unseen)
    unseen.remove(start)
    frontier = [start]
    while frontier:
        vertex = frontier.pop()
        new = {
            other for other in unseen if matrix[vertex][other] != SILVER_ZERO
        }
        unseen.difference_update(new)
        frontier.extend(new)
    return not unseen


def _subgram(matrix: SilverMatrix, vertices: tuple[int, ...]) -> SilverMatrix:
    return tuple(tuple(matrix[i][j] for j in vertices) for i in vertices)


def _is_elliptic(matrix: SilverMatrix, vertices: tuple[int, ...]) -> bool:
    return silver_inertia(_subgram(matrix, vertices), embedding=0) == Inertia(
        len(vertices), 0, 0
    )


def _spanning_wall_gram_inertias(
    wall_gram: SilverMatrix,
) -> tuple[Inertia, Inertia]:
    """Recover a large wall Gram's inertia from a spanning principal block.

    The first nine walls are the finite height-zero basis and wall 9 has
    nonzero control coordinate.  Thus the first ten wall normals span the
    rank-ten carrier as soon as their principal Gram is nonsingular.  Since
    every wall Gram factors through that carrier, its remaining dimensions
    are then exactly null.  This avoids cubic elimination on increasingly
    large, highly singular matrices while retaining an exact rank witness.
    """

    if len(wall_gram) < 10:
        raise ValueError("the silver wall Gram needs ten spanning walls")
    basis = _subgram(wall_gram, tuple(range(10)))
    identity = silver_inertia(basis, embedding=0)
    conjugate = silver_inertia(basis, embedding=1)
    if identity != Inertia(9, 1, 0) or conjugate != Inertia(10, 0, 0):
        raise AssertionError("the first ten silver walls no longer span the carrier")
    nullity = len(wall_gram) - 10
    return Inertia(9, 1, nullity), Inertia(10, 0, nullity)


def _is_connected_affine(matrix: SilverMatrix, vertices: tuple[int, ...]) -> bool:
    size = len(vertices)
    return (
        size >= 2
        and _connected(matrix, vertices)
        and silver_inertia(_subgram(matrix, vertices), embedding=0)
        == Inertia(size - 1, 0, 1)
        and all(
            _is_elliptic(matrix, tuple(node for node in vertices if node != omitted))
            for omitted in vertices
        )
    )


def _pairwise_psd_possible(
    matrix: SilverMatrix, vertices: tuple[int, ...], *, positive_definite: bool
) -> bool:
    """Apply the exact principal-two-minor prerequisite for PSD/PD subsets.

    Every wall in this module has square two and the wall Gram is acute.
    Hence a positive-semidefinite principal submatrix must have
    ``g_ij >= -2`` for every pair, while a positive-definite one must have
    ``g_ij > -2``.  This cheap necessary condition removes the overwhelmingly
    many indefinite subsets before exact inertia is computed; it does not
    replace the inertia test.
    """

    return all(
        (matrix[left][right] + 2).sign(0) > (0 if positive_definite else -1)
        for left, right in combinations(vertices, 2)
    )


def _prefix_diagram_audit(wall_gram: SilverMatrix) -> SilverPrefixDiagramAudit:
    facets = tuple(range(len(wall_gram)))
    edges = tuple(
        vertices
        for vertices in combinations(facets, 8)
        if _pairwise_psd_possible(
            wall_gram, vertices, positive_definite=True
        )
        and _is_elliptic(wall_gram, vertices)
    )
    vertices = tuple(
        nodes
        for nodes in combinations(facets, 9)
        if _pairwise_psd_possible(wall_gram, nodes, positive_definite=True)
        and _is_elliptic(wall_gram, nodes)
    )
    histogram: dict[int, int] = {}
    for edge in edges:
        edge_set = set(edge)
        endpoint_count = sum(edge_set.issubset(vertex) for vertex in vertices)
        histogram[endpoint_count] = histogram.get(endpoint_count, 0) + 1
    affine = tuple(
        nodes
        for size in range(2, len(facets) + 1)
        for nodes in combinations(facets, size)
        if _pairwise_psd_possible(wall_gram, nodes, positive_definite=False)
        and _is_connected_affine(wall_gram, nodes)
    )
    open_edges = sum(count for endpoints, count in histogram.items() if endpoints != 2)
    return SilverPrefixDiagramAudit(
        elliptic_rank_eight_count=len(edges),
        elliptic_rank_nine_count=len(vertices),
        finite_endpoint_histogram=tuple(sorted(histogram.items())),
        open_edge_count=open_edges,
        connected_affine_subdiagrams=affine,
        compact_stopping_available=not open_edges and not affine,
        # There are no affine ideal endpoints in this prefix, so the same
        # open edges also prevent a finite-volume stopping argument.
        finite_volume_stopping_available=not open_edges and not affine,
    )


@dataclass(frozen=True, slots=True)
class SilverSaturationWallAudit:
    vector: Vector
    field_norm: Silver
    trace_square: int
    reflection_coefficient_on_e0: Silver
    reflection_is_integral: bool

    @property
    def passes(self) -> bool:
        return (
            self.field_norm == Silver(3, -1)
            and self.trace_square == -2
            and self.reflection_coefficient_on_e0
            == Silver(Fraction(-4, 7), Fraction(-6, 7))
            and not self.reflection_is_integral
        )


def audit_silver_saturation_wall() -> SilverSaturationWallAudit:
    gram = silver_k3_alexeev_field_gram()
    vector = tuple(
        Silver(1 if index == 4 else 0) for index in range(10)
    )
    norm = _silver_quadratic(vector, gram)
    e0 = tuple(Silver(1 if index == 0 else 0) for index in range(10))
    coefficient = 2 * _silver_pair(e0, gram, vector) / norm
    audit = SilverSaturationWallAudit(
        vector=vector,
        field_norm=norm,
        trace_square=int(silver_trace_gram(((norm,),))[0][0]),
        reflection_coefficient_on_e0=coefficient,
        reflection_is_integral=coefficient.is_integral,
    )
    if not audit.passes:
        raise AssertionError("the saturation-wall audit failed")
    return audit


@dataclass(frozen=True, slots=True)
class SilverK3VinbergPrefix:
    identity_bound: int
    field_gram: SilverMatrix
    field_determinant: Silver
    field_identity_inertia: Inertia
    field_conjugate_inertia: Inertia
    control_vector: Vector
    control_norm: Silver
    control_covector: Vector
    height_zero: SilverHeightZeroAudit
    cone_roots: tuple[Vector, ...]
    candidate_coordinates: tuple[Silver, ...]
    levels: tuple[SilverK3HeightLevel, ...]
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    diagram: SilverPrefixDiagramAudit
    saturation_wall: SilverSaturationWallAudit
    stopping_status: str

    @property
    def positive_height_roots(self) -> tuple[Vector, ...]:
        return tuple(root for level in self.levels for root in level.roots)

    @property
    def passes(self) -> bool:
        return (
            self.identity_bound == 2
            and self.field_identity_inertia == Inertia(9, 1, 0)
            and self.field_conjugate_inertia == Inertia(10, 0, 0)
            and self.control_norm == 1 - SQRT2
            and self.control_covector == (-Silver(1), *(SILVER_ZERO for _ in range(9)))
            and self.height_zero.passes
            and self.candidate_coordinates == (Silver(1), SQRT2, Silver(2))
            and tuple(len(level.roots) for level in self.levels) == (1, 1, 1)
            and tuple(level.exact_weight_solutions for level in self.levels) == (1, 5, 7)
            and tuple(level.nonintegral_reconstructions for level in self.levels) == (0, 3, 2)
            and tuple(level.previous_wall_rejections for level in self.levels) == (0, 1, 4)
            and len(self.roots) == 12
            and self.wall_identity_inertia == Inertia(9, 1, 2)
            and self.wall_conjugate_inertia == Inertia(10, 0, 2)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.diagram.elliptic_rank_eight_count == 69
            and self.diagram.elliptic_rank_nine_count == 10
            and self.diagram.finite_endpoint_histogram == ((1, 48), (2, 21))
            and self.diagram.open_edge_count == 48
            and self.diagram.connected_affine_subdiagrams == ()
            and not self.diagram.compact_stopping_available
            and not self.diagram.finite_volume_stopping_available
            and self.saturation_wall.passes
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergExtension:
    """The first exact shell beyond the frozen ``x<=2`` prefix."""

    identity_bound: int
    completed_through_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    diagram: SilverPrefixDiagramAudit
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.identity_bound == 3
            and self.completed_through_control_coordinate == 1 + SQRT2
            and self.candidate_coordinates
            == (Silver(1), SQRT2, Silver(2), 1 + SQRT2)
            and self.prior_wall_count == 12
            and self.new_level.control_coordinate == 1 + SQRT2
            and len(self.new_level.alphabet.values) == 8
            and self.new_level.visited_by_depth
            == (1, 8, 64, 104, 120, 112, 456, 904, 512, 184)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 51, 89, 106, 55, 343, 840, 489, 103)
            and self.new_level.exact_weight_solutions == 6
            and self.new_level.nonintegral_reconstructions == 0
            and self.new_level.previous_wall_rejections == 2
            and len(self.new_level.roots) == 4
            and len(self.roots) == 16
            and self.wall_identity_inertia == Inertia(9, 1, 6)
            and self.wall_conjugate_inertia == Inertia(10, 0, 6)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.diagram.elliptic_rank_eight_count == 267
            and self.diagram.elliptic_rank_nine_count == 41
            and self.diagram.finite_endpoint_histogram
            == ((0, 1), (1, 163), (2, 103))
            and self.diagram.open_edge_count == 164
            and self.diagram.connected_affine_subdiagrams == ()
            and not self.diagram.compact_stopping_available
            and not self.diagram.finite_volume_stopping_available
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverEdgeEndpointAudit:
    """Finite pairing census and first endpoint for one open rank-eight edge."""

    edge_indices: tuple[int, ...]
    finite_type: str
    edge_gram: SilverMatrix
    edge_determinant: Silver
    pairing_alphabet: tuple[Silver, ...]
    allowable_pairing_patterns: tuple[Vector, ...]
    new_endpoint_root_indices: tuple[int, ...]
    endpoint_pairing: Vector
    endpoint_gram_determinant: Silver
    endpoint_identity_inertia: Inertia
    endpoint_conjugate_inertia: Inertia

    @property
    def passes(self) -> bool:
        return (
            self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_type == "A2 + F4 + B2"
            and self.edge_determinant == Silver(6)
            and self.pairing_alphabet == (SILVER_ZERO, Silver(-1), -SQRT2)
            and len(self.allowable_pairing_patterns) == 11
            and self.new_endpoint_root_indices == (33,)
            and self.endpoint_pairing
            == (-SQRT2, *(SILVER_ZERO for _ in range(7)))
            and self.endpoint_gram_determinant == Silver(4)
            and self.endpoint_identity_inertia == Inertia(9, 0, 0)
            and self.endpoint_conjugate_inertia == Inertia(9, 0, 0)
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergEdgeExtension:
    """The ``x=2+sqrt(2)`` shell and its targeted edge endpoint."""

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    edge_endpoint: SilverEdgeEndpointAudit
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 2 + SQRT2
            and self.next_admissible_control_coordinate == 1 + 2 * SQRT2
            and self.candidate_coordinates
            == (Silver(1), SQRT2, Silver(2), 1 + SQRT2, 2 + SQRT2)
            and self.prior_wall_count == 16
            and self.new_level.control_coordinate == 2 + SQRT2
            and len(self.new_level.alphabet.values) == 11
            and self.new_level.visited_by_depth
            == (1, 11, 121, 264, 374, 396, 1936, 4532, 3047, 1199)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 97, 230, 338, 220, 1524, 4255, 2938, 798)
            and self.new_level.exact_weight_solutions == 72
            and self.new_level.nonintegral_reconstructions == 43
            and self.new_level.previous_wall_rejections == 11
            and len(self.new_level.roots) == 18
            and len(self.roots) == 34
            and self.wall_identity_inertia == Inertia(9, 1, 24)
            and self.wall_conjugate_inertia == Inertia(10, 0, 24)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.edge_endpoint.passes
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergOnePlusTwoSqrt2Extension:
    """The exact empty simple-root shell ``x=1+2sqrt(2)``.

    The targeted rank-eight edge retains its one known finite endpoint.  An
    ideal-endpoint audit is intentionally deferred until the simple-root list
    is complete, because Vinberg's criterion requires the full parabolic
    zero-neighbourhood rather than a bounded prefix fragment.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 1 + 2 * SQRT2
            and self.next_admissible_control_coordinate == 3 + SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
            )
            and self.prior_wall_count == 34
            and self.new_level.control_coordinate == 1 + 2 * SQRT2
            and self.new_level.target == 19 + 13 * SQRT2
            and len(self.new_level.alphabet.values) == 7
            and self.new_level.visited_by_depth
            == (1, 7, 49, 77, 84, 70, 266, 469, 273, 112)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 38, 65, 74, 32, 199, 430, 257, 58)
            and self.new_level.exact_weight_solutions == 4
            and self.new_level.nonintegral_reconstructions == 0
            and self.new_level.previous_wall_rejections == 4
            and len(self.new_level.previous_wall_rejection_witnesses) == 4
            and all(
                blockers
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            and self.new_level.roots == ()
            and len(self.roots) == 34
            and self.wall_identity_inertia == Inertia(9, 1, 24)
            and self.wall_conjugate_inertia == Inertia(10, 0, 24)
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33,)
            and self.new_finite_endpoint_root_indices == ()
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status == "not_found_through_this_shell"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergThreePlusSqrt2Extension:
    """The exact three-wall shell ``x=3+sqrt(2)``.

    "Simple" here is relative to the complete lower-height prefix.  The
    three same-height candidates have pairwise nonpositive inner products,
    so none cuts off another.  This bounded computation does not assert that
    the resulting 37-wall list is the complete chamber.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        expected_new_roots = (
            (
                3 + SQRT2,
                4 + 2 * SQRT2,
                4 + 4 * SQRT2,
                8 + 6 * SQRT2,
                8 + 6 * SQRT2,
                2 + SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
                2 + 2 * SQRT2,
                1 + SQRT2,
            ),
            (
                3 + SQRT2,
                4 + 2 * SQRT2,
                4 + 4 * SQRT2,
                8 + 6 * SQRT2,
                8 + 6 * SQRT2,
                3 + 2 * SQRT2,
                6 + 4 * SQRT2,
                6 + 4 * SQRT2,
                3 + 2 * SQRT2,
                SILVER_ZERO,
            ),
            (
                3 + SQRT2,
                4 + 2 * SQRT2,
                4 + 4 * SQRT2,
                8 + 6 * SQRT2,
                8 + 6 * SQRT2,
                2 + 2 * SQRT2,
                3 + 3 * SQRT2,
                4 + 2 * SQRT2,
                2 + SQRT2,
                1 + SQRT2,
            ),
        )
        expected_pairings = (
            (
                SILVER_ZERO,
                SILVER_ZERO,
                SILVER_ZERO,
                SILVER_ZERO,
                SILVER_ZERO,
                -1 - SQRT2,
                SILVER_ZERO,
                -2 - SQRT2,
            ),
            (
                SILVER_ZERO,
                SILVER_ZERO,
                SILVER_ZERO,
                Silver(-1),
                Silver(-1),
                SILVER_ZERO,
                SILVER_ZERO,
                -8 - 5 * SQRT2,
            ),
            (
                SILVER_ZERO,
                SILVER_ZERO,
                -1 - SQRT2,
                SILVER_ZERO,
                SILVER_ZERO,
                SILVER_ZERO,
                SILVER_ZERO,
                -2 - SQRT2,
            ),
        )
        return (
            self.completed_through_control_coordinate == 3 + SQRT2
            and self.next_admissible_control_coordinate == 2 + 2 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
            )
            and self.prior_wall_count == 34
            and self.new_level.control_coordinate == 3 + SQRT2
            and self.new_level.target == 25 + 17 * SQRT2
            and len(self.new_level.alphabet.values) == 10
            and self.new_level.visited_by_depth
            == (1, 10, 100, 200, 250, 250, 1170, 2540, 1550, 570)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 80, 175, 225, 133, 916, 2385, 1493, 345)
            and self.new_level.exact_weight_solutions == 12
            and self.new_level.nonintegral_reconstructions == 0
            and self.new_level.previous_wall_rejections == 9
            and tuple(
                blockers
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            == (
                (10,),
                (9,),
                (9, 30),
                (9, 31),
                (9, 32),
                (12,),
                (13,),
                (14,),
                (9,),
            )
            and self.new_level.roots == expected_new_roots
            and len(self.roots) == 37
            and self.wall_identity_inertia == Inertia(9, 1, 27)
            and self.wall_conjugate_inertia == Inertia(10, 0, 27)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33,)
            and self.new_finite_endpoint_root_indices == ()
            and self.new_root_edge_pairings == expected_pairings
            and self.new_root_edge_extension_determinants
            == (-60 - 48 * SQRT2, -744 - 528 * SQRT2, -60 - 48 * SQRT2)
            and self.new_root_edge_extension_identity_inertias
            == (Inertia(8, 1, 0),) * 3
            and self.new_root_edge_extension_conjugate_inertias
            == (Inertia(9, 0, 0),) * 3
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status == "not_found_through_this_shell"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergTwoPlusTwoSqrt2Extension:
    """The exact 56-wall shell ``x=2+2sqrt(2)``.

    The 56 accepted same-height roots have pairwise nonpositive identity
    pairings, so they enter simultaneously and independently of candidate
    order.  The resulting 93-wall list remains only a bounded norm-two
    prefix; this object records no chamber-completeness or volume claim.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 2 + 2 * SQRT2
            and self.next_admissible_control_coordinate == 3 + 2 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
            )
            and self.prior_wall_count == 37
            and self.new_level.control_coordinate == 2 + 2 * SQRT2
            and self.new_level.target == 30 + 20 * SQRT2
            and self.new_level.alphabet.values
            == (
                SILVER_ZERO,
                Silver(1),
                SQRT2,
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
                4 + 3 * SQRT2,
                5 + 3 * SQRT2,
                4 + 4 * SQRT2,
                6 + 3 * SQRT2,
                5 + 4 * SQRT2,
            )
            and self.new_level.visited_by_depth
            == (1, 16, 256, 656, 1072, 1392, 7984, 22016, 15840, 6368)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 215, 589, 985, 893, 6608, 21026, 15442, 4706)
            and self.new_level.exact_weight_solutions == 277
            and self.new_level.nonintegral_reconstructions == 129
            and self.new_level.previous_wall_rejections == 92
            and len(self.new_level.previous_wall_rejection_witnesses) == 92
            and all(
                blockers
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            and len(self.new_level.roots) == 56
            and len(self.roots) == 93
            and self.wall_identity_inertia == Inertia(9, 1, 83)
            and self.wall_conjugate_inertia == Inertia(10, 0, 83)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 1482
            and self.same_height_zero_pair_count == 58
            and self.same_height_negative_pair_count
            + self.same_height_zero_pair_count
            == 56 * 55 // 2
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33,)
            and self.new_finite_endpoint_root_indices == ()
            and len(self.new_root_edge_pairings) == 56
            and len(self.new_root_edge_extension_determinants) == 56
            and all(
                determinant.sign(0) < 0
                for determinant in self.new_root_edge_extension_determinants
            )
            and self.new_root_edge_extension_identity_inertias
            == (Inertia(8, 1, 0),) * 56
            and self.new_root_edge_extension_conjugate_inertias
            == (Inertia(9, 0, 0),) * 56
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status == "not_found_through_this_shell"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergThreePlusTwoSqrt2Extension:
    """The exact 60-wall shell ``x=3+2sqrt(2)``.

    The 60 accepted same-height roots have pairwise nonpositive identity
    pairings, so they enter simultaneously and independently of candidate
    order.  The resulting 153-wall list remains only a bounded norm-two
    prefix; this object records no chamber-completeness or volume claim.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 3 + 2 * SQRT2
            and self.next_admissible_control_coordinate == 4 + 2 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
            )
            and self.prior_wall_count == 93
            and self.new_level.control_coordinate == 3 + 2 * SQRT2
            and self.new_level.target == 43 + 29 * SQRT2
            and self.new_level.alphabet.values
            == (
                SILVER_ZERO,
                Silver(1),
                SQRT2,
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
                4 + 3 * SQRT2,
                5 + 3 * SQRT2,
                4 + 4 * SQRT2,
                6 + 3 * SQRT2,
                5 + 4 * SQRT2,
                6 + 4 * SQRT2,
                7 + 4 * SQRT2,
            )
            and self.new_level.visited_by_depth
            == (1, 18, 324, 1080, 2016, 2772, 20142, 66078, 47934, 16290)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 264, 968, 1862, 1653, 16471, 63415, 47029, 11339)
            and self.new_level.exact_weight_solutions == 138
            and self.new_level.nonintegral_reconstructions == 0
            and self.new_level.previous_wall_rejections == 78
            and len(self.new_level.previous_wall_rejection_witnesses) == 78
            and all(
                blockers
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            and len(self.new_level.roots) == 60
            and len(self.roots) == 153
            and self.wall_identity_inertia == Inertia(9, 1, 143)
            and self.wall_conjugate_inertia == Inertia(10, 0, 143)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 1702
            and self.same_height_zero_pair_count == 68
            and self.same_height_negative_pair_count
            + self.same_height_zero_pair_count
            == 60 * 59 // 2
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33,)
            and self.new_finite_endpoint_root_indices == ()
            and len(self.new_root_edge_pairings) == 60
            and len(self.new_root_edge_extension_determinants) == 60
            and all(
                determinant.sign(0) < 0
                for determinant in self.new_root_edge_extension_determinants
            )
            and self.new_root_edge_extension_identity_inertias
            == (Inertia(8, 1, 0),) * 60
            and self.new_root_edge_extension_conjugate_inertias
            == (Inertia(9, 0, 0),) * 60
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status == "not_found_through_this_shell"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergFourPlusTwoSqrt2Extension:
    """The exact 121-wall shell ``x=4+2sqrt(2)``.

    The accepted same-height roots have pairwise nonpositive identity
    pairings, so they enter simultaneously and independently of candidate
    order.  The resulting wall list remains only a bounded norm-two prefix;
    this object records no chamber-completeness or volume claim.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 4 + 2 * SQRT2
            and self.next_admissible_control_coordinate == 3 + 3 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
            )
            and self.prior_wall_count == 153
            and self.new_level.control_coordinate == 4 + 2 * SQRT2
            and self.new_level.target == 58 + 40 * SQRT2
            and len(self.new_level.alphabet.values) == 19
            and self.new_level.visited_by_depth
            == (1, 19, 361, 1178, 2261, 3211, 22363, 73568, 56012, 21299)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 299, 1059, 2092, 2034, 18491, 70620, 54891, 15783)
            and self.new_level.exact_weight_solutions == 653
            and self.new_level.nonintegral_reconstructions == 330
            and self.new_level.previous_wall_rejections == 202
            and len(self.new_level.previous_wall_rejection_witnesses) == 202
            and all(
                blockers
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            and len(self.new_level.roots) == 121
            and len(self.roots) == 274
            and self.wall_identity_inertia == Inertia(9, 1, 264)
            and self.wall_conjugate_inertia == Inertia(10, 0, 264)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 7180
            and self.same_height_zero_pair_count == 80
            and self.same_height_negative_pair_count
            + self.same_height_zero_pair_count
            == 121 * 120 // 2
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33,)
            and self.new_finite_endpoint_root_indices == ()
            and len(self.new_root_edge_pairings) == 121
            and len(self.new_root_edge_extension_determinants) == 121
            and all(
                determinant.sign(0) < 0
                for determinant in self.new_root_edge_extension_determinants
            )
            and self.new_root_edge_extension_identity_inertias
            == (Inertia(8, 1, 0),) * 121
            and self.new_root_edge_extension_conjugate_inertias
            == (Inertia(9, 0, 0),) * 121
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status == "not_found_through_this_shell"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergThreePlusThreeSqrt2Extension:
    """The exact 59-wall shell ``x=3+3sqrt(2)``.

    The accepted same-height roots have pairwise nonpositive identity
    pairings, so they enter simultaneously and independently of candidate
    order.  The resulting wall list remains only a bounded norm-two prefix;
    this object records no chamber-completeness or volume claim.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 3 + 3 * SQRT2
            and self.next_admissible_control_coordinate == 5 + 2 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
            )
            and self.prior_wall_count == 274
            and self.new_level.control_coordinate == 3 + 3 * SQRT2
            and self.new_level.target == 65 + 45 * SQRT2
            and len(self.new_level.alphabet.values) == 19
            and self.new_level.visited_by_depth
            == (1, 19, 361, 1197, 2318, 3154, 23218, 77463, 56506, 18943)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 298, 1075, 2152, 1932, 19141, 74489, 55509, 13412)
            and self.new_level.exact_weight_solutions == 154
            and self.new_level.nonintegral_reconstructions == 0
            and self.new_level.previous_wall_rejections == 95
            and len(self.new_level.previous_wall_rejection_witnesses) == 95
            and all(
                blockers
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            and len(self.new_level.roots) == 59
            and len(self.roots) == 333
            and self.wall_identity_inertia == Inertia(9, 1, 323)
            and self.wall_conjugate_inertia == Inertia(10, 0, 323)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 1689
            and self.same_height_zero_pair_count == 22
            and self.same_height_negative_pair_count
            + self.same_height_zero_pair_count
            == 59 * 58 // 2
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33,)
            and self.new_finite_endpoint_root_indices == ()
            and len(self.new_root_edge_pairings) == 59
            and len(self.new_root_edge_extension_determinants) == 59
            and all(
                determinant.sign(0) < 0
                for determinant in self.new_root_edge_extension_determinants
            )
            and self.new_root_edge_extension_identity_inertias
            == (Inertia(8, 1, 0),) * 59
            and self.new_root_edge_extension_conjugate_inertias
            == (Inertia(9, 0, 0),) * 59
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status == "not_found_through_this_shell"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergFivePlusTwoSqrt2Extension:
    """The exact empty shell ``x=5+2sqrt(2)``.

    Its unique integral norm-equation solution is blocked by the complete
    certified 333-wall lower prefix, so the shell contributes no wall.  The
    resulting wall list remains only a bounded norm-two prefix; this object
    records no chamber-completeness or volume claim.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 5 + 2 * SQRT2
            and self.next_admissible_control_coordinate == 4 + 3 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
                5 + 2 * SQRT2,
            )
            and self.prior_wall_count == 333
            and self.new_level.control_coordinate == 5 + 2 * SQRT2
            and self.new_level.target == 75 + 53 * SQRT2
            and self.new_level.alphabet
            == SilverCoordinateAlphabet(
                values=(
                    SILVER_ZERO,
                    3 + 2 * SQRT2,
                    4 + 3 * SQRT2,
                    7 + 5 * SQRT2,
                ),
                identity_upper=18,
                conjugate_upper=1,
                coefficient_bound=19,
            )
            and self.new_level.visited_by_depth
            == (1, 4, 16, 16, 16, 12, 32, 44, 24, 12)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 12, 12, 13, 4, 21, 38, 21, 4)
            and self.new_level.exact_weight_solutions == 1
            and self.new_level.nonintegral_reconstructions == 0
            and self.new_level.previous_wall_rejections == 1
            and len(self.new_level.previous_wall_rejection_witnesses) == 1
            and self.new_level.previous_wall_rejection_witnesses[0][1] == (9,)
            and self.new_level.roots == ()
            and len(self.roots) == 333
            and self.wall_identity_inertia == Inertia(9, 1, 323)
            and self.wall_conjugate_inertia == Inertia(10, 0, 323)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 0
            and self.same_height_zero_pair_count == 0
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33,)
            and self.new_finite_endpoint_root_indices == ()
            and self.new_root_edge_pairings == ()
            and self.new_root_edge_extension_determinants == ()
            and self.new_root_edge_extension_identity_inertias == ()
            and self.new_root_edge_extension_conjugate_inertias == ()
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status == "not_found_through_this_shell"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergFourPlusThreeSqrt2Extension:
    """The exact 393-wall shell ``x=4+3sqrt(2)``.

    Every integral reconstruction is tested against the certified 333-wall
    lower prefix, and every surviving root is tested against the complete
    same-height cohort.  This remains a bounded norm-two prefix certificate,
    not a chamber-completeness or volume theorem.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        new_count = len(self.new_level.roots)
        return (
            self.completed_through_control_coordinate == 4 + 3 * SQRT2
            and self.next_admissible_control_coordinate == 5 + 3 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
                5 + 2 * SQRT2,
                4 + 3 * SQRT2,
            )
            and self.prior_wall_count == 333
            and self.new_level.control_coordinate == 4 + 3 * SQRT2
            and self.new_level.target == 84 + 58 * SQRT2
            and len(self.new_level.alphabet.values) == 26
            and self.new_level.alphabet.identity_upper == 19
            and self.new_level.alphabet.conjugate_upper == 2
            and self.new_level.alphabet.coefficient_bound == 21
            and self.new_level.visited_by_depth
            == (1, 26, 676, 3016, 7462, 12688, 117078, 493064, 425750, 156000)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 560, 2729, 6974, 8185, 98114, 476689, 419750, 117374)
            and self.new_level.exact_weight_solutions == 2507
            and self.new_level.nonintegral_reconstructions == 1544
            and self.new_level.previous_wall_rejections == 570
            and len(self.new_level.previous_wall_rejection_witnesses) == 570
            and all(
                blockers
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            and new_count == 393
            and len(self.roots) == 726
            and self.wall_identity_inertia == Inertia(9, 1, 716)
            and self.wall_conjugate_inertia == Inertia(10, 0, 716)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 76453
            and self.same_height_zero_pair_count == 575
            and self.same_height_negative_pair_count
            + self.same_height_zero_pair_count
            == 393 * 392 // 2
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33, 622)
            and self.new_finite_endpoint_root_indices == (622,)
            and len(self.new_root_edge_pairings) == 393
            and len(self.new_root_edge_extension_determinants) == 393
            and sum(
                determinant.sign(0) < 0
                for determinant in self.new_root_edge_extension_determinants
            )
            == 392
            and sum(
                determinant.sign(0) > 0
                for determinant in self.new_root_edge_extension_determinants
            )
            == 1
            and self.new_root_edge_extension_identity_inertias.count(
                Inertia(8, 1, 0)
            )
            == 392
            and self.new_root_edge_extension_identity_inertias.count(
                Inertia(9, 0, 0)
            )
            == 1
            and self.new_root_edge_extension_conjugate_inertias
            == (Inertia(9, 0, 0),) * 393
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status == "finite_endpoint_found_at_this_shell"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergFivePlusThreeSqrt2Extension:
    """The exact 350-wall shell ``x=5+3sqrt(2)``.

    Every integral reconstruction is tested against the certified 726-wall
    lower prefix, and the complete same-height cohort is acute.  The resulting
    1,076-wall list remains only a bounded norm-two prefix: targeted-edge
    closure does not imply a global Vinberg stopping or volume theorem.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        new_count = len(self.new_level.roots)
        return (
            self.completed_through_control_coordinate == 5 + 3 * SQRT2
            and self.next_admissible_control_coordinate == 4 + 4 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
                5 + 2 * SQRT2,
                4 + 3 * SQRT2,
                5 + 3 * SQRT2,
            )
            and self.prior_wall_count == 726
            and self.new_level.control_coordinate == 5 + 3 * SQRT2
            and self.new_level.target == 105 + 73 * SQRT2
            and len(self.new_level.alphabet.values) == 27
            and self.new_level.alphabet.identity_upper == 21
            and self.new_level.alphabet.conjugate_upper == 2
            and self.new_level.alphabet.coefficient_bound == 23
            and self.new_level.visited_by_depth
            == (
                1,
                27,
                729,
                3456,
                8775,
                15120,
                146556,
                639441,
                553014,
                191916,
            )
            and self.new_level.schur_pruned_by_depth
            == (
                0,
                0,
                601,
                3131,
                8215,
                9692,
                122873,
                618959,
                545906,
                142182,
            )
            and self.new_level.exact_weight_solutions == 861
            and self.new_level.nonintegral_reconstructions == 0
            and self.new_level.previous_wall_rejections == 511
            and len(self.new_level.previous_wall_rejection_witnesses) == 511
            and sum(
                len(blockers) == 1
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            == 323
            and sum(
                len(blockers) == 2
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            == 188
            and all(
                blockers
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            and new_count == 350
            and len(self.roots) == 1076
            and self.wall_identity_inertia == Inertia(9, 1, 1066)
            and self.wall_conjugate_inertia == Inertia(10, 0, 1066)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 60759
            and self.same_height_zero_pair_count == 316
            and self.same_height_negative_pair_count
            + self.same_height_zero_pair_count
            == 350 * 349 // 2
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33, 622)
            and self.new_finite_endpoint_root_indices == ()
            and len(self.new_root_edge_pairings) == 350
            and len(self.new_root_edge_extension_determinants) == 350
            and all(
                determinant.sign(0) < 0
                for determinant in self.new_root_edge_extension_determinants
            )
            and self.new_root_edge_extension_identity_inertias
            == (Inertia(8, 1, 0),) * 350
            and self.new_root_edge_extension_conjugate_inertias
            == (Inertia(9, 0, 0),) * 350
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status
            == "already_found_at_lower_shell_no_new_endpoint"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergFourPlusFourSqrt2Extension:
    """The exact 142-wall shell ``x=4+4sqrt(2)``.

    Every norm-two candidate is exhausted against the certified 1,076-wall
    lower prefix.  The resulting 1,218-wall acute list remains incomplete;
    this object records neither global Vinberg stopping nor a volume claim.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    edge_indices: tuple[int, ...]
    finite_endpoint_root_indices: tuple[int, ...]
    new_finite_endpoint_root_indices: tuple[int, ...]
    new_root_edge_pairings: tuple[Vector, ...]
    new_root_edge_extension_determinants: tuple[Silver, ...]
    new_root_edge_extension_identity_inertias: tuple[Inertia, ...]
    new_root_edge_extension_conjugate_inertias: tuple[Inertia, ...]
    ideal_endpoint_audit_status: str
    second_endpoint_status: str
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 4 + 4 * SQRT2
            and self.next_admissible_control_coordinate == 6 + 3 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
                5 + 2 * SQRT2,
                4 + 3 * SQRT2,
                5 + 3 * SQRT2,
                4 + 4 * SQRT2,
            )
            and self.prior_wall_count == 1076
            and self.new_level.control_coordinate == 4 + 4 * SQRT2
            and self.new_level.target == 114 + 80 * SQRT2
            and len(self.new_level.alphabet.values) == 21
            and self.new_level.alphabet.identity_upper == 22
            and self.new_level.alphabet.conjugate_upper == 2
            and self.new_level.alphabet.coefficient_bound == 24
            and self.new_level.visited_by_depth
            == (1, 21, 441, 1533, 3192, 4746, 35049, 122157, 95718, 36141)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 368, 1381, 2966, 3077, 29232, 117599, 93997, 27184)
            and self.new_level.exact_weight_solutions == 957
            and self.new_level.nonintegral_reconstructions == 491
            and self.new_level.previous_wall_rejections == 324
            and len(self.new_level.previous_wall_rejection_witnesses) == 324
            and sum(
                len(blockers) == 1
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            == 255
            and sum(
                len(blockers) == 2
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            == 69
            and len(self.new_level.roots) == 142
            and len(self.roots) == 1218
            and self.wall_identity_inertia == Inertia(9, 1, 1208)
            and self.wall_conjugate_inertia == Inertia(10, 0, 1208)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 10011
            and self.same_height_zero_pair_count == 0
            and self.edge_indices == (1, 2, 4, 5, 6, 7, 9, 12)
            and self.finite_endpoint_root_indices == (33, 622)
            and self.new_finite_endpoint_root_indices == ()
            and len(self.new_root_edge_pairings) == 142
            and len(self.new_root_edge_extension_determinants) == 142
            and all(
                determinant.sign(0) < 0
                for determinant in self.new_root_edge_extension_determinants
            )
            and self.new_root_edge_extension_identity_inertias
            == (Inertia(8, 1, 0),) * 142
            and self.new_root_edge_extension_conjugate_inertias
            == (Inertia(9, 0, 0),) * 142
            and self.ideal_endpoint_audit_status
            == "not_run_until_complete_simple_roots"
            and self.second_endpoint_status
            == "already_found_at_lower_shell_no_new_endpoint"
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@dataclass(frozen=True, slots=True)
class SilverK3VinbergSixPlusThreeSqrt2Extension:
    """The exact 73-wall shell ``x=6+3sqrt(2)``.

    Every norm-two candidate is exhausted against the certified 1,218-wall
    lower prefix.  The resulting 1,291-wall acute list remains incomplete;
    global Vinberg stopping is audited separately from this shell object.
    """

    completed_through_control_coordinate: Silver
    next_admissible_control_coordinate: Silver
    candidate_coordinates: tuple[Silver, ...]
    prior_wall_count: int
    new_level: SilverK3HeightLevel
    roots: tuple[Vector, ...]
    wall_gram: SilverMatrix
    wall_identity_inertia: Inertia
    wall_conjugate_inertia: Inertia
    same_height_negative_pair_count: int
    same_height_zero_pair_count: int
    stopping_status: str

    @property
    def passes(self) -> bool:
        return (
            self.completed_through_control_coordinate == 6 + 3 * SQRT2
            and self.next_admissible_control_coordinate == 5 + 4 * SQRT2
            and self.candidate_coordinates
            == (
                Silver(1),
                SQRT2,
                Silver(2),
                1 + SQRT2,
                2 + SQRT2,
                1 + 2 * SQRT2,
                3 + SQRT2,
                2 + 2 * SQRT2,
                3 + 2 * SQRT2,
                4 + 2 * SQRT2,
                3 + 3 * SQRT2,
                5 + 2 * SQRT2,
                4 + 3 * SQRT2,
                5 + 3 * SQRT2,
                4 + 4 * SQRT2,
                6 + 3 * SQRT2,
            )
            and self.prior_wall_count == 1218
            and self.new_level.control_coordinate == 6 + 3 * SQRT2
            and self.new_level.target == 128 + 90 * SQRT2
            and len(self.new_level.alphabet.values) == 20
            and self.new_level.alphabet.identity_upper == 23
            and self.new_level.alphabet.conjugate_upper == 2
            and self.new_level.alphabet.coefficient_bound == 25
            and self.new_level.visited_by_depth
            == (1, 20, 400, 1380, 2780, 3920, 28820, 98760, 76020, 27780)
            and self.new_level.schur_pruned_by_depth
            == (0, 0, 331, 1241, 2584, 2479, 23882, 94959, 74631, 20465)
            and self.new_level.exact_weight_solutions == 697
            and self.new_level.nonintegral_reconstructions == 425
            and self.new_level.previous_wall_rejections == 199
            and len(self.new_level.previous_wall_rejection_witnesses) == 199
            and sum(
                len(blockers) == 1
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            == 127
            and sum(
                len(blockers) == 2
                for _root, blockers in self.new_level.previous_wall_rejection_witnesses
            )
            == 72
            and len(self.new_level.roots) == 73
            and len(self.roots) == 1291
            and self.wall_identity_inertia == Inertia(9, 1, 1281)
            and self.wall_conjugate_inertia == Inertia(10, 0, 1281)
            and all(
                self.wall_gram[i][j].sign(0) <= 0
                for i in range(len(self.roots))
                for j in range(i)
            )
            and self.same_height_negative_pair_count == 2628
            and self.same_height_zero_pair_count == 0
            and self.stopping_status == "undecidable_incomplete_simple_roots"
        )


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_prefix(identity_bound: int = 2) -> SilverK3VinbergPrefix:
    """Exhaust the intrinsic norm-two Vinberg prefix through ``x<=2``."""

    if identity_bound != 2:
        raise ValueError("the certified silver K3 prefix has identity bound two")
    gram = silver_k3_alexeev_field_gram()
    control = silver_k3_control_vector()
    covector = _silver_mat_vec(gram, control)
    cone = silver_k3_height_zero_walls()
    candidates = admissible_silver_control_coordinates(identity_bound)
    accepted = list(cone)
    levels = []
    for coordinate in candidates:
        level = _enumerate_height_level(coordinate, accepted)
        levels.append(level)
        accepted.extend(level.roots)
    roots = tuple(accepted)
    wall_gram = tuple(
        tuple(_silver_pair(left, gram, right) for right in roots) for left in roots
    )
    prefix = SilverK3VinbergPrefix(
        identity_bound=identity_bound,
        field_gram=gram,
        field_determinant=silver_determinant(gram),
        field_identity_inertia=silver_inertia(gram, embedding=0),
        field_conjugate_inertia=silver_inertia(gram, embedding=1),
        control_vector=control,
        control_norm=_silver_quadratic(control, gram),
        control_covector=covector,
        height_zero=audit_silver_k3_height_zero(),
        cone_roots=cone,
        candidate_coordinates=candidates,
        levels=tuple(levels),
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=silver_inertia(wall_gram, embedding=0),
        wall_conjugate_inertia=silver_inertia(wall_gram, embedding=1),
        diagram=_prefix_diagram_audit(wall_gram),
        saturation_wall=audit_silver_saturation_wall(),
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not prefix.passes:
        raise AssertionError("the silver K3 Vinberg prefix audit failed")
    return prefix


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_extension(
    identity_bound: int = 3,
) -> SilverK3VinbergExtension:
    """Exhaust the first admissible level beyond ``x=2``.

    The conjugate norm bound leaves only ``x=1+sqrt(2)`` between the old
    prefix and identity bound three.  The return value remains a bounded
    prefix: failure of its stopping diagnostic says nothing by itself about
    the volume of the eventual chamber.
    """

    if identity_bound != 3:
        raise ValueError("the certified silver K3 extension has identity bound three")
    prior = run_silver_k3_vinberg_prefix()
    candidates = admissible_silver_control_coordinates(identity_bound)
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    if unseen != (1 + SQRT2,):
        raise AssertionError("the first unseen silver control level changed")
    level = _enumerate_height_level(unseen[0], prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram = tuple(
        tuple(_silver_pair(left, gram, right) for right in roots)
        for left in roots
    )
    extension = SilverK3VinbergExtension(
        identity_bound=identity_bound,
        completed_through_control_coordinate=unseen[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=silver_inertia(wall_gram, embedding=0),
        wall_conjugate_inertia=silver_inertia(wall_gram, embedding=1),
        diagram=_prefix_diagram_audit(wall_gram),
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError("the silver K3 Vinberg extension audit failed")
    return extension


def _audit_silver_open_edge_endpoint(
    roots: tuple[Vector, ...],
    wall_gram: SilverMatrix,
    *,
    first_new_index: int,
) -> SilverEdgeEndpointAudit:
    edge = (1, 2, 4, 5, 6, 7, 9, 12)
    edge_gram = _subgram(wall_gram, edge)
    edge_inverse = _silver_inverse(edge_gram)

    # If a new norm-two root extends this spherical edge, each of its eight
    # pairings lies strictly between -2 and 2 at both embeddings.  Acuteness
    # fixes the identity sign.  The exact weight alphabet for target two
    # therefore leaves only 0, -1 and -sqrt(2).
    pairing_alphabet = tuple(
        -value
        for value in _weight_alphabet(Silver(2)).values
        if (Silver(4) - value * value).sign(0) > 0
        and (Silver(4) - value * value).sign(1) > 0
    )
    patterns = tuple(
        pattern
        for pattern in product(pairing_alphabet, repeat=len(edge))
        if (Silver(2) - _silver_quadratic(pattern, edge_inverse)).sign(0) > 0
        and (Silver(2) - _silver_quadratic(pattern, edge_inverse)).sign(1) > 0
    )
    new_endpoints = tuple(
        index
        for index in range(first_new_index, len(roots))
        if _is_elliptic(wall_gram, (*edge, index))
    )
    if len(new_endpoints) != 1:
        raise AssertionError("the targeted silver edge endpoint census changed")
    endpoint = new_endpoints[0]
    endpoint_nodes = (*edge, endpoint)
    audit = SilverEdgeEndpointAudit(
        edge_indices=edge,
        finite_type="A2 + F4 + B2",
        edge_gram=edge_gram,
        edge_determinant=silver_determinant(edge_gram),
        pairing_alphabet=pairing_alphabet,
        allowable_pairing_patterns=patterns,
        new_endpoint_root_indices=new_endpoints,
        endpoint_pairing=tuple(wall_gram[endpoint][node] for node in edge),
        endpoint_gram_determinant=silver_determinant(
            _subgram(wall_gram, endpoint_nodes)
        ),
        endpoint_identity_inertia=silver_inertia(
            _subgram(wall_gram, endpoint_nodes), embedding=0
        ),
        endpoint_conjugate_inertia=silver_inertia(
            _subgram(wall_gram, endpoint_nodes), embedding=1
        ),
    )
    if not audit.passes:
        raise AssertionError("the silver open-edge endpoint audit failed")
    return audit


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_edge_extension() -> SilverK3VinbergEdgeExtension:
    """Exhaust ``x=2+sqrt(2)`` and certify one targeted finite endpoint.

    This deliberately omits a global all-subset Vinberg stopping audit on the
    resulting 34-wall prefix.  Completeness and volume therefore remain open.
    """

    prior = run_silver_k3_vinberg_extension()
    search_box = admissible_silver_control_coordinates(4)
    upper = 2 + SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 1 + 2 * SQRT2:
        raise AssertionError("the second silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram = tuple(
        tuple(_silver_pair(left, gram, right) for right in roots)
        for left in roots
    )
    extension = SilverK3VinbergEdgeExtension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=silver_inertia(wall_gram, embedding=0),
        wall_conjugate_inertia=silver_inertia(wall_gram, embedding=1),
        edge_endpoint=_audit_silver_open_edge_endpoint(
            roots, wall_gram, first_new_index=len(prior.roots)
        ),
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError("the silver K3 Vinberg edge extension audit failed")
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_one_plus_two_sqrt2_extension(
) -> SilverK3VinbergOnePlusTwoSqrt2Extension:
    """Exhaust the next admissible shell ``x=1+2sqrt(2)`` exactly.

    The shell has four integral norm-equation solutions, all cut off by
    previously accepted simple roots, and therefore contributes no wall.  A
    missing finite endpoint of the targeted edge would have to be an accepted
    wall whose rank-nine extension is positive definite.  A missing ideal
    endpoint remains a separate full-parabolic-neighbourhood obligation once
    the complete simple-root list is known.
    """

    prior = run_silver_k3_vinberg_edge_extension()
    search_box = admissible_silver_control_coordinates(5)
    upper = 1 + 2 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 3 + SQRT2:
        raise AssertionError("the third silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram = tuple(
        tuple(_silver_pair(left, gram, right) for right in roots)
        for left in roots
    )
    edge = prior.edge_endpoint.edge_indices
    finite_endpoints = tuple(
        index
        for index in range(len(roots))
        if index not in edge and _is_elliptic(wall_gram, (*edge, index))
    )
    new_endpoints = tuple(
        index for index in finite_endpoints if index >= len(prior.roots)
    )
    extension = SilverK3VinbergOnePlusTwoSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=silver_inertia(wall_gram, embedding=0),
        wall_conjugate_inertia=silver_inertia(wall_gram, embedding=1),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=new_endpoints,
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status="not_found_through_this_shell",
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=1+2sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_three_plus_sqrt2_extension(
) -> SilverK3VinbergThreePlusSqrt2Extension:
    """Exhaust ``x=3+sqrt(2)`` against the complete lower-height prefix.

    Of twelve integral norm-equation solutions, nine are cut off by old
    walls.  The remaining three are mutually acute and therefore all enter
    the current Vinberg prefix.  None extends the targeted elliptic
    ``A2 + F4 + B2`` edge to an elliptic rank-nine diagram.
    """

    prior = run_silver_k3_vinberg_one_plus_two_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(6)
    upper = 3 + SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 2 + 2 * SQRT2:
        raise AssertionError("the fourth silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram = tuple(
        tuple(_silver_pair(left, gram, right) for right in roots)
        for left in roots
    )
    edge = prior.edge_indices
    finite_endpoints = tuple(
        index
        for index in range(len(roots))
        if index not in edge and _is_elliptic(wall_gram, (*edge, index))
    )
    new_indices = tuple(range(len(prior.roots), len(roots)))
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    extension = SilverK3VinbergThreePlusSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=silver_inertia(wall_gram, embedding=0),
        wall_conjugate_inertia=silver_inertia(wall_gram, embedding=1),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=tuple(
            index for index in finite_endpoints if index in new_indices
        ),
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status="not_found_through_this_shell",
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=3+sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_two_plus_two_sqrt2_extension(
) -> SilverK3VinbergTwoPlusTwoSqrt2Extension:
    """Exhaust ``x=2+2sqrt(2)`` against the frozen 37-wall prefix.

    Exactly 56 of the 277 weight solutions reconstruct to integral roots not
    cut off by a lower-height wall.  All 1,540 mutual pairings are
    nonpositive, so the entire cohort enters at once.  None extends the
    targeted ``A2 + F4 + B2`` edge to an elliptic rank-nine diagram.
    """

    prior = run_silver_k3_vinberg_three_plus_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(8)
    upper = 2 + 2 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 3 + 2 * SQRT2:
        raise AssertionError("the fifth silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram = tuple(
        tuple(_silver_pair(left, gram, right) for right in roots)
        for left in roots
    )
    edge = prior.edge_indices
    finite_endpoints = tuple(
        index
        for index in range(len(roots))
        if index not in edge and _is_elliptic(wall_gram, (*edge, index))
    )
    new_indices = tuple(range(len(prior.roots), len(roots)))
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergTwoPlusTwoSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=tuple(
            index for index in finite_endpoints if index in new_indices
        ),
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status="not_found_through_this_shell",
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=2+2sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_three_plus_two_sqrt2_extension(
) -> SilverK3VinbergThreePlusTwoSqrt2Extension:
    """Exhaust ``x=3+2sqrt(2)`` against the frozen 93-wall prefix.

    Exactly 60 of the 138 weight solutions are not cut off by a lower-height
    wall.  All 1,770 mutual pairings are nonpositive, so the entire cohort
    enters at once.  None extends the targeted ``A2 + F4 + B2`` edge to an
    elliptic rank-nine diagram.
    """

    prior = run_silver_k3_vinberg_two_plus_two_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(10)
    upper = 3 + 2 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 4 + 2 * SQRT2:
        raise AssertionError("the sixth silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram = tuple(
        tuple(_silver_pair(left, gram, right) for right in roots)
        for left in roots
    )
    edge = prior.edge_indices
    finite_endpoints = tuple(
        index
        for index in range(len(roots))
        if index not in edge and _is_elliptic(wall_gram, (*edge, index))
    )
    new_indices = tuple(range(len(prior.roots), len(roots)))
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergThreePlusTwoSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=tuple(
            index for index in finite_endpoints if index in new_indices
        ),
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status="not_found_through_this_shell",
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=3+2sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_four_plus_two_sqrt2_extension(
) -> SilverK3VinbergFourPlusTwoSqrt2Extension:
    """Exhaust ``x=4+2sqrt(2)`` against the frozen 153-wall prefix.

    Exactly 121 of the 323 integral norm-two roots are not cut off by a
    lower-height wall.  Their complete same-height pairing census determines
    whether the cohort enters simultaneously.  The targeted
    ``A2 + F4 + B2`` edge is audited only for finite endpoints here; an ideal
    endpoint still requires the complete simple-root zero-neighbourhood.
    """

    prior = run_silver_k3_vinberg_three_plus_two_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(10)
    upper = 4 + 2 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 3 + 3 * SQRT2:
        raise AssertionError("the seventh silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    # Reuse the frozen 153-wall Gram and evaluate each new symmetric entry
    # once.  This is mathematically identical to a rectangular rebuild but
    # materially shortens certificate reproduction at 274 walls.
    wall_gram_work = [
        [SILVER_ZERO for _ in range(len(roots))] for _ in range(len(roots))
    ]
    for row in range(len(prior.roots)):
        wall_gram_work[row][: len(prior.roots)] = prior.wall_gram[row]
    for row in range(len(prior.roots), len(roots)):
        gram_root = _silver_mat_vec(gram, roots[row])
        for column in range(row + 1):
            pairing = _silver_dot(roots[column], gram_root)
            wall_gram_work[row][column] = pairing
            wall_gram_work[column][row] = pairing
    wall_gram = tuple(tuple(row) for row in wall_gram_work)
    edge = prior.edge_indices
    finite_endpoints = tuple(
        index
        for index in range(len(roots))
        if index not in edge and _is_elliptic(wall_gram, (*edge, index))
    )
    new_indices = tuple(range(len(prior.roots), len(roots)))
    new_endpoints = tuple(
        index for index in finite_endpoints if index in new_indices
    )
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergFourPlusTwoSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=new_endpoints,
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status=(
            "finite_endpoint_found_at_this_shell"
            if new_endpoints
            else "not_found_through_this_shell"
        ),
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=4+2sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_three_plus_three_sqrt2_extension(
) -> SilverK3VinbergThreePlusThreeSqrt2Extension:
    """Exhaust ``x=3+3sqrt(2)`` against the frozen 274-wall prefix.

    Exactly 59 of the 154 integral norm-two roots are not cut off by a
    lower-height wall.  Their complete same-height pairing census determines
    whether the cohort enters simultaneously.  The targeted
    ``A2 + F4 + B2`` edge is audited only for finite endpoints here; an ideal
    endpoint still requires the complete simple-root zero-neighbourhood.
    """

    prior = run_silver_k3_vinberg_four_plus_two_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(10)
    upper = 3 + 3 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 5 + 2 * SQRT2:
        raise AssertionError("the eighth silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram_work = [
        [SILVER_ZERO for _ in range(len(roots))] for _ in range(len(roots))
    ]
    for row in range(len(prior.roots)):
        wall_gram_work[row][: len(prior.roots)] = prior.wall_gram[row]
    for row in range(len(prior.roots), len(roots)):
        gram_root = _silver_mat_vec(gram, roots[row])
        for column in range(row + 1):
            pairing = _silver_dot(roots[column], gram_root)
            wall_gram_work[row][column] = pairing
            wall_gram_work[column][row] = pairing
    wall_gram = tuple(tuple(row) for row in wall_gram_work)
    edge = prior.edge_indices
    finite_endpoints = tuple(
        index
        for index in range(len(roots))
        if index not in edge and _is_elliptic(wall_gram, (*edge, index))
    )
    new_indices = tuple(range(len(prior.roots), len(roots)))
    new_endpoints = tuple(
        index for index in finite_endpoints if index in new_indices
    )
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergThreePlusThreeSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=new_endpoints,
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status=(
            "finite_endpoint_found_at_this_shell"
            if new_endpoints
            else "not_found_through_this_shell"
        ),
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=3+3sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_five_plus_two_sqrt2_extension(
) -> SilverK3VinbergFivePlusTwoSqrt2Extension:
    """Exhaust ``x=5+2sqrt(2)`` against the certified 333-wall prefix.

    The unique reconstructed integral root is tested against all lower-height
    walls and is blocked.  The targeted ``A2 + F4 + B2`` edge is checked only
    for finite endpoints; ideal endpoints require the eventual complete
    simple-root list.
    """

    prior = run_silver_k3_vinberg_three_plus_three_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(10)
    upper = 5 + 2 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 4 + 3 * SQRT2:
        raise AssertionError("the ninth silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram_work = [
        [SILVER_ZERO for _ in range(len(roots))] for _ in range(len(roots))
    ]
    for row in range(len(prior.roots)):
        wall_gram_work[row][: len(prior.roots)] = prior.wall_gram[row]
    for row in range(len(prior.roots), len(roots)):
        gram_root = _silver_mat_vec(gram, roots[row])
        for column in range(row + 1):
            pairing = _silver_dot(roots[column], gram_root)
            wall_gram_work[row][column] = pairing
            wall_gram_work[column][row] = pairing
    wall_gram = tuple(tuple(row) for row in wall_gram_work)
    edge = prior.edge_indices
    finite_endpoints = tuple(
        index
        for index in range(len(roots))
        if index not in edge and _is_elliptic(wall_gram, (*edge, index))
    )
    new_indices = tuple(range(len(prior.roots), len(roots)))
    new_endpoints = tuple(
        index for index in finite_endpoints if index in new_indices
    )
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergFivePlusTwoSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=new_endpoints,
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status=(
            "finite_endpoint_found_at_this_shell"
            if new_endpoints
            else "not_found_through_this_shell"
        ),
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=5+2sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_four_plus_three_sqrt2_extension(
) -> SilverK3VinbergFourPlusThreeSqrt2Extension:
    """Exhaust ``x=4+3sqrt(2)`` against the certified 333-wall prefix."""

    prior = run_silver_k3_vinberg_five_plus_two_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(10)
    upper = 4 + 3 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 5 + 3 * SQRT2:
        raise AssertionError("the tenth silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram_work = [
        [SILVER_ZERO for _ in range(len(roots))] for _ in range(len(roots))
    ]
    for row in range(len(prior.roots)):
        wall_gram_work[row][: len(prior.roots)] = prior.wall_gram[row]
    for row in range(len(prior.roots), len(roots)):
        gram_root = _silver_mat_vec(gram, roots[row])
        for column in range(row + 1):
            pairing = _silver_dot(roots[column], gram_root)
            wall_gram_work[row][column] = pairing
            wall_gram_work[column][row] = pairing
    wall_gram = tuple(tuple(row) for row in wall_gram_work)
    edge = prior.edge_indices
    finite_endpoints = tuple(
        index
        for index in range(len(roots))
        if index not in edge and _is_elliptic(wall_gram, (*edge, index))
    )
    new_indices = tuple(range(len(prior.roots), len(roots)))
    new_endpoints = tuple(
        index for index in finite_endpoints if index in new_indices
    )
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergFourPlusThreeSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=new_endpoints,
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status=(
            "finite_endpoint_found_at_this_shell"
            if new_endpoints
            else "not_found_through_this_shell"
        ),
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=4+3sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_five_plus_three_sqrt2_extension(
) -> SilverK3VinbergFivePlusThreeSqrt2Extension:
    """Exhaust ``x=5+3sqrt(2)`` against the certified 726-wall prefix."""

    prior = run_silver_k3_vinberg_four_plus_three_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(12)
    upper = 5 + 3 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 4 + 4 * SQRT2:
        raise AssertionError("the eleventh silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram_work = [
        [SILVER_ZERO for _ in range(len(roots))] for _ in range(len(roots))
    ]
    for row in range(len(prior.roots)):
        wall_gram_work[row][: len(prior.roots)] = prior.wall_gram[row]
    for row in range(len(prior.roots), len(roots)):
        gram_root = _silver_mat_vec(gram, roots[row])
        for column in range(row + 1):
            pairing = _silver_dot(roots[column], gram_root)
            wall_gram_work[row][column] = pairing
            wall_gram_work[column][row] = pairing
    wall_gram = tuple(tuple(row) for row in wall_gram_work)
    edge = prior.edge_indices
    new_indices = tuple(range(len(prior.roots), len(roots)))
    new_endpoints = tuple(
        index
        for index in new_indices
        if _is_elliptic(wall_gram, (*edge, index))
    )
    finite_endpoints = (*prior.finite_endpoint_root_indices, *new_endpoints)
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergFivePlusThreeSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=new_endpoints,
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status="already_found_at_lower_shell_no_new_endpoint",
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=5+3sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_four_plus_four_sqrt2_extension(
) -> SilverK3VinbergFourPlusFourSqrt2Extension:
    """Exhaust ``x=4+4sqrt(2)`` against the certified 1,076-wall prefix."""

    prior = run_silver_k3_vinberg_five_plus_three_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(14)
    upper = 4 + 4 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 6 + 3 * SQRT2:
        raise AssertionError("the twelfth silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram_work = [
        [SILVER_ZERO for _ in range(len(roots))] for _ in range(len(roots))
    ]
    for row in range(len(prior.roots)):
        wall_gram_work[row][: len(prior.roots)] = prior.wall_gram[row]
    for row in range(len(prior.roots), len(roots)):
        gram_root = _silver_mat_vec(gram, roots[row])
        for column in range(row + 1):
            pairing = _silver_dot(roots[column], gram_root)
            wall_gram_work[row][column] = pairing
            wall_gram_work[column][row] = pairing
    wall_gram = tuple(tuple(row) for row in wall_gram_work)
    edge = prior.edge_indices
    new_indices = tuple(range(len(prior.roots), len(roots)))
    new_endpoints = tuple(
        index
        for index in new_indices
        if _is_elliptic(wall_gram, (*edge, index))
    )
    finite_endpoints = (*prior.finite_endpoint_root_indices, *new_endpoints)
    edge_extensions = tuple(
        _subgram(wall_gram, (*edge, index)) for index in new_indices
    )
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergFourPlusFourSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        edge_indices=edge,
        finite_endpoint_root_indices=finite_endpoints,
        new_finite_endpoint_root_indices=new_endpoints,
        new_root_edge_pairings=tuple(
            tuple(wall_gram[index][vertex] for vertex in edge)
            for index in new_indices
        ),
        new_root_edge_extension_determinants=tuple(
            silver_determinant(matrix) for matrix in edge_extensions
        ),
        new_root_edge_extension_identity_inertias=tuple(
            silver_inertia(matrix, embedding=0) for matrix in edge_extensions
        ),
        new_root_edge_extension_conjugate_inertias=tuple(
            silver_inertia(matrix, embedding=1) for matrix in edge_extensions
        ),
        ideal_endpoint_audit_status="not_run_until_complete_simple_roots",
        second_endpoint_status="already_found_at_lower_shell_no_new_endpoint",
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=4+4sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


@lru_cache(maxsize=1)
def run_silver_k3_vinberg_six_plus_three_sqrt2_extension(
) -> SilverK3VinbergSixPlusThreeSqrt2Extension:
    """Exhaust ``x=6+3sqrt(2)`` against the certified 1,218-wall prefix."""

    prior = run_silver_k3_vinberg_four_plus_four_sqrt2_extension()
    search_box = admissible_silver_control_coordinates(14)
    upper = 6 + 3 * SQRT2
    candidates = tuple(
        coordinate
        for coordinate in search_box
        if (upper - coordinate).sign(0) >= 0
    )
    unseen = tuple(
        coordinate
        for coordinate in candidates
        if coordinate not in prior.candidate_coordinates
    )
    after = tuple(
        coordinate
        for coordinate in search_box
        if (coordinate - upper).sign(0) > 0
    )
    if unseen != (upper,) or not after or after[0] != 5 + 4 * SQRT2:
        raise AssertionError("the thirteenth silver control-coordinate gate changed")

    level = _enumerate_height_level(upper, prior.roots)
    roots = (*prior.roots, *level.roots)
    gram = silver_k3_alexeev_field_gram()
    wall_gram_work = [
        [SILVER_ZERO for _ in range(len(roots))] for _ in range(len(roots))
    ]
    for row in range(len(prior.roots)):
        wall_gram_work[row][: len(prior.roots)] = prior.wall_gram[row]
    for row in range(len(prior.roots), len(roots)):
        gram_root = _silver_mat_vec(gram, roots[row])
        for column in range(row + 1):
            pairing = _silver_dot(roots[column], gram_root)
            wall_gram_work[row][column] = pairing
            wall_gram_work[column][row] = pairing
    wall_gram = tuple(tuple(row) for row in wall_gram_work)
    new_indices = tuple(range(len(prior.roots), len(roots)))
    distinct_same_height_pairings = tuple(
        wall_gram[index][other]
        for index in new_indices
        for other in new_indices
        if other < index
    )
    wall_identity_inertia, wall_conjugate_inertia = (
        _spanning_wall_gram_inertias(wall_gram)
    )
    extension = SilverK3VinbergSixPlusThreeSqrt2Extension(
        completed_through_control_coordinate=upper,
        next_admissible_control_coordinate=after[0],
        candidate_coordinates=candidates,
        prior_wall_count=len(prior.roots),
        new_level=level,
        roots=roots,
        wall_gram=wall_gram,
        wall_identity_inertia=wall_identity_inertia,
        wall_conjugate_inertia=wall_conjugate_inertia,
        same_height_negative_pair_count=sum(
            pairing.sign(0) < 0 for pairing in distinct_same_height_pairings
        ),
        same_height_zero_pair_count=sum(
            pairing.sign(0) == 0 for pairing in distinct_same_height_pairings
        ),
        stopping_status="undecidable_incomplete_simple_roots",
    )
    if not extension.passes:
        raise AssertionError(
            "the x=6+3sqrt(2) silver K3 Vinberg extension audit failed"
        )
    return extension


__all__ = [
    "IntegerShortVectorAudit",
    "SilverCoordinateAlphabet",
    "SilverEdgeEndpointAudit",
    "SilverHeightZeroAudit",
    "SilverK3HeightLevel",
    "SilverK3VinbergExtension",
    "SilverK3VinbergEdgeExtension",
    "SilverK3VinbergFivePlusTwoSqrt2Extension",
    "SilverK3VinbergFivePlusThreeSqrt2Extension",
    "SilverK3VinbergFourPlusThreeSqrt2Extension",
    "SilverK3VinbergFourPlusFourSqrt2Extension",
    "SilverK3VinbergFourPlusTwoSqrt2Extension",
    "SilverK3VinbergOnePlusTwoSqrt2Extension",
    "SilverK3VinbergThreePlusSqrt2Extension",
    "SilverK3VinbergThreePlusThreeSqrt2Extension",
    "SilverK3VinbergThreePlusTwoSqrt2Extension",
    "SilverK3VinbergTwoPlusTwoSqrt2Extension",
    "SilverK3VinbergSixPlusThreeSqrt2Extension",
    "SilverK3VinbergPrefix",
    "SilverPrefixDiagramAudit",
    "SilverSaturationWallAudit",
    "admissible_silver_control_coordinates",
    "audit_silver_k3_height_zero",
    "audit_silver_saturation_wall",
    "run_silver_k3_vinberg_prefix",
    "run_silver_k3_vinberg_extension",
    "run_silver_k3_vinberg_edge_extension",
    "run_silver_k3_vinberg_five_plus_two_sqrt2_extension",
    "run_silver_k3_vinberg_five_plus_three_sqrt2_extension",
    "run_silver_k3_vinberg_four_plus_three_sqrt2_extension",
    "run_silver_k3_vinberg_four_plus_four_sqrt2_extension",
    "run_silver_k3_vinberg_four_plus_two_sqrt2_extension",
    "run_silver_k3_vinberg_one_plus_two_sqrt2_extension",
    "run_silver_k3_vinberg_three_plus_sqrt2_extension",
    "run_silver_k3_vinberg_three_plus_three_sqrt2_extension",
    "run_silver_k3_vinberg_three_plus_two_sqrt2_extension",
    "run_silver_k3_vinberg_two_plus_two_sqrt2_extension",
    "run_silver_k3_vinberg_six_plus_three_sqrt2_extension",
    "silver_k3_alexeev_field_gram",
    "silver_k3_control_vector",
    "silver_k3_height_zero_walls",
]
