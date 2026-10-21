"""Crystallographic trace lifts for the silver K3 eigenlattice track.

For ``F = Q(sqrt(2))`` use norm-two roots and the K3 trace scalar
``alpha=-(2+sqrt(2))/4``.  Restriction of scalars then gives

``I2(4) -> D4(-1)``, ``B3 -> D6(-1)``, and ``F4 -> E8(-1)``.

The stored unimodular matrices are literal integral isometries.  In
particular, the identifications do not rely only on rank and determinant.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations, product

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
from coxeter_polytopes.lattices import e8_cartan_gram


IntegerMatrix = tuple[tuple[int, ...], ...]
SilverMatrix = tuple[tuple[Silver, ...], ...]
LabelledEdge = tuple[int, int, int]


def _path_gram(labels: tuple[int, ...]) -> SilverMatrix:
    size = len(labels) + 1
    gram = [
        [Silver(2) if i == j else SILVER_ZERO for j in range(size)]
        for i in range(size)
    ]
    products = {3: Silver(-1), 4: -SQRT2}
    for index, label in enumerate(labels):
        try:
            value = products[label]
        except KeyError as error:
            raise ValueError("a silver finite path uses only labels 3 and 4") from error
        gram[index][index + 1] = gram[index + 1][index] = value
    return tuple(tuple(row) for row in gram)


def _a_cartan(rank: int) -> list[list[int]]:
    gram = [[2 if i == j else 0 for j in range(rank)] for i in range(rank)]
    for index in range(rank - 1):
        gram[index][index + 1] = gram[index + 1][index] = -1
    return gram


def _d_cartan(rank: int) -> list[list[int]]:
    if rank < 4:
        raise ValueError("D rank must be at least four")
    gram = [[2 if i == j else 0 for j in range(rank)] for i in range(rank)]
    chain_end = rank - 3
    for index in range(chain_end):
        gram[index][index + 1] = gram[index + 1][index] = -1
    branch = chain_end
    for leaf in (rank - 2, rank - 1):
        gram[branch][leaf] = gram[leaf][branch] = -1
    return gram


def _negative(matrix: list[list[int]]) -> IntegerMatrix:
    return tuple(tuple(-entry for entry in row) for row in matrix)


_FIELD_GRAMS: dict[str, SilverMatrix] = {
    "I2(4)": _path_gram((4,)),
    "B3": _path_gram((4, 3)),
    "F4": _path_gram((3, 4, 3)),
}

_TRACE_TYPES = {"I2(4)": "D4(-1)", "B3": "D6(-1)", "F4": "E8(-1)"}

_STANDARD_TRACE_GRAMS: dict[str, IntegerMatrix] = {
    "I2(4)": _negative(_d_cartan(4)),
    "B3": _negative(_d_cartan(6)),
    "F4": _negative(e8_cartan_gram()),
}

# Columns are trace-lattice vectors.  For each field type,
# T^t Q_trace T = Q_standard and det(T)=+/-1.
_TRACE_ISOMETRIES: dict[str, IntegerMatrix] = {
    "I2(4)": (
        (1, -2, 1, 1),
        (0, 1, -1, -1),
        (0, 1, -2, 0),
        (0, -1, 1, 0),
    ),
    "B3": (
        (1, -2, 0, 1, 0, 0),
        (0, 1, 0, -1, 0, 1),
        (0, 1, 0, -2, 0, 2),
        (0, -1, 0, 1, 0, 0),
        (0, 0, 1, -2, 1, 1),
        (0, 0, -1, 1, 0, 0),
    ),
    "F4": (
        (1, -2, 0, 0, 0, 1, 0, 1),
        (0, 1, -1, 1, 0, -1, 0, 0),
        (0, -2, -1, 0, 1, 1, 0, 2),
        (0, 1, -1, 2, -1, -1, 0, 0),
        (0, 1, -1, 2, -1, -2, 0, 1),
        (0, -1, -1, 0, 1, 1, 0, 1),
        (0, 0, 0, 1, 0, -2, 1, 0),
        (0, 0, -1, 0, 0, 1, 0, 1),
    ),
}

_LANNER_TRACE_TYPES = {
    4: "U + A6(-1)",
    5: "E8(-1) + U(2)",
}

_LANNER_TRACE_ISOMETRIES: dict[int, IntegerMatrix] = {
    4: (
        (1, 0, 0, 0, -1, -1, 1, 0),
        (0, 1, 0, 0, 0, 1, 0, 0),
        (1, 1, 0, 0, 0, 1, 1, 0),
        (0, 0, 0, 0, 0, -1, 0, 0),
        (1, 0, -1, 0, 0, 1, 0, 0),
        (0, 1, 0, 0, 0, -1, 1, 0),
        (1, 0, -1, 1, -1, 0, 0, 1),
        (0, 1, 0, 0, 0, 0, 1, -1),
    ),
    5: (
        (1, -2, 0, 0, 0, 1, 0, 1, 0, -2),
        (0, 1, -1, 1, 0, -1, 0, 0, -1, 0),
        (0, -2, -1, 0, 1, 1, 0, 2, -1, -2),
        (0, 1, -1, 2, -1, -1, 0, 0, -1, -1),
        (0, 1, -1, 2, -1, -2, 0, 1, -1, -2),
        (0, -1, -1, 0, 1, 1, 0, 1, -1, -1),
        (0, 0, 0, 1, 0, -2, 1, 0, 0, -2),
        (0, 0, -1, 0, 0, 1, 0, 1, -1, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 1, -2),
        (0, 0, 0, 0, 0, 0, 0, 0, -1, 1),
    ),
}


def _integer_determinant(matrix: IntegerMatrix) -> int:
    value = determinant(matrix)
    if value.b or value.a.denominator != 1:
        raise AssertionError("an integer matrix had a nonintegral determinant")
    return value.a.numerator


def _congruence(change: IntegerMatrix, gram: IntegerMatrix) -> IntegerMatrix:
    size = len(change)
    return tuple(
        tuple(
            sum(
                change[k][i] * gram[k][ell] * change[ell][j]
                for k in range(size)
                for ell in range(size)
            )
            for j in range(size)
        )
        for i in range(size)
    )


def _block_diagonal(*blocks: IntegerMatrix) -> IntegerMatrix:
    size = sum(len(block) for block in blocks)
    result = [[0 for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                result[offset + i][offset + j] = entry
        offset += len(block)
    return tuple(tuple(row) for row in result)


def _lanner_standard_trace_gram(order: int) -> IntegerMatrix:
    if order == 4:
        return _block_diagonal(((0, 1), (1, 0)), _negative(_a_cartan(6)))
    if order == 5:
        return _block_diagonal(
            _negative(e8_cartan_gram()), ((0, 2), (2, 0))
        )
    raise ValueError("a silver Lanner block has order four or five")


@dataclass(frozen=True, slots=True)
class SilverFiniteTraceLift:
    """Exact restriction-of-scalars certificate for a silver finite type."""

    field_type: str
    trace_type: str
    field_gram: SilverMatrix
    trace_gram: IntegerMatrix
    standard_trace_gram: IntegerMatrix
    trace_isometry: IntegerMatrix
    field_determinant: Silver
    trace_determinant: int
    field_identity_inertia: Inertia
    field_conjugate_inertia: Inertia
    trace_inertia: Inertia
    isometry_determinant: int
    isometry_verified: bool

    @property
    def field_rank(self) -> int:
        return len(self.field_gram)

    @property
    def trace_rank(self) -> int:
        return len(self.trace_gram)

    @property
    def orthogonally_splits_every_integral_overlattice(self) -> bool:
        return abs(self.trace_determinant) == 1

    @property
    def passes(self) -> bool:
        rank = self.field_rank
        return (
            self.field_identity_inertia == Inertia(rank, 0, 0)
            and self.field_conjugate_inertia == Inertia(rank, 0, 0)
            and self.trace_inertia == Inertia(0, 2 * rank, 0)
            and abs(self.isometry_determinant) == 1
            and self.isometry_verified
        )


def silver_finite_trace_lift(field_type: str) -> SilverFiniteTraceLift:
    """Return an exact crystallographic silver trace-lift certificate."""

    try:
        field_gram = _FIELD_GRAMS[field_type]
        standard = _STANDARD_TRACE_GRAMS[field_type]
        change = _TRACE_ISOMETRIES[field_type]
    except KeyError as error:
        raise ValueError("field_type must be one of I2(4), B3, F4") from error
    traced = tuple(tuple(row) for row in integral_gram(silver_trace_gram(field_gram)))
    return SilverFiniteTraceLift(
        field_type=field_type,
        trace_type=_TRACE_TYPES[field_type],
        field_gram=field_gram,
        trace_gram=traced,
        standard_trace_gram=standard,
        trace_isometry=change,
        field_determinant=silver_determinant(field_gram),
        trace_determinant=_integer_determinant(traced),
        field_identity_inertia=silver_inertia(field_gram, embedding=0),
        field_conjugate_inertia=silver_inertia(field_gram, embedding=1),
        trace_inertia=inertia(traced),
        isometry_determinant=_integer_determinant(change),
        isometry_verified=standard == _congruence(change, traced),
    )


@dataclass(frozen=True, slots=True)
class ForcedSilverCore:
    """Finite divisor core forced inside a silver Lanner missing face."""

    missing_face_size: int
    field_type: str
    trace_type: str
    trace_rank: int
    trace_determinant: int
    orthogonally_splits: bool


def forced_silver_core(missing_face_size: int) -> ForcedSilverCore:
    """Return the core forced by an order-four or order-five missing face.

    The standard missing-face reduction gives ``H3`` or ``B3`` in order
    four and ``H4`` or ``F4`` in order five.  A norm-two form over
    ``Q(sqrt(2))`` cannot contain label five, so the silver alternatives are
    respectively ``B3`` and ``F4``.
    """

    field_type_by_size = {4: "B3", 5: "F4"}
    try:
        lift = silver_finite_trace_lift(field_type_by_size[missing_face_size])
    except KeyError as error:
        raise ValueError(
            "the forced-core theorem applies to missing faces of size 4 or 5"
        ) from error
    return ForcedSilverCore(
        missing_face_size=missing_face_size,
        field_type=lift.field_type,
        trace_type=lift.trace_type,
        trace_rank=lift.trace_rank,
        trace_determinant=lift.trace_determinant,
        orthogonally_splits=lift.orthogonally_splits_every_integral_overlattice,
    )


def _silver_lanner_gram(
    order: int, edges: tuple[LabelledEdge, ...]
) -> SilverMatrix:
    gram = [
        [Silver(2) if i == j else SILVER_ZERO for j in range(order)]
        for i in range(order)
    ]
    products = {3: Silver(-1), 4: -SQRT2}
    for left, right, label in edges:
        gram[left][right] = gram[right][left] = products[label]
    return tuple(tuple(row) for row in gram)


def _silver_label_matrix(field_gram: SilverMatrix) -> tuple[tuple[int, ...], ...]:
    labels = [[2 for _ in row] for row in field_gram]
    for index in range(len(field_gram)):
        labels[index][index] = 1
        for other in range(index):
            value = field_gram[index][other]
            if value == SILVER_ZERO:
                label = 2
            elif value == -1:
                label = 3
            elif value == -SQRT2:
                label = 4
            else:
                raise ValueError("a silver Lanner diagram uses labels 2, 3, and 4")
            labels[index][other] = labels[other][index] = label
    return tuple(tuple(row) for row in labels)


def _canonical_label_key(labels: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    order = len(labels)
    pairs = tuple(combinations(range(order), 2))
    return min(
        tuple(labels[permutation[i]][permutation[j]] for i, j in pairs)
        for permutation in permutations(range(order))
    )


def _connected_edges(order: int, edges: tuple[tuple[int, int], ...]) -> bool:
    adjacency = [set() for _ in range(order)]
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    reached = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        new = adjacency[vertex] - reached
        reached.update(new)
        frontier.extend(new)
    return len(reached) == order


def _deleted_silver_grams(field_gram: SilverMatrix) -> tuple[SilverMatrix, ...]:
    order = len(field_gram)
    return tuple(
        tuple(
            tuple(field_gram[i][j] for j in range(order) if j != omitted)
            for i in range(order)
            if i != omitted
        )
        for omitted in range(order)
    )


def _silver_graph_prefilter(
    order: int, active_edges: tuple[tuple[int, int], ...]
) -> bool:
    """Forget 3/4 decorations and require every deletion to be ADE."""

    if not _connected_edges(order, active_edges):
        return False
    adjacency = {tuple(sorted(edge)) for edge in active_edges}
    cartan = [
        [
            2
            if i == j
            else (-1 if tuple(sorted((i, j))) in adjacency else 0)
            for j in range(order)
        ]
        for i in range(order)
    ]
    return all(
        inertia(
            [
                [cartan[i][j] for j in range(order) if j != omitted]
                for i in range(order)
                if i != omitted
            ]
        )
        == Inertia(order - 1, 0, 0)
        for omitted in range(order)
    )


@dataclass(frozen=True, slots=True)
class SilverLannerAudit:
    """One exact arithmetic order-four/five silver Lanner block."""

    order: int
    canonical_label_key: tuple[int, ...]
    field_gram: SilverMatrix
    field_determinant: Silver
    field_identity_inertia: Inertia
    field_conjugate_inertia: Inertia
    proper_identity_inertias: tuple[Inertia, ...]
    trace_gram: IntegerMatrix
    trace_type: str
    trace_determinant: int
    trace_inertia: Inertia
    trace_isometry: IntegerMatrix
    trace_isometry_determinant: int
    trace_isometry_verified: bool
    curve_gram: IntegerMatrix
    curve_inertia: Inertia
    fiber_square: int
    fiber_sqrt2_pairing: int
    sqrt2_fiber_square: int
    primitive_u2_gram: tuple[tuple[int, int], tuple[int, int]]

    @property
    def passes(self) -> bool:
        return (
            self.field_identity_inertia == Inertia(self.order - 1, 1, 0)
            and self.field_conjugate_inertia == Inertia(self.order, 0, 0)
            and all(
                value == Inertia(self.order - 1, 0, 0)
                for value in self.proper_identity_inertias
            )
            and self.trace_inertia == Inertia(1, 2 * self.order - 1, 0)
            and abs(self.trace_isometry_determinant) == 1
            and self.trace_isometry_verified
            and self.curve_inertia == Inertia(0, self.order - 1, 1)
            and self.primitive_u2_gram == ((0, 2), (2, 0))
        )


@dataclass(frozen=True, slots=True)
class SilverLannerEnumerationStage:
    order: int
    graph_masks: int
    graph_prefilter_survivors: int
    decorated_candidates: int
    labelled_exact_survivors: int
    isomorphism_types: int


@dataclass(frozen=True, slots=True)
class SilverLannerCensus:
    stages: tuple[SilverLannerEnumerationStage, ...]
    blocks: tuple[SilverLannerAudit, ...]

    @property
    def passes(self) -> bool:
        return (
            self.stages
            == (
                SilverLannerEnumerationStage(4, 64, 19, 176, 12, 1),
                SilverLannerEnumerationStage(5, 1024, 137, 2384, 60, 1),
            )
            and tuple(block.order for block in self.blocks) == (4, 5)
            and all(block.passes for block in self.blocks)
        )


def _audit_silver_lanner(field_gram: SilverMatrix) -> SilverLannerAudit:
    order = len(field_gram)
    traced = tuple(tuple(row) for row in integral_gram(silver_trace_gram(field_gram)))
    curve = tuple(
        tuple(traced[2 * i][2 * j] for j in range(order))
        for i in range(order)
    )
    fiber = tuple(1 for _ in range(order))
    sqrt2_fiber = tuple(
        1 if index % 2 else 0 for index in range(2 * order)
    )
    fiber_trace = tuple(
        1 if index % 2 == 0 else 0 for index in range(2 * order)
    )

    def traced_pair(left: tuple[int, ...], right: tuple[int, ...]) -> int:
        return sum(
            left[i] * traced[i][j] * right[j]
            for i in range(2 * order)
            for j in range(2 * order)
        )

    trace_determinant = _integer_determinant(traced)
    trace_isometry = _LANNER_TRACE_ISOMETRIES[order]
    return SilverLannerAudit(
        order=order,
        canonical_label_key=_canonical_label_key(_silver_label_matrix(field_gram)),
        field_gram=field_gram,
        field_determinant=silver_determinant(field_gram),
        field_identity_inertia=silver_inertia(field_gram, embedding=0),
        field_conjugate_inertia=silver_inertia(field_gram, embedding=1),
        proper_identity_inertias=tuple(
            silver_inertia(deleted, embedding=0)
            for deleted in _deleted_silver_grams(field_gram)
        ),
        trace_gram=traced,
        trace_type=_LANNER_TRACE_TYPES[order],
        trace_determinant=trace_determinant,
        trace_inertia=inertia(traced),
        trace_isometry=trace_isometry,
        trace_isometry_determinant=_integer_determinant(trace_isometry),
        trace_isometry_verified=(
            _lanner_standard_trace_gram(order)
            == _congruence(trace_isometry, traced)
        ),
        curve_gram=curve,
        curve_inertia=inertia(curve),
        fiber_square=sum(
            fiber[i] * curve[i][j] * fiber[j]
            for i in range(order)
            for j in range(order)
        ),
        fiber_sqrt2_pairing=traced_pair(fiber_trace, sqrt2_fiber),
        sqrt2_fiber_square=traced_pair(sqrt2_fiber, sqrt2_fiber),
        primitive_u2_gram=(
            (traced_pair(fiber_trace, fiber_trace), traced_pair(fiber_trace, sqrt2_fiber)),
            (traced_pair(sqrt2_fiber, fiber_trace), traced_pair(sqrt2_fiber, sqrt2_fiber)),
        ),
    )


@lru_cache(maxsize=1)
def silver_lanner_census() -> SilverLannerCensus:
    """Exhaust the arithmetic Lanner blocks over ``Z[sqrt(2)]``.

    Every graph on four or five labelled nodes and every independent 3/4
    decoration is considered.  Exact inertia at both real embeddings and on
    every one-node deletion leaves one isomorphism type in each order: a
    cycle with exactly one edge labelled four.
    """

    representatives: dict[tuple[int, ...], SilverMatrix] = {}
    stages: list[SilverLannerEnumerationStage] = []
    for order in (4, 5):
        pairs = tuple(combinations(range(order), 2))
        prefilter_survivors = 0
        decorated_candidates = 0
        labelled_survivors = 0
        keys_before = set(representatives)
        for mask in range(1 << len(pairs)):
            active = tuple(
                edge for index, edge in enumerate(pairs) if mask & (1 << index)
            )
            if not _silver_graph_prefilter(order, active):
                continue
            prefilter_survivors += 1
            for decorations in product((3, 4), repeat=len(active)):
                decorated_candidates += 1
                gram = _silver_lanner_gram(
                    order,
                    tuple(
                        (left, right, label)
                        for (left, right), label in zip(
                            active, decorations, strict=True
                        )
                    ),
                )
                if silver_inertia(gram, embedding=0) != Inertia(
                    order - 1, 1, 0
                ):
                    continue
                if silver_inertia(gram, embedding=1) != Inertia(order, 0, 0):
                    continue
                if any(
                    silver_inertia(deleted, embedding=0)
                    != Inertia(order - 1, 0, 0)
                    for deleted in _deleted_silver_grams(gram)
                ):
                    continue
                labelled_survivors += 1
                key = _canonical_label_key(_silver_label_matrix(gram))
                representatives.setdefault(key, gram)
        stages.append(
            SilverLannerEnumerationStage(
                order=order,
                graph_masks=1 << len(pairs),
                graph_prefilter_survivors=prefilter_survivors,
                decorated_candidates=decorated_candidates,
                labelled_exact_survivors=labelled_survivors,
                isomorphism_types=len(set(representatives) - keys_before),
            )
        )
    standard_grams = (
        _silver_lanner_gram(
            4,
            ((0, 1, 4), (1, 2, 3), (2, 3, 3), (3, 0, 3)),
        ),
        _silver_lanner_gram(
            5,
            (
                (0, 1, 3),
                (1, 2, 4),
                (2, 3, 3),
                (3, 4, 3),
                (4, 0, 3),
            ),
        ),
    )
    expected = {
        _canonical_label_key(_silver_label_matrix(gram)) for gram in standard_grams
    }
    if set(representatives) != expected:
        raise AssertionError("the exhaustive silver Lanner types changed")
    audits = tuple(_audit_silver_lanner(gram) for gram in standard_grams)
    if tuple(audit.order for audit in audits) != (4, 5):
        raise AssertionError("the silver Lanner census is not one block per order")
    if any(audit.canonical_label_key.count(4) != 1 for audit in audits):
        raise AssertionError("a silver Lanner block is not a one-4 cycle")
    if not all(audit.passes for audit in audits):
        raise AssertionError("a silver Lanner audit failed its exact checks")
    census = SilverLannerCensus(tuple(stages), audits)
    if not census.passes:
        raise AssertionError("the silver Lanner enumeration counts changed")
    return census


def silver_lanner_blocks() -> tuple[SilverLannerAudit, ...]:
    """Return the two exact blocks from :func:`silver_lanner_census`."""

    return silver_lanner_census().blocks
