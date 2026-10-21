"""Exact Lanner and global-Gram constraints for the golden K3 track.

For norm-two roots over ``Z[phi]`` the only finite labels relevant to the
9/14 problem are 2, 3, and 5.  This module does three related jobs:

* exhaust the order-four and order-five arithmetic Lanner diagrams using
  only those labels;
* audit the rank-one Schur complement forced by a 14-wall Gram matrix of
  rank ten at any finite nine-wall vertex; and
* provide an exact 14-wall counterexample showing that all local elliptic,
  ADE, arithmetic, and missing-face tests can pass while the global Gram
  inertia is wrong.

The last object is deliberately a counterexample to *sufficiency*, not a
hyperbolic polytope.  It isolates the global rank/inertia equations that an
inverse construction still has to solve.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations, product
from typing import Sequence

from sympy import Matrix as SympyMatrix
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.golden import Golden, PHI, ZERO, RationalLike
from coxeter_polytopes.exact.matrix import Inertia, determinant, inertia, rank
from coxeter_polytopes.exact.trace_lattice import integral_gram, trace_gram

from .k3 import K3_TRACE_SCALAR
from .k3_design import trace_curve_gram


GoldenMatrix = tuple[tuple[Golden, ...], ...]
IntegerMatrix = tuple[tuple[int, ...], ...]
LabelledEdge = tuple[int, int, int]


def _coerce_norm_two_gram(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
) -> GoldenMatrix:
    gram = tuple(tuple(Golden.coerce(entry) for entry in row) for row in field_gram)
    size = len(gram)
    if size not in (4, 5, 14):
        raise ValueError("the golden Lanner/Schur audit expects order 4, 5, or 14")
    if any(len(row) != size for row in gram):
        raise ValueError("the field Gram matrix must be square")
    if any(gram[i][j] != gram[j][i] for i in range(size) for j in range(i)):
        raise ValueError("the field Gram matrix must be symmetric")
    if any(gram[i][i] != 2 for i in range(size)):
        raise ValueError("the audit requires norm-two roots")
    return gram


def _integer_matrix(matrix: Sequence[Sequence[int]]) -> IntegerMatrix:
    return tuple(tuple(int(entry) for entry in row) for row in matrix)


def _integer_determinant(matrix: Sequence[Sequence[int]]) -> int:
    value = determinant(matrix)
    if value.b or value.a.denominator != 1:
        raise ArithmeticError("an integral matrix had a nonintegral determinant")
    return value.a.numerator


def _discriminant_invariants(matrix: IntegerMatrix) -> tuple[int, ...]:
    normal = smith_normal_form(SympyMatrix(matrix), domain=ZZ)
    return tuple(
        abs(int(normal[i, i]))
        for i in range(len(matrix))
        if abs(int(normal[i, i])) > 1
    )


def _field_gram(order: int, edges: Sequence[LabelledEdge]) -> GoldenMatrix:
    gram = [
        [Golden(2) if i == j else ZERO for j in range(order)]
        for i in range(order)
    ]
    products = {3: Golden(-1), 5: -PHI}
    for left, right, label in edges:
        if left == right or not (0 <= left < order and 0 <= right < order):
            raise ValueError("a labelled edge has invalid endpoints")
        try:
            value = products[label]
        except KeyError as error:
            raise ValueError("non-right golden edges must have label 3 or 5") from error
        if gram[left][right]:
            raise ValueError("a labelled edge was repeated")
        gram[left][right] = gram[right][left] = value
    return tuple(tuple(row) for row in gram)


def _labels(field_gram: GoldenMatrix) -> tuple[tuple[int, ...], ...]:
    labels = [[2 for _ in row] for row in field_gram]
    for i in range(len(field_gram)):
        labels[i][i] = 1
        for j in range(i):
            value = field_gram[i][j]
            if value == ZERO:
                label = 2
            elif value == -1:
                label = 3
            elif value == -PHI:
                label = 5
            else:
                raise ValueError(
                    "a golden Lanner diagram uses only products 0, -1, and -phi"
                )
            labels[i][j] = labels[j][i] = label
    return tuple(tuple(row) for row in labels)


def _canonical_label_key(labels: Sequence[Sequence[int]]) -> tuple[int, ...]:
    size = len(labels)
    pairs = tuple(combinations(range(size), 2))
    return min(
        tuple(labels[permutation[i]][permutation[j]] for i, j in pairs)
        for permutation in permutations(range(size))
    )


def _connected(labels: Sequence[Sequence[int]]) -> bool:
    size = len(labels)
    reached = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for neighbour in range(size):
            if labels[vertex][neighbour] != 2 and neighbour not in reached:
                reached.add(neighbour)
                frontier.append(neighbour)
    return len(reached) == size


def _deleted_grams(field_gram: GoldenMatrix) -> tuple[GoldenMatrix, ...]:
    size = len(field_gram)
    return tuple(
        tuple(
            tuple(field_gram[i][j] for j in range(size) if j != omitted)
            for i in range(size)
            if i != omitted
        )
        for omitted in range(size)
    )


@dataclass(frozen=True, slots=True)
class _BlockSpec:
    name: str
    order: int
    edges: tuple[LabelledEdge, ...]
    curve_type: str


_BLOCK_SPECS = (
    _BlockSpec("A4(3-5-3)", 4, ((0, 1, 3), (1, 2, 5), (2, 3, 3)), "A4(-1)"),
    _BlockSpec("A4(5-3-5)", 4, ((0, 1, 5), (1, 2, 3), (2, 3, 5)), "A4(-1)"),
    _BlockSpec("D4(5,3,3)", 4, ((0, 1, 5), (0, 2, 3), (0, 3, 3)), "D4(-1)"),
    _BlockSpec(
        "C4(one 5)",
        4,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (3, 0, 3)),
        "affine A3(-1)",
    ),
    _BlockSpec(
        "C4(opposite 5s)",
        4,
        ((0, 1, 5), (1, 2, 3), (2, 3, 5), (3, 0, 3)),
        "affine A3(-1)",
    ),
    _BlockSpec(
        "A5(5-3-3-3)",
        5,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (3, 4, 3)),
        "A5(-1)",
    ),
    _BlockSpec(
        "A5(5-3-3-5)",
        5,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (3, 4, 5)),
        "A5(-1)",
    ),
    _BlockSpec(
        "D5(terminal 5)",
        5,
        ((0, 1, 5), (1, 2, 3), (2, 3, 3), (2, 4, 3)),
        "D5(-1)",
    ),
)


_SPEC_BY_KEY = {
    _canonical_label_key(_labels(_field_gram(spec.order, spec.edges))): spec
    for spec in _BLOCK_SPECS
}
if len(_SPEC_BY_KEY) != len(_BLOCK_SPECS):
    raise AssertionError("two named golden Lanner blocks became isomorphic")


@dataclass(frozen=True, slots=True)
class GoldenLannerAudit:
    """Exact arithmetic and trace certificate for one normalized diagram."""

    name: str | None
    order: int
    field_gram: GoldenMatrix
    canonical_label_key: tuple[int, ...]
    connected: bool
    proper_identity_inertias: tuple[Inertia, ...]
    field_identity_inertia: Inertia
    field_conjugate_inertia: Inertia
    field_determinant: Golden
    curve_gram: IntegerMatrix
    curve_type: str | None
    curve_inertia: Inertia
    restriction_trace_gram: IntegerMatrix
    restriction_trace_inertia: Inertia
    restriction_trace_determinant: int
    discriminant_group_invariants: tuple[int, ...]

    @property
    def proper_subdiagrams_elliptic(self) -> bool:
        return all(
            value == Inertia(self.order - 1, 0, 0)
            for value in self.proper_identity_inertias
        )

    @property
    def is_arithmetic_lanner(self) -> bool:
        return (
            self.connected
            and self.proper_subdiagrams_elliptic
            and self.field_identity_inertia == Inertia(self.order - 1, 1, 0)
            and self.field_conjugate_inertia == Inertia(self.order, 0, 0)
            and self.restriction_trace_inertia
            == Inertia(1, 2 * self.order - 1, 0)
        )


def audit_golden_lanner(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
) -> GoldenLannerAudit:
    """Audit an order-four or order-five norm-two golden diagram exactly."""

    gram = _coerce_norm_two_gram(field_gram)
    if len(gram) not in (4, 5):
        raise ValueError("a classified golden Lanner block has order four or five")
    labels = _labels(gram)
    key = _canonical_label_key(labels)
    spec = _SPEC_BY_KEY.get(key)
    curve_fraction = trace_curve_gram(gram)
    if any(entry.denominator != 1 for row in curve_fraction for entry in row):
        raise ArithmeticError("a golden root diagram had nonintegral curve trace")
    curve = tuple(
        tuple(entry.numerator for entry in row) for row in curve_fraction
    )
    traced = _integer_matrix(integral_gram(trace_gram(gram, K3_TRACE_SCALAR)))
    return GoldenLannerAudit(
        name=spec.name if spec else None,
        order=len(gram),
        field_gram=gram,
        canonical_label_key=key,
        connected=_connected(labels),
        proper_identity_inertias=tuple(
            inertia(deleted, embedding=0) for deleted in _deleted_grams(gram)
        ),
        field_identity_inertia=inertia(gram, embedding=0),
        field_conjugate_inertia=inertia(gram, embedding=1),
        field_determinant=determinant(gram),
        curve_gram=curve,
        curve_type=spec.curve_type if spec else None,
        curve_inertia=inertia(curve),
        restriction_trace_gram=traced,
        restriction_trace_inertia=inertia(traced),
        restriction_trace_determinant=_integer_determinant(traced),
        discriminant_group_invariants=_discriminant_invariants(traced),
    )


def _graph_prefilter(order: int, active_edges: tuple[tuple[int, int], ...]) -> bool:
    labels = [[2 for _ in range(order)] for _ in range(order)]
    for i in range(order):
        labels[i][i] = 1
    for left, right in active_edges:
        labels[left][right] = labels[right][left] = 3
    if not _connected(labels):
        return False
    curve = [
        [
            -2 if i == j else (1 if labels[i][j] != 2 else 0)
            for j in range(order)
        ]
        for i in range(order)
    ]
    # Every deletion is a face and hence must already be ADE after forgetting
    # the 3/5 decoration.  This exact integral prefilter leaves only the
    # A/D, cycle, and affine-star graph shapes relevant to the proof.
    return all(
        inertia(
            [
                [curve[i][j] for j in range(order) if j != omitted]
                for i in range(order)
                if i != omitted
            ]
        )
        == Inertia(0, order - 1, 0)
        for omitted in range(order)
    )


@lru_cache(maxsize=1)
def golden_lanner_blocks() -> tuple[GoldenLannerAudit, ...]:
    """Exhaust all order-four/five arithmetic Lanner blocks over ``Z[phi]``.

    The enumeration starts from every graph on four or five labelled nodes,
    decorates every non-right edge independently by 3 or 5, and applies only
    exact inertia tests.  Isomorphic labelled diagrams are then identified by
    the lexicographically least complete edge-label tuple.
    """

    representatives: dict[tuple[int, ...], GoldenMatrix] = {}
    for order in (4, 5):
        pairs = tuple(combinations(range(order), 2))
        for mask in range(1 << len(pairs)):
            active = tuple(
                pair for index, pair in enumerate(pairs) if mask & (1 << index)
            )
            if not _graph_prefilter(order, active):
                continue
            for decorations in product((3, 5), repeat=len(active)):
                gram = _field_gram(
                    order,
                    tuple(
                        (left, right, label)
                        for (left, right), label in zip(active, decorations)
                    ),
                )
                if inertia(gram, embedding=0) != Inertia(order - 1, 1, 0):
                    continue
                if inertia(gram, embedding=1) != Inertia(order, 0, 0):
                    continue
                if any(
                    inertia(deleted, embedding=0) != Inertia(order - 1, 0, 0)
                    for deleted in _deleted_grams(gram)
                ):
                    continue
                key = _canonical_label_key(_labels(gram))
                representatives.setdefault(key, gram)

    if set(representatives) != set(_SPEC_BY_KEY):
        missing = set(_SPEC_BY_KEY) - set(representatives)
        unexpected = set(representatives) - set(_SPEC_BY_KEY)
        raise AssertionError(
            f"golden Lanner census changed: missing={missing}, unexpected={unexpected}"
        )
    audited = {
        key: audit_golden_lanner(gram) for key, gram in representatives.items()
    }
    return tuple(
        audited[
            _canonical_label_key(_labels(_field_gram(spec.order, spec.edges)))
        ]
        for spec in _BLOCK_SPECS
    )


def _inverse_golden(matrix: GoldenMatrix) -> list[list[Golden]]:
    size = len(matrix)
    work = [
        [*matrix[i], *(Golden(int(i == j)) for j in range(size))]
        for i in range(size)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot is None:
            raise ZeroDivisionError("the vertex Gram matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[column][index]
                for index in range(2 * size)
            ]
    return [row[size:] for row in work]


def _all_two_by_two_minors_zero(matrix: GoldenMatrix) -> bool:
    size = len(matrix)
    return all(
        matrix[i][j] * matrix[k][ell] - matrix[i][ell] * matrix[k][j] == ZERO
        for i, k in combinations(range(size), 2)
        for j, ell in combinations(range(size), 2)
    )


@dataclass(frozen=True, slots=True)
class VertexComplementSchurAudit:
    """The rank-one complement forced at a finite simple vertex."""

    vertex: tuple[int, ...]
    complement: tuple[int, ...]
    full_identity_inertia: Inertia
    full_conjugate_inertia: Inertia
    vertex_identity_inertia: Inertia
    vertex_conjugate_inertia: Inertia
    schur_complement: GoldenMatrix
    schur_rank: int
    schur_identity_inertia: Inertia
    schur_conjugate_inertia: Inertia
    all_two_by_two_minors_zero: bool

    @property
    def expected_corank(self) -> int:
        return len(self.complement) - 1

    @property
    def passes(self) -> bool:
        dimension = len(self.vertex)
        corank = self.expected_corank
        return (
            self.full_identity_inertia == Inertia(dimension, 1, corank)
            and self.full_conjugate_inertia == Inertia(dimension + 1, 0, corank)
            and self.vertex_identity_inertia == Inertia(dimension, 0, 0)
            and self.vertex_conjugate_inertia == Inertia(dimension, 0, 0)
            and self.schur_rank == 1
            and self.schur_identity_inertia == Inertia(0, 1, corank)
            and self.schur_conjugate_inertia == Inertia(1, 0, corank)
            and self.all_two_by_two_minors_zero
        )


def audit_vertex_complement_schur(
    field_gram: Sequence[Sequence[Golden | RationalLike]],
    vertex: Sequence[int],
) -> VertexComplementSchurAudit:
    """Audit the complement Schur matrix at one finite vertex.

    For 14 walls in dimension nine, the vertex has size nine and its
    complement has size five.  A rank-ten Lorentzian Gram matrix forces the
    resulting five-by-five Schur matrix to have rank one, identity-place
    inertia ``(0,1,4)``, and conjugate-place inertia ``(1,0,4)``.
    """

    gram = tuple(tuple(Golden.coerce(entry) for entry in row) for row in field_gram)
    size = len(gram)
    if any(len(row) != size for row in gram):
        raise ValueError("the field Gram matrix must be square")
    if any(gram[i][j] != gram[j][i] for i in range(size) for j in range(i)):
        raise ValueError("the field Gram matrix must be symmetric")
    vertex_tuple = tuple(vertex)
    if len(vertex_tuple) != len(set(vertex_tuple)):
        raise ValueError("a vertex facet list contains a repeated index")
    if any(index < 0 or index >= size for index in vertex_tuple):
        raise IndexError("a vertex facet index is out of range")
    complement = tuple(index for index in range(size) if index not in vertex_tuple)
    if not complement:
        raise ValueError("the vertex must omit at least one facet")
    vertex_gram = tuple(
        tuple(gram[i][j] for j in vertex_tuple) for i in vertex_tuple
    )
    inverse = _inverse_golden(vertex_gram)
    schur = tuple(
        tuple(
            gram[left][right]
            - sum(
                (
                    gram[left][vertex_tuple[i]]
                    * inverse[i][j]
                    * gram[vertex_tuple[j]][right]
                    for i in range(len(vertex_tuple))
                    for j in range(len(vertex_tuple))
                ),
                start=ZERO,
            )
            for right in complement
        )
        for left in complement
    )
    return VertexComplementSchurAudit(
        vertex=vertex_tuple,
        complement=complement,
        full_identity_inertia=inertia(gram, embedding=0),
        full_conjugate_inertia=inertia(gram, embedding=1),
        vertex_identity_inertia=inertia(vertex_gram, embedding=0),
        vertex_conjugate_inertia=inertia(vertex_gram, embedding=1),
        schur_complement=schur,
        schur_rank=rank(schur),
        schur_identity_inertia=inertia(schur, embedding=0),
        schur_conjugate_inertia=inertia(schur, embedding=1),
        all_two_by_two_minors_zero=_all_two_by_two_minors_zero(schur),
    )


@dataclass(frozen=True, slots=True)
class GoldenBlockProductCounterexample:
    """Exact local-but-not-global 9/14 diagram certificate."""

    field_gram: GoldenMatrix
    curve_gram: IntegerMatrix
    vertices: tuple[tuple[int, ...], ...]
    missing_faces: tuple[tuple[int, ...], ...]
    global_identity_inertia: Inertia
    global_conjugate_inertia: Inertia
    global_trace_inertia: Inertia
    every_missing_face_lanner: bool
    every_vertex_identity_elliptic: bool
    every_vertex_conjugate_elliptic: bool
    every_vertex_curve_ade: bool
    effective_curve_pairings: bool
    every_vertex_schur_rank_five: bool
    sample_schur_audit: VertexComplementSchurAudit

    @property
    def has_no_missing_face_of_order_three(self) -> bool:
        return all(len(face) != 3 for face in self.missing_faces)

    @property
    def passes_local_constraints(self) -> bool:
        return (
            len(self.field_gram) == 14
            and len(self.vertices) == 128
            and self.has_no_missing_face_of_order_three
            and self.every_missing_face_lanner
            and self.every_vertex_identity_elliptic
            and self.every_vertex_conjugate_elliptic
            and self.every_vertex_curve_ade
            and self.effective_curve_pairings
        )

    @property
    def passes_global_h9_gram(self) -> bool:
        return (
            self.global_identity_inertia == Inertia(9, 1, 4)
            and self.global_conjugate_inertia == Inertia(10, 0, 4)
            and self.global_trace_inertia == Inertia(1, 19, 8)
            and self.sample_schur_audit.passes
        )


def _is_lanner(field_gram: GoldenMatrix) -> bool:
    order = len(field_gram)
    return (
        inertia(field_gram, embedding=0) == Inertia(order - 1, 1, 0)
        and inertia(field_gram, embedding=1) == Inertia(order, 0, 0)
        and all(
            inertia(deleted, embedding=0) == Inertia(order - 1, 0, 0)
            for deleted in _deleted_grams(field_gram)
        )
    )


@lru_cache(maxsize=1)
def golden_block_product_counterexample() -> GoldenBlockProductCounterexample:
    """Return the exact ``Delta3^2 x Delta1^3`` local countercertificate.

    Its five factor blocks have sizes ``4,4,2,2,2``.  Each four-block is the
    Lanner path ``3-5-3``; each two-block has dotted product ``-phi**2``; and
    different blocks are mutually orthogonal.  A vertex omits one wall from
    each four-block and chooses one wall from each two-block.
    """

    order_four = _field_gram(4, _BLOCK_SPECS[0].edges)
    gram = [[Golden(2) if i == j else ZERO for j in range(14)] for i in range(14)]
    for offset in (0, 4):
        for i in range(4):
            for j in range(4):
                gram[offset + i][offset + j] = order_four[i][j]
    dotted_product = -(PHI * PHI)
    for left, right in ((8, 9), (10, 11), (12, 13)):
        gram[left][right] = gram[right][left] = dotted_product
    field = tuple(tuple(row) for row in gram)

    missing_faces = (
        (0, 1, 2, 3),
        (4, 5, 6, 7),
        (8, 9),
        (10, 11),
        (12, 13),
    )
    vertices: list[tuple[int, ...]] = []
    for omitted_first in range(4):
        for omitted_second in range(4, 8):
            for choices in product((0, 1), repeat=3):
                vertex = tuple(
                    sorted(
                        [i for i in range(4) if i != omitted_first]
                        + [i for i in range(4, 8) if i != omitted_second]
                        + [8 + 2 * block + choice for block, choice in enumerate(choices)]
                    )
                )
                vertices.append(vertex)
    vertex_tuple = tuple(vertices)

    curve_fraction = trace_curve_gram(field)
    curve = tuple(
        tuple(entry.numerator for entry in row) for row in curve_fraction
    )
    traced = _integer_matrix(integral_gram(trace_gram(field, K3_TRACE_SCALAR)))
    vertex_field_grams = tuple(
        tuple(tuple(field[i][j] for j in vertex) for i in vertex)
        for vertex in vertex_tuple
    )
    vertex_curve_grams = tuple(
        tuple(tuple(curve[i][j] for j in vertex) for i in vertex)
        for vertex in vertex_tuple
    )
    schur_audits = tuple(
        audit_vertex_complement_schur(field, vertex) for vertex in vertex_tuple
    )
    return GoldenBlockProductCounterexample(
        field_gram=field,
        curve_gram=curve,
        vertices=vertex_tuple,
        missing_faces=missing_faces,
        global_identity_inertia=inertia(field, embedding=0),
        global_conjugate_inertia=inertia(field, embedding=1),
        global_trace_inertia=inertia(traced),
        every_missing_face_lanner=all(
            _is_lanner(
                tuple(tuple(field[i][j] for j in face) for i in face)
            )
            for face in missing_faces
        ),
        every_vertex_identity_elliptic=all(
            inertia(vertex_gram, embedding=0) == Inertia(9, 0, 0)
            for vertex_gram in vertex_field_grams
        ),
        every_vertex_conjugate_elliptic=all(
            inertia(vertex_gram, embedding=1) == Inertia(9, 0, 0)
            for vertex_gram in vertex_field_grams
        ),
        every_vertex_curve_ade=all(
            inertia(vertex_gram) == Inertia(0, 9, 0)
            for vertex_gram in vertex_curve_grams
        ),
        effective_curve_pairings=all(
            curve[i][j] >= 0 for i in range(14) for j in range(i)
        ),
        every_vertex_schur_rank_five=all(
            audit.schur_rank == 5
            and audit.schur_identity_inertia == Inertia(0, 5, 0)
            for audit in schur_audits
        ),
        sample_schur_audit=schur_audits[0],
    )

