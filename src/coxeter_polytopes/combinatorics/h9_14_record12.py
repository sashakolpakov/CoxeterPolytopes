"""Exact rational ordinary-polytopal completion of H9/14 record 12.

Start with the cyclic simplicial 5-polytope ``C(8,5)``.  Stellar-subdivide
the triangle ``{0,2,4}`` by vertex 8, then ``{0,1,8}`` by vertex 9, and
apply one simplicial wedge at each of vertices 1, 2, 3, and 4.  A fixed
relabeling gives record 12.  Exact rational coordinates certify both stellar
subdivisions, while minimal-nonface expansion certifies the four wedges.

This proves ordinary convex polytopality only.  It supplies neither Coxeter
labels nor a Lorentzian Gram matrix.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 12
TARGET_GRAPH6 = "M???????????o?o??"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
TARGET_MISSING_EDGES = (
    (0, 12),
    (1, 12),
    (0, 13),
    (1, 13),
)
TARGET_MISSING_FOURS = (
    (2, 3, 8, 12),
    (5, 7, 10, 13),
)
TARGET_MISSING_FIVES = (
    (1, 2, 4, 7, 9),
    (1, 2, 5, 7, 10),
    (0, 3, 5, 8, 10),
    (3, 5, 7, 8, 10),
    (0, 3, 6, 8, 11),
    (0, 4, 6, 9, 11),
    (2, 4, 6, 9, 11),
    (4, 6, 9, 11, 13),
)
TARGET_F_VECTOR = (14, 87, 320, 778, 1308, 1526, 1192, 567, 126)

BASE_DIMENSION = 5
CYCLIC_LABELS = tuple(range(8))
BASE_LABELS = tuple(range(10))
FIRST_STELLAR_FACE = (0, 2, 4)
FIRST_STELLAR_VERTEX = 8
SECOND_STELLAR_FACE = (0, 1, 8)
SECOND_STELLAR_VERTEX = 9
FIRST_EPSILON = Fraction(1, 8)
SECOND_EPSILON = Fraction(1, 32)
BASE_F_VECTOR = (10, 41, 84, 85, 34)
BASE_MINIMAL_NONFACES = (
    (6, 8),
    (6, 9),
    (7, 8),
    (7, 9),
    (0, 1, 8),
    (0, 2, 4),
    (1, 3, 5),
    (1, 3, 6),
    (1, 4, 6),
    (2, 4, 6),
    (2, 4, 9),
    (3, 5, 9),
    (0, 2, 5, 7),
    (0, 3, 5, 7),
)

# One new clone for each carrier.  The old carrier and its new clone form the
# two-vertex block used in minimal-nonface expansion.
CLONE_BLOCKS: Mapping[int, tuple[int, int]] = {
    1: (1, 10),
    2: (2, 11),
    3: (3, 12),
    4: (4, 13),
}

# Relabel the expanded fourteen-vertex complex onto the frozen record-12
# graph6 labels.  In particular {6,7}|{8,9} becomes {0,1}|{12,13}.
TARGET_RELABEL: Mapping[int, int] = {
    6: 0,
    7: 1,
    0: 2,
    1: 3,
    2: 4,
    3: 5,
    4: 6,
    5: 7,
    10: 8,
    11: 9,
    12: 10,
    13: 11,
    8: 12,
    9: 13,
}


def _moment_point(parameter: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(parameter**power) for power in range(1, 6))


BASE_COORDINATES: Mapping[int, tuple[Fraction, ...]] = {
    **{label: _moment_point(label) for label in CYCLIC_LABELS},
    8: (
        Fraction(29, 16),
        Fraction(85, 16),
        Fraction(59, 4),
        Fraction(463, 16),
        Fraction(-229, 4),
    ),
    9: (
        Fraction(221, 256),
        Fraction(1279, 768),
        Fraction(169, 64),
        Fraction(-4643, 768),
        Fraction(-7671, 64),
    ),
}


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex))


def _determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    rows = [list(map(Fraction, row)) for row in matrix]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("record-12 determinant input must be square")
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next(
            (row for row in range(column, len(rows)) if rows[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            result = -result
        value = rows[column][column]
        result *= value
        for row in range(column + 1, len(rows)):
            scale = rows[row][column] / value
            for later in range(column + 1, len(rows)):
                rows[row][later] -= scale * rows[column][later]
    return result


def supporting_facets(
    labels: Sequence[int], coordinates: Mapping[int, tuple[Fraction, ...]]
) -> tuple[frozenset[int], ...]:
    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, BASE_DIMENSION):
        outside = tuple(label for label in labels if label not in chosen)
        values = tuple(
            _determinant(
                tuple(
                    (Fraction(1), *coordinates[label])
                    for label in (*chosen, outside_label)
                )
            )
            for outside_label in outside
        )
        nonzero = tuple(value for value in values if value)
        if not nonzero or not (
            all(value > 0 for value in nonzero)
            or all(value < 0 for value in nonzero)
        ):
            continue
        supports.add(
            frozenset(chosen)
            | frozenset(
                label for label, value in zip(outside, values) if not value
            )
        )
    maximal = tuple(
        sorted(
            (face for face in supports if not any(face < other for other in supports)),
            key=lambda face: tuple(sorted(face)),
        )
    )
    if any(len(facet) != BASE_DIMENSION for facet in maximal):
        raise AssertionError("record-12 base realization is not simplicial")
    return maximal


def cyclic_facets() -> tuple[frozenset[int], ...]:
    result = []
    for positions in combinations(range(8), 5):
        omitted = tuple(position for position in range(8) if position not in positions)
        if all(
            sum(first < selected < second for selected in positions) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            result.append(frozenset(positions))
    return tuple(sorted(result, key=lambda face: tuple(sorted(face))))


def stellar_subdivision(
    facets: Iterable[frozenset[int]], carrier: Iterable[int], new_vertex: int
) -> tuple[frozenset[int], ...]:
    source = set(facets)
    face = frozenset(carrier)
    star = {facet for facet in source if face <= facet}
    if not star or any(new_vertex in facet for facet in source):
        raise ValueError("invalid record-12 stellar subdivision")
    target = source - star
    target.update(
        (facet - {removed}) | {new_vertex}
        for facet in star
        for removed in face
    )
    return tuple(sorted(target, key=lambda item: tuple(sorted(item))))


def radial_beyond_point(
    coordinates: Mapping[int, tuple[Fraction, ...]],
    labels: Sequence[int],
    carrier: Sequence[int],
    epsilon: Fraction,
) -> tuple[Fraction, ...]:
    barycentre = tuple(
        sum((coordinates[vertex][axis] for vertex in carrier), Fraction(0))
        / len(carrier)
        for axis in range(BASE_DIMENSION)
    )
    interior = tuple(
        sum((coordinates[vertex][axis] for vertex in labels), Fraction(0))
        / len(labels)
        for axis in range(BASE_DIMENSION)
    )
    return tuple(
        barycentre[axis] + epsilon * (barycentre[axis] - interior[axis])
        for axis in range(BASE_DIMENSION)
    )


def _all_faces(facets: Iterable[frozenset[int]]) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        for size in range(1, BASE_DIMENSION + 1):
            result.update(frozenset(value) for value in combinations(facet, size))
    return result


def minimal_nonfaces(
    labels: Sequence[int], facets: Iterable[frozenset[int]]
) -> tuple[frozenset[int], ...]:
    faces = _all_faces(facets)
    result = []
    for size in range(1, len(labels) + 1):
        for chosen in combinations(labels, size):
            candidate = frozenset(chosen)
            if candidate not in faces and all(
                candidate - {vertex} in faces for vertex in candidate
            ):
                result.append(candidate)
    return tuple(sorted(result, key=lambda face: (len(face), tuple(sorted(face)))))


def _face_vector(facets: Iterable[frozenset[int]]) -> tuple[int, ...]:
    faces = _all_faces(facets)
    return tuple(sum(len(face) == size for face in faces) for size in range(1, 6))


def expand_and_relabel(nonface: Iterable[int]) -> frozenset[int]:
    expanded: set[int] = set()
    for vertex in nonface:
        expanded.update(CLONE_BLOCKS.get(vertex, (vertex,)))
    return frozenset(TARGET_RELABEL[vertex] for vertex in expanded)


def target_clutter_masks() -> tuple[int, ...]:
    return tuple(
        mask(value)
        for value in TARGET_MISSING_EDGES
        + TARGET_MISSING_FOURS
        + TARGET_MISSING_FIVES
    )


def target_facet_complements() -> tuple[int, ...]:
    clutter = target_clutter_masks()

    def hits(candidate: int) -> bool:
        return all(candidate & nonface for nonface in clutter)

    return tuple(
        candidate
        for candidate in range(1, 1 << TARGET_VERTEX_COUNT)
        if hits(candidate)
        and all(
            not hits(candidate ^ (1 << vertex)) for vertex in vertices(candidate)
        )
    )


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True)
class Record12PolytopalityAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    cyclic_facet_count: int
    first_stellar_facet_count: int
    second_stellar_facet_count: int
    first_radial_formula_verified: bool
    second_radial_formula_verified: bool
    first_coordinate_support_match: bool
    second_coordinate_support_match: bool
    base_face_vector: tuple[int, ...]
    base_minimal_nonfaces: tuple[frozenset[int], ...]
    expanded_minimal_nonfaces: tuple[frozenset[int], ...]
    target_minimal_nonfaces: tuple[frozenset[int], ...]
    blocker_order_histogram: tuple[tuple[int, int], ...]

    @property
    def verified(self) -> bool:
        return (
            self.decoded_edges == TARGET_MISSING_EDGES
            and (
                self.cyclic_facet_count,
                self.first_stellar_facet_count,
                self.second_stellar_facet_count,
            )
            == (20, 26, 34)
            and self.first_radial_formula_verified
            and self.second_radial_formula_verified
            and self.first_coordinate_support_match
            and self.second_coordinate_support_match
            and self.base_face_vector == BASE_F_VECTOR
            and self.base_minimal_nonfaces
            == tuple(
                sorted(
                    (frozenset(value) for value in BASE_MINIMAL_NONFACES),
                    key=lambda face: (len(face), tuple(sorted(face))),
                )
            )
            and self.expanded_minimal_nonfaces == self.target_minimal_nonfaces
            and self.blocker_order_histogram == ((5, 126),)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "cyclic_order": list(CYCLIC_LABELS),
            "cyclic_facet_count": self.cyclic_facet_count,
            "stellar_subdivisions": [
                {
                    "carrier": list(FIRST_STELLAR_FACE),
                    "new_vertex": FIRST_STELLAR_VERTEX,
                    "epsilon": _fraction_text(FIRST_EPSILON),
                    "facet_count": self.first_stellar_facet_count,
                    "radial_formula_verified": self.first_radial_formula_verified,
                    "coordinate_support_match": self.first_coordinate_support_match,
                },
                {
                    "carrier": list(SECOND_STELLAR_FACE),
                    "new_vertex": SECOND_STELLAR_VERTEX,
                    "epsilon": _fraction_text(SECOND_EPSILON),
                    "facet_count": self.second_stellar_facet_count,
                    "radial_formula_verified": self.second_radial_formula_verified,
                    "coordinate_support_match": self.second_coordinate_support_match,
                },
            ],
            "base_coordinates": {
                str(label): [_fraction_text(value) for value in BASE_COORDINATES[label]]
                for label in BASE_LABELS
            },
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces": [
                list(sorted(face)) for face in self.base_minimal_nonfaces
            ],
            "clone_blocks": {
                str(vertex): list(block) for vertex, block in CLONE_BLOCKS.items()
            },
            "target_relabel": {
                str(source): target for source, target in TARGET_RELABEL.items()
            },
            "expanded_minimal_nonfaces": [
                list(sorted(face)) for face in self.expanded_minimal_nonfaces
            ],
            "target_minimal_nonfaces": [
                list(sorted(face)) for face in self.target_minimal_nonfaces
            ],
            "blocker_order_histogram": {
                str(order): count for order, count in self.blocker_order_histogram
            },
            "simplicial_wedge_count": 4,
            "target_dimension": TARGET_DIMENSION,
            "target_vertex_count": TARGET_VERTEX_COUNT,
            "ordinary_rational_polytopal": self.verified,
        }


def audit_record12_polytopality() -> Record12PolytopalityAudit:
    vertex_count, edges = decode_graph6(TARGET_GRAPH6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise AssertionError("record-12 graph6 vertex count changed")
    cyclic = cyclic_facets()
    first = stellar_subdivision(cyclic, FIRST_STELLAR_FACE, FIRST_STELLAR_VERTEX)
    second = stellar_subdivision(first, SECOND_STELLAR_FACE, SECOND_STELLAR_VERTEX)
    first_coordinates = {label: BASE_COORDINATES[label] for label in range(9)}
    base_nonfaces = minimal_nonfaces(BASE_LABELS, second)
    expanded = tuple(
        sorted(
            (expand_and_relabel(nonface) for nonface in base_nonfaces),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    target = tuple(
        sorted(
            (
                frozenset(value)
                for value in TARGET_MISSING_EDGES
                + TARGET_MISSING_FOURS
                + TARGET_MISSING_FIVES
            ),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    blocker = target_facet_complements()
    audit = Record12PolytopalityAudit(
        decoded_edges=tuple(edges),
        cyclic_facet_count=len(cyclic),
        first_stellar_facet_count=len(first),
        second_stellar_facet_count=len(second),
        first_radial_formula_verified=(
            BASE_COORDINATES[8]
            == radial_beyond_point(
                {label: BASE_COORDINATES[label] for label in CYCLIC_LABELS},
                CYCLIC_LABELS,
                FIRST_STELLAR_FACE,
                FIRST_EPSILON,
            )
        ),
        second_radial_formula_verified=(
            BASE_COORDINATES[9]
            == radial_beyond_point(
                first_coordinates,
                tuple(range(9)),
                SECOND_STELLAR_FACE,
                SECOND_EPSILON,
            )
        ),
        first_coordinate_support_match=(
            set(supporting_facets(tuple(range(9)), first_coordinates)) == set(first)
        ),
        second_coordinate_support_match=(
            set(supporting_facets(BASE_LABELS, BASE_COORDINATES)) == set(second)
        ),
        base_face_vector=_face_vector(second),
        base_minimal_nonfaces=base_nonfaces,
        expanded_minimal_nonfaces=expanded,
        target_minimal_nonfaces=target,
        blocker_order_histogram=tuple(
            sorted(Counter(value.bit_count() for value in blocker).items())
        ),
    )
    if not audit.verified:
        raise AssertionError("record-12 rational polytopality audit failed")
    return audit

