"""Exact rational polytopality proof for the selected H9/14 record 17 nerve.

Four clone pairs collapse the fourteen-vertex nerve to a ten-vertex
simplicial 4-sphere.  The quotient is obtained from the cyclic 5-polytope
``C(8,5)`` by two stellar subdivisions at triangular faces.  Exact rational
coordinates verify all supporting facets.  Four simplicial wedges then
recover the selected record-17 nerve.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Iterable, Mapping, Sequence


TARGET_RECORD_INDEX = 17
TARGET_GRAPH6 = "M???????????o?g?_"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
TARGET_MISSING_EDGES = (
    (0, 12),
    (1, 12),
    (0, 13),
    (2, 13),
    (12, 13),
)
TARGET_MISSING_FIVES = (
    (0, 3, 4, 5, 6),
    (1, 3, 4, 5, 6),
    (2, 3, 4, 5, 6),
    (0, 5, 6, 7, 8),
    (0, 7, 8, 9, 10),
    (1, 2, 3, 4, 11),
    (1, 2, 9, 10, 11),
    (7, 8, 9, 10, 11),
)

# Fresh labels 14,...,17 are the inverse-multiwedge carriers.
CLONE_BLOCKS: Mapping[int, tuple[int, ...]] = {
    14: (3, 4),
    15: (5, 6),
    16: (7, 8),
    17: (9, 10),
}
BASE_VERTEX_LABELS = (0, 1, 2, 11, 12, 13, 14, 15, 16, 17)
BASE_MISSING_EDGES = TARGET_MISSING_EDGES
BASE_MISSING_TRIPLES = (
    (0, 14, 15),
    (1, 14, 15),
    (2, 14, 15),
    (0, 15, 16),
    (0, 16, 17),
    (11, 16, 17),
)
BASE_MISSING_FOURS = (
    (1, 2, 11, 14),
    (1, 2, 11, 17),
)
BASE_F_VECTOR = (10, 40, 80, 80, 32)

# At t=0,...,7 the moment-curve labels occur in this order.
CYCLIC_ORDER = (1, 0, 11, 15, 17, 14, 16, 2)
FIRST_STELLAR_FACE = (1, 14, 15)
FIRST_STELLAR_VERTEX = 13
SECOND_STELLAR_FACE = (2, 14, 15)
SECOND_STELLAR_VERTEX = 12

# The first eight rows are moment-curve points.  The last two are exact
# beyond/beneath placements for the stated stellar subdivisions.
BASE_COORDINATES: Mapping[int, tuple[Fraction, ...]] = {
    0: (Fraction(1), Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    1: (Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0)),
    2: (Fraction(7), Fraction(49), Fraction(343), Fraction(2401), Fraction(16807)),
    11: (Fraction(2), Fraction(4), Fraction(8), Fraction(16), Fraction(32)),
    12: (
        Fraction(39419, 7614),
        Fraction(219847, 7614),
        Fraction(658768, 3807),
        Fraction(8300773, 7614),
        Fraction(27033112, 3807),
    ),
    13: (
        Fraction(249, 94),
        Fraction(1053, 94),
        Fraction(2334, 47),
        Fraction(21423, 94),
        Fraction(50262, 47),
    ),
    14: (Fraction(5), Fraction(25), Fraction(125), Fraction(625), Fraction(3125)),
    15: (Fraction(3), Fraction(9), Fraction(27), Fraction(81), Fraction(243)),
    16: (Fraction(6), Fraction(36), Fraction(216), Fraction(1296), Fraction(7776)),
    17: (Fraction(4), Fraction(16), Fraction(64), Fraction(256), Fraction(1024)),
}


def _determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    rows = [list(map(Fraction, row)) for row in matrix]
    order = len(rows)
    if not order or any(len(row) != order for row in rows):
        raise ValueError("determinant input must be square")
    result = Fraction(1)
    for column in range(order):
        pivot = next(
            (row for row in range(column, order) if rows[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            result = -result
        pivot_value = rows[column][column]
        result *= pivot_value
        for row in range(column + 1, order):
            quotient = rows[row][column] / pivot_value
            for later in range(column + 1, order):
                rows[row][later] -= quotient * rows[column][later]
    return result


def _homogeneous_determinant(
    labels: Sequence[int], coordinates: Mapping[int, tuple[Fraction, ...]]
) -> Fraction:
    return _determinant(
        tuple((Fraction(1),) + coordinates[label] for label in labels)
    )


def supporting_facets(
    labels: Sequence[int], coordinates: Mapping[int, tuple[Fraction, ...]]
) -> tuple[frozenset[int], ...]:
    """Enumerate every supporting hyperplane and reject nonsimplicial ones."""

    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, 5):
        outside = tuple(label for label in labels if label not in chosen)
        values = tuple(
            _homogeneous_determinant((*chosen, label), coordinates)
            for label in outside
        )
        nonzero = tuple(value for value in values if value)
        if not nonzero:
            continue
        if not (all(value > 0 for value in nonzero) or all(value < 0 for value in nonzero)):
            continue
        support = frozenset(chosen) | frozenset(
            label for label, value in zip(outside, values) if not value
        )
        supports.add(support)
    maximal = tuple(
        sorted(
            (support for support in supports if not any(support < other for other in supports)),
            key=lambda face: tuple(sorted(face)),
        )
    )
    if any(len(facet) != 5 for facet in maximal):
        raise AssertionError("record-17 base realization has a nonsimplicial facet")
    return maximal


def cyclic_facets() -> tuple[frozenset[int], ...]:
    """Return the facets of C(8,5) by Gale's evenness criterion."""

    facets = []
    for positions in combinations(range(8), 5):
        omitted = tuple(position for position in range(8) if position not in positions)
        if all(
            sum(first < selected < second for selected in positions) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            facets.append(frozenset(CYCLIC_ORDER[position] for position in positions))
    return tuple(sorted(facets, key=lambda face: tuple(sorted(face))))


def stellar_subdivision(
    facets: Iterable[frozenset[int]], face: Iterable[int], new_vertex: int
) -> tuple[frozenset[int], ...]:
    source = set(facets)
    carrier = frozenset(face)
    star = {facet for facet in source if carrier <= facet}
    if not star or any(new_vertex in facet for facet in source):
        raise ValueError("invalid record-17 stellar subdivision")
    target = source - star
    target.update(
        (facet - {removed}) | {new_vertex}
        for facet in star
        for removed in carrier
    )
    return tuple(sorted(target, key=lambda item: tuple(sorted(item))))


def _all_faces(facets: Iterable[frozenset[int]]) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        ordered = tuple(sorted(facet))
        for size in range(1, len(ordered) + 1):
            result.update(frozenset(face) for face in combinations(ordered, size))
    return result


def _minimal_nonfaces(
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


def _target_base_facets() -> tuple[frozenset[int], ...]:
    nonfaces = tuple(
        frozenset(face)
        for face in BASE_MISSING_EDGES + BASE_MISSING_TRIPLES + BASE_MISSING_FOURS
    )
    return tuple(
        frozenset(chosen)
        for chosen in combinations(BASE_VERTEX_LABELS, 5)
        if not any(nonface <= frozenset(chosen) for nonface in nonfaces)
    )


def _expand_nonface(nonface: frozenset[int]) -> frozenset[int]:
    result: set[int] = set()
    for vertex in nonface:
        result.update(CLONE_BLOCKS.get(vertex, (vertex,)))
    return frozenset(result)


def _fraction_record(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True)
class Record17PolytopalityAudit:
    cyclic_facet_count: int
    first_stellar_facet_count: int
    second_stellar_facet_count: int
    coordinate_support_facet_count: int
    exact_coordinate_support_match: bool
    exact_target_base_match: bool
    base_face_vector: tuple[int, ...]
    base_minimal_nonfaces: tuple[frozenset[int], ...]
    expanded_minimal_nonfaces: tuple[frozenset[int], ...]
    target_minimal_nonfaces: tuple[frozenset[int], ...]
    exact_multiwedge_nonface_match: bool

    @property
    def verified(self) -> bool:
        return (
            (self.cyclic_facet_count, self.first_stellar_facet_count, self.second_stellar_facet_count)
            == (20, 26, 32)
            and self.coordinate_support_facet_count == 32
            and self.exact_coordinate_support_match
            and self.exact_target_base_match
            and self.base_face_vector == BASE_F_VECTOR
            and self.exact_multiwedge_nonface_match
        )

    def to_record(self) -> dict[str, object]:
        return {
            "cyclic_order": list(CYCLIC_ORDER),
            "cyclic_facet_count": self.cyclic_facet_count,
            "stellar_subdivisions": [
                {"new_vertex": FIRST_STELLAR_VERTEX, "face": list(FIRST_STELLAR_FACE), "facet_count": self.first_stellar_facet_count},
                {"new_vertex": SECOND_STELLAR_VERTEX, "face": list(SECOND_STELLAR_FACE), "facet_count": self.second_stellar_facet_count},
            ],
            "base_coordinates": {
                str(label): [_fraction_record(value) for value in BASE_COORDINATES[label]]
                for label in BASE_VERTEX_LABELS
            },
            "coordinate_support_facet_count": self.coordinate_support_facet_count,
            "exact_coordinate_support_match": self.exact_coordinate_support_match,
            "exact_target_base_match": self.exact_target_base_match,
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces": [list(sorted(face)) for face in self.base_minimal_nonfaces],
            "clone_blocks": {str(carrier): list(block) for carrier, block in CLONE_BLOCKS.items()},
            "expanded_minimal_nonfaces": [list(sorted(face)) for face in self.expanded_minimal_nonfaces],
            "target_minimal_nonfaces": [list(sorted(face)) for face in self.target_minimal_nonfaces],
            "exact_multiwedge_nonface_match": self.exact_multiwedge_nonface_match,
            "base_dimension": 5,
            "simplicial_wedge_count": 4,
            "target_dimension": 9,
            "target_vertex_count": 14,
            "ordinary_rational_polytopal": self.verified,
        }


def audit_record17_polytopality() -> Record17PolytopalityAudit:
    cyclic = cyclic_facets()
    first = stellar_subdivision(cyclic, FIRST_STELLAR_FACE, FIRST_STELLAR_VERTEX)
    second = stellar_subdivision(first, SECOND_STELLAR_FACE, SECOND_STELLAR_VERTEX)
    target_base = _target_base_facets()
    geometric = supporting_facets(BASE_VERTEX_LABELS, BASE_COORDINATES)
    base_nonfaces = _minimal_nonfaces(BASE_VERTEX_LABELS, second)
    expected_base_nonfaces = tuple(
        sorted(
            (frozenset(face) for face in BASE_MISSING_EDGES + BASE_MISSING_TRIPLES + BASE_MISSING_FOURS),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    if base_nonfaces != expected_base_nonfaces:
        raise AssertionError("record-17 collapsed minimal nonfaces changed")
    expanded = tuple(
        sorted(
            (_expand_nonface(nonface) for nonface in base_nonfaces),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    target = tuple(
        sorted(
            (frozenset(face) for face in TARGET_MISSING_EDGES + TARGET_MISSING_FIVES),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    audit = Record17PolytopalityAudit(
        cyclic_facet_count=len(cyclic),
        first_stellar_facet_count=len(first),
        second_stellar_facet_count=len(second),
        coordinate_support_facet_count=len(geometric),
        exact_coordinate_support_match=set(geometric) == set(second),
        exact_target_base_match=set(second) == set(target_base),
        base_face_vector=_face_vector(second),
        base_minimal_nonfaces=base_nonfaces,
        expanded_minimal_nonfaces=expanded,
        target_minimal_nonfaces=target,
        exact_multiwedge_nonface_match=expanded == target,
    )
    if not audit.verified:
        raise AssertionError("record-17 rational polytopality audit failed")
    return audit

