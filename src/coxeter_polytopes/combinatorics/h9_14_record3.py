"""Exact ordinary-polytopal completion of H9/14 graph-shadow record 3.

The nerve is a four-fold simplicial multiwedge of a rational simplicial
5-polytope obtained from ``C(8,5)`` by two edge stellar subdivisions.  Exact
rational supporting-hyperplane determinants certify the base realization.
This supplies ordinary convex polytopality only, not Coxeter geometry.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import generated_face_table
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 3
TARGET_GRAPH6 = "M???????????_?O??"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
TARGET_F_VECTOR = (14, 89, 340, 866, 1528, 1862, 1504, 729, 162)
TARGET_MISSING_EDGES = ((0, 12), (1, 13))
TARGET_MISSING_FOURS = (
    (1, 5, 10, 12),
    (3, 4, 9, 12),
    (3, 4, 9, 13),
    (5, 7, 10, 11),
)
TARGET_MISSING_FIVES = (
    (0, 1, 2, 6, 8),
    (0, 2, 4, 6, 8),
    (0, 3, 4, 6, 9),
    (1, 2, 5, 8, 10),
    (2, 4, 6, 8, 13),
    (2, 6, 7, 8, 11),
    (3, 4, 7, 9, 11),
    (3, 5, 9, 10, 12),
    (3, 5, 9, 10, 13),
)

BASE_DIMENSION = 5
BASE_VERTEX_LABELS = tuple(range(10))
CYCLIC_ORDER = tuple(range(8))
FIRST_STELLAR_FACE = (0, 1)
FIRST_STELLAR_VERTEX = 8
FIRST_RADIAL_EPSILON = Fraction(1, 16)
SECOND_STELLAR_FACE = (4, 8)
SECOND_STELLAR_VERTEX = 9
SECOND_RADIAL_EPSILON = Fraction(1, 64)
BASE_F_VECTOR = (10, 42, 88, 90, 36)
BASE_MINIMAL_NONFACES = (
    (0, 1),
    (4, 8),
    (6, 9),
    (1, 3, 5),
    (1, 3, 6),
    (1, 4, 6),
    (2, 4, 6),
    (2, 7, 9),
    (3, 5, 8),
    (3, 5, 9),
    (3, 6, 8),
    (0, 2, 4, 7),
    (0, 2, 5, 7),
    (0, 3, 5, 7),
    (2, 5, 7, 8),
)

# Four clone pairs, including the construction-to-target relabeling.
BASE_BLOCKS: Mapping[int, tuple[int, ...]] = {
    0: (0,),
    1: (12,),
    2: (2, 8),
    3: (3, 9),
    4: (1,),
    5: (4,),
    6: (5, 10),
    7: (6,),
    8: (13,),
    9: (7, 11),
}


def _moment_point(parameter: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(parameter**power) for power in range(1, 6))


BASE_COORDINATES: Mapping[int, tuple[Fraction, ...]] = {
    **{label: _moment_point(label) for label in CYCLIC_ORDER},
    8: (
        Fraction(5, 16),
        Fraction(-9, 16),
        Fraction(-179, 32),
        Fraction(-36),
        Fraction(-7235, 32),
    ),
    9: (
        Fraction(13153, 6144),
        Fraction(140033, 18432),
        Fraction(347849, 12288),
        Fraction(29855, 288),
        Fraction(4364921, 12288),
    ),
}


def mask(items: Iterable[int]) -> int:
    result = 0
    for vertex in items:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-3 target vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if value & (1 << vertex))


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    rows = [list(map(Fraction, row)) for row in matrix]
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
    for chosen in combinations(labels, 5):
        outside = tuple(label for label in labels if label not in chosen)
        values = tuple(
            determinant(
                tuple(
                    (Fraction(1),) + coordinates[item]
                    for item in (*chosen, label)
                )
            )
            for label in outside
        )
        nonzero = tuple(value for value in values if value)
        if nonzero and (
            all(value > 0 for value in nonzero)
            or all(value < 0 for value in nonzero)
        ):
            supports.add(
                frozenset(chosen)
                | frozenset(
                    label for label, value in zip(outside, values) if not value
                )
            )
    maximal = tuple(
        sorted(
            (
                support
                for support in supports
                if not any(support < other for other in supports)
            ),
            key=lambda face: tuple(sorted(face)),
        )
    )
    if any(len(facet) != 5 for facet in maximal):
        raise AssertionError("record-3 base realization is not simplicial")
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
    facets: Iterable[frozenset[int]], face: Iterable[int], new_vertex: int
) -> tuple[frozenset[int], ...]:
    source = set(facets)
    carrier = frozenset(face)
    star = {facet for facet in source if carrier <= facet}
    if not star:
        raise ValueError("record-3 stellar carrier is absent")
    return tuple(
        sorted(
            (source - star)
            | {
                (facet - {removed}) | {new_vertex}
                for facet in star
                for removed in carrier
            },
            key=lambda item: tuple(sorted(item)),
        )
    )


def radial_beyond(
    coordinates: Mapping[int, tuple[Fraction, ...]],
    labels: Sequence[int],
    face: Sequence[int],
    epsilon: Fraction,
) -> tuple[Fraction, ...]:
    barycentre = tuple(
        sum((coordinates[vertex][axis] for vertex in face), Fraction(0))
        / len(face)
        for axis in range(5)
    )
    interior = tuple(
        sum((coordinates[vertex][axis] for vertex in labels), Fraction(0))
        / len(labels)
        for axis in range(5)
    )
    return tuple(
        barycentre[axis] + epsilon * (barycentre[axis] - interior[axis])
        for axis in range(5)
    )


def all_faces(facets: Iterable[frozenset[int]]) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        for size in range(1, 6):
            result.update(frozenset(face) for face in combinations(facet, size))
    return result


def minimal_nonfaces(
    labels: Sequence[int], facets: Iterable[frozenset[int]]
) -> tuple[frozenset[int], ...]:
    faces = all_faces(facets)
    result = []
    for size in range(1, len(labels) + 1):
        for chosen in combinations(labels, size):
            candidate = frozenset(chosen)
            if candidate not in faces and all(
                candidate - {vertex} in faces for vertex in candidate
            ):
                result.append(candidate)
    return tuple(sorted(result, key=lambda face: (len(face), tuple(sorted(face)))))


def target_clutter_masks() -> tuple[int, ...]:
    return tuple(
        mask(face)
        for face in TARGET_MISSING_EDGES + TARGET_MISSING_FOURS + TARGET_MISSING_FIVES
    )


def target_facet_complements() -> tuple[int, ...]:
    clutter = target_clutter_masks()

    def hits(candidate: int) -> bool:
        return all(candidate & row for row in clutter)

    return tuple(
        candidate
        for candidate in range(1, 1 << 14)
        if hits(candidate)
        and all(not hits(candidate ^ (1 << vertex)) for vertex in vertices(candidate))
    )


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True)
class Record3PolytopalityAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    stellar_facet_counts: tuple[int, int, int]
    radial_formulas_verified: tuple[bool, bool]
    coordinate_support_facet_count: int
    exact_coordinate_support_match: bool
    base_face_vector: tuple[int, ...]
    base_minimal_nonfaces: tuple[frozenset[int], ...]
    expanded_minimal_nonfaces: tuple[frozenset[int], ...]
    target_minimal_nonfaces: tuple[frozenset[int], ...]
    blocker_order_histogram: tuple[tuple[int, int], ...]
    target_face_vector: tuple[int, ...]

    @property
    def verified(self) -> bool:
        return (
            self.decoded_edges == TARGET_MISSING_EDGES
            and self.stellar_facet_counts == (20, 28, 36)
            and self.radial_formulas_verified == (True, True)
            and self.coordinate_support_facet_count == 36
            and self.exact_coordinate_support_match
            and self.base_face_vector == BASE_F_VECTOR
            and self.base_minimal_nonfaces == tuple(frozenset(face) for face in BASE_MINIMAL_NONFACES)
            and self.expanded_minimal_nonfaces == self.target_minimal_nonfaces
            and self.blocker_order_histogram == ((5, 162),)
            and self.target_face_vector == TARGET_F_VECTOR
        )

    def to_record(self) -> dict[str, object]:
        return {
            "decoded_edges": [list(edge) for edge in self.decoded_edges],
            "cyclic_order": list(CYCLIC_ORDER),
            "stellar_facet_counts": list(self.stellar_facet_counts),
            "stellar_subdivisions": [
                {
                    "face": list(FIRST_STELLAR_FACE),
                    "new_vertex": FIRST_STELLAR_VERTEX,
                    "radial_beyond_epsilon": fraction_text(FIRST_RADIAL_EPSILON),
                    "radial_formula_verified": self.radial_formulas_verified[0],
                    "facet_count": self.stellar_facet_counts[1],
                },
                {
                    "face": list(SECOND_STELLAR_FACE),
                    "new_vertex": SECOND_STELLAR_VERTEX,
                    "radial_beyond_epsilon": fraction_text(SECOND_RADIAL_EPSILON),
                    "radial_formula_verified": self.radial_formulas_verified[1],
                    "facet_count": self.stellar_facet_counts[2],
                },
            ],
            "base_coordinates": {
                str(label): [fraction_text(value) for value in row]
                for label, row in sorted(BASE_COORDINATES.items())
            },
            "coordinate_support_facet_count": self.coordinate_support_facet_count,
            "exact_coordinate_support_match": self.exact_coordinate_support_match,
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces": [list(sorted(face)) for face in self.base_minimal_nonfaces],
            "base_blocks": {str(key): list(value) for key, value in sorted(BASE_BLOCKS.items())},
            "expanded_minimal_nonfaces": [list(sorted(face)) for face in self.expanded_minimal_nonfaces],
            "target_minimal_nonfaces": [list(sorted(face)) for face in self.target_minimal_nonfaces],
            "exact_multiwedge_nonface_match": self.expanded_minimal_nonfaces == self.target_minimal_nonfaces,
            "simplicial_wedge_count": 4,
            "blocker_order_histogram": [
                {"order": order, "count": count} for order, count in self.blocker_order_histogram
            ],
            "target_face_vector": list(self.target_face_vector),
            "ordinary_rational_polytopal": self.verified,
        }


def audit_record3_polytopality() -> Record3PolytopalityAudit:
    order, edges = decode_graph6(TARGET_GRAPH6)
    if order != 14:
        raise AssertionError("record-3 graph order changed")
    initial = cyclic_facets()
    first = stellar_subdivision(initial, FIRST_STELLAR_FACE, FIRST_STELLAR_VERTEX)
    second = stellar_subdivision(first, SECOND_STELLAR_FACE, SECOND_STELLAR_VERTEX)
    geometric = supporting_facets(BASE_VERTEX_LABELS, BASE_COORDINATES)
    base_nonfaces = minimal_nonfaces(BASE_VERTEX_LABELS, second)
    expanded = tuple(
        sorted(
            (
                frozenset(vertex for base in face for vertex in BASE_BLOCKS[base])
                for face in base_nonfaces
            ),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    target = tuple(
        sorted(
            (
                frozenset(face)
                for face in TARGET_MISSING_EDGES + TARGET_MISSING_FOURS + TARGET_MISSING_FIVES
            ),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    facets = target_facet_complements()
    table = generated_face_table(facets)
    base_face_set = all_faces(second)
    audit = Record3PolytopalityAudit(
        decoded_edges=tuple(edges),
        stellar_facet_counts=(len(initial), len(first), len(second)),
        radial_formulas_verified=(
            radial_beyond(BASE_COORDINATES, CYCLIC_ORDER, FIRST_STELLAR_FACE, FIRST_RADIAL_EPSILON)
            == BASE_COORDINATES[8],
            radial_beyond(BASE_COORDINATES, (*CYCLIC_ORDER, 8), SECOND_STELLAR_FACE, SECOND_RADIAL_EPSILON)
            == BASE_COORDINATES[9],
        ),
        coordinate_support_facet_count=len(geometric),
        exact_coordinate_support_match=(geometric == second),
        base_face_vector=tuple(sum(len(face) == size for face in base_face_set) for size in range(1, 6)),
        base_minimal_nonfaces=base_nonfaces,
        expanded_minimal_nonfaces=expanded,
        target_minimal_nonfaces=target,
        blocker_order_histogram=tuple(sorted(Counter(map(int.bit_count, facets)).items())),
        target_face_vector=tuple(
            sum(table[value] for value in range(1 << 14) if value.bit_count() == size)
            for size in range(1, 10)
        ),
    )
    if not audit.verified:
        raise AssertionError("record-3 polytopality audit did not reproduce")
    return audit


__all__ = [
    "BASE_F_VECTOR", "TARGET_F_VECTOR", "TARGET_GRAPH6", "TARGET_MISSING_EDGES",
    "TARGET_MISSING_FOURS", "TARGET_MISSING_FIVES", "TARGET_RECORD_INDEX",
    "audit_record3_polytopality", "target_facet_complements",
]
