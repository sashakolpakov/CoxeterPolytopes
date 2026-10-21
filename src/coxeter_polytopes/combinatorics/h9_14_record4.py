"""Exact rational ordinary-polytopal completion of H9/14 record 4.

The construction takes ``C(8,5)``, stellar-subdivides the edges ``{0,1}``
and ``{0,7}``, and then applies simplicial wedges at vertices 3, 4, 5, 6.
Exact rational coordinates certify the stellar subdivisions.  This proves
ordinary convex polytopality only, not a hyperbolic Coxeter realization.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 4
TARGET_GRAPH6 = "M???????????_?o??"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
TARGET_MISSING_EDGES = ((0, 12), (0, 13), (1, 13))
TARGET_MISSING_FOURS = ((1, 2, 4, 9), (1, 2, 5, 10))
TARGET_MISSING_FIVES = (
    (1, 3, 5, 8, 10),
    (2, 4, 6, 9, 11),
    (2, 4, 7, 9, 12),
    (2, 5, 7, 10, 12),
    (3, 5, 7, 8, 10),
    (3, 5, 8, 10, 13),
    (3, 6, 7, 8, 11),
    (3, 6, 8, 11, 13),
    (4, 6, 7, 9, 11),
    (4, 6, 9, 11, 13),
)
TARGET_F_VECTOR = (14, 88, 330, 823, 1424, 1708, 1364, 657, 146)

BASE_DIMENSION = 5
CYCLIC_LABELS = tuple(range(8))
BASE_LABELS = tuple(range(10))
FIRST_STELLAR_FACE = (0, 1)
FIRST_STELLAR_VERTEX = 8
SECOND_STELLAR_FACE = (0, 7)
SECOND_STELLAR_VERTEX = 9
FIRST_EPSILON = SECOND_EPSILON = Fraction(1, 16)
BASE_F_VECTOR = (10, 42, 88, 90, 36)
BASE_MINIMAL_NONFACES = (
    (0, 1),
    (0, 7),
    (1, 9),
    (1, 3, 5),
    (1, 3, 6),
    (1, 4, 6),
    (2, 4, 6),
    (3, 5, 8),
    (3, 6, 8),
    (4, 6, 8),
    (2, 4, 9),
    (2, 5, 9),
    (3, 5, 9),
    (2, 4, 7, 8),
    (2, 5, 7, 8),
)
CLONE_BLOCKS: Mapping[int, tuple[int, int]] = {
    3: (3, 10),
    4: (4, 11),
    5: (5, 12),
    6: (6, 13),
}
TARGET_RELABEL: Mapping[int, int] = {
    0: 0,
    1: 13,
    2: 2,
    3: 3,
    4: 4,
    5: 5,
    6: 6,
    7: 12,
    8: 7,
    9: 1,
    10: 8,
    11: 9,
    12: 10,
    13: 11,
}


def _moment_point(parameter: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(parameter**power) for power in range(1, 6))


BASE_COORDINATES: Mapping[int, tuple[Fraction, ...]] = {
    **{label: _moment_point(label) for label in CYCLIC_LABELS},
    8: (
        Fraction(5, 16),
        Fraction(-9, 16),
        Fraction(-179, 32),
        Fraction(-36),
        Fraction(-7235, 32),
    ),
    9: (
        Fraction(2705, 768),
        Fraction(57745, 2304),
        Fraction(271585, 1536),
        Fraction(358073, 288),
        Fraction(13407505, 1536),
    ),
}


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex))


def _determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
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
            (face for face in supports if not any(face < other for other in supports)),
            key=lambda face: tuple(sorted(face)),
        )
    )
    if any(len(face) != BASE_DIMENSION for face in maximal):
        raise AssertionError("record-4 rational hull is not simplicial")
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
        raise ValueError("invalid record-4 stellar subdivision")
    return tuple(
        sorted(
            (source - star)
            | {
                (facet - {removed}) | {new_vertex}
                for facet in star
                for removed in face
            },
            key=lambda item: tuple(sorted(item)),
        )
    )


def radial_beyond_point(
    coordinates: Mapping[int, tuple[Fraction, ...]],
    labels: Sequence[int],
    carrier: Sequence[int],
    epsilon: Fraction,
) -> tuple[Fraction, ...]:
    barycentre = tuple(
        sum((coordinates[vertex][axis] for vertex in carrier), Fraction())
        / len(carrier)
        for axis in range(BASE_DIMENSION)
    )
    interior = tuple(
        sum((coordinates[vertex][axis] for vertex in labels), Fraction())
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
        and all(not hits(candidate ^ (1 << vertex)) for vertex in vertices(candidate))
    )


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True)
class Record4PolytopalityAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    stellar_facet_counts: tuple[int, int, int]
    radial_formulas_verified: tuple[bool, bool]
    coordinate_support_matches: tuple[bool, bool]
    base_face_vector: tuple[int, ...]
    base_minimal_nonfaces: tuple[frozenset[int], ...]
    expanded_minimal_nonfaces: tuple[frozenset[int], ...]
    target_minimal_nonfaces: tuple[frozenset[int], ...]
    blocker_order_histogram: tuple[tuple[int, int], ...]

    @property
    def verified(self) -> bool:
        expected_base = tuple(
            sorted(
                (frozenset(value) for value in BASE_MINIMAL_NONFACES),
                key=lambda face: (len(face), tuple(sorted(face))),
            )
        )
        return (
            self.decoded_edges == TARGET_MISSING_EDGES
            and self.stellar_facet_counts == (20, 28, 36)
            and all(self.radial_formulas_verified)
            and all(self.coordinate_support_matches)
            and self.base_face_vector == BASE_F_VECTOR
            and self.base_minimal_nonfaces == expected_base
            and self.expanded_minimal_nonfaces == self.target_minimal_nonfaces
            and self.blocker_order_histogram == ((5, 146),)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "cyclic_order": list(CYCLIC_LABELS),
            "stellar_subdivisions": [
                {
                    "carrier": list(FIRST_STELLAR_FACE),
                    "new_vertex": FIRST_STELLAR_VERTEX,
                    "epsilon": _fraction_text(FIRST_EPSILON),
                    "facet_count": self.stellar_facet_counts[1],
                    "radial_formula_verified": self.radial_formulas_verified[0],
                    "coordinate_support_match": self.coordinate_support_matches[0],
                },
                {
                    "carrier": list(SECOND_STELLAR_FACE),
                    "new_vertex": SECOND_STELLAR_VERTEX,
                    "epsilon": _fraction_text(SECOND_EPSILON),
                    "facet_count": self.stellar_facet_counts[2],
                    "radial_formula_verified": self.radial_formulas_verified[1],
                    "coordinate_support_match": self.coordinate_support_matches[1],
                },
            ],
            "cyclic_facet_count": self.stellar_facet_counts[0],
            "base_coordinates": {
                str(label): [_fraction_text(value) for value in BASE_COORDINATES[label]]
                for label in BASE_LABELS
            },
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces": [list(sorted(face)) for face in self.base_minimal_nonfaces],
            "clone_blocks": {str(vertex): list(block) for vertex, block in CLONE_BLOCKS.items()},
            "target_relabel": {str(source): target for source, target in TARGET_RELABEL.items()},
            "expanded_minimal_nonfaces": [list(sorted(face)) for face in self.expanded_minimal_nonfaces],
            "target_minimal_nonfaces": [list(sorted(face)) for face in self.target_minimal_nonfaces],
            "blocker_order_histogram": {
                str(order): count for order, count in self.blocker_order_histogram
            },
            "simplicial_wedge_count": 4,
            "target_dimension": TARGET_DIMENSION,
            "target_vertex_count": TARGET_VERTEX_COUNT,
            "ordinary_rational_polytopal": self.verified,
        }


def audit_record4_polytopality() -> Record4PolytopalityAudit:
    vertex_count, edges = decode_graph6(TARGET_GRAPH6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise AssertionError("record-4 graph6 order changed")
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
    audit = Record4PolytopalityAudit(
        decoded_edges=tuple(edges),
        stellar_facet_counts=(len(cyclic), len(first), len(second)),
        radial_formulas_verified=(
            BASE_COORDINATES[8]
            == radial_beyond_point(
                {label: BASE_COORDINATES[label] for label in CYCLIC_LABELS},
                CYCLIC_LABELS,
                FIRST_STELLAR_FACE,
                FIRST_EPSILON,
            ),
            BASE_COORDINATES[9]
            == radial_beyond_point(
                first_coordinates,
                tuple(range(9)),
                SECOND_STELLAR_FACE,
                SECOND_EPSILON,
            ),
        ),
        coordinate_support_matches=(
            set(supporting_facets(tuple(range(9)), first_coordinates)) == set(first),
            set(supporting_facets(BASE_LABELS, BASE_COORDINATES)) == set(second),
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
        raise AssertionError("record-4 rational polytopality audit failed")
    return audit

