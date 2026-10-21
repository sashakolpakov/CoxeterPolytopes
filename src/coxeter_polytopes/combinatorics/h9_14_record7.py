"""Exact ordinary-polytopal completion of H9/14 graph-shadow record 7.

The selected nerve is a four-fold simplicial multiwedge of a rational
ten-vertex simplicial 5-polytope.  The base is obtained from ``C(8,5)`` by
stellar subdivision at an edge and then a triangle.  Exact rational
radial-beyond coordinates certify both subdivisions.

This is an ordinary convex-polytopality result.  It supplies no Coxeter
labels, Lorentzian Gram matrix, compactness, or finite hyperbolic volume.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import generated_face_table
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 7
TARGET_GRAPH6 = "M???????????_?w??"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
TARGET_F_VECTOR = (14, 87, 320, 780, 1320, 1554, 1224, 585, 130)

TARGET_MISSING_EDGES = (
    (0, 12),
    (0, 13),
    (1, 13),
    (2, 13),
)
TARGET_MISSING_FIVES = (
    (0, 4, 6, 8, 10),
    (0, 4, 7, 8, 11),
    (0, 5, 7, 9, 11),
    (1, 2, 3, 5, 9),
    (1, 2, 3, 6, 10),
    (1, 3, 5, 9, 12),
    (1, 3, 6, 10, 12),
    (2, 4, 6, 8, 10),
    (2, 4, 7, 8, 11),
    (2, 5, 7, 9, 11),
    (3, 5, 7, 9, 11),
    (4, 6, 8, 10, 12),
)

BASE_DIMENSION = 5
BASE_VERTEX_LABELS = tuple(range(10))
CYCLIC_ORDER = tuple(range(8))
FIRST_STELLAR_FACE = (0, 1)
FIRST_STELLAR_VERTEX = 8
FIRST_RADIAL_EPSILON = Fraction(1, 16)
SECOND_STELLAR_FACE = (0, 3, 5)
SECOND_STELLAR_VERTEX = 9
SECOND_RADIAL_EPSILON = Fraction(1, 128)
BASE_F_VECTOR = (10, 41, 84, 85, 34)
BASE_MINIMAL_NONFACES = (
    (0, 1),
    (1, 9),
    (7, 9),
    (8, 9),
    (0, 3, 5),
    (1, 3, 5),
    (1, 3, 6),
    (1, 4, 6),
    (2, 4, 6),
    (3, 5, 8),
    (3, 6, 8),
    (4, 6, 8),
    (0, 2, 4, 7),
    (0, 2, 5, 7),
    (2, 4, 7, 8),
    (2, 5, 7, 8),
)

# The four size-two blocks implement four simplicial wedges.  The singleton
# blocks incorporate the construction-to-archive relabeling.
BASE_BLOCKS: Mapping[int, tuple[int, ...]] = {
    0: (12,),
    1: (0,),
    2: (3,),
    3: (4, 8),
    4: (5, 9),
    5: (6, 10),
    6: (7, 11),
    7: (1,),
    8: (2,),
    9: (13,),
}


def _moment_point(parameter: int) -> tuple[Fraction, ...]:
    return tuple(
        Fraction(parameter**power) for power in range(1, BASE_DIMENSION + 1)
    )


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
        Fraction(16361, 6144),
        Fraction(208297, 18432),
        Fraction(619153, 12288),
        Fraction(134291, 576),
        Fraction(13596097, 12288),
    ),
}


def mask(items: Iterable[int]) -> int:
    result = 0
    for vertex in items:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-7 target vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    if value < 0 or value & ~FULL_MASK:
        raise ValueError("record-7 mask lies outside [14]")
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex)
    )


def _determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    rows = [list(map(Fraction, row)) for row in matrix]
    order = len(rows)
    if not order or any(len(row) != order for row in rows):
        raise ValueError("record-7 determinant input must be square")
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


def supporting_facets(
    labels: Sequence[int], coordinates: Mapping[int, tuple[Fraction, ...]]
) -> tuple[frozenset[int], ...]:
    """Enumerate every exact supporting facet of the realized base."""

    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, BASE_DIMENSION):
        outside = tuple(label for label in labels if label not in chosen)
        values = tuple(
            _determinant(
                tuple(
                    (Fraction(1),) + coordinates[item]
                    for item in (*chosen, label)
                )
            )
            for label in outside
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
            (
                support
                for support in supports
                if not any(support < other for other in supports)
            ),
            key=lambda face: tuple(sorted(face)),
        )
    )
    if any(len(facet) != BASE_DIMENSION for facet in maximal):
        raise AssertionError("record-7 base realization is not simplicial")
    return maximal


def cyclic_facets() -> tuple[frozenset[int], ...]:
    result = []
    for positions in combinations(range(len(CYCLIC_ORDER)), BASE_DIMENSION):
        omitted = tuple(
            position
            for position in range(len(CYCLIC_ORDER))
            if position not in positions
        )
        if all(
            sum(first < selected < second for selected in positions) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            result.append(
                frozenset(CYCLIC_ORDER[position] for position in positions)
            )
    return tuple(sorted(result, key=lambda face: tuple(sorted(face))))


def stellar_subdivision(
    facets: Iterable[frozenset[int]], face: Iterable[int], new_vertex: int
) -> tuple[frozenset[int], ...]:
    source = set(facets)
    carrier = frozenset(face)
    star = {facet for facet in source if carrier <= facet}
    if not star or any(new_vertex in facet for facet in source):
        raise ValueError("invalid record-7 stellar subdivision")
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


def radial_beyond_point(
    coordinates: Mapping[int, tuple[Fraction, ...]],
    current_labels: Sequence[int],
    face: Sequence[int],
    epsilon: Fraction,
) -> tuple[Fraction, ...]:
    barycentre = tuple(
        sum((coordinates[vertex][axis] for vertex in face), Fraction(0))
        / len(face)
        for axis in range(BASE_DIMENSION)
    )
    interior = tuple(
        sum(
            (coordinates[vertex][axis] for vertex in current_labels),
            Fraction(0),
        )
        / len(current_labels)
        for axis in range(BASE_DIMENSION)
    )
    return tuple(
        barycentre[axis] + epsilon * (barycentre[axis] - interior[axis])
        for axis in range(BASE_DIMENSION)
    )


def _all_faces(facets: Iterable[frozenset[int]]) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        ordered = tuple(sorted(facet))
        for size in range(1, len(ordered) + 1):
            result.update(frozenset(face) for face in combinations(ordered, size))
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
    return tuple(
        sorted(result, key=lambda face: (len(face), tuple(sorted(face))))
    )


def expand_nonface(nonface: Iterable[int]) -> frozenset[int]:
    result: set[int] = set()
    for vertex in nonface:
        result.update(BASE_BLOCKS[vertex])
    return frozenset(result)


def target_clutter_masks() -> tuple[int, ...]:
    return tuple(
        mask(face) for face in TARGET_MISSING_EDGES + TARGET_MISSING_FIVES
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
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


@dataclass(frozen=True)
class Record7PolytopalityAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    cyclic_facet_count: int
    first_stellar_facet_count: int
    second_stellar_facet_count: int
    first_radial_formula_verified: bool
    second_radial_formula_verified: bool
    coordinate_support_facet_count: int
    exact_coordinate_support_match: bool
    base_face_vector: tuple[int, ...]
    base_minimal_nonfaces: tuple[frozenset[int], ...]
    expanded_minimal_nonfaces: tuple[frozenset[int], ...]
    target_minimal_nonfaces: tuple[frozenset[int], ...]
    exact_multiwedge_nonface_match: bool
    blocker_order_histogram: tuple[tuple[int, int], ...]
    target_face_vector: tuple[int, ...]

    @property
    def verified(self) -> bool:
        return (
            self.decoded_edges == TARGET_MISSING_EDGES
            and (
                self.cyclic_facet_count,
                self.first_stellar_facet_count,
                self.second_stellar_facet_count,
            )
            == (20, 28, 34)
            and self.first_radial_formula_verified
            and self.second_radial_formula_verified
            and self.coordinate_support_facet_count == 34
            and self.exact_coordinate_support_match
            and self.base_face_vector == BASE_F_VECTOR
            and self.base_minimal_nonfaces
            == tuple(frozenset(face) for face in BASE_MINIMAL_NONFACES)
            and self.exact_multiwedge_nonface_match
            and self.blocker_order_histogram == ((5, 130),)
            and self.target_face_vector == TARGET_F_VECTOR
        )

    def to_record(self) -> dict[str, object]:
        return {
            "decoded_edges": [list(edge) for edge in self.decoded_edges],
            "cyclic_order": list(CYCLIC_ORDER),
            "cyclic_facet_count": self.cyclic_facet_count,
            "stellar_subdivisions": [
                {
                    "face": list(FIRST_STELLAR_FACE),
                    "new_vertex": FIRST_STELLAR_VERTEX,
                    "radial_beyond_epsilon": _fraction_text(
                        FIRST_RADIAL_EPSILON
                    ),
                    "radial_formula_verified": self.first_radial_formula_verified,
                    "facet_count": self.first_stellar_facet_count,
                },
                {
                    "face": list(SECOND_STELLAR_FACE),
                    "new_vertex": SECOND_STELLAR_VERTEX,
                    "radial_beyond_epsilon": _fraction_text(
                        SECOND_RADIAL_EPSILON
                    ),
                    "radial_formula_verified": self.second_radial_formula_verified,
                    "facet_count": self.second_stellar_facet_count,
                },
            ],
            "base_coordinates": {
                str(label): [_fraction_text(value) for value in row]
                for label, row in sorted(BASE_COORDINATES.items())
            },
            "coordinate_support_facet_count": self.coordinate_support_facet_count,
            "exact_coordinate_support_match": self.exact_coordinate_support_match,
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces": [
                list(sorted(face)) for face in self.base_minimal_nonfaces
            ],
            "base_blocks": {
                str(carrier): list(block)
                for carrier, block in sorted(BASE_BLOCKS.items())
            },
            "expanded_minimal_nonfaces": [
                list(sorted(face)) for face in self.expanded_minimal_nonfaces
            ],
            "target_minimal_nonfaces": [
                list(sorted(face)) for face in self.target_minimal_nonfaces
            ],
            "exact_multiwedge_nonface_match": self.exact_multiwedge_nonface_match,
            "simplicial_wedge_count": sum(
                len(block) - 1 for block in BASE_BLOCKS.values()
            ),
            "blocker_order_histogram": [
                {"order": order, "count": count}
                for order, count in self.blocker_order_histogram
            ],
            "target_face_vector": list(self.target_face_vector),
            "ordinary_rational_polytopal": self.verified,
        }


def audit_record7_polytopality() -> Record7PolytopalityAudit:
    vertex_count, decoded_edges = decode_graph6(TARGET_GRAPH6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise AssertionError("record-7 graph order changed")
    cyclic = cyclic_facets()
    first = stellar_subdivision(cyclic, FIRST_STELLAR_FACE, FIRST_STELLAR_VERTEX)
    second = stellar_subdivision(
        first, SECOND_STELLAR_FACE, SECOND_STELLAR_VERTEX
    )
    first_radial = radial_beyond_point(
        BASE_COORDINATES,
        CYCLIC_ORDER,
        FIRST_STELLAR_FACE,
        FIRST_RADIAL_EPSILON,
    )
    second_radial = radial_beyond_point(
        BASE_COORDINATES,
        (*CYCLIC_ORDER, FIRST_STELLAR_VERTEX),
        SECOND_STELLAR_FACE,
        SECOND_RADIAL_EPSILON,
    )
    geometric = supporting_facets(BASE_VERTEX_LABELS, BASE_COORDINATES)
    base_faces = _all_faces(second)
    base_vector = tuple(
        sum(len(face) == size for face in base_faces)
        for size in range(1, BASE_DIMENSION + 1)
    )
    base_nonfaces = minimal_nonfaces(BASE_VERTEX_LABELS, second)
    expanded = tuple(
        sorted(
            (expand_nonface(nonface) for nonface in base_nonfaces),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    target = tuple(
        sorted(
            (
                frozenset(face)
                for face in TARGET_MISSING_EDGES + TARGET_MISSING_FIVES
            ),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    facets = target_facet_complements()
    blocker_histogram = tuple(sorted(Counter(map(int.bit_count, facets)).items()))
    table = generated_face_table(facets)
    target_vector = tuple(
        sum(
            table[subset]
            for subset in range(1 << TARGET_VERTEX_COUNT)
            if subset.bit_count() == size
        )
        for size in range(1, TARGET_DIMENSION + 1)
    )
    audit = Record7PolytopalityAudit(
        decoded_edges=tuple(decoded_edges),
        cyclic_facet_count=len(cyclic),
        first_stellar_facet_count=len(first),
        second_stellar_facet_count=len(second),
        first_radial_formula_verified=(
            first_radial == BASE_COORDINATES[FIRST_STELLAR_VERTEX]
        ),
        second_radial_formula_verified=(
            second_radial == BASE_COORDINATES[SECOND_STELLAR_VERTEX]
        ),
        coordinate_support_facet_count=len(geometric),
        exact_coordinate_support_match=(geometric == second),
        base_face_vector=base_vector,
        base_minimal_nonfaces=base_nonfaces,
        expanded_minimal_nonfaces=expanded,
        target_minimal_nonfaces=target,
        exact_multiwedge_nonface_match=(expanded == target),
        blocker_order_histogram=blocker_histogram,
        target_face_vector=target_vector,
    )
    if not audit.verified:
        raise AssertionError("record-7 polytopality audit did not reproduce")
    return audit


__all__ = [
    "BASE_BLOCKS",
    "BASE_COORDINATES",
    "BASE_F_VECTOR",
    "BASE_MINIMAL_NONFACES",
    "CYCLIC_ORDER",
    "FIRST_RADIAL_EPSILON",
    "FIRST_STELLAR_FACE",
    "FIRST_STELLAR_VERTEX",
    "SECOND_RADIAL_EPSILON",
    "SECOND_STELLAR_FACE",
    "SECOND_STELLAR_VERTEX",
    "TARGET_F_VECTOR",
    "TARGET_GRAPH6",
    "TARGET_MISSING_EDGES",
    "TARGET_MISSING_FIVES",
    "TARGET_RECORD_INDEX",
    "audit_record7_polytopality",
    "target_facet_complements",
]
