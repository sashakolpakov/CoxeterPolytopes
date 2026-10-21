"""Exact ordinary-polytopal completion of H9/14 graph-shadow record 31.

The selected nerve is obtained from the cyclic simplicial 5-polytope
``C(8,5)`` by stellar subdivision at two triangular boundary faces and then
by four simplicial wedges.  Exact rational coordinates certify the two
stellar subdivisions.  Expanding the base minimal nonfaces through the four
clone blocks gives exactly the record-31 missing-edge path and nine missing
five-faces.

This is an ordinary convex-polytopality result.  It does not supply Coxeter
labels, a Lorentzian Gram matrix, compactness, or finite hyperbolic volume.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 31
TARGET_GRAPH6 = "M?????????O?O?o??"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
TARGET_MISSING_EDGES = (
    (0, 11),
    (1, 12),
    (0, 13),
    (1, 13),
)
TARGET_MISSING_FIVES = (
    (0, 5, 6, 9, 10),
    (2, 3, 4, 7, 8),
    (2, 5, 6, 9, 10),
    (2, 5, 6, 11, 12),
    (2, 7, 8, 11, 12),
    (3, 4, 7, 8, 11),
    (3, 4, 7, 8, 13),
    (3, 4, 9, 10, 13),
    (5, 6, 9, 10, 13),
)
TARGET_F_VECTOR = (14, 87, 319, 771, 1287, 1491, 1158, 549, 122)

BASE_DIMENSION = 5
BASE_VERTEX_LABELS = (0, 1, 2, 11, 12, 13, 14, 15, 16, 17)
BASE_F_VECTOR = (10, 41, 84, 85, 34)
CYCLIC_ORDER = (2, 14, 15, 16, 17, 11, 13, 12)
FIRST_STELLAR_FACE = (2, 15, 17)
FIRST_STELLAR_VERTEX = 1
SECOND_STELLAR_FACE = (2, 14, 16)
SECOND_STELLAR_VERTEX = 0
RADIAL_BEYOND_EPSILON = Fraction(1, 32)
CLONE_BLOCKS: Mapping[int, tuple[int, ...]] = {
    14: (3, 4),
    15: (5, 6),
    16: (7, 8),
    17: (9, 10),
}
BASE_MINIMAL_NONFACES = (
    (0, 11),
    (0, 13),
    (1, 12),
    (1, 13),
    (0, 15, 17),
    (2, 14, 16),
    (2, 15, 17),
    (11, 14, 16),
    (13, 14, 16),
    (13, 14, 17),
    (13, 15, 17),
    (2, 11, 12, 15),
    (2, 11, 12, 16),
)


def _moment_point(parameter: int) -> tuple[Fraction, ...]:
    return tuple(
        Fraction(parameter**power) for power in range(1, BASE_DIMENSION + 1)
    )


BASE_COORDINATES: Mapping[int, tuple[Fraction, ...]] = {
    **{
        label: _moment_point(parameter)
        for parameter, label in enumerate(CYCLIC_ORDER)
    },
    1: (
        Fraction(125, 64),
        Fraction(405, 64),
        Fraction(347, 16),
        Fraction(4815, 64),
        Fraction(3995, 16),
    ),
    0: (
        Fraction(2603, 2048),
        Fraction(53995, 18432),
        Fraction(10487, 1536),
        Fraction(215473, 18432),
        Fraction(-27209, 1536),
    ),
}


def mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-31 target vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex)
    )


def _determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    rows = [list(map(Fraction, row)) for row in matrix]
    order = len(rows)
    if not order or any(len(row) != order for row in rows):
        raise ValueError("record-31 determinant input must be square")
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
    """Enumerate all exact supporting facets of the realized base."""

    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, BASE_DIMENSION):
        outside = tuple(label for label in labels if label not in chosen)
        values = tuple(
            _homogeneous_determinant((*chosen, label), coordinates)
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
        raise AssertionError("record-31 base realization is not simplicial")
    return maximal


def cyclic_facets() -> tuple[frozenset[int], ...]:
    """Return the twenty facets of ``C(8,5)`` by Gale evenness."""

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
        raise ValueError("invalid record-31 stellar subdivision")
    target = source - star
    target.update(
        (facet - {removed}) | {new_vertex}
        for facet in star
        for removed in carrier
    )
    return tuple(sorted(target, key=lambda item: tuple(sorted(item))))


def radial_beyond_point(
    coordinates: Mapping[int, tuple[Fraction, ...]],
    current_labels: Sequence[int],
    face: Sequence[int],
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
        barycentre[axis]
        + RADIAL_BEYOND_EPSILON * (barycentre[axis] - interior[axis])
        for axis in range(BASE_DIMENSION)
    )


def _all_faces(
    facets: Iterable[frozenset[int]],
) -> set[frozenset[int]]:
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
        result.update(CLONE_BLOCKS.get(vertex, (vertex,)))
    return frozenset(result)


def target_clutter_masks() -> tuple[int, ...]:
    return tuple(mask(face) for face in TARGET_MISSING_EDGES + TARGET_MISSING_FIVES)


def target_facet_complements() -> tuple[int, ...]:
    """Enumerate the complete blocker of the target minimal nonfaces."""

    clutter = target_clutter_masks()

    def hits(candidate: int) -> bool:
        return all(candidate & nonface for nonface in clutter)

    result = []
    for candidate in range(1, 1 << TARGET_VERTEX_COUNT):
        if hits(candidate) and all(
            not hits(candidate ^ (1 << vertex)) for vertex in vertices(candidate)
        ):
            result.append(candidate)
    return tuple(result)


def _fraction_text(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


@dataclass(frozen=True)
class Record31PolytopalityAudit:
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

    @property
    def verified(self) -> bool:
        expected_base = tuple(
            sorted(
                (frozenset(face) for face in BASE_MINIMAL_NONFACES),
                key=lambda face: (len(face), tuple(sorted(face))),
            )
        )
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
            and self.coordinate_support_facet_count == 34
            and self.exact_coordinate_support_match
            and self.base_face_vector == BASE_F_VECTOR
            and self.base_minimal_nonfaces == expected_base
            and self.exact_multiwedge_nonface_match
            and self.blocker_order_histogram == ((5, 122),)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "decoded_missing_edges": [list(edge) for edge in self.decoded_edges],
            "cyclic_order": list(CYCLIC_ORDER),
            "cyclic_facet_count": self.cyclic_facet_count,
            "stellar_subdivisions": [
                {
                    "new_vertex": FIRST_STELLAR_VERTEX,
                    "face": list(FIRST_STELLAR_FACE),
                    "radial_beyond_epsilon": str(RADIAL_BEYOND_EPSILON),
                    "facet_count": self.first_stellar_facet_count,
                },
                {
                    "new_vertex": SECOND_STELLAR_VERTEX,
                    "face": list(SECOND_STELLAR_FACE),
                    "radial_beyond_epsilon": str(RADIAL_BEYOND_EPSILON),
                    "facet_count": self.second_stellar_facet_count,
                },
            ],
            "base_coordinates": {
                str(label): [_fraction_text(value) for value in BASE_COORDINATES[label]]
                for label in BASE_VERTEX_LABELS
            },
            "first_radial_formula_verified": self.first_radial_formula_verified,
            "second_radial_formula_verified": self.second_radial_formula_verified,
            "coordinate_support_facet_count": self.coordinate_support_facet_count,
            "exact_coordinate_support_match": self.exact_coordinate_support_match,
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces": [
                list(sorted(face)) for face in self.base_minimal_nonfaces
            ],
            "clone_blocks": {
                str(carrier): list(block) for carrier, block in CLONE_BLOCKS.items()
            },
            "expanded_minimal_nonfaces": [
                list(sorted(face)) for face in self.expanded_minimal_nonfaces
            ],
            "target_minimal_nonfaces": [
                list(sorted(face)) for face in self.target_minimal_nonfaces
            ],
            "exact_multiwedge_nonface_match": self.exact_multiwedge_nonface_match,
            "blocker_order_histogram": {
                str(order): count for order, count in self.blocker_order_histogram
            },
            "base_dimension": BASE_DIMENSION,
            "simplicial_wedge_count": len(CLONE_BLOCKS),
            "target_dimension": TARGET_DIMENSION,
            "target_vertex_count": TARGET_VERTEX_COUNT,
            "ordinary_rational_polytopal": self.verified,
        }


def audit_record31_polytopality() -> Record31PolytopalityAudit:
    order, decoded = decode_graph6(TARGET_GRAPH6)
    if order != TARGET_VERTEX_COUNT:
        raise AssertionError("record-31 graph order changed")

    cyclic = cyclic_facets()
    first = stellar_subdivision(
        cyclic, FIRST_STELLAR_FACE, FIRST_STELLAR_VERTEX
    )
    second = stellar_subdivision(
        first, SECOND_STELLAR_FACE, SECOND_STELLAR_VERTEX
    )
    geometric = supporting_facets(BASE_VERTEX_LABELS, BASE_COORDINATES)
    faces = _all_faces(second)
    base_vector = tuple(
        sum(len(face) == size for face in faces)
        for size in range(1, BASE_DIMENSION + 1)
    )
    base_nonfaces = minimal_nonfaces(BASE_VERTEX_LABELS, second)
    expanded = tuple(
        sorted(
            (expand_nonface(face) for face in base_nonfaces),
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
    blockers = target_facet_complements()
    audit = Record31PolytopalityAudit(
        decoded_edges=decoded,
        cyclic_facet_count=len(cyclic),
        first_stellar_facet_count=len(first),
        second_stellar_facet_count=len(second),
        first_radial_formula_verified=(
            radial_beyond_point(
                BASE_COORDINATES, CYCLIC_ORDER, FIRST_STELLAR_FACE
            )
            == BASE_COORDINATES[FIRST_STELLAR_VERTEX]
        ),
        second_radial_formula_verified=(
            radial_beyond_point(
                BASE_COORDINATES,
                (*CYCLIC_ORDER, FIRST_STELLAR_VERTEX),
                SECOND_STELLAR_FACE,
            )
            == BASE_COORDINATES[SECOND_STELLAR_VERTEX]
        ),
        coordinate_support_facet_count=len(geometric),
        exact_coordinate_support_match=(geometric == second),
        base_face_vector=base_vector,
        base_minimal_nonfaces=base_nonfaces,
        expanded_minimal_nonfaces=expanded,
        target_minimal_nonfaces=target,
        exact_multiwedge_nonface_match=(expanded == target),
        blocker_order_histogram=tuple(
            sorted(Counter(value.bit_count() for value in blockers).items())
        ),
    )
    if not audit.verified:
        raise AssertionError(f"record-31 polytopality audit failed: {audit!r}")
    return audit

