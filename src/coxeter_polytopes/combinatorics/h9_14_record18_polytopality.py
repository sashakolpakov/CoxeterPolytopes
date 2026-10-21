"""Exact rational polytopality proof for the selected H9/14 record 18 nerve.

Two clone pairs collapse the fourteen-vertex nerve to a twelve-vertex
simplicial 6-sphere.  The quotient is obtained from the cyclic 7-polytope
``C(10,7)`` by two boundary stellar subdivisions.  Exact rational
coordinates verify all supporting facets.  Two simplicial wedges then
recover the selected record-18 nerve.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Iterable, Mapping, Sequence


TARGET_RECORD_INDEX = 18
TARGET_GRAPH6 = "M???????????o?K?_"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
TARGET_MISSING_EDGES = (
    (0, 12),
    (1, 12),
    (2, 13),
    (3, 13),
    (12, 13),
)
TARGET_MISSING_FIVES = (
    (0, 2, 4, 5, 6),
    (0, 5, 6, 7, 8),
    (0, 3, 7, 8, 9),
    (3, 7, 8, 9, 10),
    (2, 4, 5, 6, 11),
    (1, 2, 4, 10, 11),
    (1, 3, 9, 10, 11),
    (1, 4, 9, 10, 11),
)

# A simplicial wedge replaces a carrier by its clone block.  We retain the
# first original label as the carrier in the collapsed complex.
CLONE_BLOCKS: Mapping[int, tuple[int, ...]] = {
    5: (5, 6),
    7: (7, 8),
}
BASE_VERTEX_LABELS = (0, 1, 2, 3, 4, 5, 7, 9, 10, 11, 12, 13)
BASE_DIMENSION = 7
BASE_F_VECTOR = (12, 61, 175, 310, 341, 217, 62)
BASE_MISSING_EDGES = TARGET_MISSING_EDGES
BASE_MISSING_TRIPLES = ((0, 5, 7),)
BASE_MISSING_FOURS = (
    (0, 2, 4, 5),
    (0, 3, 7, 9),
    (3, 7, 9, 10),
    (2, 4, 5, 11),
)
BASE_MISSING_FIVES = (
    (1, 2, 4, 10, 11),
    (1, 3, 9, 10, 11),
    (1, 4, 9, 10, 11),
)

# At t=0,...,9 the moment-curve labels occur in this order.
CYCLIC_ORDER = (1, 0, 10, 5, 9, 2, 3, 4, 7, 11)
FIRST_STELLAR_FACE = (2, 4, 5, 11)
FIRST_STELLAR_VERTEX = 12
SECOND_STELLAR_FACE = (0, 5, 7)
SECOND_STELLAR_VERTEX = 13
RADIAL_BEYOND_EPSILON = Fraction(1, 32)

# The first ten rows are moment-curve points (t,...,t^7).  Vertices 12 and
# 13 are exact radial-beyond placements.  If y is the barycentre of the face
# being subdivided and p is the barycentre of all current vertices, the new
# point is y + (y-p)/32.  This lies beyond precisely the facets containing
# that face and beneath all other facets, as the exhaustive determinant audit
# independently confirms.
BASE_COORDINATES: Mapping[int, tuple[Fraction, ...]] = {
    0: (Fraction(1), Fraction(1), Fraction(1), Fraction(1), Fraction(1), Fraction(1), Fraction(1)),
    1: (Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0), Fraction(0)),
    2: (Fraction(5), Fraction(25), Fraction(125), Fraction(625), Fraction(3125), Fraction(15625), Fraction(78125)),
    3: (Fraction(6), Fraction(36), Fraction(216), Fraction(1296), Fraction(7776), Fraction(46656), Fraction(279936)),
    4: (Fraction(7), Fraction(49), Fraction(343), Fraction(2401), Fraction(16807), Fraction(117649), Fraction(823543)),
    5: (Fraction(3), Fraction(9), Fraction(27), Fraction(81), Fraction(243), Fraction(729), Fraction(2187)),
    7: (Fraction(8), Fraction(64), Fraction(512), Fraction(4096), Fraction(32768), Fraction(262144), Fraction(2097152)),
    9: (Fraction(4), Fraction(16), Fraction(64), Fraction(256), Fraction(1024), Fraction(4096), Fraction(16384)),
    10: (Fraction(2), Fraction(4), Fraction(8), Fraction(16), Fraction(32), Fraction(64), Fraction(128)),
    11: (Fraction(9), Fraction(81), Fraction(729), Fraction(6561), Fraction(59049), Fraction(531441), Fraction(4782969)),
    12: (
        Fraction(387, 64),
        Fraction(2649, 64),
        Fraction(19791, 64),
        Fraction(782277, 320),
        Fraction(1283031, 64),
        Fraction(10784145, 64),
        Fraction(92216511, 64),
    ),
    13: (
        Fraction(8151, 2048),
        Fraction(50197, 2048),
        Fraction(366579, 2048),
        Fraction(14189393, 10240),
        Fraction(22420827, 2048),
        Fraction(178390381, 2048),
        Fraction(1422538659, 2048),
    ),
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
    """Enumerate every exact supporting hyperplane of a 7-polytope."""

    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, BASE_DIMENSION):
        outside = tuple(label for label in labels if label not in chosen)
        values = tuple(
            _homogeneous_determinant((*chosen, label), coordinates)
            for label in outside
        )
        nonzero = tuple(value for value in values if value)
        if not nonzero:
            continue
        if not (
            all(value > 0 for value in nonzero)
            or all(value < 0 for value in nonzero)
        ):
            continue
        support = frozenset(chosen) | frozenset(
            label for label, value in zip(outside, values) if not value
        )
        supports.add(support)
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
        raise AssertionError("record-18 base realization has a nonsimplicial facet")
    return maximal


def cyclic_facets() -> tuple[frozenset[int], ...]:
    """Return the facets of ``C(10,7)`` by Gale's evenness criterion."""

    facets = []
    for positions in combinations(range(10), BASE_DIMENSION):
        omitted = tuple(position for position in range(10) if position not in positions)
        if all(
            sum(first < selected < second for selected in positions) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            facets.append(
                frozenset(CYCLIC_ORDER[position] for position in positions)
            )
    return tuple(sorted(facets, key=lambda face: tuple(sorted(face))))


def stellar_subdivision(
    facets: Iterable[frozenset[int]], face: Iterable[int], new_vertex: int
) -> tuple[frozenset[int], ...]:
    source = set(facets)
    carrier = frozenset(face)
    star = {facet for facet in source if carrier <= facet}
    if not star or any(new_vertex in facet for facet in source):
        raise ValueError("invalid record-18 stellar subdivision")
    target = source - star
    target.update(
        (facet - {removed}) | {new_vertex}
        for facet in star
        for removed in carrier
    )
    return tuple(sorted(target, key=lambda item: tuple(sorted(item))))


def _radial_beyond_point(
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
    return tuple(
        sum(len(face) == size for face in faces)
        for size in range(1, BASE_DIMENSION + 1)
    )


def _target_base_facets() -> tuple[frozenset[int], ...]:
    nonfaces = tuple(
        frozenset(face)
        for face in (
            BASE_MISSING_EDGES
            + BASE_MISSING_TRIPLES
            + BASE_MISSING_FOURS
            + BASE_MISSING_FIVES
        )
    )
    return tuple(
        frozenset(chosen)
        for chosen in combinations(BASE_VERTEX_LABELS, BASE_DIMENSION)
        if not any(nonface <= frozenset(chosen) for nonface in nonfaces)
    )


def _expand_nonface(nonface: frozenset[int]) -> frozenset[int]:
    result: set[int] = set()
    for vertex in nonface:
        result.update(CLONE_BLOCKS.get(vertex, (vertex,)))
    return frozenset(result)


def _fraction_record(value: Fraction) -> str:
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"{value.numerator}/{value.denominator}"
    )


@dataclass(frozen=True)
class Record18PolytopalityAudit:
    cyclic_facet_count: int
    first_stellar_facet_count: int
    second_stellar_facet_count: int
    first_radial_formula_verified: bool
    second_radial_formula_verified: bool
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
            (
                self.cyclic_facet_count,
                self.first_stellar_facet_count,
                self.second_stellar_facet_count,
            )
            == (40, 52, 62)
            and self.first_radial_formula_verified
            and self.second_radial_formula_verified
            and self.coordinate_support_facet_count == 62
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
                {
                    "new_vertex": FIRST_STELLAR_VERTEX,
                    "face": list(FIRST_STELLAR_FACE),
                    "facet_count": self.first_stellar_facet_count,
                    "radial_beyond_epsilon": _fraction_record(
                        RADIAL_BEYOND_EPSILON
                    ),
                    "radial_formula_verified": self.first_radial_formula_verified,
                },
                {
                    "new_vertex": SECOND_STELLAR_VERTEX,
                    "face": list(SECOND_STELLAR_FACE),
                    "facet_count": self.second_stellar_facet_count,
                    "radial_beyond_epsilon": _fraction_record(
                        RADIAL_BEYOND_EPSILON
                    ),
                    "radial_formula_verified": self.second_radial_formula_verified,
                },
            ],
            "base_coordinates": {
                str(label): [
                    _fraction_record(value) for value in BASE_COORDINATES[label]
                ]
                for label in BASE_VERTEX_LABELS
            },
            "coordinate_support_facet_count": self.coordinate_support_facet_count,
            "exact_coordinate_support_match": self.exact_coordinate_support_match,
            "exact_target_base_match": self.exact_target_base_match,
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces": [
                list(sorted(face)) for face in self.base_minimal_nonfaces
            ],
            "clone_blocks": {
                str(carrier): list(block)
                for carrier, block in CLONE_BLOCKS.items()
            },
            "expanded_minimal_nonfaces": [
                list(sorted(face)) for face in self.expanded_minimal_nonfaces
            ],
            "target_minimal_nonfaces": [
                list(sorted(face)) for face in self.target_minimal_nonfaces
            ],
            "exact_multiwedge_nonface_match": (
                self.exact_multiwedge_nonface_match
            ),
            "base_dimension": BASE_DIMENSION,
            "simplicial_wedge_count": 2,
            "target_dimension": TARGET_DIMENSION,
            "target_vertex_count": TARGET_VERTEX_COUNT,
            "ordinary_rational_polytopal": self.verified,
        }


def audit_record18_polytopality() -> Record18PolytopalityAudit:
    cyclic = cyclic_facets()
    first = stellar_subdivision(
        cyclic, FIRST_STELLAR_FACE, FIRST_STELLAR_VERTEX
    )
    second = stellar_subdivision(
        first, SECOND_STELLAR_FACE, SECOND_STELLAR_VERTEX
    )
    target_base = _target_base_facets()
    geometric = supporting_facets(BASE_VERTEX_LABELS, BASE_COORDINATES)

    initial_labels = CYCLIC_ORDER
    first_coordinates = {
        label: BASE_COORDINATES[label]
        for label in (*initial_labels, FIRST_STELLAR_VERTEX)
    }
    first_radial = _radial_beyond_point(
        BASE_COORDINATES, initial_labels, FIRST_STELLAR_FACE
    )
    second_radial = _radial_beyond_point(
        first_coordinates,
        (*initial_labels, FIRST_STELLAR_VERTEX),
        SECOND_STELLAR_FACE,
    )

    base_nonfaces = _minimal_nonfaces(BASE_VERTEX_LABELS, second)
    expected_base_nonfaces = tuple(
        sorted(
            (
                frozenset(face)
                for face in (
                    BASE_MISSING_EDGES
                    + BASE_MISSING_TRIPLES
                    + BASE_MISSING_FOURS
                    + BASE_MISSING_FIVES
                )
            ),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    if base_nonfaces != expected_base_nonfaces:
        raise AssertionError("record-18 collapsed minimal nonfaces changed")
    expanded = tuple(
        sorted(
            (_expand_nonface(nonface) for nonface in base_nonfaces),
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
    audit = Record18PolytopalityAudit(
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
        exact_coordinate_support_match=set(geometric) == set(second),
        exact_target_base_match=set(second) == set(target_base),
        base_face_vector=_face_vector(second),
        base_minimal_nonfaces=base_nonfaces,
        expanded_minimal_nonfaces=expanded,
        target_minimal_nonfaces=target,
        exact_multiwedge_nonface_match=expanded == target,
    )
    if not audit.verified:
        raise AssertionError("record-18 rational polytopality audit failed")
    return audit


__all__ = [
    "BASE_COORDINATES",
    "BASE_F_VECTOR",
    "BASE_VERTEX_LABELS",
    "CLONE_BLOCKS",
    "CYCLIC_ORDER",
    "Record18PolytopalityAudit",
    "audit_record18_polytopality",
    "cyclic_facets",
    "stellar_subdivision",
    "supporting_facets",
]
