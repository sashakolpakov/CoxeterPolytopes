"""Exact incidence and rational-polytopality certificate for H9/14 record 35.

The fourteen-vertex nerve has two clone blocks.  Collapsing them gives a
nine-vertex simplicial 3-sphere.  An explicit integral realization of that
base as the boundary of a simplicial 4-polytope is checked by exact
determinants.  Simplicial multiwedges of multiplicities three and four then
recover the original nerve and prove ordinary convex polytopality.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from typing import Iterable, Sequence

from coxeter_polytopes.combinatorics.facet_completion import generated_face_table
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 35
TARGET_GRAPH6 = "M?????????O?O?w??"
TARGET_VERTEX_COUNT = 14
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
TARGET_MISSING_EDGES = (
    (0, 11),
    (1, 12),
    (0, 13),
    (1, 13),
    (2, 13),
)
TARGET_MISSING_FIVES = (
    (2, 3, 4, 6, 8),
    (0, 5, 7, 9, 10),
    (1, 5, 7, 9, 10),
    (4, 5, 7, 9, 10),
    (3, 6, 8, 11, 12),
)
TARGET_F_VECTOR = (14, 86, 309, 726, 1171, 1309, 986, 459, 102)

# The two blocks are the fibres of the inverse simplicial multiwedge.  The
# fresh labels 14 and 15 name their collapsed carrier vertices.
FIRST_CLONE_BLOCK = (3, 6, 8)
SECOND_CLONE_BLOCK = (5, 7, 9, 10)
FIRST_CARRIER = 14
SECOND_CARRIER = 15
BASE_VERTEX_LABELS = (0, 1, 2, 4, 11, 12, 13, FIRST_CARRIER, SECOND_CARRIER)
BASE_MISSING_EDGES = (
    (0, 11),
    (1, 12),
    (0, 13),
    (1, 13),
    (2, 13),
    (0, SECOND_CARRIER),
    (1, SECOND_CARRIER),
    (4, SECOND_CARRIER),
)
BASE_MISSING_TRIPLES = (
    (2, 4, FIRST_CARRIER),
    (11, 12, FIRST_CARRIER),
)

# Coordinates are in BASE_VERTEX_LABELS order.  They realize the base as a
# simplicial integral 4-polytope.  Translation is irrelevant to the support
# determinant audit.
BASE_COORDINATES = (
    (0, 3, 0, -5),
    (2, 5, -1, -1),
    (2, 4, 3, -3),
    (0, 5, -4, 0),
    (5, -2, 0, 3),
    (-4, -4, -2, 1),
    (3, -5, -1, 0),
    (2, -1, 0, -6),
    (3, -5, 1, 1),
)


def mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-35 vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    if value < 0 or value & ~FULL_MASK:
        raise ValueError("record-35 mask lies outside [14]")
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex)
    )


def record35_clutter() -> tuple[int, ...]:
    return tuple(mask(edge) for edge in TARGET_MISSING_EDGES) + tuple(
        mask(face) for face in TARGET_MISSING_FIVES
    )


def _is_transversal(candidate: int, clutter: tuple[int, ...]) -> bool:
    return all(candidate & nonface for nonface in clutter)


def _is_minimal_transversal(candidate: int, clutter: tuple[int, ...]) -> bool:
    if not _is_transversal(candidate, clutter):
        return False
    remaining = candidate
    while remaining:
        bit = remaining & -remaining
        if _is_transversal(candidate ^ bit, clutter):
            return False
        remaining ^= bit
    return True


def record35_minimal_transversals() -> tuple[int, ...]:
    clutter = record35_clutter()
    return tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if _is_minimal_transversal(candidate, clutter)
    )


def _packed_face_table(face_table: tuple[bool, ...]) -> bytes:
    if len(face_table) != 1 << TARGET_VERTEX_COUNT:
        raise ValueError("record-35 face table must have 2^14 entries")
    return bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )


@dataclass(frozen=True)
class Record35IncidenceAudit:
    minimal_transversal_size_histogram: tuple[tuple[int, int], ...]
    facet_complement_masks: tuple[int, ...]
    facet_complement_sha256: str
    face_table_true_count: int
    packed_face_table_sha256: str

    @property
    def verified(self) -> bool:
        return (
            self.minimal_transversal_size_histogram == ((5, 102),)
            and len(self.facet_complement_masks) == TARGET_F_VECTOR[-1]
        )

    def to_record(self) -> dict[str, object]:
        return {
            "missing_edges": [list(edge) for edge in TARGET_MISSING_EDGES],
            "missing_fives": [list(face) for face in TARGET_MISSING_FIVES],
            "minimal_transversal_size_histogram": {
                str(size): count
                for size, count in self.minimal_transversal_size_histogram
            },
            "facet_complement_masks": list(self.facet_complement_masks),
            "facet_complement_sha256": self.facet_complement_sha256,
            "face_table_true_count": self.face_table_true_count,
            "packed_face_table_sha256": self.packed_face_table_sha256,
        }


def audit_record35_incidence(
    graph6: str = TARGET_GRAPH6,
) -> Record35IncidenceAudit:
    order, edges = decode_graph6(graph6)
    if order != TARGET_VERTEX_COUNT or edges != TARGET_MISSING_EDGES:
        raise ValueError("record-35 missing-edge graph changed")
    transversals = tuple(sorted(record35_minimal_transversals()))
    histogram = tuple(sorted(Counter(item.bit_count() for item in transversals).items()))
    payload = b"".join(item.to_bytes(2, "little") for item in transversals)
    face_table = generated_face_table(transversals)
    audit = Record35IncidenceAudit(
        minimal_transversal_size_histogram=histogram,
        facet_complement_masks=transversals,
        facet_complement_sha256=sha256(payload).hexdigest(),
        face_table_true_count=sum(face_table),
        packed_face_table_sha256=sha256(_packed_face_table(face_table)).hexdigest(),
    )
    if not audit.verified:
        raise AssertionError("record-35 blocker is not pure with 102 members")
    return audit


def _determinant(matrix: Sequence[Sequence[int]]) -> int:
    """Return an exact integer determinant by fraction-free elimination."""

    rows = [list(map(int, row)) for row in matrix]
    order = len(rows)
    if order == 0 or any(len(row) != order for row in rows):
        raise ValueError("determinant input must be a nonempty square matrix")
    sign = 1
    denominator = 1
    for column in range(order - 1):
        pivot = next(
            (row for row in range(column, order) if rows[row][column]), None
        )
        if pivot is None:
            return 0
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            sign *= -1
        pivot_value = rows[column][column]
        for row in range(column + 1, order):
            for later in range(column + 1, order):
                numerator = (
                    rows[row][later] * pivot_value
                    - rows[row][column] * rows[column][later]
                )
                if numerator % denominator:
                    raise ArithmeticError("Bareiss division ceased to be exact")
                rows[row][later] = numerator // denominator
            rows[row][column] = 0
        denominator = pivot_value
    return sign * rows[-1][-1]


def _all_faces(facets: Iterable[frozenset[int]]) -> set[frozenset[int]]:
    faces = {frozenset()}
    for facet in facets:
        ordered = tuple(sorted(facet))
        for size in range(1, len(ordered) + 1):
            faces.update(frozenset(face) for face in combinations(ordered, size))
    return faces


def _minimal_nonfaces(
    labels: tuple[int, ...], facets: Iterable[frozenset[int]]
) -> tuple[frozenset[int], ...]:
    face_table = _all_faces(facets)
    result = []
    for size in range(1, len(labels) + 1):
        for subset in combinations(labels, size):
            candidate = frozenset(subset)
            if candidate in face_table:
                continue
            if all(candidate - {vertex} in face_table for vertex in candidate):
                result.append(candidate)
    return tuple(
        sorted(result, key=lambda item: (len(item), tuple(sorted(item))))
    )


def _target_base_facets() -> tuple[frozenset[int], ...]:
    nonfaces = tuple(frozenset(edge) for edge in BASE_MISSING_EDGES) + tuple(
        frozenset(face) for face in BASE_MISSING_TRIPLES
    )
    return tuple(
        frozenset(candidate)
        for candidate in combinations(BASE_VERTEX_LABELS, 4)
        if not any(nonface <= frozenset(candidate) for nonface in nonfaces)
    )


@dataclass(frozen=True)
class Record35BaseRealizationAudit:
    support_facets: tuple[tuple[int, ...], ...]
    target_facets: tuple[tuple[int, ...], ...]
    realized_minimal_nonfaces: tuple[tuple[int, ...], ...]
    facet_side_determinant_evaluation_count: int
    five_point_determinant_count: int
    zero_five_point_determinant_count: int
    minimum_absolute_determinant: int
    minimum_facet_side_absolute_determinant: int

    @property
    def verified(self) -> bool:
        expected_nonfaces = tuple(
            sorted(
                (*BASE_MISSING_EDGES, *BASE_MISSING_TRIPLES),
                key=lambda item: (len(item), item),
            )
        )
        return (
            self.support_facets == self.target_facets
            and len(self.support_facets) == 19
            and self.realized_minimal_nonfaces == expected_nonfaces
            and self.facet_side_determinant_evaluation_count == 630
            and self.five_point_determinant_count == 126
            and self.zero_five_point_determinant_count == 0
            and self.minimum_absolute_determinant == 9
            and self.minimum_facet_side_absolute_determinant == 21
        )

    def to_record(self) -> dict[str, object]:
        return {
            "dimension": 4,
            "vertex_labels_in_coordinate_order": list(BASE_VERTEX_LABELS),
            "integral_coordinates": [list(point) for point in BASE_COORDINATES],
            "support_facets": [list(facet) for facet in self.support_facets],
            "target_facets": [list(facet) for facet in self.target_facets],
            "realized_minimal_nonfaces": [
                list(face) for face in self.realized_minimal_nonfaces
            ],
            "facet_side_determinant_evaluation_count": (
                self.facet_side_determinant_evaluation_count
            ),
            "homogeneous_five_point_determinant_count": (
                self.five_point_determinant_count
            ),
            "zero_five_point_determinant_count": (
                self.zero_five_point_determinant_count
            ),
            "minimum_absolute_determinant": self.minimum_absolute_determinant,
            "minimum_facet_side_absolute_determinant": (
                self.minimum_facet_side_absolute_determinant
            ),
            "exact_support_facet_match": self.support_facets == self.target_facets,
            "verified_integral_simplicial_4_polytope": self.verified,
        }


def audit_record35_base_realization() -> Record35BaseRealizationAudit:
    if len(BASE_COORDINATES) != len(BASE_VERTEX_LABELS):
        raise AssertionError("base coordinate/label count changed")
    facet_side_evaluations: list[int] = []
    facet_side_determinants: list[int] = []
    support_facets = []
    for indices in combinations(range(len(BASE_VERTEX_LABELS)), 4):
        outside_values = []
        for outside in range(len(BASE_VERTEX_LABELS)):
            if outside in indices:
                continue
            value = _determinant(
                tuple(
                    (1, *BASE_COORDINATES[index])
                    for index in (*indices, outside)
                )
            )
            facet_side_evaluations.append(value)
            outside_values.append(value)
        if all(value > 0 for value in outside_values) or all(
            value < 0 for value in outside_values
        ):
            facet = tuple(BASE_VERTEX_LABELS[index] for index in indices)
            support_facets.append(facet)
            facet_side_determinants.extend(outside_values)
    support = tuple(sorted(support_facets))
    target = tuple(
        sorted(tuple(sorted(facet)) for facet in _target_base_facets())
    )
    realized_nonfaces = tuple(
        tuple(sorted(face))
        for face in _minimal_nonfaces(
            BASE_VERTEX_LABELS,
            tuple(frozenset(facet) for facet in support),
        )
    )
    five_point_determinants = tuple(
        _determinant(
            tuple((1, *BASE_COORDINATES[index]) for index in indices)
        )
        for indices in combinations(range(len(BASE_VERTEX_LABELS)), 5)
    )
    audit = Record35BaseRealizationAudit(
        support_facets=support,
        target_facets=target,
        realized_minimal_nonfaces=realized_nonfaces,
        facet_side_determinant_evaluation_count=len(facet_side_evaluations),
        five_point_determinant_count=len(five_point_determinants),
        zero_five_point_determinant_count=sum(
            value == 0 for value in five_point_determinants
        ),
        minimum_absolute_determinant=min(map(abs, five_point_determinants)),
        minimum_facet_side_absolute_determinant=min(
            map(abs, facet_side_determinants)
        ),
    )
    if not audit.verified:
        raise AssertionError("record-35 integral base realization failed")
    return audit


def _expand_base_nonface(nonface: frozenset[int]) -> frozenset[int]:
    expanded: set[int] = set()
    for vertex in nonface:
        if vertex == FIRST_CARRIER:
            expanded.update(FIRST_CLONE_BLOCK)
        elif vertex == SECOND_CARRIER:
            expanded.update(SECOND_CLONE_BLOCK)
        else:
            expanded.add(vertex)
    return frozenset(expanded)


@dataclass(frozen=True)
class Record35MultiwedgeAudit:
    base_minimal_nonfaces: tuple[tuple[int, ...], ...]
    expanded_minimal_nonfaces: tuple[tuple[int, ...], ...]
    target_minimal_nonfaces: tuple[tuple[int, ...], ...]
    base_polytope_dimension: int
    wedge_count: int
    lifted_polytope_dimension: int
    lifted_facet_count: int

    @property
    def verified(self) -> bool:
        return (
            self.expanded_minimal_nonfaces == self.target_minimal_nonfaces
            and self.base_polytope_dimension == 4
            and self.wedge_count == 5
            and self.lifted_polytope_dimension == 9
            and self.lifted_facet_count == 14
        )

    def to_record(self) -> dict[str, object]:
        return {
            "base_minimal_nonfaces": [
                list(face) for face in self.base_minimal_nonfaces
            ],
            "first_carrier": FIRST_CARRIER,
            "first_clone_block": list(FIRST_CLONE_BLOCK),
            "second_carrier": SECOND_CARRIER,
            "second_clone_block": list(SECOND_CLONE_BLOCK),
            "expanded_minimal_nonfaces": [
                list(face) for face in self.expanded_minimal_nonfaces
            ],
            "target_minimal_nonfaces": [
                list(face) for face in self.target_minimal_nonfaces
            ],
            "exact_minimal_nonface_match": (
                self.expanded_minimal_nonfaces == self.target_minimal_nonfaces
            ),
            "base_polytope_dimension": self.base_polytope_dimension,
            "simplicial_wedge_count": self.wedge_count,
            "lifted_polytope_dimension": self.lifted_polytope_dimension,
            "lifted_facet_count": self.lifted_facet_count,
            "ordinary_polytopal": self.verified,
        }


def audit_record35_multiwedge() -> Record35MultiwedgeAudit:
    base = tuple(frozenset(edge) for edge in BASE_MISSING_EDGES) + tuple(
        frozenset(face) for face in BASE_MISSING_TRIPLES
    )
    expanded = tuple(_expand_base_nonface(nonface) for nonface in base)
    target = tuple(frozenset(edge) for edge in TARGET_MISSING_EDGES) + tuple(
        frozenset(face) for face in TARGET_MISSING_FIVES
    )
    key = lambda face: (len(face), tuple(sorted(face)))
    audit = Record35MultiwedgeAudit(
        base_minimal_nonfaces=tuple(
            tuple(sorted(face)) for face in sorted(base, key=key)
        ),
        expanded_minimal_nonfaces=tuple(
            tuple(sorted(face)) for face in sorted(expanded, key=key)
        ),
        target_minimal_nonfaces=tuple(
            tuple(sorted(face)) for face in sorted(target, key=key)
        ),
        base_polytope_dimension=4,
        wedge_count=(len(FIRST_CLONE_BLOCK) - 1)
        + (len(SECOND_CLONE_BLOCK) - 1),
        lifted_polytope_dimension=9,
        lifted_facet_count=14,
    )
    if not audit.verified:
        raise AssertionError("record-35 inverse multiwedge failed")
    return audit


__all__ = [
    "BASE_COORDINATES",
    "BASE_VERTEX_LABELS",
    "Record35BaseRealizationAudit",
    "Record35IncidenceAudit",
    "Record35MultiwedgeAudit",
    "TARGET_F_VECTOR",
    "TARGET_GRAPH6",
    "TARGET_MISSING_EDGES",
    "TARGET_MISSING_FIVES",
    "TARGET_RECORD_INDEX",
    "audit_record35_base_realization",
    "audit_record35_incidence",
    "audit_record35_multiwedge",
    "record35_clutter",
    "record35_minimal_transversals",
]
