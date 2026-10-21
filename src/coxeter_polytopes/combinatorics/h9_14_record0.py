"""Exact rational ordinary-polytopal completion of H9/14 record 0.

The construction starts from ``C(8,5)``, stellar-subdivides edges ``{0,1}``
and ``{6,8}``, then applies simplicial wedges at vertices 0, 1, 6, 8.
It proves ordinary rational convex polytopality only, not a Coxeter or
hyperbolic realization.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Mapping

from coxeter_polytopes.combinatorics.h9_14_record12 import (
    _face_vector,
    cyclic_facets,
    minimal_nonfaces,
    radial_beyond_point,
    stellar_subdivision,
    supporting_facets,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 0
TARGET_GRAPH6 = "M?????????????o??"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
TARGET_MISSING_EDGES = ((0, 13), (1, 13))
TARGET_MISSING_FOURS = (
    (0, 3, 5, 10),
    (0, 5, 8, 12),
    (1, 4, 6, 11),
    (2, 3, 9, 10),
    (4, 5, 7, 13),
    (6, 8, 11, 12),
)
TARGET_MISSING_FIVES = (
    (0, 2, 5, 7, 9),
    (0, 3, 6, 10, 11),
    (1, 2, 4, 7, 9),
    (1, 3, 6, 10, 11),
    (1, 4, 7, 8, 12),
    (2, 4, 5, 7, 9),
    (4, 5, 7, 8, 12),
)
TARGET_F_VECTOR = (14, 89, 341, 874, 1555, 1911, 1554, 756, 168)

CYCLIC_LABELS = tuple(range(8))
BASE_LABELS = tuple(range(10))
FIRST_STELLAR_FACE = (0, 1)
FIRST_STELLAR_VERTEX = 8
SECOND_STELLAR_FACE = (6, 8)
SECOND_STELLAR_VERTEX = 9
FIRST_EPSILON = SECOND_EPSILON = Fraction(1, 16)
BASE_F_VECTOR = (10, 41, 84, 85, 34)
BASE_MINIMAL_NONFACES = (
    (0, 1),
    (6, 8),
    (3, 9),
    (4, 9),
    (1, 3, 5),
    (1, 3, 6),
    (1, 4, 6),
    (2, 4, 6),
    (3, 5, 8),
    (0, 2, 4, 7),
    (0, 2, 5, 7),
    (0, 3, 5, 7),
    (2, 4, 7, 8),
    (2, 5, 7, 8),
    (2, 5, 7, 9),
)
CLONE_BLOCKS: Mapping[int, tuple[int, int]] = {
    0: (0, 10),
    1: (1, 11),
    6: (6, 12),
    8: (8, 13),
}
TARGET_RELABEL: Mapping[int, int] = {
    0: 2,
    1: 3,
    2: 4,
    3: 0,
    4: 1,
    5: 5,
    6: 6,
    7: 7,
    8: 8,
    9: 13,
    10: 9,
    11: 10,
    12: 11,
    13: 12,
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
        Fraction(4849, 1536),
        Fraction(82289, 4608),
        Fraction(326777, 3072),
        Fraction(45875, 72),
        Fraction(11707433, 3072),
    ),
}


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex))


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
class Record0PolytopalityAudit:
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
            and self.stellar_facet_counts == (20, 28, 34)
            and all(self.radial_formulas_verified)
            and all(self.coordinate_support_matches)
            and self.base_face_vector == BASE_F_VECTOR
            and self.base_minimal_nonfaces == expected_base
            and self.expanded_minimal_nonfaces == self.target_minimal_nonfaces
            and self.blocker_order_histogram == ((5, 168),)
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


def audit_record0_polytopality() -> Record0PolytopalityAudit:
    vertex_count, edges = decode_graph6(TARGET_GRAPH6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise AssertionError("record-0 graph6 order changed")
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
    audit = Record0PolytopalityAudit(
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
        raise AssertionError("record-0 rational polytopality audit failed")
    return audit

