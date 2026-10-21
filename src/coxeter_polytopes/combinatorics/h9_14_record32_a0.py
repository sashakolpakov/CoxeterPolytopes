"""Exact ``a=0`` completion and compact-Coxeter gate for H9/14 record 32.

This is the second explicit completion currently frozen over the record-32
missing-edge graph.  Its six higher minimal nonfaces all have order five.
The blocker is pure: its 114 minimal transversals are the complements of the
rank-nine vertices of an abstract homology 8-sphere candidate.

No ordinary-polytopality assertion is made here.  For compact Coxeter
realizability it is enough to use three of the six missing faces.  They share
the four-node core ``(5, 7, 8, 9)``.  Exhausting the classical order-five
Lanner placements on those supports leaves no labeling whose restrictions
to all actual vertices are elliptic.  Adding any other non-right edge only
increases the nonnegative cosine-adjacency matrix, so Perron--Frobenius
monotonicity makes the obstruction unrestricted in the remaining labels.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from typing import Iterable

from coxeter_polytopes.combinatorics.facet_completion import generated_face_table
from coxeter_polytopes.combinatorics.h9_14_record35_coxeter_labels import (
    CommonFourCoreLannerShadowAudit,
    audit_common_four_core_lanner_shadow,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 32
TARGET_SCHEDULE_INDEX = 0
TARGET_GRAPH6 = "M?????????O?O?g??"
TARGET_VERTEX_COUNT = 14
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
TARGET_MISSING_EDGES = (
    (0, 11),
    (1, 12),
    (0, 13),
    (2, 13),
)
TARGET_MISSING_FIVES = (
    (4, 5, 7, 8, 9),
    (1, 2, 3, 6, 10),
    (2, 3, 4, 6, 10),
    (3, 4, 6, 10, 11),
    (5, 7, 8, 9, 12),
    (5, 7, 8, 9, 13),
)
TARGET_F_VECTOR = (14, 87, 318, 762, 1254, 1428, 1092, 513, 114)

COMMON_FOUR_CORE = (5, 7, 8, 9)
COMMON_FOUR_OUTERS = (4, 12, 13)
COMMON_FOUR_SUPPORTS = tuple(
    tuple(sorted((*COMMON_FOUR_CORE, outer))) for outer in COMMON_FOUR_OUTERS
)


def mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-32 a=0 vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    if value < 0 or value & ~FULL_MASK:
        raise ValueError("record-32 a=0 mask lies outside [14]")
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex)
    )


def record32_a0_clutter() -> tuple[int, ...]:
    return tuple(mask(face) for face in (*TARGET_MISSING_EDGES, *TARGET_MISSING_FIVES))


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


def record32_a0_minimal_transversals() -> tuple[int, ...]:
    clutter = record32_a0_clutter()
    return tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if _is_minimal_transversal(candidate, clutter)
    )


def record32_a0_facet_complements() -> tuple[int, ...]:
    transversals = record32_a0_minimal_transversals()
    if Counter(value.bit_count() for value in transversals) != {5: 114}:
        raise AssertionError("record-32 a=0 blocker is not pure of order five")
    return transversals


def record32_a0_vertex_facets() -> tuple[tuple[int, ...], ...]:
    """Return all 114 rank-nine vertex diagrams deterministically."""

    return tuple(
        sorted(
            vertices(FULL_MASK ^ complement)
            for complement in record32_a0_facet_complements()
        )
    )


def _packed_face_table(face_table: tuple[bool, ...]) -> bytes:
    if len(face_table) != 1 << TARGET_VERTEX_COUNT:
        raise ValueError("record-32 a=0 face table must have 2^14 entries")
    return bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )


@dataclass(frozen=True, slots=True)
class Record32A0IncidenceAudit:
    missing_edges: tuple[tuple[int, int], ...]
    missing_fives: tuple[tuple[int, ...], ...]
    minimal_transversal_size_histogram: tuple[tuple[int, int], ...]
    facet_complement_masks: tuple[int, ...]
    facet_complement_sha256: str
    face_table_true_count: int
    packed_face_table_sha256: str

    @property
    def verified(self) -> bool:
        return (
            self.missing_edges == TARGET_MISSING_EDGES
            and self.missing_fives == TARGET_MISSING_FIVES
            and self.minimal_transversal_size_histogram == ((5, 114),)
            and len(self.facet_complement_masks) == TARGET_F_VECTOR[-1]
        )

    def to_record(self) -> dict[str, object]:
        return {
            "missing_edges": [list(face) for face in self.missing_edges],
            "missing_fours": [],
            "missing_fives": [list(face) for face in self.missing_fives],
            "minimal_transversal_size_histogram": {
                str(size): count
                for size, count in self.minimal_transversal_size_histogram
            },
            "facet_complement_masks": list(self.facet_complement_masks),
            "facet_complement_vertices": [
                list(vertices(value)) for value in self.facet_complement_masks
            ],
            "facet_complement_sha256": self.facet_complement_sha256,
            "face_table_encoding": (
                "mask order 0..2^14-1; one bit per mask; least-significant "
                "bit first in each byte"
            ),
            "face_table_true_count": self.face_table_true_count,
            "packed_face_table_sha256": self.packed_face_table_sha256,
        }


def audit_record32_a0_incidence(
    graph6: str = TARGET_GRAPH6,
) -> Record32A0IncidenceAudit:
    order, edges = decode_graph6(graph6)
    if order != TARGET_VERTEX_COUNT or edges != TARGET_MISSING_EDGES:
        raise ValueError("record-32 a=0 missing-edge graph changed")
    facets = tuple(sorted(record32_a0_minimal_transversals()))
    histogram = tuple(
        sorted(Counter(value.bit_count() for value in facets).items())
    )
    facet_payload = b"".join(value.to_bytes(2, "little") for value in facets)
    face_table = generated_face_table(facets)
    audit = Record32A0IncidenceAudit(
        missing_edges=edges,
        missing_fives=TARGET_MISSING_FIVES,
        minimal_transversal_size_histogram=histogram,
        facet_complement_masks=facets,
        facet_complement_sha256=sha256(facet_payload).hexdigest(),
        face_table_true_count=sum(face_table),
        packed_face_table_sha256=sha256(_packed_face_table(face_table)).hexdigest(),
    )
    if not audit.verified:
        raise AssertionError("record-32 a=0 incidence audit failed")
    return audit


def audit_record32_a0_common_four_core() -> CommonFourCoreLannerShadowAudit:
    audit = audit_common_four_core_lanner_shadow(
        COMMON_FOUR_OUTERS,
        COMMON_FOUR_CORE,
        record32_a0_vertex_facets(),
    )
    expected_stages = ((4, 13, 4620), (5, 12, 96), (6, 4, 0))
    observed_stages = tuple(
        (stage.order, len(stage.restrictions), stage.passing_state_count)
        for stage in audit.stages
    )
    if audit.supports != COMMON_FOUR_SUPPORTS:
        raise AssertionError("record-32 a=0 common-core supports changed")
    if audit.join_counts != (420, 1332, 4620):
        raise AssertionError("record-32 a=0 common-core join census changed")
    if audit.core_extension_count_histogram != (
        (1, 36),
        (2, 30),
        (3, 40),
        (4, 51),
    ):
        raise AssertionError("record-32 a=0 core-extension census changed")
    if observed_stages != expected_stages:
        raise AssertionError("record-32 a=0 restriction census changed")
    if audit.sparse_shadow_passing_state_count != 0:
        raise AssertionError("record-32 a=0 common-core obstruction disappeared")
    if not audit.unrestricted_external_label_no_go:
        raise AssertionError("record-32 a=0 monotone no-go was not certified")
    return audit


__all__ = [
    "COMMON_FOUR_CORE",
    "COMMON_FOUR_OUTERS",
    "COMMON_FOUR_SUPPORTS",
    "Record32A0IncidenceAudit",
    "TARGET_F_VECTOR",
    "TARGET_GRAPH6",
    "TARGET_MISSING_EDGES",
    "TARGET_MISSING_FIVES",
    "TARGET_RECORD_INDEX",
    "TARGET_SCHEDULE_INDEX",
    "audit_record32_a0_common_four_core",
    "audit_record32_a0_incidence",
    "record32_a0_clutter",
    "record32_a0_facet_complements",
    "record32_a0_minimal_transversals",
    "record32_a0_vertex_facets",
]
