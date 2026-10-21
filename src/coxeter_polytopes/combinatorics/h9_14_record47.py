"""Exact incidence certificate for the positive H9/14 record-47 completion."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from typing import Iterable

from coxeter_polytopes.combinatorics.facet_completion import generated_face_table
from coxeter_polytopes.combinatorics.h9_14_record54 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 47
TARGET_GRAPH6 = "M?????????O?o?g?_"
TARGET_MISSING_EDGES = (
    (0, 11),
    (0, 12),
    (1, 12),
    (0, 13),
    (2, 13),
    (12, 13),
)
TARGET_MISSING_FIVES = (
    (3, 4, 5, 6, 7),
    (1, 2, 8, 9, 10),
    (1, 8, 9, 10, 11),
    (2, 8, 9, 10, 11),
)
TARGET_F_VECTOR = (14, 85, 300, 690, 1088, 1190, 880, 405, 90)
TARGET_BASE_LINK_G = (1, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 1, 1)
TARGET_VERTEX_COUNT = 14
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1


def mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-47 vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    if value < 0 or value & ~FULL_MASK:
        raise ValueError("record-47 mask lies outside [14]")
    return tuple(
        vertex
        for vertex in range(TARGET_VERTEX_COUNT)
        if value & (1 << vertex)
    )


def record47_clutter() -> tuple[int, ...]:
    return tuple(mask(edge) for edge in TARGET_MISSING_EDGES) + tuple(
        mask(face) for face in TARGET_MISSING_FIVES
    )


def _is_transversal(candidate: int, hyperedges: tuple[int, ...]) -> bool:
    return all(candidate & edge for edge in hyperedges)


def _is_minimal_transversal(candidate: int, hyperedges: tuple[int, ...]) -> bool:
    if not _is_transversal(candidate, hyperedges):
        return False
    remaining = candidate
    while remaining:
        bit = remaining & -remaining
        if _is_transversal(candidate ^ bit, hyperedges):
            return False
        remaining ^= bit
    return True


def record47_minimal_transversals() -> tuple[int, ...]:
    clutter = record47_clutter()
    return tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if _is_minimal_transversal(candidate, clutter)
    )


def record47_facet_complements() -> tuple[int, ...]:
    transversals = record47_minimal_transversals()
    if Counter(value.bit_count() for value in transversals) != {5: 90}:
        raise AssertionError("record-47 blocker is no longer pure of order five")
    return transversals


def _packed_face_table(face_table: tuple[bool, ...]) -> bytes:
    if len(face_table) != 1 << TARGET_VERTEX_COUNT:
        raise ValueError("record-47 face table must have 2^14 entries")
    return bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )


@dataclass(frozen=True)
class Record47IncidenceAudit:
    missing_edges: tuple[tuple[int, int], ...]
    missing_fives: tuple[tuple[int, ...], ...]
    minimal_transversal_size_histogram: tuple[tuple[int, int], ...]
    facet_complement_masks: tuple[int, ...]
    facet_complement_sha256: str
    face_table_true_count: int
    packed_face_table_sha256: str
    base_vertex_link_g3: tuple[int, ...]
    base_vertex_link_g4: tuple[int, ...]

    @property
    def verified(self) -> bool:
        return (
            self.missing_edges == TARGET_MISSING_EDGES
            and self.missing_fives == TARGET_MISSING_FIVES
            and self.minimal_transversal_size_histogram == ((5, 90),)
            and len(self.facet_complement_masks) == 90
            and self.base_vertex_link_g3 == TARGET_BASE_LINK_G
            and self.base_vertex_link_g4 == TARGET_BASE_LINK_G
        )

    def to_record(self) -> dict[str, object]:
        return {
            "missing_edges": [list(edge) for edge in self.missing_edges],
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
            "base_vertex_link_g3": list(self.base_vertex_link_g3),
            "base_vertex_link_g4": list(self.base_vertex_link_g4),
        }


def audit_record47_incidence(
    graph6: str = TARGET_GRAPH6,
) -> Record47IncidenceAudit:
    order, edges = decode_graph6(graph6)
    if order != TARGET_VERTEX_COUNT or edges != TARGET_MISSING_EDGES:
        raise ValueError("record-47 missing-edge graph changed")
    transversals = record47_minimal_transversals()
    histogram = tuple(
        sorted(Counter(value.bit_count() for value in transversals).items())
    )
    if histogram != ((5, 90),):
        raise AssertionError("record-47 blocker has a non-five minimal transversal")
    facets = tuple(sorted(transversals))
    facet_payload = b"".join(value.to_bytes(2, "little") for value in facets)
    face_table = generated_face_table(facets)
    packed = _packed_face_table(face_table)
    base_g3, base_g4 = _base_link_g3_g4(graph6)
    audit = Record47IncidenceAudit(
        missing_edges=edges,
        missing_fives=TARGET_MISSING_FIVES,
        minimal_transversal_size_histogram=histogram,
        facet_complement_masks=facets,
        facet_complement_sha256=sha256(facet_payload).hexdigest(),
        face_table_true_count=sum(face_table),
        packed_face_table_sha256=sha256(packed).hexdigest(),
        base_vertex_link_g3=base_g3,
        base_vertex_link_g4=base_g4,
    )
    if not audit.verified:
        raise AssertionError("record-47 incidence audit failed")
    return audit


__all__ = [
    "Record47IncidenceAudit",
    "TARGET_F_VECTOR",
    "TARGET_GRAPH6",
    "TARGET_MISSING_EDGES",
    "TARGET_MISSING_FIVES",
    "TARGET_RECORD_INDEX",
    "audit_record47_incidence",
    "record47_clutter",
    "record47_facet_complements",
    "record47_minimal_transversals",
]
