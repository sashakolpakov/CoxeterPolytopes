"""Exact incidence and convex multiwedge certificate for H9/14 record 32.

The selected ``a=2`` completion has two missing four-faces and one missing
five-face.  Its three-vertex clone block collapses to a seven-vertex base
whose missing graph is ``C5 + K2``.  The base is the boundary of a
pentagonal bipyramid.  Two simplicial wedges and a join with the boundary of
a four-simplex recover the fourteen-vertex, nine-dimensional nerve.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from typing import Iterable

from coxeter_polytopes.combinatorics.facet_completion import generated_face_table
from coxeter_polytopes.combinatorics.h9_14_record54 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 32
TARGET_GRAPH6 = "M?????????O?O?g??"
TARGET_VERTEX_COUNT = 14
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
TARGET_MISSING_EDGES = (
    (0, 11),
    (1, 12),
    (0, 13),
    (2, 13),
)
TARGET_MISSING_FOURS = (
    (2, 3, 4, 5),
    (3, 4, 5, 11),
)
TARGET_MISSING_FIVES = ((6, 7, 8, 9, 10),)
TARGET_F_VECTOR = (14, 87, 318, 760, 1242, 1400, 1060, 495, 110)
TARGET_BASE_LINK_G = (2, 3, 4, 6, 6, 6, 6, 6, 6, 6, 6, 4, 3, 2)

CLONE_BLOCK = (3, 4, 5)
CLONE_CARRIER = 14
SIMPLEX_JOIN_VERTICES = (6, 7, 8, 9, 10)
BASE_VERTICES = (0, 1, 2, 11, 12, 13, CLONE_CARRIER)
BASE_MISSING_EDGES = (
    (0, 11),
    (1, 12),
    (0, 13),
    (2, 13),
    (2, CLONE_CARRIER),
    (11, CLONE_CARRIER),
)
# The five allowed equatorial edges form this cycle; 1 and 12 are its two
# suspension vertices.  Thus the base is a pentagonal bipyramid.
BASE_EQUATOR_CYCLE = (11, 13, CLONE_CARRIER, 0, 2)
BASE_SUSPENSION_VERTICES = (1, 12)


def mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-32 vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    if value < 0 or value & ~FULL_MASK:
        raise ValueError("record-32 mask lies outside [14]")
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex)
    )


def record32_clutter() -> tuple[int, ...]:
    return (
        tuple(mask(edge) for edge in TARGET_MISSING_EDGES)
        + tuple(mask(face) for face in TARGET_MISSING_FOURS)
        + tuple(mask(face) for face in TARGET_MISSING_FIVES)
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


def record32_minimal_transversals() -> tuple[int, ...]:
    clutter = record32_clutter()
    return tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if _is_minimal_transversal(candidate, clutter)
    )


def record32_facet_complements() -> tuple[int, ...]:
    transversals = record32_minimal_transversals()
    if Counter(value.bit_count() for value in transversals) != {5: 110}:
        raise AssertionError("record-32 blocker is no longer pure of order five")
    return transversals


def _packed_face_table(face_table: tuple[bool, ...]) -> bytes:
    if len(face_table) != 1 << TARGET_VERTEX_COUNT:
        raise ValueError("record-32 face table must have 2^14 entries")
    return bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )


@dataclass(frozen=True)
class Record32IncidenceAudit:
    missing_edges: tuple[tuple[int, int], ...]
    missing_fours: tuple[tuple[int, ...], ...]
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
            and self.missing_fours == TARGET_MISSING_FOURS
            and self.missing_fives == TARGET_MISSING_FIVES
            and self.minimal_transversal_size_histogram == ((5, 110),)
            and len(self.facet_complement_masks) == TARGET_F_VECTOR[-1]
            and self.base_vertex_link_g3 == TARGET_BASE_LINK_G
            and self.base_vertex_link_g4 == TARGET_BASE_LINK_G
        )

    def to_record(self) -> dict[str, object]:
        return {
            "missing_edges": [list(face) for face in self.missing_edges],
            "missing_fours": [list(face) for face in self.missing_fours],
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


def audit_record32_incidence(
    graph6: str = TARGET_GRAPH6,
) -> Record32IncidenceAudit:
    order, edges = decode_graph6(graph6)
    if order != TARGET_VERTEX_COUNT or edges != TARGET_MISSING_EDGES:
        raise ValueError("record-32 missing-edge graph changed")
    facets = tuple(sorted(record32_minimal_transversals()))
    histogram = tuple(
        sorted(Counter(value.bit_count() for value in facets).items())
    )
    facet_payload = b"".join(value.to_bytes(2, "little") for value in facets)
    face_table = generated_face_table(facets)
    base_g3, base_g4 = _base_link_g3_g4(graph6)
    audit = Record32IncidenceAudit(
        missing_edges=edges,
        missing_fours=TARGET_MISSING_FOURS,
        missing_fives=TARGET_MISSING_FIVES,
        minimal_transversal_size_histogram=histogram,
        facet_complement_masks=facets,
        facet_complement_sha256=sha256(facet_payload).hexdigest(),
        face_table_true_count=sum(face_table),
        packed_face_table_sha256=sha256(_packed_face_table(face_table)).hexdigest(),
        base_vertex_link_g3=base_g3,
        base_vertex_link_g4=base_g4,
    )
    if not audit.verified:
        raise AssertionError("record-32 incidence audit failed")
    return audit


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
    return tuple(sorted(result, key=lambda face: (len(face), tuple(sorted(face)))))


def _base_facets() -> tuple[frozenset[int], ...]:
    facets = []
    for index, first in enumerate(BASE_EQUATOR_CYCLE):
        second = BASE_EQUATOR_CYCLE[(index + 1) % len(BASE_EQUATOR_CYCLE)]
        for apex in BASE_SUSPENSION_VERTICES:
            facets.append(frozenset((first, second, apex)))
    return tuple(sorted(facets, key=lambda face: tuple(sorted(face))))


def _expand_base_nonface(nonface: frozenset[int]) -> frozenset[int]:
    expanded: set[int] = set()
    for vertex in nonface:
        if vertex == CLONE_CARRIER:
            expanded.update(CLONE_BLOCK)
        else:
            expanded.add(vertex)
    return frozenset(expanded)


@dataclass(frozen=True)
class Record32MultiwedgeAudit:
    base_facets: tuple[tuple[int, ...], ...]
    base_face_vector: tuple[int, ...]
    base_minimal_nonfaces: tuple[tuple[int, ...], ...]
    expanded_minimal_nonfaces: tuple[tuple[int, ...], ...]
    target_minimal_nonfaces: tuple[tuple[int, ...], ...]
    base_polytope_dimension: int
    wedge_count: int
    lifted_factor_dimension: int
    simplex_join_dimension: int
    target_polytope_dimension: int
    target_vertex_count: int

    @property
    def verified(self) -> bool:
        return (
            self.base_face_vector == (7, 15, 10)
            and self.base_minimal_nonfaces
            == tuple(tuple(sorted(edge)) for edge in sorted(BASE_MISSING_EDGES))
            and self.expanded_minimal_nonfaces == self.target_minimal_nonfaces
            and self.base_polytope_dimension == 3
            and self.wedge_count == 2
            and self.lifted_factor_dimension == 5
            and self.simplex_join_dimension == 4
            and self.target_polytope_dimension == 9
            and self.target_vertex_count == 14
        )

    def to_record(self) -> dict[str, object]:
        return {
            "base_description": (
                "join of a pentagon boundary and a zero-sphere: the boundary "
                "of a convex pentagonal bipyramid"
            ),
            "base_vertices": list(BASE_VERTICES),
            "base_equator_cycle": list(BASE_EQUATOR_CYCLE),
            "base_suspension_vertices": list(BASE_SUSPENSION_VERTICES),
            "base_facets": [list(face) for face in self.base_facets],
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces": [
                list(face) for face in self.base_minimal_nonfaces
            ],
            "clone_carrier": CLONE_CARRIER,
            "clone_block": list(CLONE_BLOCK),
            "simplicial_wedge_count": self.wedge_count,
            "expanded_minimal_nonfaces": [
                list(face) for face in self.expanded_minimal_nonfaces
            ],
            "simplex_join_vertices": list(SIMPLEX_JOIN_VERTICES),
            "target_minimal_nonfaces": [
                list(face) for face in self.target_minimal_nonfaces
            ],
            "exact_minimal_nonface_match": (
                self.expanded_minimal_nonfaces == self.target_minimal_nonfaces
            ),
            "base_polytope_dimension": self.base_polytope_dimension,
            "lifted_factor_dimension": self.lifted_factor_dimension,
            "simplex_join_dimension": self.simplex_join_dimension,
            "target_polytope_dimension": self.target_polytope_dimension,
            "target_vertex_count": self.target_vertex_count,
            "ordinary_polytopal": self.verified,
        }


def audit_record32_multiwedge() -> Record32MultiwedgeAudit:
    base_facets = _base_facets()
    base_faces = _all_faces(base_facets)
    base_nonfaces = _minimal_nonfaces(BASE_VERTICES, base_facets)
    expanded = tuple(_expand_base_nonface(face) for face in base_nonfaces) + (
        frozenset(SIMPLEX_JOIN_VERTICES),
    )
    target = (
        tuple(frozenset(face) for face in TARGET_MISSING_EDGES)
        + tuple(frozenset(face) for face in TARGET_MISSING_FOURS)
        + tuple(frozenset(face) for face in TARGET_MISSING_FIVES)
    )
    key = lambda face: (len(face), tuple(sorted(face)))
    audit = Record32MultiwedgeAudit(
        base_facets=tuple(tuple(sorted(face)) for face in base_facets),
        base_face_vector=tuple(
            sum(len(face) == size for face in base_faces) for size in range(1, 4)
        ),
        base_minimal_nonfaces=tuple(
            tuple(sorted(face)) for face in sorted(base_nonfaces, key=key)
        ),
        expanded_minimal_nonfaces=tuple(
            tuple(sorted(face)) for face in sorted(expanded, key=key)
        ),
        target_minimal_nonfaces=tuple(
            tuple(sorted(face)) for face in sorted(target, key=key)
        ),
        base_polytope_dimension=3,
        wedge_count=len(CLONE_BLOCK) - 1,
        lifted_factor_dimension=3 + len(CLONE_BLOCK) - 1,
        simplex_join_dimension=4,
        target_polytope_dimension=3 + len(CLONE_BLOCK) - 1 + 4,
        target_vertex_count=(
            len(BASE_VERTICES) + len(CLONE_BLOCK) - 1 + len(SIMPLEX_JOIN_VERTICES)
        ),
    )
    if not audit.verified:
        raise AssertionError("record-32 multiwedge factorization failed")
    return audit


__all__ = [
    "Record32IncidenceAudit",
    "Record32MultiwedgeAudit",
    "TARGET_F_VECTOR",
    "TARGET_GRAPH6",
    "TARGET_MISSING_EDGES",
    "TARGET_MISSING_FIVES",
    "TARGET_MISSING_FOURS",
    "TARGET_RECORD_INDEX",
    "audit_record32_incidence",
    "audit_record32_multiwedge",
    "record32_clutter",
    "record32_facet_complements",
    "record32_minimal_transversals",
]
