"""Exact incidence data for an abstract H9/14 record-18 candidate nerve.

The five missing edges have three inclusion-minimal graph covers.  After
adjoining the eight declared missing five-faces, the complete blocker has
exactly 110 members, all of order five.  Their complements generate a pure
eight-dimensional complex on fourteen vertices.

This module certifies incidence only.  In particular, a homology-sphere or
oriented-matroid audit does not prove PL sphericity, convex polytopality, a
Coxeter labelling, a Lorentzian Gram realization, or finite hyperbolic volume.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
from typing import Iterable

from coxeter_polytopes.combinatorics.facet_completion import generated_face_table
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 18
TARGET_GRAPH6 = "M???????????o?K?_"
TARGET_VERTEX_COUNT = 14
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
TARGET_F_VECTOR = (14, 86, 310, 735, 1204, 1372, 1052, 495, 110)
TARGET_BASE_LINK_G3 = (5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1)
TARGET_BASE_LINK_G4 = (6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 1, 1)

LEFT = (0, 1)
RIGHT = (2, 3)
INTERIOR = tuple(range(4, 12))
WORK = tuple(range(12))
SPECIAL = (12, 13)
EXPECTED_MINIMAL_GRAPH_COVERS = (
    (2, 3, 12),
    (0, 1, 13),
    (12, 13),
)
EXPECTED_Z_CUT_PARTS = ((2, 4, 5, 6, 11), (3, 7, 8, 9, 10))
EXPECTED_W_CUT_PARTS = ((0, 5, 6, 7, 8), (1, 4, 9, 10, 11))
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1


def mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-18 vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    if value < 0 or value & ~FULL_MASK:
        raise ValueError("record-18 mask lies outside [14]")
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex)
    )


EDGE_MASKS = tuple(mask(edge) for edge in TARGET_MISSING_EDGES)


def record18_clutter() -> tuple[int, ...]:
    """Return all thirteen declared inclusion-minimal nonfaces as masks."""

    return EDGE_MASKS + tuple(mask(face) for face in TARGET_MISSING_FIVES)


def is_transversal(candidate: int, hyperedges: tuple[int, ...]) -> bool:
    return all(candidate & edge for edge in hyperedges)


def is_minimal_transversal(candidate: int, hyperedges: tuple[int, ...]) -> bool:
    if not is_transversal(candidate, hyperedges):
        return False
    return all(
        not is_transversal(candidate ^ (1 << vertex), hyperedges)
        for vertex in vertices(candidate)
    )


def record18_minimal_transversals() -> tuple[int, ...]:
    """Enumerate the complete blocker directly over all ``2^14`` masks."""

    clutter = record18_clutter()
    return tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if is_minimal_transversal(candidate, clutter)
    )


def record18_facet_complements() -> tuple[int, ...]:
    transversals = record18_minimal_transversals()
    if Counter(value.bit_count() for value in transversals) != {5: 110}:
        raise AssertionError("record-18 blocker is no longer pure of order five")
    return transversals


def minimal_graph_covers() -> tuple[int, ...]:
    """Enumerate the minimal transversals of the missing-edge graph."""

    return tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if is_minimal_transversal(candidate, EDGE_MASKS)
    )


def packed_face_table(face_table: tuple[bool, ...]) -> bytes:
    if len(face_table) != 1 << TARGET_VERTEX_COUNT:
        raise ValueError("record-18 face table must have 2^14 entries")
    return bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )


def _base_vertex_link_g3_g4() -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Compute link ``g3,g4`` from the missing-edge shadow alone."""

    independent: dict[int, tuple[int, ...]] = {}
    for size in range(2, 6):
        independent[size] = tuple(
            mask(chosen)
            for chosen in combinations(range(TARGET_VERTEX_COUNT), size)
            if not any(mask(chosen) & edge == edge for edge in EDGE_MASKS)
        )
    g3_values = []
    g4_values = []
    for vertex in range(TARGET_VERTEX_COUNT):
        counts = {
            size: sum(bool(face & (1 << vertex)) for face in independent[size])
            for size in independent
        }
        f0, f1, f2, f3 = (counts[size] for size in range(2, 6))
        g3_values.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4_values.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3_values), tuple(g4_values)


def _complete_bipartite_parts(
    domain: tuple[int, ...], edge_pairs: tuple[tuple[int, int], ...]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return canonical bipartition, rejecting anything other than ``K5,5``."""

    adjacency = {vertex: set() for vertex in domain}
    for first, second in edge_pairs:
        if first not in adjacency or second not in adjacency or first == second:
            raise AssertionError("cut edge lies outside its ten-vertex domain")
        adjacency[first].add(second)
        adjacency[second].add(first)
    colour = {domain[0]: 0}
    queue = deque([domain[0]])
    while queue:
        vertex = queue.popleft()
        for neighbor in adjacency[vertex]:
            expected = 1 - colour[vertex]
            if neighbor in colour and colour[neighbor] != expected:
                raise AssertionError("selected blocker-pair graph is not bipartite")
            if neighbor not in colour:
                colour[neighbor] = expected
                queue.append(neighbor)
    if len(colour) != len(domain):
        raise AssertionError("selected blocker-pair graph is disconnected")
    parts = tuple(
        tuple(sorted(vertex for vertex in domain if colour[vertex] == side))
        for side in (0, 1)
    )
    parts = tuple(sorted(parts))
    if tuple(map(len, parts)) != (5, 5):
        raise AssertionError("selected blocker-pair cut is not balanced")
    expected_edges = {
        tuple(sorted((first, second))) for first in parts[0] for second in parts[1]
    }
    if set(edge_pairs) != expected_edges:
        raise AssertionError("selected blocker-pair graph is not complete bipartite")
    return parts  # type: ignore[return-value]


@dataclass(frozen=True)
class Record18StructuralAudit:
    minimal_graph_covers: tuple[tuple[int, ...], ...]
    graph_cover_five_set_count: int
    potential_y_count: int
    potential_z_count: int
    potential_w_count: int
    nonempty_symbolic_ridge_groups: int
    two_special_ridge_groups: int
    z_ridge_groups: int
    w_ridge_groups: int
    y_ridge_incidences: int
    z_ridge_incidences: int
    w_ridge_incidences: int
    selected_y_count: int
    selected_z_count: int
    selected_w_count: int
    z_cut_parts: tuple[tuple[int, ...], tuple[int, ...]]
    w_cut_parts: tuple[tuple[int, ...], tuple[int, ...]]

    @property
    def verified(self) -> bool:
        return (
            self.minimal_graph_covers == EXPECTED_MINIMAL_GRAPH_COVERS
            and (
                self.graph_cover_five_set_count,
                self.potential_y_count,
                self.potential_z_count,
                self.potential_w_count,
            )
            == (310, 220, 45, 45)
            and (
                self.nonempty_symbolic_ridge_groups,
                self.two_special_ridge_groups,
                self.z_ridge_groups,
                self.w_ridge_groups,
            )
            == (735, 495, 120, 120)
            and (
                self.y_ridge_incidences,
                self.z_ridge_incidences,
                self.w_ridge_incidences,
            )
            == (1980, 405, 405)
            and (self.selected_y_count, self.selected_z_count, self.selected_w_count)
            == (60, 25, 25)
            and self.z_cut_parts == EXPECTED_Z_CUT_PARTS
            and self.w_cut_parts == EXPECTED_W_CUT_PARTS
        )

    def to_record(self) -> dict[str, object]:
        return {
            "minimal_graph_covers": [list(value) for value in self.minimal_graph_covers],
            "graph_cover_five_set_count": self.graph_cover_five_set_count,
            "potential_facet_type_counts": {
                "y_12_13_plus_W_triple": self.potential_y_count,
                "z_0_1_13_plus_RI_pair": self.potential_z_count,
                "w_2_3_12_plus_LI_pair": self.potential_w_count,
            },
            "symbolic_ridge_system": {
                "nonempty_groups": self.nonempty_symbolic_ridge_groups,
                "two_special_groups": self.two_special_ridge_groups,
                "z_triple_groups": self.z_ridge_groups,
                "w_triple_groups": self.w_ridge_groups,
                "y_incidences": self.y_ridge_incidences,
                "z_incidences": self.z_ridge_incidences,
                "w_incidences": self.w_ridge_incidences,
            },
            "selected_facet_type_counts": {
                "y": self.selected_y_count,
                "z": self.selected_z_count,
                "w": self.selected_w_count,
            },
            "z_complete_bipartite_cut_parts": [list(part) for part in self.z_cut_parts],
            "w_complete_bipartite_cut_parts": [list(part) for part in self.w_cut_parts],
            "verified": self.verified,
        }


def audit_record18_structure() -> Record18StructuralAudit:
    """Replay the graph-cover types, complete ridge support, and two cuts."""

    graph_covers = minimal_graph_covers()
    graph_cover_vertices = tuple(vertices(value) for value in graph_covers)
    if graph_cover_vertices != EXPECTED_MINIMAL_GRAPH_COVERS:
        raise AssertionError("record-18 minimal graph covers changed")

    special = mask(SPECIAL)
    left = mask(LEFT)
    right = mask(RIGHT)
    possible_y: list[int] = []
    possible_z: list[int] = []
    possible_w: list[int] = []
    for chosen in combinations(range(TARGET_VERTEX_COUNT), 5):
        candidate = mask(chosen)
        if not is_transversal(candidate, EDGE_MASKS):
            continue
        if candidate & special == special:
            possible_y.append(candidate)
        elif (
            candidate & left == left
            and candidate & (1 << 13)
            and not candidate & (1 << 12)
        ):
            possible_z.append(candidate)
        elif (
            candidate & right == right
            and candidate & (1 << 12)
            and not candidate & (1 << 13)
        ):
            possible_w.append(candidate)
        else:
            raise AssertionError("unclassified record-18 graph-cover five-set")

    terms = tuple((value, "y") for value in possible_y)
    terms += tuple((value, "z") for value in possible_z)
    terms += tuple((value, "w") for value in possible_w)
    group_counts = Counter()
    incidence_counts = Counter()
    nonempty = 0
    for chosen in combinations(range(TARGET_VERTEX_COUNT), 6):
        six_set = mask(chosen)
        contained = tuple(kind for value, kind in terms if value & six_set == value)
        if not contained:
            continue
        nonempty += 1
        incidence_counts.update(contained)
        if six_set & special == special:
            group_counts["two_special"] += 1
        elif six_set & mask((0, 1, 13)) == mask((0, 1, 13)):
            group_counts["z"] += 1
        elif six_set & mask((2, 3, 12)) == mask((2, 3, 12)):
            group_counts["w"] += 1
        else:
            raise AssertionError("unclassified nonempty symbolic ridge group")

    selected = record18_facet_complements()
    selected_types = Counter()
    z_pairs = []
    w_pairs = []
    for candidate in selected:
        if candidate & special == special:
            selected_types["y"] += 1
        elif candidate in possible_z:
            selected_types["z"] += 1
            z_pairs.append(tuple(v for v in vertices(candidate) if v not in (0, 1, 13)))
        elif candidate in possible_w:
            selected_types["w"] += 1
            w_pairs.append(tuple(v for v in vertices(candidate) if v not in (2, 3, 12)))
        else:
            raise AssertionError("selected record-18 facet has no blocker type")

    audit = Record18StructuralAudit(
        minimal_graph_covers=graph_cover_vertices,
        graph_cover_five_set_count=len(possible_y) + len(possible_z) + len(possible_w),
        potential_y_count=len(possible_y),
        potential_z_count=len(possible_z),
        potential_w_count=len(possible_w),
        nonempty_symbolic_ridge_groups=nonempty,
        two_special_ridge_groups=group_counts["two_special"],
        z_ridge_groups=group_counts["z"],
        w_ridge_groups=group_counts["w"],
        y_ridge_incidences=incidence_counts["y"],
        z_ridge_incidences=incidence_counts["z"],
        w_ridge_incidences=incidence_counts["w"],
        selected_y_count=selected_types["y"],
        selected_z_count=selected_types["z"],
        selected_w_count=selected_types["w"],
        z_cut_parts=_complete_bipartite_parts(tuple(range(2, 12)), tuple(z_pairs)),
        w_cut_parts=_complete_bipartite_parts((0, 1, *INTERIOR), tuple(w_pairs)),
    )
    if not audit.verified:
        raise AssertionError("record-18 structural audit failed")
    return audit


@dataclass(frozen=True)
class Record18IncidenceAudit:
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
            and self.minimal_transversal_size_histogram == ((5, 110),)
            and len(self.facet_complement_masks) == 110
            and self.base_vertex_link_g3 == TARGET_BASE_LINK_G3
            and self.base_vertex_link_g4 == TARGET_BASE_LINK_G4
        )

    def to_record(self) -> dict[str, object]:
        return {
            "missing_edges": [list(edge) for edge in self.missing_edges],
            "missing_fives": [list(face) for face in self.missing_fives],
            "minimal_transversal_size_histogram": {
                str(size): count for size, count in self.minimal_transversal_size_histogram
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


def audit_record18_incidence(graph6: str = TARGET_GRAPH6) -> Record18IncidenceAudit:
    order, edges = decode_graph6(graph6)
    if order != TARGET_VERTEX_COUNT or edges != TARGET_MISSING_EDGES:
        raise ValueError("record-18 missing-edge graph changed")
    transversals = record18_minimal_transversals()
    histogram = tuple(sorted(Counter(value.bit_count() for value in transversals).items()))
    if histogram != ((5, 110),):
        raise AssertionError("record-18 blocker has a non-five minimal transversal")
    facets = tuple(sorted(transversals))
    facet_payload = b"".join(value.to_bytes(2, "little") for value in facets)
    face_table = generated_face_table(facets)
    packed = packed_face_table(face_table)
    base_g3, base_g4 = _base_vertex_link_g3_g4()
    audit = Record18IncidenceAudit(
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
        raise AssertionError("record-18 incidence audit failed")
    return audit


__all__ = [
    "Record18IncidenceAudit",
    "Record18StructuralAudit",
    "TARGET_F_VECTOR",
    "TARGET_GRAPH6",
    "TARGET_MISSING_EDGES",
    "TARGET_MISSING_FIVES",
    "TARGET_RECORD_INDEX",
    "audit_record18_incidence",
    "audit_record18_structure",
    "record18_clutter",
    "record18_facet_complements",
    "record18_minimal_transversals",
]
