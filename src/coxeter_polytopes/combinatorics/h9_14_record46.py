"""Structural reduction for the solver-independent record-46 certificate."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    independent_masks,
    minimal_graph_vertex_covers,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_VERTEX_COUNT,
    decode_graph6,
)


TARGET_RECORD_INDEX = 46
TARGET_GRAPH6 = "M?????????O?o?o?_"
TARGET_MISSING_EDGES = (
    (0, 11),
    (0, 12),
    (1, 12),
    (0, 13),
    (1, 13),
    (12, 13),
)
ACTIVE_VERTICES = (0, 1, 11, 12, 13)
ISOLATED_VERTICES = tuple(range(2, 11))
MINIMAL_GRAPH_COVERS = (
    (0, 1, 12),
    (0, 1, 13),
    (0, 12, 13),
    (11, 12, 13),
)


def _mask(vertices: Sequence[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask & (1 << vertex))


def _base_link_g3_g4(graph6: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    independent = {size: independent_masks(graph6, size) for size in range(2, 6)}
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


def _a4_interval_feasible(face: int, local) -> bool:
    for figure in local.vertex_figures:
        incidence = int(bool(face & (1 << figure.vertex)))
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.edge_figures:
        subset = _mask(figure.edge)
        incidence = int(face & subset == subset)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.triple_figures:
        subset = _mask(figure.triple)
        if int(face & subset == subset) > figure.missing_four_incidence_max:
            return False
    return True


@dataclass(frozen=True)
class A4TypeAudit:
    active_vertices: tuple[int, ...]
    multiplicity: int
    representative_mask: int
    shadow_size: int
    nonpositive_g4_vertices: tuple[int, ...]

    @property
    def impossible(self) -> bool:
        return bool(self.nonpositive_g4_vertices)

    def to_record(self) -> dict[str, object]:
        return {
            "active_vertices": list(self.active_vertices),
            "multiplicity": self.multiplicity,
            "representative_mask": self.representative_mask,
            "representative_vertices": list(_vertices(self.representative_mask)),
            "forced_five_shadow_size": self.shadow_size,
            "nonpositive_g4_vertices": list(self.nonpositive_g4_vertices),
            "impossible": self.impossible,
        }


def audit_a4_types(graph6: str) -> tuple[A4TypeAudit, ...]:
    local = analyze_local_link_graph(graph6)
    active_mask = _mask(ACTIVE_VERTICES)
    candidates = tuple(
        face
        for face in independent_masks(graph6, 4)
        if _a4_interval_feasible(face, local)
    )
    grouped: dict[tuple[int, ...], list[int]] = {}
    for face in candidates:
        signature = _vertices(face & active_mask)
        grouped.setdefault(signature, []).append(face)
    fives = set(independent_masks(graph6, 5))
    _, base_g4 = _base_link_g3_g4(graph6)
    result = []
    for signature, faces in sorted(grouped.items()):
        representative = min(faces)
        shadow = tuple(
            face for face in fives if face & representative == representative
        )
        high = tuple(
            base_g4[vertex]
            + 6 * int(bool(representative & (1 << vertex)))
            - sum(bool(face & (1 << vertex)) for face in shadow)
            for vertex in range(14)
        )
        result.append(
            A4TypeAudit(
                active_vertices=signature,
                multiplicity=len(faces),
                representative_mask=representative,
                shadow_size=len(shadow),
                nonpositive_g4_vertices=tuple(
                    vertex for vertex, value in enumerate(high) if value <= 0
                ),
            )
        )
    return tuple(result)


@dataclass(frozen=True)
class Record46StructureAudit:
    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    minimal_graph_covers: tuple[int, ...]
    base_vertex_link_g3: tuple[int, ...]
    base_vertex_link_g4: tuple[int, ...]
    schedule_link_g4_sums: tuple[tuple[int, int], ...]
    excluded_high_a_branches: tuple[tuple[int, int, int], ...]
    a4_types: tuple[A4TypeAudit, ...]
    a0_missing_five_count: int
    a0_required_facet_count: int
    a0_forbidden_vertices: tuple[int, ...]
    a0_individually_allowed_a5_count: int

    @property
    def verified(self) -> bool:
        expected_base = (1, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 1, 1)
        return (
            self.graph6 == TARGET_GRAPH6
            and self.missing_edges == TARGET_MISSING_EDGES
            and self.minimal_graph_covers
            == tuple(sorted(_mask(cover) for cover in MINIMAL_GRAPH_COVERS))
            and self.base_vertex_link_g3 == expected_base
            and self.base_vertex_link_g4 == expected_base
            and self.schedule_link_g4_sums
            == ((0, 24), (1, 18), (2, 12), (3, 6))
            and self.excluded_high_a_branches == ((2, 12, 14), (3, 6, 14))
            and [item.active_vertices for item in self.a4_types]
            == [(), (1,), (11,)]
            and [item.multiplicity for item in self.a4_types] == [126, 84, 84]
            and [item.shadow_size for item in self.a4_types] == [10, 8, 9]
            and all(item.impossible for item in self.a4_types)
            and self.a0_missing_five_count == 4
            and self.a0_required_facet_count == 90
            and self.a0_forbidden_vertices == (0, 12, 13)
            and self.a0_individually_allowed_a5_count == 378
        )

    def to_record(self) -> dict[str, object]:
        return {
            "graph6": self.graph6,
            "missing_edges": [list(edge) for edge in self.missing_edges],
            "active_vertices": list(ACTIVE_VERTICES),
            "isolated_vertices": list(ISOLATED_VERTICES),
            "minimal_graph_cover_masks": list(self.minimal_graph_covers),
            "minimal_graph_cover_vertices": [
                list(_vertices(mask)) for mask in self.minimal_graph_covers
            ],
            "base_vertex_link_g3": list(self.base_vertex_link_g3),
            "base_vertex_link_g4": list(self.base_vertex_link_g4),
            "schedule_link_g4_sums": [
                {"missing_four_count": a, "summed_vertex_link_g4": total}
                for a, total in self.schedule_link_g4_sums
            ],
            "excluded_high_a_branches": [
                {
                    "missing_four_count": a,
                    "summed_vertex_link_g4": total,
                    "required_minimum": minimum,
                }
                for a, total, minimum in self.excluded_high_a_branches
            ],
            "a1": {
                "independent_a4_count": 690,
                "interval_feasible_a4_count": sum(
                    item.multiplicity for item in self.a4_types
                ),
                "interval_feasible_a4_orbit_count": len(self.a4_types),
                "types": [item.to_record() for item in self.a4_types],
                "shadow_strict_survivor_count": 0,
                "conclusion": (
                    "all three interval-feasible A4 types force a nonpositive "
                    "vertex-link g4 before any additional minimal A5"
                ),
            },
            "a0": {
                "minimal_a5_count": self.a0_missing_five_count,
                "required_facet_count": self.a0_required_facet_count,
                "forbidden_vertices": list(self.a0_forbidden_vertices),
                "allowed_support_vertices": [1, 11, *ISOLATED_VERTICES],
                "pair_1_11_a5_incidence_capacity": 0,
                "individually_allowed_a5_count": (
                    self.a0_individually_allowed_a5_count
                ),
            },
            "verified": self.verified,
        }


def audit_record46_structure(
    graph6: str, schedules: Sequence[Mapping[str, object]]
) -> Record46StructureAudit:
    if graph6 != TARGET_GRAPH6:
        raise ValueError("certificate is scoped to H9/14 record 46")
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != 14 or tuple(edges) != TARGET_MISSING_EDGES:
        raise ValueError("record-46 missing graph changed")
    parsed = {}
    missing_fives = {}
    for schedule in schedules:
        a = schedule.get("missing_four_count")
        vector = schedule.get("f_vector")
        offset = schedule.get("missing_five_count_offset")
        if not isinstance(a, int) or not isinstance(vector, (list, tuple)):
            raise ValueError("malformed record-46 schedule")
        if not isinstance(offset, int):
            raise ValueError("record-46 schedule lacks A5 count")
        parsed[a] = tuple(int(value) for value in vector)
        missing_fives[a] = offset
    if tuple(sorted(parsed)) != (0, 1, 2, 3):
        raise ValueError("record-46 schedules must be a=0,1,2,3")
    base_g3, base_g4 = _base_link_g3_g4(graph6)
    link_sums = tuple((a, summed_vertex_link_g4(parsed[a])) for a in sorted(parsed))
    forbidden = tuple(vertex for vertex, value in enumerate(base_g4) if value == 1)
    forbidden_mask = _mask(forbidden)
    allowed_rows = tuple(
        face
        for face in independent_masks(graph6, 5)
        if not face & forbidden_mask
        and face & _mask((1, 11)) != _mask((1, 11))
    )
    audit = Record46StructureAudit(
        graph6=graph6,
        missing_edges=tuple(edges),
        minimal_graph_covers=minimal_graph_vertex_covers(graph6),
        base_vertex_link_g3=base_g3,
        base_vertex_link_g4=base_g4,
        schedule_link_g4_sums=link_sums,
        excluded_high_a_branches=tuple(
            (a, total, TARGET_VERTEX_COUNT)
            for a, total in link_sums
            if total < TARGET_VERTEX_COUNT
        ),
        a4_types=audit_a4_types(graph6),
        a0_missing_five_count=missing_fives[0],
        a0_required_facet_count=parsed[0][8],
        a0_forbidden_vertices=forbidden,
        a0_individually_allowed_a5_count=len(allowed_rows),
    )
    if not audit.verified:
        raise AssertionError("record-46 structural reduction did not reproduce")
    return audit
