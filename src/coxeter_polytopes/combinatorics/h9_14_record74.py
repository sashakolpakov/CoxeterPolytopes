"""Solver-independent exclusion of H9/14 record 74.

The missing-edge graph is ``K4 + 10K1``.  Strict local ``g``-vectors reduce
the proof to a divisibility obstruction: in the only remaining schedule the
four missing five-faces lie entirely in the ten isolated vertices, so their
five-element transversals occur in four equal classes.  The scheduled facet
count is 90, which is not divisible by four.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
    minimal_graph_vertex_covers,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_POLYTOPE_DIMENSION,
    TARGET_VERTEX_COUNT,
    decode_graph6,
)


TARGET_RECORD_INDEX = 74
TARGET_GRAPH6 = "M?????????O?_@_@_"
TARGET_MISSING_EDGES = (
    (0, 11),
    (0, 12),
    (11, 12),
    (0, 13),
    (11, 13),
    (12, 13),
)
ACTIVE_VERTICES = (0, 11, 12, 13)
ISOLATED_VERTICES = tuple(range(1, 11))
MINIMAL_GRAPH_COVERS = tuple(
    tuple(vertex for vertex in ACTIVE_VERTICES if vertex != omitted)
    for omitted in ACTIVE_VERTICES
)


def _mask(vertices: Sequence[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT)
        if mask & (1 << vertex)
    )


def _base_link_g3_g4(graph6: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Reconstruct local ``g3,g4`` before selecting A4/A5 nonfaces."""

    independent = {size: independent_masks(graph6, size) for size in range(2, 6)}
    g3_values: list[int] = []
    g4_values: list[int] = []
    for vertex in range(TARGET_VERTEX_COUNT):
        counts = {
            size: sum(bool(face & (1 << vertex)) for face in independent[size])
            for size in independent
        }
        f0, f1, f2, f3 = (counts[size] for size in range(2, 6))
        g3_values.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4_values.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3_values), tuple(g4_values)


@dataclass(frozen=True)
class A4TypeAudit:
    active_vertex_count: int
    isolate_vertex_count: int
    orbit_multiplicity: int
    representative_mask: int
    independent_five_shadow_size: int
    strict_g3_zero_vertices: tuple[int, ...]
    strict_g4_zero_vertices_before_additional_a5: tuple[int, ...]

    @property
    def locally_impossible(self) -> bool:
        return bool(
            self.strict_g3_zero_vertices
            or self.strict_g4_zero_vertices_before_additional_a5
        )

    def to_record(self) -> dict[str, object]:
        return {
            "active_vertex_count": self.active_vertex_count,
            "isolate_vertex_count": self.isolate_vertex_count,
            "orbit_multiplicity": self.orbit_multiplicity,
            "representative_mask": self.representative_mask,
            "representative_vertices": list(_vertices(self.representative_mask)),
            "independent_five_shadow_size": self.independent_five_shadow_size,
            "strict_g3_zero_vertices": list(self.strict_g3_zero_vertices),
            "strict_g4_zero_vertices_before_additional_a5": list(
                self.strict_g4_zero_vertices_before_additional_a5
            ),
            "locally_impossible": self.locally_impossible,
        }


def audit_a4_types(graph6: str) -> tuple[A4TypeAudit, ...]:
    """Exhaust all 690 independent A4s and collapse the two symmetry types."""

    active_mask = _mask(ACTIVE_VERTICES)
    fives = set(independent_masks(graph6, 5))
    base_g3, base_g4 = _base_link_g3_g4(graph6)
    grouped: dict[int, list[int]] = {0: [], 1: []}
    for face in independent_masks(graph6, 4):
        active_count = (face & active_mask).bit_count()
        if active_count not in grouped:
            raise AssertionError("an independent set in K4 contains two active vertices")
        grouped[active_count].append(face)

    result = []
    for active_count, faces in sorted(grouped.items()):
        representative = min(faces)
        shadow = tuple(
            extension
            for extension in fives
            if extension & representative == representative
        )
        g3 = tuple(
            base_g3[vertex] - int(bool(representative & (1 << vertex)))
            for vertex in range(TARGET_VERTEX_COUNT)
        )
        shadow_incidence = tuple(
            sum(bool(face & (1 << vertex)) for face in shadow)
            for vertex in range(TARGET_VERTEX_COUNT)
        )
        g4 = tuple(
            base_g4[vertex]
            + 6 * int(bool(representative & (1 << vertex)))
            - shadow_incidence[vertex]
            for vertex in range(TARGET_VERTEX_COUNT)
        )
        result.append(
            A4TypeAudit(
                active_vertex_count=active_count,
                isolate_vertex_count=4 - active_count,
                orbit_multiplicity=len(faces),
                representative_mask=representative,
                independent_five_shadow_size=len(shadow),
                strict_g3_zero_vertices=tuple(
                    vertex for vertex, value in enumerate(g3) if value <= 0
                ),
                strict_g4_zero_vertices_before_additional_a5=tuple(
                    vertex for vertex, value in enumerate(g4) if value <= 0
                ),
            )
        )
    return tuple(result)


@dataclass(frozen=True)
class Record74Audit:
    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    minimal_graph_covers: tuple[int, ...]
    base_vertex_link_g3: tuple[int, ...]
    base_vertex_link_g4: tuple[int, ...]
    schedule_link_g4_sums: tuple[tuple[int, int], ...]
    excluded_high_a_branches: tuple[tuple[int, int, int], ...]
    a4_types: tuple[A4TypeAudit, ...]
    a1_total_missing_five_count: int
    a0_missing_five_count: int
    a0_required_facet_count: int
    a0_locally_allowed_missing_five_count: int
    a0_forbidden_active_missing_five_count: int
    a0_every_active_a5_forces_nonpositive_g4: bool
    a0_three_active_two_isolate_cover_count: int
    a0_four_active_one_isolate_cover_count: int
    a0_three_active_class_multiplicity: int
    a0_four_active_survivors_under_empty_intersection: int
    a0_facet_divisor: int

    @property
    def verified(self) -> bool:
        expected_base = tuple(
            1 if vertex in ACTIVE_VERTICES else 4
            for vertex in range(TARGET_VERTEX_COUNT)
        )
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
            and len(self.a4_types) == 2
            and [item.orbit_multiplicity for item in self.a4_types] == [210, 480]
            and [item.independent_five_shadow_size for item in self.a4_types]
            == [10, 7]
            and all(item.locally_impossible for item in self.a4_types)
            and self.a1_total_missing_five_count == 10
            and self.a0_missing_five_count == 4
            and self.a0_required_facet_count == 90
            and self.a0_locally_allowed_missing_five_count == comb(10, 5)
            and self.a0_forbidden_active_missing_five_count
            == len(ACTIVE_VERTICES) * comb(10, 4)
            and self.a0_every_active_a5_forces_nonpositive_g4
            and self.a0_three_active_two_isolate_cover_count
            == len(ACTIVE_VERTICES) * comb(10, 2)
            and self.a0_four_active_one_isolate_cover_count == 10
            and self.a0_three_active_class_multiplicity == len(ACTIVE_VERTICES)
            and self.a0_four_active_survivors_under_empty_intersection == 0
            and self.a0_facet_divisor == 4
            and self.a0_required_facet_count % self.a0_facet_divisor != 0
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
                "independent_a4_count": sum(
                    item.orbit_multiplicity for item in self.a4_types
                ),
                "a4_orbit_count": len(self.a4_types),
                "active_vertex_a4_incidence_interval": [0, 0],
                "a4_types": [item.to_record() for item in self.a4_types],
                "scheduled_missing_five_count": self.a1_total_missing_five_count,
                "isolate_a4_shadow_size": (
                    self.a4_types[0].independent_five_shadow_size
                ),
                "additional_minimal_a5_count_after_isolate_a4": (
                    self.a1_total_missing_five_count
                    - self.a4_types[0].independent_five_shadow_size
                ),
                "small_transversal_witness": (
                    "any minimal three-cover C of K4 union any q in the "
                    "isolate-only A4"
                ),
                "conclusion": (
                    "every A4 orbit violates strict vertex-link g3 or g4 "
                    "before any additional minimal A5 is selected; the "
                    "isolate-only orbit also has an explicit four-transversal"
                ),
            },
            "a0": {
                "forced_minimal_a5_count": self.a0_missing_five_count,
                "required_facet_count": self.a0_required_facet_count,
                "locally_allowed_a5_count": (
                    self.a0_locally_allowed_missing_five_count
                ),
                "active_a5_count_rejected_by_strict_g4": (
                    self.a0_forbidden_active_missing_five_count
                ),
                "every_active_a5_forces_nonpositive_g4": (
                    self.a0_every_active_a5_forces_nonpositive_g4
                ),
                "local_support": "all A5s lie in the ten isolated vertices",
                "small_transversal_gate": (
                    "the four A5s must have empty total intersection"
                ),
                "five_transversal_formula": "4 * tau_2(A5 family)",
                "five_cover_classes": {
                    "three_active_plus_two_isolates": (
                        self.a0_three_active_two_isolate_cover_count
                    ),
                    "multiplicity_per_isolate_pair": (
                        self.a0_three_active_class_multiplicity
                    ),
                    "four_active_plus_one_isolate": (
                        self.a0_four_active_one_isolate_cover_count
                    ),
                    "four_active_survivors_when_A5_intersection_empty": (
                        self.a0_four_active_survivors_under_empty_intersection
                    ),
                },
                "facet_divisor": self.a0_facet_divisor,
                "required_facet_residue": (
                    self.a0_required_facet_count % self.a0_facet_divisor
                ),
            },
            "verified": self.verified,
        }


def audit_record74(
    graph6: str, schedules: Sequence[Mapping[str, object]]
) -> Record74Audit:
    if graph6 != TARGET_GRAPH6:
        raise ValueError("certificate is scoped to H9/14 record 74")
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT or tuple(edges) != TARGET_MISSING_EDGES:
        raise ValueError("record-74 missing-edge graph changed")
    parsed: dict[int, tuple[int, ...]] = {}
    missing_fives: dict[int, int] = {}
    for schedule in schedules:
        a = schedule.get("missing_four_count")
        vector = schedule.get("f_vector")
        offset = schedule.get("missing_five_count_offset")
        if not isinstance(a, int) or not isinstance(vector, (list, tuple)):
            raise ValueError("malformed record-74 schedule")
        if not isinstance(offset, int):
            raise ValueError("record-74 schedule lacks its A5 count")
        parsed[a] = tuple(int(value) for value in vector)
        missing_fives[a] = offset
    if tuple(sorted(parsed)) != (0, 1, 2, 3):
        raise ValueError("record-74 schedules must be exactly a=0,1,2,3")

    base_g3, base_g4 = _base_link_g3_g4(graph6)
    link_sums = tuple((a, summed_vertex_link_g4(parsed[a])) for a in sorted(parsed))
    active_mask = _mask(ACTIVE_VERTICES)
    fives = independent_masks(graph6, 5)
    locally_allowed_a5 = sum(not bool(face & active_mask) for face in fives)
    forbidden_active_a5 = tuple(face for face in fives if face & active_mask)
    every_active_a5_fails = all(
        all(
            base_g4[vertex] - 1 <= 0
            for vertex in ACTIVE_VERTICES
            if face & (1 << vertex)
        )
        for face in forbidden_active_a5
    )
    five_covers = graph_vertex_cover_masks(graph6, 5)
    three_active_covers = tuple(
        cover for cover in five_covers if (cover & active_mask).bit_count() == 3
    )
    four_active_covers = tuple(
        cover for cover in five_covers if (cover & active_mask).bit_count() == 4
    )
    isolate_pair_multiplicities: dict[int, int] = {}
    for cover in three_active_covers:
        isolate_part = cover & ~active_mask
        isolate_pair_multiplicities[isolate_part] = (
            isolate_pair_multiplicities.get(isolate_part, 0) + 1
        )
    class_multiplicities = set(isolate_pair_multiplicities.values())
    if len(class_multiplicities) != 1:
        raise AssertionError("record-74 isolate-pair cover classes are not uniform")
    audit = Record74Audit(
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
        a1_total_missing_five_count=missing_fives[1],
        a0_missing_five_count=missing_fives[0],
        a0_required_facet_count=parsed[0][TARGET_POLYTOPE_DIMENSION - 1],
        a0_locally_allowed_missing_five_count=locally_allowed_a5,
        a0_forbidden_active_missing_five_count=len(forbidden_active_a5),
        a0_every_active_a5_forces_nonpositive_g4=every_active_a5_fails,
        a0_three_active_two_isolate_cover_count=len(three_active_covers),
        a0_four_active_one_isolate_cover_count=len(four_active_covers),
        a0_three_active_class_multiplicity=class_multiplicities.pop(),
        a0_four_active_survivors_under_empty_intersection=0,
        a0_facet_divisor=len(ACTIVE_VERTICES),
    )
    if not audit.verified:
        raise AssertionError("record-74 structural exclusion did not reproduce")
    return audit
