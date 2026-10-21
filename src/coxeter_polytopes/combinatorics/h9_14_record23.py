"""Exact structural reduction for the H9/14 record-23 certificate.

The missing-edge graph is ``K_{3,2} + 9K_1``.  This module performs the
small, self-contained part of the proof: Dehn--Sommerville/link ``g4``
pruning and a complete orbit audit of up to three missing four-faces.  The
larger blocker incidence censuses are implemented by the companion C++
enumerators and frozen by ``scripts/certify_h9_14_record23.py``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import permutations
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    independent_masks,
    minimal_graph_vertex_covers,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_VERTEX_COUNT,
    decode_graph6,
    macaulay_raise,
)


TARGET_RECORD_INDEX = 23
TARGET_GRAPH6 = "M???????????w?w??"
TARGET_MISSING_EDGES = (
    (0, 12),
    (1, 12),
    (2, 12),
    (0, 13),
    (1, 13),
    (2, 13),
)
LEFT_VERTICES = (0, 1, 2)
RIGHT_VERTICES = (12, 13)
ISOLATED_VERTICES = tuple(range(3, 12))
USABLE_VERTICES = LEFT_VERTICES + ISOLATED_VERTICES
VERTEX_CLASSES = (LEFT_VERTICES, RIGHT_VERTICES, ISOLATED_VERTICES)
MINIMAL_GRAPH_COVERS = (LEFT_VERTICES, RIGHT_VERTICES)


def _mask(vertices: Sequence[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(
        vertex for vertex in range(TARGET_VERTEX_COUNT)
        if mask & (1 << vertex)
    )


def _base_vertex_link_g3_g4(graph6: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
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


def _canonical_a4_key(family: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Canonicalize under ``S3 x S2 x S9`` and missing-face permutations."""

    row_count = len(family)
    best: tuple[tuple[int, ...], ...] | None = None
    for row_permutation in permutations(range(row_count)):
        candidate = tuple(
            tuple(
                sorted(
                    sum(
                        bool(family[row_permutation[row]] & (1 << vertex))
                        << row
                        for row in range(row_count)
                    )
                    for vertex in vertex_class
                )
            )
            for vertex_class in VERTEX_CLASSES
        )
        if best is None or candidate < best:
            best = candidate
    if best is None:
        return ((), (), ())
    return best


def _family_from_key(key: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    row_count = max(
        (column for vertex_class in key for column in vertex_class),
        default=0,
    ).bit_length()
    return tuple(
        sum(
            bool(column & (1 << row)) << vertex
            for vertex_class, vertices in zip(key, VERTEX_CLASSES)
            for column, vertex in zip(vertex_class, vertices)
        )
        for row in range(row_count)
    )


@dataclass(frozen=True)
class A4OrbitLayer:
    missing_four_count: int
    all_orbit_count: int
    vertex_interval_feasible_orbit_count: int
    feasible_minimal_five_shadow_histogram: tuple[tuple[tuple[int, int], int], ...]
    all_feasible_fours_meet_left_and_avoid_right: bool

    def to_record(self) -> dict[str, object]:
        return {
            "missing_four_count": self.missing_four_count,
            "all_orbit_count": self.all_orbit_count,
            "vertex_interval_feasible_orbit_count": (
                self.vertex_interval_feasible_orbit_count
            ),
            "feasible_minimal_five_shadow_histogram": [
                {
                    "minimal_missing_five_count": key[0],
                    "forced_five_shadow_size": key[1],
                    "orbit_count": multiplicity,
                }
                for key, multiplicity in self.feasible_minimal_five_shadow_histogram
            ],
            "all_feasible_fours_meet_left_and_avoid_right": (
                self.all_feasible_fours_meet_left_and_avoid_right
            ),
        }


def enumerate_a4_orbit_layers(graph6: str) -> tuple[A4OrbitLayer, ...]:
    """Audit every A4 orbit through order three under the full graph symmetry."""

    independent_fours = independent_masks(graph6, 4)
    independent_fives = independent_masks(graph6, 5)
    five_index = {face: index for index, face in enumerate(independent_fives)}
    shadows = {
        face: sum(
            1 << five_index[extension]
            for extension in independent_fives
            if extension & face == face
        )
        for face in independent_fours
    }
    base_g3, base_g4 = _base_vertex_link_g3_g4(graph6)
    levels: dict[int, set[tuple[tuple[int, ...], ...]]] = {
        0: {((), (), ())}
    }
    result: list[A4OrbitLayer] = []

    for missing_four_count in range(4):
        if missing_four_count:
            next_level: set[tuple[tuple[int, ...], ...]] = set()
            for key in levels[missing_four_count - 1]:
                family = _family_from_key(key)
                for face in independent_fours:
                    if face not in family:
                        next_level.add(_canonical_a4_key(family + (face,)))
            levels[missing_four_count] = next_level

        feasible: list[tuple[tuple[tuple[int, ...], ...], int, int]] = []
        for key in levels[missing_four_count]:
            family = _family_from_key(key)
            shadow_bits = 0
            for face in family:
                shadow_bits |= shadows[face]
            shadow_faces = tuple(
                independent_fives[index]
                for index in range(len(independent_fives))
                if shadow_bits & (1 << index)
            )
            shadow_size = len(shadow_faces)
            minimal_five_count = 7 + 6 * missing_four_count - shadow_size
            if minimal_five_count < 0:
                continue
            four_incidence = tuple(
                sum(bool(face & (1 << vertex)) for face in family)
                for vertex in range(TARGET_VERTEX_COUNT)
            )
            shadow_incidence = tuple(
                sum(bool(face & (1 << vertex)) for face in shadow_faces)
                for vertex in range(TARGET_VERTEX_COUNT)
            )
            lower: list[int] = []
            upper: list[int] = []
            possible = True
            for vertex in range(TARGET_VERTEX_COUNT):
                previous = base_g3[vertex] - four_incidence[vertex]
                if previous < 0:
                    possible = False
                    break
                before_chosen_fives = (
                    base_g4[vertex]
                    + 6 * four_incidence[vertex]
                    - shadow_incidence[vertex]
                )
                lower.append(
                    max(0, before_chosen_fives - macaulay_raise(previous, 3))
                )
                upper.append(before_chosen_fives - 1)
            if not possible or min(upper, default=0) < 0:
                continue
            if not sum(lower) <= 5 * minimal_five_count <= sum(upper):
                continue
            feasible.append((key, minimal_five_count, shadow_size))

        left_mask = _mask(LEFT_VERTICES)
        right_mask = _mask(RIGHT_VERTICES)
        correct_support = all(
            all((face & left_mask) and not (face & right_mask)
                for face in _family_from_key(key))
            for key, _, _ in feasible
        )
        histogram = Counter(
            (minimal_five_count, shadow_size)
            for _, minimal_five_count, shadow_size in feasible
        )
        result.append(
            A4OrbitLayer(
                missing_four_count=missing_four_count,
                all_orbit_count=len(levels[missing_four_count]),
                vertex_interval_feasible_orbit_count=len(feasible),
                feasible_minimal_five_shadow_histogram=tuple(
                    sorted(histogram.items())
                ),
                all_feasible_fours_meet_left_and_avoid_right=correct_support,
            )
        )
    return tuple(result)


@dataclass(frozen=True)
class Record23StructureAudit:
    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    minimal_graph_covers: tuple[int, ...]
    schedule_link_g4_sums: tuple[tuple[int, int], ...]
    excluded_high_a_branches: tuple[tuple[int, int, int], ...]
    base_vertex_link_g3: tuple[int, ...]
    base_vertex_link_g4: tuple[int, ...]
    a4_orbit_layers: tuple[A4OrbitLayer, ...]

    @property
    def verified(self) -> bool:
        layers = self.a4_orbit_layers
        return (
            self.graph6 == TARGET_GRAPH6
            and self.missing_edges == TARGET_MISSING_EDGES
            and self.minimal_graph_covers
            == tuple(_mask(cover) for cover in MINIMAL_GRAPH_COVERS)
            and self.schedule_link_g4_sums
            == ((0, 36), (1, 30), (2, 24), (3, 18), (4, 12), (5, 6))
            and self.excluded_high_a_branches == ((4, 12, 14), (5, 6, 14))
            and self.base_vertex_link_g3
            == (4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1)
            and self.base_vertex_link_g4
            == (5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1)
            and [layer.all_orbit_count for layer in layers] == [1, 6, 72, 1_357]
            and [layer.vertex_interval_feasible_orbit_count for layer in layers]
            == [1, 3, 23, 86]
            and [dict(layer.feasible_minimal_five_shadow_histogram) for layer in layers]
            == [
                {(7, 0): 1},
                {(5, 8): 3},
                {(3, 16): 16, (4, 15): 7},
                {(1, 24): 7, (2, 23): 51, (3, 22): 21, (4, 21): 7},
            ]
            and all(
                layer.all_feasible_fours_meet_left_and_avoid_right
                for layer in layers
            )
        )

    def to_record(self) -> dict[str, object]:
        return {
            "graph6": self.graph6,
            "missing_edges": [list(edge) for edge in self.missing_edges],
            "minimal_graph_cover_masks": list(self.minimal_graph_covers),
            "minimal_graph_cover_vertices": [
                list(_vertices(mask)) for mask in self.minimal_graph_covers
            ],
            "schedule_link_g4_sums": [
                {"missing_four_count": a, "summed_vertex_link_g4": total}
                for a, total in self.schedule_link_g4_sums
            ],
            "strict_vertex_link_g4_required_sum": 14,
            "excluded_high_a_branches": [
                {
                    "missing_four_count": a,
                    "summed_vertex_link_g4": total,
                    "required_minimum": minimum,
                }
                for a, total, minimum in self.excluded_high_a_branches
            ],
            "base_vertex_link_g3": list(self.base_vertex_link_g3),
            "base_vertex_link_g4": list(self.base_vertex_link_g4),
            "a4_orbit_layers": [layer.to_record() for layer in self.a4_orbit_layers],
            "verified": self.verified,
        }


def audit_record23_structure(
    graph6: str, schedules: Sequence[Mapping[str, object]]
) -> Record23StructureAudit:
    if graph6 != TARGET_GRAPH6:
        raise ValueError("record-23 certificate is scoped to one graph shadow")
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT or tuple(edges) != TARGET_MISSING_EDGES:
        raise ValueError("record-23 missing-edge graph changed")
    parsed: dict[int, tuple[int, ...]] = {}
    for schedule in schedules:
        a = schedule.get("missing_four_count")
        vector = schedule.get("f_vector")
        if not isinstance(a, int) or not isinstance(vector, (list, tuple)):
            raise ValueError("malformed record-23 schedule")
        parsed[a] = tuple(int(value) for value in vector)
    if tuple(sorted(parsed)) != tuple(range(6)):
        raise ValueError("record-23 schedules must be exactly a=0,...,5")
    link_sums = tuple((a, summed_vertex_link_g4(parsed[a])) for a in sorted(parsed))
    base_g3, base_g4 = _base_vertex_link_g3_g4(graph6)
    audit = Record23StructureAudit(
        graph6=graph6,
        missing_edges=tuple(edges),
        minimal_graph_covers=minimal_graph_vertex_covers(graph6),
        schedule_link_g4_sums=link_sums,
        excluded_high_a_branches=tuple(
            (a, total, TARGET_VERTEX_COUNT)
            for a, total in link_sums
            if a >= 4
        ),
        base_vertex_link_g3=base_g3,
        base_vertex_link_g4=base_g4,
        a4_orbit_layers=enumerate_a4_orbit_layers(graph6),
    )
    if not audit.verified:
        raise AssertionError("record-23 structural audit did not reproduce")
    return audit
