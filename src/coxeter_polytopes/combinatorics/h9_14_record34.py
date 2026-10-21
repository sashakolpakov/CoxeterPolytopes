"""Exact positive-``a`` reduction for the H9/14 record-34 shadow.

The missing-edge graph is ``K3 + K2 + 9I``.  Dehn--Sommerville leaves
``a=0,1,2,3`` after the summed strict vertex-link ``g4`` inequality.  This
module proves that every branch with a missing four-face is impossible.

No SAT/SMT result is used.  The final ``a=1`` calculation is a labelled,
finite blocker census; retaining labels only makes it a deliberately
redundant exhaustive proof.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import reduce
from itertools import combinations
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record46 import (
    _a4_interval_feasible,
    _base_link_g3_g4,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 34
TARGET_GRAPH6 = "M?????????O?O?_@?"
TARGET_MISSING_EDGES = (
    (0, 11),
    (1, 12),
    (0, 13),
    (11, 13),
)
K3_VERTICES = (0, 11, 13)
K2_VERTICES = (1, 12)
ACTIVE_VERTICES = K3_VERTICES + K2_VERTICES
ISOLATED_VERTICES = tuple(range(2, 11))
CANONICAL_A4 = (2, 3, 4, 5)
TARGET_A1_FACET_COUNT = 112


def _mask(vertices: Sequence[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask & (1 << vertex))


def _base_transversals(
    edge_masks: tuple[int, ...], four: int, size: int
) -> tuple[int, ...]:
    return tuple(
        candidate
        for candidate in (
            _mask(vertices) for vertices in combinations(range(14), size)
        )
        if candidate & four and all(candidate & edge for edge in edge_masks)
    )


@dataclass(frozen=True)
class Record34PositiveAudit:
    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    base_vertex_link_g3: tuple[int, ...]
    base_vertex_link_g4: tuple[int, ...]
    schedule_link_g4_sums: tuple[tuple[int, int], ...]
    summed_link_excluded_branches: tuple[tuple[int, int, int], ...]
    active_a4_counts: tuple[tuple[int, int], ...]
    a4_shadow_histogram: tuple[tuple[int, int, int], ...]
    active_a4_single_link_g4_histogram: tuple[tuple[int, int], ...]
    k2_double_shadow_histogram: tuple[tuple[int, int], ...]
    a2_forced_shadow_histogram: tuple[tuple[int, int], ...]
    a2_scheduled_missing_fives: int
    a3_forced_shadow_histogram: tuple[tuple[int, int], ...]
    a3_scheduled_missing_fives: int
    canonical_a4_mask: int
    canonical_a4_shadow_size: int
    a1_residual_candidate_count: int
    a1_vertex_capacity_pair_count: int
    a1_no_small_transversal_pair_count: int
    a1_facet_count_histogram: tuple[tuple[int, int], ...]
    a1_required_facet_count: int

    @property
    def verified(self) -> bool:
        return (
            self.graph6 == TARGET_GRAPH6
            and self.missing_edges == TARGET_MISSING_EDGES
            and self.base_vertex_link_g3
            == (2, 3, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 3, 2)
            and self.base_vertex_link_g4 == self.base_vertex_link_g3
            and self.schedule_link_g4_sums
            == ((0, 36), (1, 30), (2, 24), (3, 18), (4, 12), (5, 6))
            and self.summed_link_excluded_branches
            == ((4, 12, 14), (5, 6, 14))
            and self.active_a4_counts == ((0, 126), (1, 420))
            and self.a4_shadow_histogram
            == ((0, 10, 126), (1, 8, 252), (1, 9, 168))
            and self.active_a4_single_link_g4_histogram == ((0, 420),)
            and self.k2_double_shadow_histogram == ((17, 1512), (18, 5460))
            and self.a2_forced_shadow_histogram == ((19, 1260), (20, 6615))
            and self.a2_scheduled_missing_fives == 18
            and self.a3_forced_shadow_histogram
            == ((27, 1680), (28, 16380), (29, 117180), (30, 190260))
            and self.a3_scheduled_missing_fives == 24
            and self.canonical_a4_mask == _mask(CANONICAL_A4)
            and self.canonical_a4_shadow_size == 10
            and self.a1_residual_candidate_count == 371
            and self.a1_vertex_capacity_pair_count == 15_645
            and self.a1_no_small_transversal_pair_count == 1_645
            and dict(self.a1_facet_count_histogram)
            == {90: 480, 96: 600, 102: 490, 108: 75}
            and self.a1_required_facet_count == TARGET_A1_FACET_COUNT
            and TARGET_A1_FACET_COUNT not in dict(self.a1_facet_count_histogram)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "graph6": self.graph6,
            "missing_edges": [list(edge) for edge in self.missing_edges],
            "missing_graph_decomposition": {
                "k3_vertices": list(K3_VERTICES),
                "k2_vertices": list(K2_VERTICES),
                "isolated_vertices": list(ISOLATED_VERTICES),
            },
            "base_vertex_link_g3": list(self.base_vertex_link_g3),
            "base_vertex_link_g4": list(self.base_vertex_link_g4),
            "schedule_link_g4_sums": [
                {"missing_four_count": a, "summed_vertex_link_g4": total}
                for a, total in self.schedule_link_g4_sums
            ],
            "summed_link_excluded_branches": [
                {
                    "missing_four_count": a,
                    "summed_vertex_link_g4": total,
                    "required_minimum": minimum,
                }
                for a, total, minimum in self.summed_link_excluded_branches
            ],
            "active_a4_obstruction": {
                "a4_counts_by_active_vertex_count": {
                    str(count): multiplicity
                    for count, multiplicity in self.active_a4_counts
                },
                "single_a4_shadow_histogram": [
                    {
                        "active_vertex_count": count,
                        "forced_shadow_size": shadow,
                        "multiplicity": multiplicity,
                    }
                    for count, shadow, multiplicity in self.a4_shadow_histogram
                ],
                "single_active_a4_link_g4_histogram": {
                    str(value): multiplicity
                    for value, multiplicity in (
                        self.active_a4_single_link_g4_histogram
                    )
                },
                "two_k2_incident_a4_shadow_union_histogram": {
                    str(size): multiplicity
                    for size, multiplicity in self.k2_double_shadow_histogram
                },
                "argument": (
                    "Strict link g3 allows a K3 vertex in at most one A4 "
                    "and a K2 vertex in at most two. One incident A4 forces "
                    "respectively 8 or 9 incident A5s; two distinct K2-"
                    "incident A4 shadows have union at least 17. The resulting "
                    "link g4 is nonpositive in every case. Hence every A4 "
                    "lies entirely in the nine-vertex isolate block."
                ),
            },
            "a2_a3_forced_shadow_obstruction": {
                "a2_forced_shadow_histogram": {
                    str(size): multiplicity
                    for size, multiplicity in self.a2_forced_shadow_histogram
                },
                "a2_scheduled_missing_fives": self.a2_scheduled_missing_fives,
                "a3_forced_shadow_histogram": {
                    str(size): multiplicity
                    for size, multiplicity in self.a3_forced_shadow_histogram
                },
                "a3_scheduled_missing_fives": self.a3_scheduled_missing_fives,
                "argument": (
                    "Each isolate A4 has ten independent A5 extensions. "
                    "Two distinct shadows intersect in at most one A5."
                ),
            },
            "a1_labelled_blocker_census": {
                "canonical_a4_mask": self.canonical_a4_mask,
                "canonical_a4_vertices": list(CANONICAL_A4),
                "forced_shadow_size": self.canonical_a4_shadow_size,
                "scheduled_missing_fives": 12,
                "additional_minimal_a5_count": 2,
                "residual_candidate_a5_count": (
                    self.a1_residual_candidate_count
                ),
                "vertex_capacity_pair_count": (
                    self.a1_vertex_capacity_pair_count
                ),
                "no_small_transversal_pair_count": (
                    self.a1_no_small_transversal_pair_count
                ),
                "facet_count_histogram": {
                    str(count): multiplicity
                    for count, multiplicity in self.a1_facet_count_histogram
                },
                "required_facet_count": self.a1_required_facet_count,
                "target_pair_count": 0,
            },
            "conclusion": "all branches with at least one missing A4 are impossible",
            "verified": self.verified,
        }


def audit_record34_positive_branches(
    graph6: str, schedules: Sequence[Mapping[str, object]]
) -> Record34PositiveAudit:
    """Recompute the complete ``a>0`` exclusion using integer arithmetic."""

    if graph6 != TARGET_GRAPH6:
        raise ValueError("certificate is scoped to H9/14 record 34")
    order, edges = decode_graph6(graph6)
    if order != 14 or tuple(edges) != TARGET_MISSING_EDGES:
        raise ValueError("record-34 missing-edge graph changed")
    vectors: dict[int, tuple[int, ...]] = {}
    missing_fives: dict[int, int] = {}
    for schedule in schedules:
        a = schedule.get("missing_four_count")
        vector = schedule.get("f_vector")
        offset = schedule.get("missing_five_count_offset")
        if not isinstance(a, int) or not isinstance(vector, (list, tuple)):
            raise ValueError("malformed record-34 schedule")
        if not isinstance(offset, int):
            raise ValueError("record-34 schedule lacks its missing-A5 count")
        vectors[a] = tuple(int(value) for value in vector)
        missing_fives[a] = offset
    if tuple(sorted(vectors)) != tuple(range(6)):
        raise ValueError("record-34 schedules must be a=0,...,5")

    base_g3, base_g4 = _base_link_g3_g4(graph6)
    link_sums = tuple(
        (a, summed_vertex_link_g4(vectors[a])) for a in sorted(vectors)
    )
    high_a = tuple(
        (a, total, 14) for a, total in link_sums if total < 14
    )

    independent_fours = independent_masks(graph6, 4)
    active_mask = _mask(ACTIVE_VERTICES)
    active_count_histogram = Counter(
        (four & active_mask).bit_count() for four in independent_fours
    )
    # There are also 216 independent A4s with two or more active vertices,
    # but the local A4 interval already rejects them.  The 546 interval-
    # feasible rows are exactly the two classes recorded below.
    local = analyze_local_link_graph(graph6)
    interval_feasible = tuple(
        four
        for four in independent_fours
        if _a4_interval_feasible(four, local)
    )
    expected_interval_feasible = tuple(
        four
        for four in independent_fours
        if (four & active_mask).bit_count() <= 1
    )
    if interval_feasible != expected_interval_feasible:
        raise AssertionError("record-34 local A4 interval classification changed")
    if Counter((four & active_mask).bit_count() for four in interval_feasible) != {
        0: 126,
        1: 420,
    }:
        raise AssertionError("record-34 interval-feasible A4 classes changed")

    fives = independent_masks(graph6, 5)
    five_indices = {five: index for index, five in enumerate(fives)}

    def shadow_word(candidate: int) -> int:
        return sum(
            1 << five_indices[five]
            for five in fives
            if five & candidate == candidate
        )

    shadow_words = {candidate: shadow_word(candidate) for candidate in interval_feasible}
    shadow_histogram = Counter(
        ((candidate & active_mask).bit_count(), shadow_words[candidate].bit_count())
        for candidate in interval_feasible
    )

    active_single_g4 = Counter()
    by_k2_vertex: dict[int, list[int]] = {vertex: [] for vertex in K2_VERTICES}
    for candidate in interval_feasible:
        active_part = candidate & active_mask
        if active_part.bit_count() != 1:
            continue
        vertex = _vertices(active_part)[0]
        incident_forced = sum(
            bool(five & (1 << vertex))
            for index, five in enumerate(fives)
            if shadow_words[candidate] & (1 << index)
        )
        active_single_g4[base_g4[vertex] + 6 - incident_forced] += 1
        if vertex in by_k2_vertex:
            by_k2_vertex[vertex].append(shadow_words[candidate])
    k2_double_shadow = Counter(
        (first | second).bit_count()
        for words in by_k2_vertex.values()
        for first, second in combinations(words, 2)
    )

    isolate_fours = tuple(
        candidate for candidate in interval_feasible if not candidate & active_mask
    )
    isolate_words = tuple(shadow_words[candidate] for candidate in isolate_fours)
    a2_shadow_histogram = Counter(
        (first | second).bit_count()
        for first, second in combinations(isolate_words, 2)
    )
    a3_shadow_histogram = Counter(
        reduce(int.__or__, words).bit_count()
        for words in combinations(isolate_words, 3)
    )

    four = _mask(CANONICAL_A4)
    shadow = tuple(five for five in fives if five & four == four)
    if len(shadow) != 10:
        raise AssertionError("canonical isolate A4 no longer has ten extensions")

    before_additional = tuple(
        base_g4[vertex]
        + 6 * int(bool(four & (1 << vertex)))
        - sum(bool(five & (1 << vertex)) for five in shadow)
        for vertex in range(14)
    )
    capacities = tuple(value - 1 for value in before_additional)
    k3_mask = _mask(K3_VERTICES)
    residual = tuple(
        five
        for five in fives
        if five & four != four and not five & k3_mask
    )
    capacity_pairs = tuple(
        (first, second)
        for first, second in combinations(residual, 2)
        if all(
            int(bool(first & (1 << vertex)))
            + int(bool(second & (1 << vertex)))
            <= capacities[vertex]
            for vertex in range(14)
        )
    )

    edge_masks = tuple(_mask(edge) for edge in edges)
    small_base = tuple(
        candidate
        for size in range(1, 5)
        for candidate in _base_transversals(edge_masks, four, size)
    )
    five_base = _base_transversals(edge_masks, four, 5)
    hit_words = {
        face: sum(
            bool(transversal & face) << index
            for index, transversal in enumerate(five_base)
        )
        for face in residual
    }
    histogram: Counter[int] = Counter()
    no_small = 0
    for first, second in capacity_pairs:
        if any(
            transversal & first and transversal & second
            for transversal in small_base
        ):
            continue
        no_small += 1
        # With no transversal of order <=4, every order-five transversal is
        # inclusion-minimal and therefore complements one candidate facet.
        facet_count = (hit_words[first] & hit_words[second]).bit_count()
        histogram[facet_count] += 1

    audit = Record34PositiveAudit(
        graph6=graph6,
        missing_edges=tuple(edges),
        base_vertex_link_g3=base_g3,
        base_vertex_link_g4=base_g4,
        schedule_link_g4_sums=link_sums,
        summed_link_excluded_branches=high_a,
        active_a4_counts=((0, 126), (1, 420)),
        a4_shadow_histogram=tuple(
            (active_count, shadow_size, multiplicity)
            for (active_count, shadow_size), multiplicity in sorted(
                shadow_histogram.items()
            )
        ),
        active_a4_single_link_g4_histogram=tuple(sorted(active_single_g4.items())),
        k2_double_shadow_histogram=tuple(sorted(k2_double_shadow.items())),
        a2_forced_shadow_histogram=tuple(sorted(a2_shadow_histogram.items())),
        a2_scheduled_missing_fives=missing_fives[2],
        a3_forced_shadow_histogram=tuple(sorted(a3_shadow_histogram.items())),
        a3_scheduled_missing_fives=missing_fives[3],
        canonical_a4_mask=four,
        canonical_a4_shadow_size=len(shadow),
        a1_residual_candidate_count=len(residual),
        a1_vertex_capacity_pair_count=len(capacity_pairs),
        a1_no_small_transversal_pair_count=no_small,
        a1_facet_count_histogram=tuple(sorted(histogram.items())),
        a1_required_facet_count=vectors[1][8],
    )
    if active_count_histogram != {0: 126, 1: 420, 2: 216}:
        raise AssertionError("record-34 independent A4 census changed")
    if not audit.verified:
        raise AssertionError("record-34 positive-branch audit did not reproduce")
    return audit


__all__ = [
    "TARGET_GRAPH6",
    "TARGET_RECORD_INDEX",
    "Record34PositiveAudit",
    "audit_record34_positive_branches",
]
