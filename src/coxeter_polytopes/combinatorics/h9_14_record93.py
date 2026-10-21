"""Solver-independent exclusion certificate for H9/14 record 93.

This is the sibling of :mod:`h9_14_record92`.  The same vertex-link ``g_4``
identity removes the schedules with missing four-faces.  In the remaining
branch, a complete blocker enumeration gives 54, 57, 63, or 75 facet
complements, never the required 78.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    TARGET_LINK_COUNT,
    TransversalFamily,
    five_transversal_count,
    graph_vertex_cover_masks,
    has_transversal_of_size_at_most_four,
    hit_pattern,
    independent_masks,
    minimal_graph_vertex_covers,
    pattern_compatible_missing_five_triples,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_VERTEX_COUNT,
    decode_graph6,
)


TARGET_RECORD_INDEX = 93
TARGET_GRAPH6 = "M?????????W?o?o@?"
TARGET_MISSING_EDGES = (
    (0, 11),
    (1, 11),
    (0, 12),
    (1, 12),
    (0, 13),
    (1, 13),
    (11, 13),
)
MINIMAL_GRAPH_COVERS = (
    (0, 1, 11),
    (0, 1, 13),
    (11, 12, 13),
)
EXPECTED_PATTERN_COUNTS = {
    0b000: 126,
    0b011: 336,
    0b100: 126,
    0b101: 210,
    0b110: 210,
}


def _mask(vertices: Sequence[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(
        vertex
        for vertex in range(TARGET_VERTEX_COUNT)
        if mask & (1 << vertex)
    )


@dataclass(frozen=True)
class Record93Audit:
    """Complete arithmetic certificate for the three record-93 schedules."""

    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    minimal_graph_covers: tuple[int, ...]
    schedule_link_g4_sums: tuple[tuple[int, int], ...]
    positive_a_branch_contradictions: tuple[tuple[int, int, int], ...]
    a0_missing_five_count: int
    a0_required_facet_count: int
    a0_hit_pattern_counts: tuple[tuple[int, int], ...]
    a0_pairwise_pattern_compatible_triples: int
    a0_no_small_transversal_triples: int
    a0_five_transversal_histogram: tuple[tuple[int, int], ...]
    a0_zero_pattern_signature_histogram: tuple[
        tuple[tuple[int, int, int], int, int, int], ...
    ]
    a0_families: tuple[TransversalFamily, ...]

    @property
    def verified(self) -> bool:
        return (
            self.minimal_graph_covers
            == tuple(_mask(cover) for cover in MINIMAL_GRAPH_COVERS)
            and all(
                total < minimum
                for _, total, minimum in self.positive_a_branch_contradictions
            )
            and self.a0_missing_five_count == 3
            and self.a0_required_facet_count == 78
            and dict(self.a0_hit_pattern_counts) == EXPECTED_PATTERN_COUNTS
            and self.a0_pairwise_pattern_compatible_triples == 12_605_586
            and self.a0_no_small_transversal_triples == 32_256
            and dict(self.a0_five_transversal_histogram)
            == {54: 7_560, 57: 22_680, 63: 1_260, 75: 756}
            and {
                (pairwise, triple, count): multiplicity
                for pairwise, triple, count, multiplicity in (
                    self.a0_zero_pattern_signature_histogram
                )
            }
            == {
                ((0, 0, 3), 0, 63): 1_260,
                ((0, 1, 2), 0, 57): 15_120,
                ((1, 1, 1), 0, 54): 7_560,
                ((1, 1, 2), 1, 57): 7_560,
            }
            and sum(
                family.multiplicity for family in self.a0_families
            )
            == self.a0_no_small_transversal_triples
            and all(
                family.five_transversal_count != self.a0_required_facet_count
                for family in self.a0_families
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
            "strict_vertex_link_g4_lower_bound": 1,
            "vertex_link_count": TARGET_LINK_COUNT,
            "schedule_link_g4_sums": [
                {"missing_four_count": a, "summed_vertex_link_g4": total}
                for a, total in self.schedule_link_g4_sums
            ],
            "positive_a_branch_contradictions": [
                {
                    "missing_four_count": a,
                    "summed_vertex_link_g4": total,
                    "required_minimum": minimum,
                }
                for a, total, minimum in self.positive_a_branch_contradictions
            ],
            "a0": {
                "forced_missing_five_count": self.a0_missing_five_count,
                "required_facet_count": self.a0_required_facet_count,
                "hit_pattern_counts": {
                    str(pattern): count
                    for pattern, count in self.a0_hit_pattern_counts
                },
                "pairwise_pattern_compatible_triple_count": (
                    self.a0_pairwise_pattern_compatible_triples
                ),
                "no_small_transversal_triple_count": (
                    self.a0_no_small_transversal_triples
                ),
                "five_transversal_count_histogram": {
                    str(count): multiplicity
                    for count, multiplicity in self.a0_five_transversal_histogram
                },
                "zero_pattern_complement_intersection_signatures": [
                    {
                        "sorted_pairwise_intersection_sizes": list(pairwise),
                        "triple_intersection_size": triple,
                        "five_transversal_count": count,
                        "multiplicity": multiplicity,
                    }
                    for pairwise, triple, count, multiplicity in (
                        self.a0_zero_pattern_signature_histogram
                    )
                ],
                "family_count": len(self.a0_families),
                "families": [family.to_record() for family in self.a0_families],
            },
            "verified": self.verified,
        }


def audit_record93(
    graph6: str, schedules: Sequence[Mapping[str, object]]
) -> Record93Audit:
    """Recompute the complete record-93 exclusion from frozen inputs."""

    if graph6 != TARGET_GRAPH6:
        raise ValueError("certificate is scoped to H9/14 record 93")
    _, decoded_edges = decode_graph6(graph6)
    if tuple(decoded_edges) != TARGET_MISSING_EDGES:
        raise ValueError("record-93 missing-edge graph does not verify")

    parsed_schedules: dict[int, tuple[int, ...]] = {}
    for schedule in schedules:
        a = schedule.get("missing_four_count")
        f_vector = schedule.get("f_vector")
        if not isinstance(a, int) or not isinstance(f_vector, (list, tuple)):
            raise ValueError("invalid record-93 Dehn--Sommerville schedule")
        parsed_schedules[a] = tuple(int(value) for value in f_vector)
    if tuple(sorted(parsed_schedules)) != (0, 1, 2):
        raise ValueError("record-93 schedules must be exactly a=0,1,2")

    link_sums = tuple(
        (a, summed_vertex_link_g4(parsed_schedules[a]))
        for a in sorted(parsed_schedules)
    )
    positive_contradictions = tuple(
        (a, total, TARGET_LINK_COUNT)
        for a, total in link_sums
        if a > 0
    )
    if any(total >= minimum for _, total, minimum in positive_contradictions):
        raise AssertionError("positive-a branch lacks the strict-g4 contradiction")

    minimal_covers = minimal_graph_vertex_covers(graph6)
    expected_covers = tuple(_mask(cover) for cover in MINIMAL_GRAPH_COVERS)
    if minimal_covers != expected_covers:
        raise AssertionError("record-93 minimal graph covers do not verify")

    independent_fives = independent_masks(graph6, 5)
    a0_vector = parsed_schedules[0]
    missing_five_count = len(independent_fives) - a0_vector[4]
    if missing_five_count != 3:
        raise AssertionError("a=0 does not force exactly three missing fives")

    by_pattern: dict[int, list[int]] = {}
    for face in independent_fives:
        by_pattern.setdefault(hit_pattern(face, minimal_covers), []).append(face)
    pattern_counts = tuple(
        (pattern, len(by_pattern[pattern])) for pattern in sorted(by_pattern)
    )
    if dict(pattern_counts) != EXPECTED_PATTERN_COUNTS:
        raise AssertionError("record-93 A5 hit-pattern table does not verify")

    graph_cover_fives = graph_vertex_cover_masks(graph6, 5)
    compatible_count = 0
    no_small_count = 0
    transversal_histogram: Counter[int] = Counter()
    zero_pattern_signatures: Counter[
        tuple[tuple[int, int, int], int, int]
    ] = Counter()
    family_counts: Counter[tuple[tuple[int, int, int], int]] = Counter()
    family_representatives: dict[
        tuple[tuple[int, int, int], int], tuple[int, int, int]
    ] = {}
    for missing_fives in pattern_compatible_missing_five_triples(by_pattern):
        compatible_count += 1
        if has_transversal_of_size_at_most_four(
            missing_fives, minimal_covers
        ):
            continue
        no_small_count += 1
        count = five_transversal_count(graph_cover_fives, missing_fives)
        transversal_histogram[count] += 1
        pattern_multiset = tuple(
            sorted(hit_pattern(face, minimal_covers) for face in missing_fives)
        )
        if pattern_multiset == (0, 0, 0):
            isolated_mask = _mask(tuple(range(2, 11)))
            complements = tuple(
                isolated_mask ^ face for face in missing_fives
            )
            pairwise = tuple(
                sorted(
                    (complements[first] & complements[second]).bit_count()
                    for first, second in ((0, 1), (0, 2), (1, 2))
                )
            )
            triple = (
                complements[0] & complements[1] & complements[2]
            ).bit_count()
            zero_pattern_signatures[(pairwise, triple, count)] += 1
        family_key = (pattern_multiset, count)
        family_counts[family_key] += 1
        family_representatives.setdefault(family_key, tuple(missing_fives))

    families = tuple(
        TransversalFamily(
            pattern_multiset=pattern_multiset,
            multiplicity=family_counts[(pattern_multiset, count)],
            five_transversal_count=count,
            representative_masks=family_representatives[
                (pattern_multiset, count)
            ],
        )
        for pattern_multiset, count in sorted(family_counts)
    )
    audit = Record93Audit(
        graph6=graph6,
        missing_edges=tuple(decoded_edges),
        minimal_graph_covers=minimal_covers,
        schedule_link_g4_sums=link_sums,
        positive_a_branch_contradictions=positive_contradictions,
        a0_missing_five_count=missing_five_count,
        a0_required_facet_count=a0_vector[8],
        a0_hit_pattern_counts=pattern_counts,
        a0_pairwise_pattern_compatible_triples=compatible_count,
        a0_no_small_transversal_triples=no_small_count,
        a0_five_transversal_histogram=tuple(sorted(transversal_histogram.items())),
        a0_zero_pattern_signature_histogram=tuple(
            (pairwise, triple, count, multiplicity)
            for (pairwise, triple, count), multiplicity in sorted(
                zero_pattern_signatures.items()
            )
        ),
        a0_families=families,
    )
    if not audit.verified:
        raise AssertionError("record-93 arithmetic certificate did not verify")
    return audit
