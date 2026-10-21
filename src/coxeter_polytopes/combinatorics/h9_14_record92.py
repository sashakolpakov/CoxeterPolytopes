"""Solver-independent exclusion certificate for H9/14 record 92.

The three Dehn--Sommerville branches are excluded by two small arguments.
For one or two missing four-faces, double-counting the sum of ``g_4`` over
the fourteen vertex links contradicts strict positivity.  With no missing
four-face, Dehn--Sommerville forces three missing five-faces.  Minimal
transversals of the resulting nonface clutter are complements of maximal
faces.  A finite symmetry classification leaves two families, with 57 or 75
five-element transversals instead of the required 78 facets.

Everything here is integer arithmetic and finite enumeration.  The module
does not import or trust a SAT/SMT backend.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, product
from math import comb
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_POLYTOPE_DIMENSION,
    TARGET_VERTEX_COUNT,
    decode_graph6,
)


TARGET_RECORD_INDEX = 92
TARGET_GRAPH6 = "M?????????W?o?w??"
TARGET_MISSING_EDGES = (
    (0, 11),
    (1, 11),
    (0, 12),
    (1, 12),
    (0, 13),
    (1, 13),
    (2, 13),
)
VERTEX_LINK_POLYTOPE_DIMENSION = TARGET_POLYTOPE_DIMENSION - 1
TARGET_LINK_COUNT = TARGET_VERTEX_COUNT
ACTIVE_VERTICES = (0, 1, 2, 11, 12, 13)
ISOLATED_VERTICES = tuple(range(3, 11))
MINIMAL_GRAPH_COVERS = (
    (0, 1, 2),
    (0, 1, 13),
    (11, 12, 13),
)
EXPECTED_PATTERN_COUNTS = {
    0b000: 56,
    0b001: 70,
    0b011: 336,
    0b100: 196,
    0b101: 140,
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


def _edge_masks(graph6: str) -> tuple[int, ...]:
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise ValueError("record-92 certificate requires fourteen vertices")
    return tuple((1 << first) | (1 << second) for first, second in edges)


def independent_masks(graph6: str, size: int) -> tuple[int, ...]:
    """Return all graph-independent masks of one order in canonical order."""

    edge_masks = _edge_masks(graph6)
    values = []
    for vertices in combinations(range(TARGET_VERTEX_COUNT), size):
        mask = _mask(vertices)
        if all(mask & edge != edge for edge in edge_masks):
            values.append(mask)
    return tuple(values)


def graph_vertex_cover_masks(graph6: str, size: int) -> tuple[int, ...]:
    """Return all vertex covers of the missing-edge graph of one size."""

    edge_masks = _edge_masks(graph6)
    values = []
    for vertices in combinations(range(TARGET_VERTEX_COUNT), size):
        mask = _mask(vertices)
        if all(mask & edge for edge in edge_masks):
            values.append(mask)
    return tuple(values)


def minimal_graph_vertex_covers(graph6: str) -> tuple[int, ...]:
    """Enumerate every inclusion-minimal cover without a graph library."""

    edge_masks = _edge_masks(graph6)
    values = []
    for mask in range(1 << TARGET_VERTEX_COUNT):
        if not all(mask & edge for edge in edge_masks):
            continue
        if all(
            not all((mask ^ bit) & edge for edge in edge_masks)
            for bit in (1 << vertex for vertex in _vertices(mask))
        ):
            values.append(mask)
    return tuple(values)


def _h_entry_from_f(
    polytope_dimension: int,
    f_by_order: Mapping[int, int],
    index: int,
) -> int:
    return sum(
        (-1) ** (index - order)
        * comb(
            polytope_dimension - order,
            polytope_dimension - index,
        )
        * f_by_order[order]
        for order in range(index + 1)
    )


def summed_vertex_link_g4(f_vector: Sequence[int]) -> int:
    """Double-count the sum of ``g4`` over all fourteen vertex links.

    An order-``j`` face of a vertex link corresponds to an order-``j+1``
    face of the global boundary containing that vertex.  Summing over all
    vertices therefore multiplies the global count by ``j+1``.
    """

    if len(f_vector) != TARGET_POLYTOPE_DIMENSION:
        raise ValueError("the H9 polar boundary needs f0 through f8")
    summed_link_f = {0: TARGET_VERTEX_COUNT}
    for order in range(1, 5):
        summed_link_f[order] = (order + 1) * int(f_vector[order])
    h3 = _h_entry_from_f(
        VERTEX_LINK_POLYTOPE_DIMENSION, summed_link_f, 3
    )
    h4 = _h_entry_from_f(
        VERTEX_LINK_POLYTOPE_DIMENSION, summed_link_f, 4
    )
    return h4 - h3


def hit_pattern(mask: int, minimal_covers: Sequence[int]) -> int:
    """Encode which of the three minimal graph covers meet ``mask``."""

    return sum(
        bool(mask & cover) << index
        for index, cover in enumerate(minimal_covers)
    )


def has_transversal_of_size_at_most_four(
    missing_fives: Sequence[int], minimal_covers: Sequence[int]
) -> bool:
    """Test the small-transversal obstruction through the minimal covers.

    Every transversal of the missing-edge graph contains a minimal graph
    cover.  A cover ``C`` extends to a transversal using at most one extra
    vertex exactly when the missing five-faces disjoint from ``C`` have a
    common vertex (with the empty missed family meaning ``C`` already works).
    """

    for cover in minimal_covers:
        missed = tuple(face for face in missing_fives if not face & cover)
        if not missed:
            return True
        common = missed[0]
        for face in missed[1:]:
            common &= face
        if common:
            return True
    return False


def five_transversal_count(
    graph_cover_fives: Sequence[int], missing_fives: Sequence[int]
) -> int:
    """Count five-sets hitting the graph edges and all missing five-faces."""

    return sum(
        all(transversal & face for face in missing_fives)
        for transversal in graph_cover_fives
    )


def pattern_compatible_missing_five_triples(
    by_pattern: Mapping[int, Sequence[int]],
):
    """Yield all unordered triples with pairwise-disjoint hit patterns."""

    patterns = tuple(sorted(by_pattern))
    for first_index, first_pattern in enumerate(patterns):
        for second_index in range(first_index, len(patterns)):
            second_pattern = patterns[second_index]
            if first_pattern & second_pattern:
                continue
            for third_index in range(second_index, len(patterns)):
                third_pattern = patterns[third_index]
                if (
                    first_pattern & third_pattern
                    or second_pattern & third_pattern
                ):
                    continue
                first_values = by_pattern[first_pattern]
                second_values = by_pattern[second_pattern]
                third_values = by_pattern[third_pattern]
                if first_pattern == second_pattern == third_pattern:
                    iterator = combinations(first_values, 3)
                elif first_pattern == second_pattern:
                    iterator = (
                        (first, second, third)
                        for first, second in combinations(first_values, 2)
                        for third in third_values
                    )
                elif second_pattern == third_pattern:
                    iterator = (
                        (first, second, third)
                        for first in first_values
                        for second, third in combinations(second_values, 2)
                    )
                else:
                    iterator = product(
                        first_values, second_values, third_values
                    )
                yield from iterator


@dataclass(frozen=True)
class TransversalFamily:
    """One theorem-facing family of A5 triples surviving small transversals."""

    pattern_multiset: tuple[int, int, int]
    multiplicity: int
    five_transversal_count: int
    representative_masks: tuple[int, int, int]

    def to_record(self) -> dict[str, object]:
        return {
            "hit_pattern_multiset": list(self.pattern_multiset),
            "multiplicity": self.multiplicity,
            "five_transversal_count": self.five_transversal_count,
            "representative_masks": list(self.representative_masks),
            "representative_vertices": [
                list(_vertices(mask)) for mask in self.representative_masks
            ],
        }


@dataclass(frozen=True)
class Record92Audit:
    """Complete arithmetic certificate for the three record-92 schedules."""

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
            and self.a0_pairwise_pattern_compatible_triples == 6_773_256
            and self.a0_no_small_transversal_triples == 1_176
            and dict(self.a0_five_transversal_histogram) == {57: 840, 75: 336}
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
                "family_count": len(self.a0_families),
                "families": [family.to_record() for family in self.a0_families],
            },
            "verified": self.verified,
        }


def audit_record92(
    graph6: str, schedules: Sequence[Mapping[str, object]]
) -> Record92Audit:
    """Recompute the complete record-92 exclusion from frozen inputs."""

    if graph6 != TARGET_GRAPH6:
        raise ValueError("certificate is scoped to H9/14 record 92")
    _, decoded_edges = decode_graph6(graph6)
    if tuple(decoded_edges) != TARGET_MISSING_EDGES:
        raise ValueError("record-92 missing-edge graph does not verify")

    parsed_schedules: dict[int, tuple[int, ...]] = {}
    for schedule in schedules:
        a = schedule.get("missing_four_count")
        f_vector = schedule.get("f_vector")
        if not isinstance(a, int) or not isinstance(f_vector, (list, tuple)):
            raise ValueError("invalid record-92 Dehn--Sommerville schedule")
        parsed_schedules[a] = tuple(int(value) for value in f_vector)
    if tuple(sorted(parsed_schedules)) != (0, 1, 2):
        raise ValueError("record-92 schedules must be exactly a=0,1,2")

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
        raise AssertionError("record-92 minimal graph covers do not verify")

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
        raise AssertionError("record-92 A5 hit-pattern table does not verify")

    graph_cover_fives = graph_vertex_cover_masks(graph6, 5)
    compatible_count = 0
    no_small_count = 0
    transversal_histogram: Counter[int] = Counter()
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
    audit = Record92Audit(
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
        a0_families=families,
    )
    if not audit.verified:
        raise AssertionError("record-92 arithmetic certificate did not verify")
    return audit
