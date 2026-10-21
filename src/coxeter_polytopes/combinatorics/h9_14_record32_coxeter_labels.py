"""Unrestricted compact-Coxeter obstruction for H9/14 record 32.

The ordinary-polytopal record-32 nerve has two overlapping order-four
minimal nonfaces and one disjoint order-five minimal nonface.  In a compact
Coxeter realization these supports must carry Lanner diagrams, while each of
the 110 rank-nine vertex diagrams must be elliptic.

This module exhausts the classical labelled Lanner placements using exact
finite-Coxeter classification.  The two order-four supports have 546
compatible joins and 48 survive their actual-vertex restrictions.  All 420
placements on the order-five support survive its proper restrictions.  For
the cross gate, the 48 distinct full overlap states restrict to 24 first-A4
labelings, each with multiplicity two; retaining those multiplicities gives
20,160 full sparse states.  For each state, every one of the twenty pairs
between the first order-four support and the order-five support is tested at
the smallest non-right label, ``m=3``.  None passes all actual vertices.

This finite gate is unrestricted.  For normalized Gram ``G = I-A``, lowering
external labels to two and a selected cross label to three decreases the
nonnegative matrix ``A`` entrywise.  Perron--Frobenius monotonicity therefore
shows that labels four/five and arbitrary extra finite labels cannot repair a
failed rank-nine vertex.  If no cross pair is non-right, the principal Gram
on the two disjoint Lanner supports is their orthogonal direct sum and has
negative index two.  Interlacing then excludes ambient inertia ``(9,1,4)``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_disjoint_lanner_gate import (
    audit_disjoint_lanner_cross_gate,
)
from coxeter_polytopes.combinatorics.h9_14_record32 import (
    FULL_MASK,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    record32_facet_complements,
    vertices,
)


FACET_COUNT = 14
ALL_PAIRS = tuple(combinations(range(FACET_COUNT), 2))
DOTTED_PAIRS = frozenset(pair(*edge) for edge in TARGET_MISSING_EDGES)
MEETING_PAIRS = tuple(edge for edge in ALL_PAIRS if edge not in DOTTED_PAIRS)

FIRST_A4 = TARGET_MISSING_FOURS[0]
SECOND_A4 = TARGET_MISSING_FOURS[1]
A4_UNION = tuple(sorted(set(FIRST_A4) | set(SECOND_A4)))
A5_SUPPORT = TARGET_MISSING_FIVES[0]
SPARSE_UNION = tuple(sorted(set(A4_UNION) | set(A5_SUPPORT)))
CROSS_PAIRS = tuple(pair(left, right) for left in FIRST_A4 for right in A5_SUPPORT)

LabelState = tuple[tuple[Pair, int], ...]


def record32_vertex_facets() -> tuple[tuple[int, ...], ...]:
    """Return the complete deterministic list of 110 rank-nine vertices."""

    return tuple(
        sorted(
            vertices(FULL_MASK ^ complement)
            for complement in record32_facet_complements()
        )
    )


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


def _merge(left: LabelState, right: LabelState) -> LabelState | None:
    labels = dict(left)
    for edge, label in right:
        previous = labels.get(edge)
        if previous is not None and previous != label:
            return None
        labels[edge] = label
    return tuple(sorted(labels.items()))


def _full_labels(
    support: Sequence[int], state: Mapping[Pair, int] | LabelState
) -> dict[Pair, int]:
    labels = {
        pair(left, right): 2
        for left, right in combinations(tuple(sorted(support)), 2)
    }
    labels.update(dict(state))
    return labels


def _restrictions(
    support: Iterable[int], vertex_facets: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    support_set = set(support)
    return tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in support_set)
                for vertex in vertex_facets
            },
            key=lambda item: (len(item), item),
        )
    )


def _restriction_histogram(
    restrictions: Sequence[Sequence[int]],
) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(map(len, restrictions)).items()))


@dataclass(frozen=True, slots=True)
class Record32CoxeterLabelAudit:
    vertex_count: int
    meeting_pair_count: int
    order_four_placement_count: int
    order_five_placement_count: int
    a4_join_count: int
    a4_restriction_size_histogram: tuple[tuple[int, int], ...]
    a4_passing_state_count: int
    a4_first_support_projection_count: int
    a4_first_support_projection_multiplicity_histogram: tuple[tuple[int, int], ...]
    a5_restriction_size_histogram: tuple[tuple[int, int], ...]
    a5_passing_state_count: int
    sparse_locally_elliptic_state_count: int
    cross_supports: tuple[tuple[int, ...], tuple[int, ...]]
    cross_pair_count: int
    minimal_cross_label: int
    cross_label_three_test_count: int
    cross_edge_passing_state_counts: tuple[tuple[Pair, int], ...]
    first_bad_restriction_size_histogram: tuple[tuple[int, int], ...]
    cross_label_three_passing_count: int

    @property
    def unrestricted_compact_coxeter_no_go(self) -> bool:
        return (
            self.cross_label_three_passing_count == 0
            and all(count == 0 for _edge, count in self.cross_edge_passing_state_counts)
        )

    @property
    def verified(self) -> bool:
        return (
            self.vertex_count == 110
            and self.meeting_pair_count == 87
            and self.order_four_placement_count == 108
            and self.order_five_placement_count == 420
            and self.a4_join_count == 546
            and self.a4_restriction_size_histogram == ((3, 7), (4, 3))
            and self.a4_passing_state_count == 48
            and self.a4_first_support_projection_count == 24
            and self.a4_first_support_projection_multiplicity_histogram == ((2, 24),)
            and self.a5_restriction_size_histogram == ((4, 5),)
            and self.a5_passing_state_count == 420
            and self.sparse_locally_elliptic_state_count == 20_160
            and self.cross_supports == (FIRST_A4, A5_SUPPORT)
            and self.cross_pair_count == 20
            and self.minimal_cross_label == 3
            and self.cross_label_three_test_count == 403_200
            and len(self.cross_edge_passing_state_counts) == 20
            and self.cross_label_three_passing_count == 0
            and sum(count for _size, count in self.first_bad_restriction_size_histogram)
            == self.cross_label_three_test_count
            and self.unrestricted_compact_coxeter_no_go
        )

    def to_record(self) -> dict[str, object]:
        return {
            "facet_count": FACET_COUNT,
            "rank_nine_vertex_count": self.vertex_count,
            "meeting_pair_count": self.meeting_pair_count,
            "minimal_nonfaces": {
                "dotted_pairs": [list(edge) for edge in TARGET_MISSING_EDGES],
                "order_four_lanner_supports": [
                    list(face) for face in TARGET_MISSING_FOURS
                ],
                "order_five_lanner_support": list(A5_SUPPORT),
            },
            "classical_lanner_placements": {
                "order_four": self.order_four_placement_count,
                "order_five": self.order_five_placement_count,
            },
            "overlapping_order_four_gate": {
                "supports": [list(FIRST_A4), list(SECOND_A4)],
                "compatible_join_count": self.a4_join_count,
                "actual_vertex_restriction_size_histogram": {
                    str(size): count
                    for size, count in self.a4_restriction_size_histogram
                },
                "passing_state_count": self.a4_passing_state_count,
                "distinct_first_support_projection_count": (
                    self.a4_first_support_projection_count
                ),
                "first_support_projection_multiplicity_histogram": {
                    str(multiplicity): count
                    for multiplicity, count in (
                        self.a4_first_support_projection_multiplicity_histogram
                    )
                },
            },
            "order_five_gate": {
                "support": list(A5_SUPPORT),
                "actual_vertex_restriction_size_histogram": {
                    str(size): count
                    for size, count in self.a5_restriction_size_histogram
                },
                "passing_state_count": self.a5_passing_state_count,
            },
            "sparse_locally_elliptic_state_count": (
                self.sparse_locally_elliptic_state_count
            ),
            "cross_support_gate": {
                "disjoint_lanner_supports": [
                    list(support) for support in self.cross_supports
                ],
                "cross_pair_count": self.cross_pair_count,
                "tested_label": self.minimal_cross_label,
                "test_count": self.cross_label_three_test_count,
                "passing_count": self.cross_label_three_passing_count,
                "passing_states_by_cross_pair": [
                    {"pair": list(edge), "passing_states": count}
                    for edge, count in self.cross_edge_passing_state_counts
                ],
                "first_bad_restriction_size_histogram": {
                    str(size): count
                    for size, count in self.first_bad_restriction_size_histogram
                },
            },
            "unrestricted_compact_coxeter_no_go": {
                "proved": self.unrestricted_compact_coxeter_no_go,
                "cross_edge_case": (
                    "Any non-right edge between the disjoint Lanner supports "
                    "already fails an actual vertex at m=3; labels 4/5 and "
                    "arbitrary additional non-right labels only increase the "
                    "nonnegative cosine adjacency, so Perron--Frobenius "
                    "monotonicity forbids repair."
                ),
                "orthogonal_case": (
                    "If all cross pairs are right, the principal Gram on the "
                    "order-four and order-five supports is their orthogonal "
                    "direct sum and has negative index two; interlacing "
                    "excludes ambient inertia (9,1,4)."
                ),
                "dotted_variables": (
                    "The four dotted variables lie outside this decisive "
                    "principal block and cannot change either alternative."
                ),
            },
            "verified": self.verified,
        }


@lru_cache(maxsize=1)
def audit_record32_coxeter_labels() -> Record32CoxeterLabelAudit:
    """Exhaust the forced Lanner joins and certify the unrestricted no-go."""

    vertex_facets = record32_vertex_facets()
    if len(vertex_facets) != 110 or any(len(vertex) != 9 for vertex in vertex_facets):
        raise AssertionError("record-32 vertex list changed")

    first_options = _placement_options(FIRST_A4)
    second_options = _placement_options(SECOND_A4)
    a4_states = tuple(
        sorted(
            {
                merged
                for first in first_options
                for second in second_options
                if (merged := _merge(first, second)) is not None
            }
        )
    )
    a4_restrictions = _restrictions(A4_UNION, vertex_facets)
    a4_passing = tuple(
        state
        for state in a4_states
        if all(
            is_finite_coxeter_labeling(
                restriction, _full_labels(A4_UNION, state)
            )
            for restriction in a4_restrictions
        )
    )

    a5_states = _placement_options(A5_SUPPORT)
    a5_restrictions = _restrictions(A5_SUPPORT, vertex_facets)
    a5_passing = tuple(
        state
        for state in a5_states
        if all(
            is_finite_coxeter_labeling(
                restriction, _full_labels(A5_SUPPORT, state)
            )
            for restriction in a5_restrictions
        )
    )

    a4_projected_states = tuple(
        tuple(
            sorted(
                (
                    pair(left, right),
                    dict(state)[pair(left, right)],
                )
                for left, right in combinations(FIRST_A4, 2)
            )
        )
        for state in a4_passing
    )
    a4_projection_multiplicities = Counter(a4_projected_states)

    cross_gate = audit_disjoint_lanner_cross_gate(
        vertex_facets,
        FIRST_A4,
        A5_SUPPORT,
        tuple(dict(state) for state in a4_projected_states),
        tuple(dict(state) for state in a5_passing),
    )

    audit = Record32CoxeterLabelAudit(
        vertex_count=len(vertex_facets),
        meeting_pair_count=len(MEETING_PAIRS),
        order_four_placement_count=len(first_options),
        order_five_placement_count=len(a5_states),
        a4_join_count=len(a4_states),
        a4_restriction_size_histogram=_restriction_histogram(a4_restrictions),
        a4_passing_state_count=len(a4_passing),
        a4_first_support_projection_count=len(a4_projection_multiplicities),
        a4_first_support_projection_multiplicity_histogram=tuple(
            sorted(Counter(a4_projection_multiplicities.values()).items())
        ),
        a5_restriction_size_histogram=_restriction_histogram(a5_restrictions),
        a5_passing_state_count=len(a5_passing),
        sparse_locally_elliptic_state_count=len(a4_passing) * len(a5_passing),
        cross_supports=(FIRST_A4, A5_SUPPORT),
        cross_pair_count=len(CROSS_PAIRS),
        minimal_cross_label=3,
        cross_label_three_test_count=cross_gate.label_three_test_count,
        cross_edge_passing_state_counts=cross_gate.passing_count_by_cross_pair,
        first_bad_restriction_size_histogram=(
            cross_gate.first_bad_restriction_size_histogram
        ),
        cross_label_three_passing_count=len(cross_gate.passing_trials),
    )
    if not audit.verified:
        raise AssertionError("record-32 Coxeter-label audit failed")
    return audit


__all__ = [
    "A4_UNION",
    "A5_SUPPORT",
    "CROSS_PAIRS",
    "FIRST_A4",
    "Record32CoxeterLabelAudit",
    "SECOND_A4",
    "audit_record32_coxeter_labels",
    "record32_vertex_facets",
]
