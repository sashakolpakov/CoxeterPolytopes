"""Reusable exact cross-edge gate for two disjoint Lanner supports.

The input consists of the complete actual rank-nine vertex list, two disjoint
Lanner supports, and finite lists of labelled placements on those supports.
Every support pair is joined in the all-right shadow and then each possible
cross edge is tested at the smallest non-right label ``m=3``.  The test uses
only exact finite-Coxeter classification.

When the passing count is zero, Perron--Frobenius monotonicity promotes the
finite census to an unrestricted alternative: every cross edge is forbidden,
or the orthogonal direct sum of the two Lanner blocks has negative index two.
Callers remain responsible for proving that their supplied placement lists
cover all labelings, either literally or by a stated nerve symmetry.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    pair,
)


@dataclass(frozen=True, slots=True)
class DisjointLannerCrossGateAudit:
    left_support: tuple[int, ...]
    right_support: tuple[int, ...]
    left_state_count: int
    right_state_count: int
    sparse_state_count: int
    cross_pairs: tuple[Pair, ...]
    unique_actual_restriction_size_histogram: tuple[tuple[int, int], ...]
    label_three_test_count: int
    passing_trials: tuple[tuple[int, int, Pair], ...]
    passing_count_by_cross_pair: tuple[tuple[Pair, int], ...]
    first_bad_restriction_size_histogram: tuple[tuple[int, int], ...]

    @property
    def all_cross_edges_forbidden(self) -> bool:
        return not self.passing_trials

    def to_record(self) -> dict[str, object]:
        return {
            "left_support": list(self.left_support),
            "right_support": list(self.right_support),
            "left_state_count": self.left_state_count,
            "right_state_count": self.right_state_count,
            "sparse_state_count": self.sparse_state_count,
            "cross_pairs": [list(edge) for edge in self.cross_pairs],
            "cross_pair_count": len(self.cross_pairs),
            "unique_actual_restriction_size_histogram": {
                str(size): count
                for size, count in self.unique_actual_restriction_size_histogram
            },
            "tested_label": 3,
            "label_three_test_count": self.label_three_test_count,
            "passing_trials": [
                {
                    "left_state_index": left,
                    "right_state_index": right,
                    "cross_pair": list(edge),
                }
                for left, right, edge in self.passing_trials
            ],
            "passing_count": len(self.passing_trials),
            "passing_count_by_cross_pair": [
                {"pair": list(edge), "passing_count": count}
                for edge, count in self.passing_count_by_cross_pair
            ],
            "first_bad_restriction_size_histogram": {
                str(size): count
                for size, count in self.first_bad_restriction_size_histogram
            },
            "all_cross_edges_forbidden": self.all_cross_edges_forbidden,
        }


def audit_disjoint_lanner_cross_gate(
    vertex_facets: Sequence[Sequence[int]],
    left_support: Sequence[int],
    right_support: Sequence[int],
    left_states: Sequence[Mapping[Pair, int]],
    right_states: Sequence[Mapping[Pair, int]],
) -> DisjointLannerCrossGateAudit:
    """Test every supplied sparse state and cross pair at label three."""

    vertices = tuple(tuple(sorted(set(vertex))) for vertex in vertex_facets)
    if not vertices or any(len(vertex) != 9 for vertex in vertices):
        raise ValueError("the disjoint-Lanner gate needs all rank-nine vertices")
    left = tuple(sorted(set(left_support)))
    right = tuple(sorted(set(right_support)))
    if len(left) not in (4, 5) or len(right) not in (4, 5) or set(left) & set(right):
        raise ValueError("the gate needs two disjoint order-four/five supports")
    if not left_states or not right_states:
        raise ValueError("the gate needs nonempty placement lists")
    union = tuple(sorted((*left, *right)))
    cross_pairs = tuple(pair(a, b) for a in left for b in right)
    template = {pair(a, b): 2 for a, b in combinations(union, 2)}

    left_pairs = {pair(a, b) for a, b in combinations(left, 2)}
    right_pairs = {pair(a, b) for a, b in combinations(right, 2)}
    normalized_left = tuple(dict(state) for state in left_states)
    normalized_right = tuple(dict(state) for state in right_states)
    if any(set(state) != left_pairs for state in normalized_left):
        raise ValueError("a left placement does not label every internal pair")
    if any(set(state) != right_pairs for state in normalized_right):
        raise ValueError("a right placement does not label every internal pair")

    restrictions = tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in union)
                for vertex in vertices
            },
            key=lambda item: (len(item), item),
        )
    )
    restrictions_by_edge = {
        edge: tuple(
            restriction
            for restriction in restrictions
            if edge[0] in restriction and edge[1] in restriction
        )
        for edge in cross_pairs
    }
    if any(not values for values in restrictions_by_edge.values()):
        raise ValueError("a cross pair occurs in no actual vertex")

    passing: list[tuple[int, int, Pair]] = []
    pass_by_edge: Counter[Pair] = Counter({edge: 0 for edge in cross_pairs})
    bad_sizes: Counter[int] = Counter()
    tests = 0
    for left_index, left_state in enumerate(normalized_left):
        for right_index, right_state in enumerate(normalized_right):
            sparse = dict(template)
            sparse.update(left_state)
            sparse.update(right_state)
            for edge in cross_pairs:
                tests += 1
                labels = dict(sparse)
                labels[edge] = 3
                bad = next(
                    (
                        restriction
                        for restriction in restrictions_by_edge[edge]
                        if not is_finite_coxeter_labeling(restriction, labels)
                    ),
                    None,
                )
                if bad is None:
                    passing.append((left_index, right_index, edge))
                    pass_by_edge[edge] += 1
                else:
                    bad_sizes[len(bad)] += 1

    return DisjointLannerCrossGateAudit(
        left_support=left,
        right_support=right,
        left_state_count=len(normalized_left),
        right_state_count=len(normalized_right),
        sparse_state_count=len(normalized_left) * len(normalized_right),
        cross_pairs=cross_pairs,
        unique_actual_restriction_size_histogram=tuple(
            sorted(Counter(map(len, restrictions)).items())
        ),
        label_three_test_count=tests,
        passing_trials=tuple(passing),
        passing_count_by_cross_pair=tuple(
            (edge, pass_by_edge[edge]) for edge in cross_pairs
        ),
        first_bad_restriction_size_histogram=tuple(sorted(bad_sizes.items())),
    )


__all__ = [
    "DisjointLannerCrossGateAudit",
    "audit_disjoint_lanner_cross_gate",
]
