#!/usr/bin/env python3
"""Unrestricted cross-edge scouts for the two ridge-valid ``a=1`` targets.

Each target contains two disjoint Lanner supports.  If every meeting pair
between those supports is right-angled, their principal Gram is an orthogonal
sum with negative index two.  Otherwise, lower one non-right cross label to
three and every other unforced finite label to two.  Perron--Frobenius
monotonicity makes the resulting single-edge actual-vertex test necessary.

This file exhausts those tests for every sparse locally elliptic Lanner state.
A zero single-edge survivor count is therefore an unrestricted compact
obstruction.  Positive trials are reported honestly for stronger-label,
multi-edge, and exact-inertia follow-up.  No finite-volume inference is made,
and results remain scouts until independently replayed.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import Pair, pair
from scripts.explore_h9_14_record32_a0_target_gate import (
    EDGE_MASKS,
    LabelState,
    next_support,
    placement_options,
    rank_nine_vertices,
    restriction_finite,
    vertices,
)
from scripts.explore_h9_14_record32_a2_additions import mask


TARGETS = (
    {
        "name": "a1_orbit0_unique_ridge_target",
        "minimal_faces": (
            (3, 4, 5, 6),
            (2, 7, 8, 9, 10),
            (7, 8, 9, 10, 11),
        ),
        "disjoint_cross_support_indices": (0, 1),
        "expected_sparse_state_count": 10_368,
    },
    {
        "name": "a1_orbit1_unique_ridge_target",
        "minimal_faces": (
            (2, 3, 4, 5),
            (6, 7, 8, 9, 10),
            (3, 4, 5, 6, 11),
            (7, 8, 9, 10, 13),
        ),
        "disjoint_cross_support_indices": (0, 1),
        "expected_sparse_state_count": 864,
    },
)


def _cluster_restrictions(
    actual_vertices: Sequence[Sequence[int]], cluster: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    cluster_set = set(cluster)
    return tuple(
        sorted(
            {
                tuple(vertex for vertex in actual if vertex in cluster_set)
                for actual in actual_vertices
            },
            key=lambda row: (len(row), row),
        )
    )


def enumerate_sparse_states(
    minimal_faces: Sequence[Sequence[int]],
) -> tuple[tuple[LabelState, ...], tuple[tuple[int, ...], ...], list[dict[str, int]]]:
    masks = tuple(mask(face) for face in minimal_faces)
    supports = tuple(tuple(sorted(face)) for face in minimal_faces)
    actual_vertices = rank_nine_vertices(masks)
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    while states and remaining:
        support_index = next_support(supports, remaining, union)
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = _cluster_restrictions(actual_vertices, cluster)
        options = placement_options(support)
        passing: set[LabelState] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in options:
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
            }
        )
        remaining.remove(support_index)
    return states, actual_vertices, stages


def audit_target(target: Mapping[str, object]) -> dict[str, object]:
    minimal_faces = tuple(tuple(face) for face in target["minimal_faces"])
    states, actual_vertices, stages = enumerate_sparse_states(minimal_faces)
    expected = int(target["expected_sparse_state_count"])
    if len(states) != expected:
        raise AssertionError("record-32 a=1 sparse state count changed")
    left_index, right_index = target["disjoint_cross_support_indices"]
    left = minimal_faces[int(left_index)]
    right = minimal_faces[int(right_index)]
    if set(left) & set(right):
        raise AssertionError("cross-gate supports are not disjoint")
    dotted = set(EDGE_MASKS)
    cross_pairs = tuple(
        pair(left_vertex, right_vertex)
        for left_vertex in left
        for right_vertex in right
        if mask((left_vertex, right_vertex)) not in dotted
    )
    if len(cross_pairs) != len(left) * len(right):
        raise AssertionError("a dotted pair entered the disjoint-support cross gate")

    cluster = tuple(sorted(set().union(*map(set, minimal_faces))))
    restrictions = _cluster_restrictions(actual_vertices, cluster)
    restrictions_by_edge = {
        edge: tuple(row for row in restrictions if set(edge) <= set(row))
        for edge in cross_pairs
    }
    if any(not rows for rows in restrictions_by_edge.values()):
        raise AssertionError("a cross edge occurs in no actual vertex")
    template = {
        pair(left_vertex, right_vertex): 2
        for left_vertex, right_vertex in combinations(cluster, 2)
    }
    passing: list[tuple[int, Pair]] = []
    passing_by_state: Counter[int] = Counter()
    passing_by_edge: Counter[Pair] = Counter()
    first_bad_sizes: Counter[int] = Counter()
    for state_index, state in enumerate(states):
        base = dict(template)
        base.update(dict(state))
        for edge in cross_pairs:
            labels = dict(base)
            labels[edge] = 3
            bad = next(
                (
                    restriction
                    for restriction in restrictions_by_edge[edge]
                    if not restriction_finite(restriction, labels)
                ),
                None,
            )
            if bad is None:
                passing.append((state_index, edge))
                passing_by_state[state_index] += 1
                passing_by_edge[edge] += 1
            else:
                first_bad_sizes[len(bad)] += 1

    stronger = []
    for label in (4, 5):
        count = 0
        for state_index, edge in passing:
            labels = dict(template)
            labels.update(dict(states[state_index]))
            labels[edge] = label
            count += all(
                restriction_finite(restriction, labels)
                for restriction in restrictions_by_edge[edge]
            )
        stronger.append((label, count))

    unrestricted_no_go = not passing
    return {
        "name": target["name"],
        "minimal_faces": [list(face) for face in minimal_faces],
        "rank_nine_vertex_count": len(actual_vertices),
        "sparse_join_stages": stages,
        "sparse_locally_elliptic_state_count": len(states),
        "disjoint_lanner_supports": [list(left), list(right)],
        "cross_pair_count": len(cross_pairs),
        "minimum_nonright_cross_edge_count_from_inertia": 1,
        "single_m3_test_count": len(states) * len(cross_pairs),
        "single_m3_passing_count": len(passing),
        "passing_state_count": len(passing_by_state),
        "passing_edges_per_state_histogram": dict(
            sorted(Counter(passing_by_state.values()).items())
        ),
        "passing_count_by_cross_pair": [
            {"pair": list(edge), "count": passing_by_edge[edge]}
            for edge in cross_pairs
        ],
        "first_bad_restriction_size_histogram": dict(sorted(first_bad_sizes.items())),
        "stronger_single_edge_survivors": dict(stronger),
        "unrestricted_compact_label_no_go": unrestricted_no_go,
        "residual_scope": (
            "closed by PF monotonicity and the orthogonal negative-index-two "
            "alternative"
            if unrestricted_no_go
            else "positive single-edge trials require multi-edge and exact-inertia audit"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=int, choices=range(len(TARGETS)), action="append")
    args = parser.parse_args()
    indices = tuple(args.target) if args.target else tuple(range(len(TARGETS)))
    audits = tuple(audit_target(TARGETS[index]) for index in indices)
    print(
        json.dumps(
            {
                "kind": "h9_14_record32_a1_unrestricted_cross_edge_scout",
                "theorem_certificate": False,
                "target_audits": audits,
                "scope": (
                    "Compact Lanner, actual-vertex, PF monotonicity, and principal "
                    "negative-index reasoning only. Independent replay required; "
                    "no finite-volume claim."
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
