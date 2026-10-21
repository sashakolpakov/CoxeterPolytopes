#!/usr/bin/env python3
"""Exact exploratory label census for record 35's residual a=1 nerve."""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go import (
    _a5_witnesses,
    _canonical_lanner_states,
)
from coxeter_polytopes.combinatorics.h9_14_record35_all_completions import (
    _one_bridge_class_key,
)


EDGES = ((0, 11), (1, 12), (0, 13), (1, 13), (2, 13))
Q = (2, 3, 4, 5)
R = (3, 4, 5, 11, 12)
B = (6, 7, 8, 9, 10)
LEFT = tuple(sorted(set(Q) | set(R)))
UNION = tuple(sorted(set(LEFT) | set(B)))
INTERNAL_FREE = tuple(pair(a, b) for a in set(Q) - set(R) for b in set(R) - set(Q))
CROSS = tuple(pair(a, b) for a in LEFT for b in B)


@lru_cache(maxsize=None)
def placements(support):
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(dict(zip(edges, code, strict=True)) for code in lanner_placement_codes(len(support)))


def merge(left, right):
    answer = dict(left)
    for edge, label in right.items():
        if edge in answer and answer[edge] != label:
            return None
        answer[edge] = label
    return answer


def actual_vertices():
    nonfaces = tuple(frozenset(face) for face in (*EDGES, Q, R, B))
    return tuple(
        chosen
        for chosen in combinations(range(14), 9)
        if not any(face <= set(chosen) for face in nonfaces)
    )


def restrictions(vertices, union):
    support = set(union)
    return tuple(sorted({tuple(v for v in vertex if v in support) for vertex in vertices}, key=lambda row: (len(row), row)))


@lru_cache(maxsize=262_144)
def finite_code(order, code):
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def finite_row(row, labels):
    code = tuple(labels[pair(a, b)] for a, b in combinations(row, 2))
    return finite_code(len(row), code)


def finite_on(labels, rows):
    return all(finite_row(row, labels) for row in rows)


def full_template():
    return {pair(a, b): 2 for a, b in combinations(UNION, 2)}


def main():
    vertices = actual_vertices()
    all_q = placements(Q)
    all_r = placements(R)
    all_b = placements(B)
    left_rows = restrictions(vertices, LEFT)
    left_states = []
    compatible = 0
    for q in all_q:
        for r in all_r:
            joined = merge(q, r)
            if joined is None:
                continue
            compatible += 1
            labels = {pair(a, b): 2 for a, b in combinations(LEFT, 2)}
            labels.update(joined)
            if finite_on(labels, left_rows):
                left_states.append(joined)
    left_states = tuple(left_states)
    rows = restrictions(vertices, UNION)
    rows_by_edge = {
        edge: tuple(row for row in rows if set(edge) <= set(row))
        for edge in (*CROSS, *INTERNAL_FREE)
    }
    template = full_template()
    sparse_bad = 0
    passing = []
    by_state = defaultdict(list)
    bad_sizes = Counter()
    tests = 0
    for li, left in enumerate(left_states):
        for ri, right in enumerate(all_b):
            base = dict(template)
            base.update(left)
            base.update(right)
            if not finite_on(base, rows):
                sparse_bad += 1
            for edge in CROSS:
                tests += 1
                labels = dict(base)
                labels[edge] = 3
                bad = next((row for row in rows_by_edge[edge] if not finite_row(row, labels)), None)
                if bad is None:
                    passing.append((li, ri, edge))
                    by_state[(li, ri)].append(edge)
                else:
                    bad_sizes[len(bad)] += 1
    print("VERTICES", len(vertices), "RESTRICTIONS", tuple(sorted(Counter(map(len, rows)).items())))
    print("PLACEMENTS", len(all_q), len(all_r), len(all_b), "COMPAT", compatible, "LEFT", len(left_states), "INTERNAL", INTERNAL_FREE)
    print("SPARSE_BAD", sparse_bad, "TESTS", tests, "PASS", len(passing), "PASS_EDGE_HIST", tuple(sorted(Counter(len(v) for v in by_state.values()).items())), "BAD_SIZES", tuple(sorted(bad_sizes.items())))
    stronger = Counter()
    single_records = []
    for li, ri, edge in passing:
        base = dict(template)
        base.update(left_states[li])
        base.update(all_b[ri])
        for label in (3, 4, 5):
            labels = dict(base)
            labels[edge] = label
            if finite_on(labels, rows_by_edge[edge]):
                stronger[label] += 1
                single_records.append((li, ri, edge, label))
    print("STRONGER", tuple(sorted(stronger.items())))

    pair_trials = pair_passing = 0
    cross_internal_trials = cross_internal_passing = 0
    for (li, ri), cross_edges in by_state.items():
        base = dict(template)
        base.update(left_states[li])
        base.update(all_b[ri])
        for first, second in combinations(sorted(set(cross_edges)), 2):
            pair_trials += 1
            labels = dict(base)
            labels[first] = labels[second] = 3
            affected = tuple(row for row in rows if set(first) <= set(row) or set(second) <= set(row))
            pair_passing += finite_on(labels, affected)
        for cross_edge in sorted(set(cross_edges)):
            for internal in INTERNAL_FREE:
                cross_internal_trials += 1
                labels = dict(base)
                labels[cross_edge] = labels[internal] = 3
                affected = tuple(row for row in rows if set(cross_edge) <= set(row) or set(internal) <= set(row))
                cross_internal_passing += finite_on(labels, affected)
    print("DOUBLE_CROSS", pair_trials, pair_passing, "CROSS_INTERNAL", cross_internal_trials, cross_internal_passing)
    print("SINGLE_RECORDS", len(single_records))
    canonical_r = _canonical_lanner_states(R)
    canonical_b = _canonical_lanner_states(B)
    references = (
        _one_bridge_class_key(R, B, canonical_r[1], canonical_b[1], pair(R[0], B[0])),
        _one_bridge_class_key(R, B, canonical_r[4], canonical_b[1], pair(R[4], B[0])),
    )
    witnesses = _a5_witnesses(R, B)
    classes = Counter()
    endpoint_histogram = Counter()
    for li, ri, edge, label in single_records:
        r_endpoints = tuple(vertex for vertex in edge if vertex in R)
        b_endpoints = tuple(vertex for vertex in edge if vertex in B)
        assert label == 3 and len(r_endpoints) == len(b_endpoints) == 1
        graph = _one_bridge_class_key(R, B, left_states[li], all_b[ri], edge)
        matches = tuple(index for index, reference in enumerate(references) if graph == reference)
        if len(matches) != 1:
            raise AssertionError((li, ri, edge, matches))
        classes[witnesses[matches[0]].name] += 1
        endpoint_histogram[r_endpoints[0]] += 1
    print("CLASSES", tuple(sorted(classes.items())), "ENDPOINTS", tuple(sorted(endpoint_histogram.items())), "INERTIAS", tuple(w.inertia for w in witnesses))


if __name__ == "__main__":
    main()
