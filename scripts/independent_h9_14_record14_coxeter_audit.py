#!/usr/bin/env python3
"""Independent hostile replay of the record-14 compact-Coxeter no-go.

This script does not import the production record-14 Coxeter audit or the
generic disjoint-Lanner gate.  It reconstructs the symmetry quotient, local
finite-diagram census, exceptional weighted blocks, and exact inertias from
the low-level Lanner tables and the frozen record-14 nerve.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from typing import Mapping, Sequence

import networkx as nx
import sympy as sp

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_isomorphism_codes,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record14 import (
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    record14_facet_complements,
)


FULL_MASK = (1 << 14) - 1
LEFT = TARGET_MISSING_FIVES[2]
RIGHT = TARGET_MISSING_FIVES[1]
BLOCKS = ((3, 4, 5), (6, 7), (9, 10, 11))
EXPECTED_SUMMARY_SHA256 = (
    "b754de10cfd74822fb6e8708546b872d3b3e3cf9892a31e8af320db955567eec"
)


def _states(
    support: Sequence[int], codes: Sequence[tuple[int, ...]]
) -> tuple[dict[Pair, int], ...]:
    edges = tuple(combinations(sorted(support), 2))
    return tuple(dict(zip(edges, code, strict=True)) for code in codes)


def _all_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    return _states(support, lanner_placement_codes(len(tuple(support))))


def _canonical_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    return _states(support, lanner_isomorphism_codes(len(tuple(support))))


def _group() -> tuple[dict[int, int], ...]:
    result = []
    for images in product(*(tuple(permutations(block)) for block in BLOCKS)):
        move = {vertex: vertex for vertex in range(14)}
        for block, image in zip(BLOCKS, images, strict=True):
            move.update(zip(block, image, strict=True))
        result.append(move)
    return tuple(result)


def _state_key(state: Mapping[Pair, int], support: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        state[pair(a, b)] for a, b in combinations(sorted(support), 2)
    )


def _move_state(
    state: Mapping[Pair, int], move: Mapping[int, int]
) -> dict[Pair, int]:
    return {
        pair(move[a], move[b]): label for (a, b), label in state.items()
    }


def _quotient(
    support: Sequence[int], group: Sequence[Mapping[int, int]]
) -> tuple[tuple[dict[Pair, int], ...], tuple[tuple[int, int], ...]]:
    states = _all_states(support)
    seen = set()
    representatives = []
    sizes = []
    for state in states:
        if _state_key(state, support) in seen:
            continue
        orbit = {
            _state_key(_move_state(state, move), support) for move in group
        }
        seen.update(orbit)
        representatives.append(state)
        sizes.append(len(orbit))
    if len(seen) != 420:
        raise AssertionError("independent Lanner quotient is incomplete")
    return tuple(representatives), tuple(sorted(Counter(sizes).items()))


def _vertices() -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(
            tuple(
                vertex
                for vertex in range(14)
                if (FULL_MASK ^ complement) >> vertex & 1
            )
            for complement in record14_facet_complements()
        )
    )


def _restrictions(
    vertices: Sequence[Sequence[int]], union: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    selected = set(union)
    return tuple(
        sorted(
            {
                tuple(vertex for vertex in actual if vertex in selected)
                for actual in vertices
            },
            key=lambda item: (len(item), item),
        )
    )


def _labels(
    union: Sequence[int],
    left: Mapping[Pair, int],
    right: Mapping[Pair, int],
    cross: Mapping[Pair, int] | None = None,
) -> dict[Pair, int]:
    result = {pair(a, b): 2 for a, b in combinations(union, 2)}
    result.update(left)
    result.update(right)
    if cross:
        result.update(cross)
    return result


def _affected_finite(
    restrictions: Sequence[Sequence[int]],
    labels: Mapping[Pair, int],
    edges: Sequence[Pair],
) -> bool:
    return all(
        is_finite_coxeter_labeling(restriction, labels)
        for restriction in restrictions
        if any(set(edge) <= set(restriction) for edge in edges)
    )


def _weighted_graph(
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    edge: Pair,
) -> nx.Graph:
    union = tuple((*LEFT, *RIGHT))
    labels = _labels(union, left_state, right_state, {edge: 3})
    graph = nx.Graph()
    graph.add_nodes_from(union)
    for a, b in combinations(union, 2):
        graph.add_edge(a, b, label=labels[pair(a, b)])
    return graph


def _isomorphic(first: nx.Graph, second: nx.Graph) -> bool:
    return nx.is_isomorphic(
        first,
        second,
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("label", None),
    )


def _twice_gram(order: Sequence[int], labels: Mapping[Pair, int]) -> sp.Matrix:
    values = {
        2: sp.Integer(0),
        3: -sp.Integer(1),
        4: -sp.sqrt(2),
        5: -(1 + sp.sqrt(5)) / 2,
    }
    return sp.Matrix(
        [
            [
                2 if a == b else values[labels[pair(a, b)]]
                for b in order
            ]
            for a in order
        ]
    )


def _exact_inertia(
    order: Sequence[int],
    labels: Mapping[Pair, int],
    expected: Sequence[sp.Expr],
) -> dict[str, object]:
    matrix = _twice_gram(order, labels)
    observed = tuple(
        sp.factor(matrix[:size, :size].det()) for size in range(1, len(order) + 1)
    )
    expected = tuple(sp.sympify(value) for value in expected)
    if any(
        sp.simplify(actual - wanted) != 0
        for actual, wanted in zip(observed, expected, strict=True)
    ):
        raise AssertionError("independent record-14 leading minors changed")
    signs = []
    for value in observed:
        if value.is_positive is True:
            signs.append("+")
        elif value.is_negative is True:
            signs.append("-")
        else:
            raise AssertionError(f"independent sign undecided: {value}")
    augmented = "+" + "".join(signs)
    negative = sum(
        left != right
        for left, right in zip(augmented[:-1], augmented[1:], strict=True)
    )
    return {
        "leading_minors": [sp.sstr(value) for value in expected],
        "signs": "".join(signs),
        "inertia": [len(order) - negative, negative, 0],
    }


def replay() -> dict[str, object]:
    vertices = _vertices()
    group = _group()
    nonfaces = {
        frozenset(row)
        for row in (*TARGET_MISSING_EDGES, *TARGET_MISSING_FOURS, *TARGET_MISSING_FIVES)
    }
    vertex_set = {frozenset(vertex) for vertex in vertices}
    symmetry_ok = all(
        {
            frozenset(move[vertex] for vertex in row) for row in nonfaces
        }
        == nonfaces
        and {
            frozenset(move[vertex] for vertex in actual) for actual in vertices
        }
        == vertex_set
        for move in group
    )
    left_states, left_orbits = _quotient(LEFT, group)
    right_states, right_orbits = _quotient(RIGHT, group)
    if (
        len(vertices),
        len(group),
        symmetry_ok,
        len(left_states),
        len(right_states),
        left_orbits,
        right_orbits,
    ) != (
        106,
        72,
        True,
        39,
        73,
        ((6, 8), (12, 31)),
        ((3, 6), (6, 67)),
    ):
        raise AssertionError("independent record-14 symmetry census changed")

    union = tuple(sorted((*LEFT, *RIGHT)))
    restrictions = _restrictions(vertices, union)
    cross_pairs = tuple(pair(a, b) for a in LEFT for b in RIGHT)
    baseline_failures = 0
    bad_sizes: Counter[int] = Counter()
    passing = []
    pass_by_edge: Counter[Pair] = Counter({edge: 0 for edge in cross_pairs})
    for left_index, left_state in enumerate(left_states):
        for right_index, right_state in enumerate(right_states):
            baseline = _labels(union, left_state, right_state)
            if not all(
                is_finite_coxeter_labeling(restriction, baseline)
                for restriction in restrictions
            ):
                baseline_failures += 1
            for edge in cross_pairs:
                labels = dict(baseline)
                labels[edge] = 3
                bad = next(
                    (
                        restriction
                        for restriction in restrictions
                        if set(edge) <= set(restriction)
                        and not is_finite_coxeter_labeling(restriction, labels)
                    ),
                    None,
                )
                if bad is None:
                    passing.append((left_index, right_index, edge))
                    pass_by_edge[edge] += 1
                else:
                    bad_sizes[len(bad)] += 1

    stronger = Counter()
    for label in (4, 5):
        stronger[label] = sum(
            _affected_finite(
                restrictions,
                _labels(
                    union,
                    left_states[left_index],
                    right_states[right_index],
                    {edge: label},
                ),
                (edge,),
            )
            for left_index, right_index, edge in passing
        )
    by_state: defaultdict[tuple[int, int], list[Pair]] = defaultdict(list)
    for left_index, right_index, edge in passing:
        by_state[(left_index, right_index)].append(edge)
    residual_double = sum(
        len(tuple(combinations(set(edges), 2))) for edges in by_state.values()
    )

    triple_orbits = {}
    for left_index, right_index, edge in passing:
        orbit = {
            (
                _state_key(_move_state(left_states[left_index], move), LEFT),
                _state_key(_move_state(right_states[right_index], move), RIGHT),
                pair(move[edge[0]], move[edge[1]]),
            )
            for move in group
        }
        triple_orbits.setdefault(
            min(orbit), (left_index, right_index, edge, len(orbit))
        )

    canonical_left = _canonical_states(LEFT)
    canonical_right = _canonical_states(RIGHT)
    references = (
        _weighted_graph(
            canonical_left[1], canonical_right[1], pair(LEFT[0], RIGHT[0])
        ),
        _weighted_graph(
            canonical_left[4], canonical_right[1], pair(LEFT[4], RIGHT[0])
        ),
    )
    orbit_classes: Counter[str] = Counter()
    full_classes: Counter[str] = Counter()
    names = ("A5_A5_type_1_1", "A5_A5_type_4_1")
    for left_index, right_index, edge, orbit_size in triple_orbits.values():
        graph = _weighted_graph(
            left_states[left_index], right_states[right_index], edge
        )
        matches = tuple(
            index for index, reference in enumerate(references) if _isomorphic(graph, reference)
        )
        if len(matches) != 1:
            raise AssertionError("independent exceptional class changed")
        orbit_classes[names[matches[0]]] += 1
        full_classes[names[matches[0]]] += orbit_size

    labels_a = _labels(
        tuple((*LEFT, *RIGHT)),
        canonical_left[1],
        canonical_right[1],
        {pair(LEFT[0], RIGHT[0]): 3},
    )
    inertia_a = _exact_inertia(
        tuple((*LEFT, *RIGHT)),
        labels_a,
        (
            2,
            4,
            8,
            6 - 2 * sp.sqrt(5),
            4 - 2 * sp.sqrt(5),
            (9 - 5 * sp.sqrt(5)) / 2,
            9 - 5 * sp.sqrt(5),
            18 - 10 * sp.sqrt(5),
            26 - 12 * sp.sqrt(5),
            (25 - 11 * sp.sqrt(5)) / 2,
        ),
    )
    labels_b = _labels(
        tuple((*LEFT, *RIGHT)),
        canonical_left[4],
        canonical_right[1],
        {pair(LEFT[4], RIGHT[0]): 3},
    )
    inertia_b = _exact_inertia(
        tuple(reversed(tuple((*LEFT, *RIGHT)))),
        labels_b,
        (
            2,
            4,
            5 - sp.sqrt(5),
            (7 - 3 * sp.sqrt(5)) / 2,
            4 - 2 * sp.sqrt(5),
            (9 - 5 * sp.sqrt(5)) / 2,
            9 - 5 * sp.sqrt(5),
            9 - 5 * sp.sqrt(5),
            1 - sp.sqrt(5),
            (-23 - 16 * sp.sqrt(2) + 11 * sp.sqrt(5) + 8 * sp.sqrt(10)) / 2,
        ),
    )

    pass_histogram = tuple(sorted(Counter(pass_by_edge.values()).items()))
    triple_orbit_histogram = tuple(
        sorted(Counter(value[3] for value in triple_orbits.values()).items())
    )
    full_passing = sum(value[3] for value in triple_orbits.values())
    if (
        tuple(sorted(Counter(map(len, restrictions)).items()))
        != ((7, 47), (8, 25))
        or baseline_failures != 0
        or len(left_states) * len(right_states) * len(cross_pairs) != 71_175
        or len(passing) != 420
        or bad_sizes != {7: 70_454, 8: 301}
        or stronger != {4: 0, 5: 0}
        or residual_double != 0
        or len(triple_orbits) != 420
        or triple_orbit_histogram != ((36, 40), (72, 380))
        or full_passing != 28_800
        or orbit_classes != {names[0]: 200, names[1]: 220}
        or full_classes != {names[0]: 14_400, names[1]: 14_400}
        or inertia_a["inertia"] != [8, 2, 0]
        or inertia_b["inertia"] != [8, 2, 0]
    ):
        raise AssertionError("independent record-14 compact gate changed")

    result: dict[str, object] = {
        "actual_vertex_count": len(vertices),
        "symmetry_group_order": len(group),
        "symmetry_verified": symmetry_ok,
        "full_lanner_state_counts": [420, 420],
        "quotient_state_counts": [len(left_states), len(right_states)],
        "state_orbit_histograms": [dict(left_orbits), dict(right_orbits)],
        "restriction_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, restrictions)).items())
        },
        "quotient_single_edge_test_count": 71_175,
        "quotient_single_edge_passing_count": len(passing),
        "quotient_pass_per_cross_pair_histogram": {
            str(count): multiplicity for count, multiplicity in pass_histogram
        },
        "first_bad_restriction_size_histogram": {
            str(size): count for size, count in sorted(bad_sizes.items())
        },
        "full_single_edge_test_count": 4_410_000,
        "full_single_edge_passing_count": full_passing,
        "passing_triple_orbit_size_histogram": dict(triple_orbit_histogram),
        "stronger_label_survivors": dict(stronger),
        "full_raw_double_edge_trial_count": 52_920_000,
        "quotient_residual_double_edge_trial_count": residual_double,
        "exceptional_orbit_class_histogram": dict(sorted(orbit_classes.items())),
        "exceptional_full_class_histogram": dict(sorted(full_classes.items())),
        "exact_exceptional_inertias": {
            names[0]: inertia_a,
            names[1]: inertia_b,
        },
        "unrestricted_compact_coxeter_no_go": True,
        "selected_nerve_only": True,
        "finite_volume_noncompact_checked": False,
    }
    digest = sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    result["summary_sha256"] = digest
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("independent record-14 compact-gate digest changed")
    return result


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
