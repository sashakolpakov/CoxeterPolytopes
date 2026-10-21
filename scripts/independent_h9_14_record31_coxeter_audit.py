#!/usr/bin/env python3
"""Independent replay of the selected record-31 compact-Coxeter no-go.

This verifier imports neither the production record-31 Coxeter audit nor the
generic disjoint-Lanner gate.  It reads the frozen selected nerve, rebuilds
the clone-symmetry quotient, and independently runs the actual-vertex,
rooted-block, residual-double-edge, and seven-support Lanner-CSP checks.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Mapping, Sequence

import sympy as sp

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_isomorphism_codes,
    lanner_placement_codes,
    pair,
)


ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / "data/enumerations/h9_14_record31_polytopality.json"
EXPECTED_SELECTION_CERTIFICATE_SHA256 = (
    "12b01aec7f3de0bb4e71040291753e8362787b9c7063704877bd24693884a799"
)
EXPECTED_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EXPECTED_FIVES = (
    (0, 5, 6, 9, 10),
    (2, 3, 4, 7, 8),
    (2, 5, 6, 9, 10),
    (2, 5, 6, 11, 12),
    (2, 7, 8, 11, 12),
    (3, 4, 7, 8, 11),
    (3, 4, 7, 8, 13),
    (3, 4, 9, 10, 13),
    (5, 6, 9, 10, 13),
)
FULL_MASK = (1 << 14) - 1
EXPECTED_SUMMARY_SHA256 = (
    "4cab6a9ed8ad1ae6e921a91c41e5de00e65495b627c627f7529e6444af5b8c05"
)


def _selection_data() -> tuple[
    tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...], tuple[int, ...]
]:
    document = json.loads(SELECTION.read_text(encoding="utf-8"))
    if document.get("certificate_sha256") != EXPECTED_SELECTION_CERTIFICATE_SHA256:
        raise AssertionError("independent record-31 selection source changed")
    target = document["target"]
    edges = tuple(map(tuple, target["missing_edges"]))
    fives = tuple(map(tuple, target["missing_fives"]))
    facets = tuple(target["facet_complement_masks"])
    if edges != EXPECTED_EDGES or fives != EXPECTED_FIVES or len(facets) != 122:
        raise AssertionError("independent record-31 selected nerve changed")
    return edges, fives, facets


def _states(
    support: Sequence[int], codes: Sequence[tuple[int, ...]]
) -> tuple[dict[Pair, int], ...]:
    edges = tuple(combinations(sorted(support), 2))
    return tuple(dict(zip(edges, code, strict=True)) for code in codes)


def _all_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    return _states(support, lanner_placement_codes(5))


def _canonical_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    return _states(support, lanner_isomorphism_codes(5))


def _state_key(
    state: Mapping[Pair, int], support: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        state[pair(a, b)] for a, b in combinations(sorted(support), 2)
    )


def _move_state(
    state: Mapping[Pair, int], move: Mapping[int, int]
) -> dict[Pair, int]:
    return {pair(move[a], move[b]): label for (a, b), label in state.items()}


def _clone_group() -> tuple[dict[int, int], ...]:
    swaps = ((3, 4), (5, 6), (7, 8), (9, 10))
    result = []
    for bits in range(16):
        move = {vertex: vertex for vertex in range(14)}
        for index, (first, second) in enumerate(swaps):
            if bits & (1 << index):
                move[first], move[second] = second, first
        result.append(move)
    return tuple(result)


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
        raise AssertionError("independent record-31 quotient is incomplete")
    return tuple(representatives), tuple(sorted(Counter(sizes).items()))


def _actual_vertices(facets: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(
            tuple(
                vertex
                for vertex in range(14)
                if (FULL_MASK ^ complement) & (1 << vertex)
            )
            for complement in facets
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
            key=lambda value: (len(value), value),
        )
    )


def _labels(
    union: Sequence[int],
    left: Mapping[Pair, int],
    right: Mapping[Pair, int],
    cross: Mapping[Pair, int] | None = None,
) -> dict[Pair, int]:
    labels = {pair(a, b): 2 for a, b in combinations(union, 2)}
    labels.update(left)
    labels.update(right)
    if cross:
        labels.update(cross)
    return labels


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


def _rooted_code(
    state: Mapping[Pair, int], support: Sequence[int], root: int
) -> tuple[int, ...]:
    remainder = tuple(vertex for vertex in support if vertex != root)
    return min(
        tuple(
            state[pair(order[a], order[b])]
            for a, b in combinations(range(5), 2)
        )
        for image in permutations(remainder)
        for order in ((root, *image),)
    )


def _block_code(
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    edge: Pair,
    left: Sequence[int],
    right: Sequence[int],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    left_root = edge[0] if edge[0] in left else edge[1]
    right_root = edge[1] if edge[1] in right else edge[0]
    return tuple(
        sorted(
            (
                _rooted_code(left_state, left, left_root),
                _rooted_code(right_state, right, right_root),
            )
        )
    )  # type: ignore[return-value]


def _reference_codes(
    left: Sequence[int], right: Sequence[int]
) -> tuple[tuple[str, tuple[tuple[int, ...], tuple[int, ...]]], ...]:
    left_states = _canonical_states(left)
    right_states = _canonical_states(right)
    specifications = (
        ("A5_A5_type_0_1", 0, 0, 1, 0),
        ("A5_A5_type_1_1", 1, 0, 1, 0),
        ("A5_A5_type_1_1_opposite_roots", 1, 4, 1, 0),
        ("A5_A5_type_4_1", 4, 4, 1, 0),
    )
    return tuple(
        (
            name,
            _block_code(
                left_states[left_type],
                right_states[right_type],
                pair(left[left_root], right[right_root]),
                left,
                right,
            ),
        )
        for name, left_type, left_root, right_type, right_root in specifications
    )


def _csp(
    fixed: dict[Pair, int],
    domains: Sequence[Sequence[Mapping[Pair, int]]],
) -> tuple[bool, int]:
    nodes = 0

    def compatible(state: Mapping[Pair, int]) -> bool:
        return all(
            edge not in fixed or fixed[edge] == label
            for edge, label in state.items()
        )

    def visit(remaining: tuple[int, ...]) -> bool:
        nonlocal nodes
        nodes += 1
        if not remaining:
            return True
        candidates = []
        for index in remaining:
            domain = tuple(state for state in domains[index] if compatible(state))
            if not domain:
                return False
            candidates.append((len(domain), index, domain))
        _size, chosen, domain = min(candidates, key=lambda value: (value[0], value[1]))
        rest = tuple(index for index in remaining if index != chosen)
        for state in domain:
            added = {edge: label for edge, label in state.items() if edge not in fixed}
            fixed.update(added)
            if visit(rest):
                return True
            for edge in added:
                fixed.pop(edge)
        return False

    return visit(tuple(range(len(domains)))), nodes


def _twice_gram(order: Sequence[int], labels: Mapping[Pair, int]) -> sp.Matrix:
    value = {
        2: sp.Integer(0),
        3: -sp.Integer(1),
        4: -sp.sqrt(2),
        5: -(1 + sp.sqrt(5)) / 2,
    }
    return sp.Matrix(
        [
            [2 if a == b else value[labels[pair(a, b)]] for b in order]
            for a in order
        ]
    )


def _exact_inertia(
    order: Sequence[int],
    labels: Mapping[Pair, int],
    expected: Sequence[sp.Expr],
) -> tuple[str, tuple[int, int, int]]:
    matrix = _twice_gram(order, labels)
    observed = tuple(
        sp.factor(matrix[:size, :size].det()) for size in range(1, len(order) + 1)
    )
    expected = tuple(sp.sympify(value) for value in expected)
    if any(
        sp.simplify(actual - wanted) != 0
        for actual, wanted in zip(observed, expected, strict=True)
    ):
        raise AssertionError("independent record-31 leading minors changed")
    signs = []
    for value in observed:
        if value.is_positive is True:
            signs.append("+")
        elif value.is_negative is True:
            signs.append("-")
        else:
            raise AssertionError(f"independent record-31 sign undecided: {value}")
    word = "".join(signs)
    augmented = "+" + word
    negative = sum(a != b for a, b in zip(augmented[:-1], augmented[1:], strict=True))
    return word, (len(order) - negative, negative, 0)


@lru_cache(maxsize=1)
def replay() -> dict[str, object]:
    edges, fives, facets = _selection_data()
    left, right = fives[:2]
    residual = fives[2:]
    vertices = _actual_vertices(facets)
    group = _clone_group()
    nonfaces = {frozenset(face) for face in (*edges, *fives)}
    vertex_set = {frozenset(vertex) for vertex in vertices}
    symmetry_ok = all(
        {
            frozenset(move[vertex] for vertex in face) for face in nonfaces
        }
        == nonfaces
        and {
            frozenset(move[vertex] for vertex in actual) for actual in vertex_set
        }
        == vertex_set
        and {move[vertex] for vertex in left} == set(left)
        and {move[vertex] for vertex in right} == set(right)
        for move in group
    )
    left_states, left_orbits = _quotient(left, group)
    right_states, right_orbits = _quotient(right, group)
    if (
        len(vertices),
        len(group),
        symmetry_ok,
        len(left_states),
        len(right_states),
        left_orbits,
        right_orbits,
    ) != (122, 16, True, 110, 110, ((2, 10), (4, 100)), ((2, 10), (4, 100))):
        raise AssertionError("independent record-31 symmetry census changed")

    union = tuple(sorted((*left, *right)))
    restrictions = _restrictions(vertices, union)
    cross_pairs = tuple(pair(a, b) for a in left for b in right)
    bad_sizes: Counter[int] = Counter()
    pass_by_edge: Counter[Pair] = Counter({edge: 0 for edge in cross_pairs})
    passing = []
    for left_index, left_state in enumerate(left_states):
        for right_index, right_state in enumerate(right_states):
            baseline = _labels(union, left_state, right_state)
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

    by_state: defaultdict[tuple[int, int], set[Pair]] = defaultdict(set)
    for left_index, right_index, edge in passing:
        by_state[(left_index, right_index)].add(edge)
    residual_pairs = tuple(
        (left_index, right_index, edges)
        for (left_index, right_index), values in sorted(by_state.items())
        for edges in combinations(sorted(values), 2)
    )
    double_passing = sum(
        _affected_finite(
            restrictions,
            _labels(
                union,
                left_states[left_index],
                right_states[right_index],
                {edge: 3 for edge in edges_},
            ),
            edges_,
        )
        for left_index, right_index, edges_ in residual_pairs
    )

    references = _reference_codes(left, right)
    classes: Counter[str] = Counter()
    classified = []
    for left_index, right_index, edge in passing:
        code = _block_code(
            left_states[left_index], right_states[right_index], edge, left, right
        )
        matches = tuple(name for name, reference in references if code == reference)
        if len(matches) != 1:
            raise AssertionError("independent record-31 rooted class changed")
        classes[matches[0]] += 1
        classified.append((left_index, right_index, edge, matches[0]))

    residual_names = {"A5_A5_type_0_1", "A5_A5_type_1_1_opposite_roots"}
    domains = tuple(_all_states(support) for support in residual)
    csp_nodes = Counter()
    csp_completing = 0
    csp_trials = 0
    for left_index, right_index, edge, name in classified:
        if name not in residual_names:
            continue
        csp_trials += 1
        fixed = _labels(
            union,
            left_states[left_index],
            right_states[right_index],
            {edge: 3},
        )
        completing, nodes = _csp(fixed, domains)
        csp_nodes[nodes] += 1
        csp_completing += completing

    canonical_left = _canonical_states(left)
    canonical_right = _canonical_states(right)
    labels_a = _labels(
        (*left, *right),
        canonical_left[1],
        canonical_right[1],
        {pair(left[0], right[0]): 3},
    )
    signs_a, inertia_a = _exact_inertia(
        (*left, *right),
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
        (*left, *right),
        canonical_left[4],
        canonical_right[1],
        {pair(left[4], right[0]): 3},
    )
    reverse_order = tuple(reversed((*left, *right)))
    signs_b, inertia_b = _exact_inertia(
        reverse_order,
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

    expected_classes = {
        "A5_A5_type_0_1": 36,
        "A5_A5_type_1_1": 900,
        "A5_A5_type_1_1_opposite_roots": 36,
        "A5_A5_type_4_1": 960,
    }
    if (
        tuple(sorted(Counter(map(len, restrictions)).items()))
        != ((6, 28), (7, 20), (8, 21))
        or len(passing) != 1_932
        or bad_sizes != {6: 214_866, 7: 17_754, 8: 67_948}
        or stronger != {4: 0, 5: 0}
        or len(residual_pairs) != 36
        or double_passing != 0
        or classes != expected_classes
        or csp_trials != 72
        or csp_completing != 0
        or csp_nodes != {1: 50, 4: 12, 23: 2, 38: 2, 63: 2, 68: 2, 113: 2}
        or signs_a != "++++-----+"
        or signs_b != "++++-----+"
        or inertia_a != (8, 2, 0)
        or inertia_b != (8, 2, 0)
    ):
        raise AssertionError("independent record-31 compact gate changed")

    result: dict[str, object] = {
        "actual_vertex_count": len(vertices),
        "disjoint_lanner_supports": [list(left), list(right)],
        "symmetry_group_order": len(group),
        "symmetry_verified": symmetry_ok,
        "full_lanner_state_counts": [420, 420],
        "orbit_representative_counts": [len(left_states), len(right_states)],
        "state_orbit_size_histograms": [dict(left_orbits), dict(right_orbits)],
        "restriction_size_histogram": dict(Counter(map(len, restrictions))),
        "quotient_single_edge_trial_count": 110 * 110 * 25,
        "quotient_single_edge_passing_count": len(passing),
        "passing_count_by_cross_pair": [
            {"pair": list(edge), "count": pass_by_edge[edge]}
            for edge in cross_pairs
        ],
        "first_bad_restriction_size_histogram": dict(bad_sizes),
        "stronger_label_survivors": dict(stronger),
        "residual_double_edge_trials": len(residual_pairs),
        "residual_double_edge_passing": double_passing,
        "exceptional_class_histogram": dict(classes),
        "index_two_exceptional_count": (
            classes["A5_A5_type_1_1"] + classes["A5_A5_type_4_1"]
        ),
        "residual_lanner_csp_trials": csp_trials,
        "residual_lanner_csp_completing": csp_completing,
        "residual_lanner_csp_node_histogram": dict(csp_nodes),
        "exact_index_two_sign_words": [signs_a, signs_b],
        "exact_index_two_inertias": [list(inertia_a), list(inertia_b)],
        "unrestricted_compact_coxeter_no_go": True,
        "selected_nerve_only": True,
        "finite_volume_noncompact_checked": False,
    }
    payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    result["summary_sha256"] = sha256(payload).hexdigest()
    return result


def main() -> None:
    result = replay()
    if (
        EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN"
        and result["summary_sha256"] != EXPECTED_SUMMARY_SHA256
    ):
        raise AssertionError("independent record-31 summary is not frozen")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
