#!/usr/bin/env python3
"""Independent hostile replay of the selected H9/14 compact no-gos.

This verifier intentionally does not import the aggregate audit module or the
generic disjoint-Lanner gate.  It rebuilds the placement joins, actual-vertex
restrictions, one/two-cross-edge reductions, symmetry quotients, weighted
block classes, and exact leading-minor inertias from lower-level tables.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
import json
from pathlib import Path
from time import perf_counter
from typing import Iterable, Mapping, Sequence

import networkx as nx
import sympy as sp

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_isomorphism_codes,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record32 import (
    record32_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record43 import (
    TARGET_MISSING_EDGES as RECORD43_EDGES,
    TARGET_MISSING_FIVES as RECORD43_FIVES,
    record43_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record44 import (
    TARGET_MISSING_EDGES as RECORD44_EDGES,
    TARGET_MISSING_FIVES as RECORD44_FIVES,
    record44_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record47 import (
    TARGET_MISSING_EDGES as RECORD47_EDGES,
    TARGET_MISSING_FIVES as RECORD47_FIVES,
    record47_facet_complements,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_selected_nerve_coxeter_no_go.json"
RECORD35_ARCHIVE = ROOT / "data/enumerations/h9_14_record35_coxeter_labels.json"
FULL_MASK = (1 << 14) - 1
RECORD30_NONFACES = (
    (0, 11),
    (1, 12),
    (2, 13),
    (3, 4, 5, 6),
    (7, 8, 9, 10),
)
RECORD42_NONFACES = (
    (0, 11),
    (0, 12),
    (1, 12),
    (1, 13),
    (11, 13),
    (2, 3, 4, 5),
    (6, 7, 8, 9, 10),
)
RECORD32_A4S = ((2, 3, 4, 5), (3, 4, 5, 11))
RECORD32_A5 = (6, 7, 8, 9, 10)


@dataclass(frozen=True, slots=True)
class GateReplay:
    restriction_size_histogram: tuple[tuple[int, int], ...]
    single_test_count: int
    passing_trials: tuple[tuple[int, int, Pair], ...]
    passing_per_cross_pair_histogram: tuple[tuple[int, int], ...]
    first_bad_restriction_size_histogram: tuple[tuple[int, int], ...]
    stronger_survivors: tuple[tuple[int, int], ...]
    raw_double_trials: int
    residual_double_trials: int
    double_survivors: int
    baseline_failure_count: int


def _vertices_from_nonfaces(
    nonfaces: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(face) for face in nonfaces)
    return tuple(
        subset
        for subset in combinations(range(14), 9)
        if not any(face <= set(subset) for face in forbidden)
    )


def _vertices_from_complements(complements: Iterable[int]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(
            tuple(vertex for vertex in range(14) if (FULL_MASK ^ value) >> vertex & 1)
            for value in complements
        )
    )


def _states(
    support: Sequence[int], codes: Sequence[tuple[int, ...]]
) -> tuple[dict[Pair, int], ...]:
    ordered = tuple(sorted(support))
    edges = tuple(combinations(ordered, 2))
    return tuple(dict(zip(edges, code, strict=True)) for code in codes)


def _canonical_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    return _states(support, lanner_isomorphism_codes(len(tuple(support))))


def _all_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    return _states(support, lanner_placement_codes(len(tuple(support))))


def _restriction_family(
    vertices: Sequence[Sequence[int]], union: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    union_set = set(union)
    return tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in union_set)
                for vertex in vertices
            },
            key=lambda item: (len(item), item),
        )
    )


def _labels(
    union: Sequence[int],
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    cross: Mapping[Pair, int] = {},
) -> dict[Pair, int]:
    result = {pair(a, b): 2 for a, b in combinations(tuple(union), 2)}
    result.update(left_state)
    result.update(right_state)
    result.update(cross)
    return result


def _affected_finite(
    restrictions: Sequence[Sequence[int]],
    labels: Mapping[Pair, int],
    affected_edges: Sequence[Pair],
) -> bool:
    return all(
        is_finite_coxeter_labeling(restriction, labels)
        for restriction in restrictions
        if any(set(edge) <= set(restriction) for edge in affected_edges)
    )


def _replay_gate(
    vertices: Sequence[Sequence[int]],
    left: Sequence[int],
    right: Sequence[int],
    left_states: Sequence[Mapping[Pair, int]],
    right_states: Sequence[Mapping[Pair, int]],
) -> GateReplay:
    left = tuple(sorted(left))
    right = tuple(sorted(right))
    union = tuple(sorted((*left, *right)))
    restrictions = _restriction_family(vertices, union)
    cross_pairs = tuple(pair(a, b) for a in left for b in right)

    baseline_failures = 0
    passing: list[tuple[int, int, Pair]] = []
    pass_by_edge: Counter[Pair] = Counter({edge: 0 for edge in cross_pairs})
    bad_sizes: Counter[int] = Counter()
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

    stronger = []
    for label in (4, 5):
        count = 0
        for left_index, right_index, edge in passing:
            labels = _labels(
                union,
                left_states[left_index],
                right_states[right_index],
                {edge: label},
            )
            count += _affected_finite(restrictions, labels, (edge,))
        stronger.append((label, count))

    by_state: defaultdict[tuple[int, int], list[Pair]] = defaultdict(list)
    for left_index, right_index, edge in passing:
        by_state[(left_index, right_index)].append(edge)
    residual = double_passing = 0
    for (left_index, right_index), edges in by_state.items():
        for first, second in combinations(sorted(set(edges)), 2):
            residual += 1
            labels = _labels(
                union,
                left_states[left_index],
                right_states[right_index],
                {first: 3, second: 3},
            )
            double_passing += _affected_finite(
                restrictions, labels, (first, second)
            )

    return GateReplay(
        restriction_size_histogram=tuple(
            sorted(Counter(map(len, restrictions)).items())
        ),
        single_test_count=len(left_states) * len(right_states) * len(cross_pairs),
        passing_trials=tuple(passing),
        passing_per_cross_pair_histogram=tuple(
            sorted(Counter(pass_by_edge.values()).items())
        ),
        first_bad_restriction_size_histogram=tuple(sorted(bad_sizes.items())),
        stronger_survivors=tuple(stronger),
        raw_double_trials=(
            len(left_states)
            * len(right_states)
            * len(tuple(combinations(cross_pairs, 2)))
        ),
        residual_double_trials=residual,
        double_survivors=double_passing,
        baseline_failure_count=baseline_failures,
    )


def _permuted_code(
    order: int, code: tuple[int, ...], permutation: tuple[int, ...]
) -> tuple[int, ...]:
    local_pairs = tuple(combinations(range(order), 2))
    source = dict(zip(local_pairs, code, strict=True))
    moved = {
        pair(permutation[a], permutation[b]): label
        for (a, b), label in source.items()
    }
    return tuple(moved[edge] for edge in local_pairs)


def _audit_lanner_quotient(order: int) -> dict[str, object]:
    full = lanner_placement_codes(order)
    canonical = lanner_isomorphism_codes(order)
    orbit_counts = Counter(
        min(_permuted_code(order, code, image) for image in permutations(range(order)))
        for code in full
    )
    if tuple(sorted(orbit_counts)) != canonical:
        raise AssertionError(f"independent order-{order} Lanner quotient changed")
    return {
        "placement_count": len(full),
        "isomorphism_type_count": len(canonical),
        "orbit_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(orbit_counts.values()).items())
        },
    }


def _permutation_preserves(
    nonfaces: Sequence[Sequence[int]],
    vertices: Sequence[Sequence[int]],
    block: Sequence[int],
) -> bool:
    target_nonfaces = {frozenset(face) for face in nonfaces}
    target_vertices = {frozenset(vertex) for vertex in vertices}
    for index in range(len(block) - 1):
        first, second = block[index : index + 2]
        image = {vertex: vertex for vertex in range(14)}
        image[first], image[second] = second, first
        if {
            frozenset(image[vertex] for vertex in face) for face in target_nonfaces
        } != target_nonfaces:
            return False
        if {
            frozenset(image[vertex] for vertex in face) for face in target_vertices
        } != target_vertices:
            return False
    return True


def _weighted_graph(
    left: Sequence[int],
    right: Sequence[int],
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    edge: Pair,
) -> nx.Graph:
    union = tuple((*left, *right))
    labels = _labels(union, left_state, right_state, {edge: 3})
    graph = nx.Graph()
    graph.add_nodes_from(union)
    for first, second in combinations(union, 2):
        graph.add_edge(first, second, label=labels[pair(first, second)])
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
                2 if left == right else values[labels[pair(left, right)]]
                for right in order
            ]
            for left in order
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
        raise AssertionError("independent exact leading-minor formula changed")
    signs = []
    for value in observed:
        if value.is_positive is True:
            signs.append("+")
        elif value.is_negative is True:
            signs.append("-")
        else:
            raise AssertionError(f"independent exact sign undecided: {value}")
    augmented = "+" + "".join(signs)
    negative = sum(a != b for a, b in zip(augmented[:-1], augmented[1:], strict=True))
    inertia = (len(order) - negative, negative, 0)
    return {
        "leading_minors": [sp.sstr(value) for value in observed],
        "signs": "".join(signs),
        "inertia": list(inertia),
    }


def _record32_replay() -> dict[str, object]:
    vertices = _vertices_from_complements(record32_facet_complements())
    first, second = RECORD32_A4S
    first_states = _all_states(first)
    second_states = _all_states(second)
    merged = set()
    for left_state in first_states:
        for right_state in second_states:
            labels = dict(left_state)
            compatible = True
            for edge, label in right_state.items():
                if edge in labels and labels[edge] != label:
                    compatible = False
                    break
                labels[edge] = label
            if compatible:
                merged.add(tuple(sorted(labels.items())))
    union = tuple(sorted(set(first) | set(second)))
    restrictions = _restriction_family(vertices, union)
    passing = []
    for state in sorted(merged):
        labels = {pair(a, b): 2 for a, b in combinations(union, 2)}
        labels.update(dict(state))
        if all(is_finite_coxeter_labeling(face, labels) for face in restrictions):
            passing.append(dict(state))
    projected = tuple(
        {
            pair(a, b): labels[pair(a, b)]
            for a, b in combinations(first, 2)
        }
        for labels in passing
    )
    projection_multiplicities = Counter(
        tuple(sorted(state.items())) for state in projected
    )
    # The 48 objects are genuinely distinct full joins on the two overlapping
    # A4 supports.  Forgetting the second support gives 24 first-A4 states,
    # each with two distinct extensions.  The production 403,200 count keeps
    # this extension multiplicity; it is redundant but exhaustive.
    if (
        len({tuple(sorted(state)) for state in merged}) != 546
        or Counter(projection_multiplicities.values()) != {2: 24}
    ):
        raise AssertionError("record-32 A4 projection multiplicities changed")
    a5_states = _all_states(RECORD32_A5)
    gate = _replay_gate(vertices, first, RECORD32_A5, projected, a5_states)
    if (
        len(vertices),
        len(merged),
        len(passing),
        len(a5_states),
        gate.single_test_count,
        len(gate.passing_trials),
    ) != (110, 546, 48, 420, 403_200, 0):
        raise AssertionError("independent record-32 census changed")
    return {
        "vertex_count": len(vertices),
        "A4_compatible_join_count": len(merged),
        "A4_passing_count": len(passing),
        "distinct_FIRST_A4_projection_count": len(projection_multiplicities),
        "FIRST_A4_projection_multiplicity_histogram": {
            str(size): count
            for size, count in sorted(Counter(projection_multiplicities.values()).items())
        },
        "A5_state_count": len(a5_states),
        "gate": _gate_record(gate),
        "unrestricted_compact_coxeter_no_go": True,
    }


def _gate_record(gate: GateReplay) -> dict[str, object]:
    return {
        "restriction_size_histogram": {
            str(size): count for size, count in gate.restriction_size_histogram
        },
        "single_test_count": gate.single_test_count,
        "single_passing_count": len(gate.passing_trials),
        "passing_per_cross_pair_histogram": {
            str(count): multiplicity
            for count, multiplicity in gate.passing_per_cross_pair_histogram
        },
        "first_bad_restriction_size_histogram": {
            str(size): count
            for size, count in gate.first_bad_restriction_size_histogram
        },
        "stronger_survivors": {
            str(label): count for label, count in gate.stronger_survivors
        },
        "raw_double_trials": gate.raw_double_trials,
        "residual_double_trials": gate.residual_double_trials,
        "double_survivors": gate.double_survivors,
        "baseline_failure_count": gate.baseline_failure_count,
    }


def replay() -> dict[str, object]:
    started = perf_counter()
    quotient = {str(order): _audit_lanner_quotient(order) for order in (4, 5)}

    vertices30 = _vertices_from_nonfaces(RECORD30_NONFACES)
    left30, right30 = RECORD30_NONFACES[-2:]
    if not all(
        _permutation_preserves(RECORD30_NONFACES, vertices30, block)
        for block in (left30, right30)
    ):
        raise AssertionError("record-30 quotient symmetry failed")
    gate30 = _replay_gate(
        vertices30,
        left30,
        right30,
        _canonical_states(left30),
        _canonical_states(right30),
    )
    if (
        len(vertices30),
        gate30.restriction_size_histogram,
        gate30.single_test_count,
        len(gate30.passing_trials),
        gate30.baseline_failure_count,
    ) != (128, ((6, 16),), 1296, 0, 0):
        raise AssertionError("independent record-30 gate changed")

    vertices42 = _vertices_from_nonfaces(RECORD42_NONFACES)
    left42, right42 = RECORD42_NONFACES[-2:]
    if not all(
        _permutation_preserves(RECORD42_NONFACES, vertices42, block)
        for block in (left42, right42)
    ):
        raise AssertionError("record-42 quotient symmetry failed")
    left42_states = _canonical_states(left42)
    right42_states = _canonical_states(right42)
    gate42 = _replay_gate(
        vertices42, left42, right42, left42_states, right42_states
    )
    expected42 = ((4, 1, pair(2, 6)), (4, 1, pair(4, 6)))
    if (
        len(vertices42),
        gate42.single_test_count,
        gate42.passing_trials,
        gate42.stronger_survivors,
        gate42.raw_double_trials,
        gate42.residual_double_trials,
        gate42.double_survivors,
        gate42.baseline_failure_count,
    ) != (100, 900, expected42, ((4, 0), (5, 0)), 8550, 1, 0, 0):
        raise AssertionError("independent record-42 gate changed")
    labels42 = _labels(
        tuple((*left42, *right42)),
        left42_states[4],
        right42_states[1],
        {pair(2, 6): 3},
    )
    inertia42 = _exact_inertia(
        (9, 3, 7, 4, 10, 5, 6, 8, 2),
        labels42,
        (
            2,
            4,
            6,
            9,
            12,
            8,
            10,
            8 - 4 * sp.sqrt(5),
            -8 * sp.sqrt(2) - 11 + 5 * sp.sqrt(5) + 4 * sp.sqrt(10),
        ),
    )
    if inertia42["inertia"] != [7, 2, 0]:
        raise AssertionError("independent record-42 inertia changed")

    multi_specs = (
        (43, RECORD43_EDGES, RECORD43_FIVES, record43_facet_complements()),
        (44, RECORD44_EDGES, RECORD44_FIVES, record44_facet_complements()),
        (47, RECORD47_EDGES, RECORD47_FIVES, record47_facet_complements()),
    )
    reference_left, reference_right = RECORD44_FIVES[:2]
    reference_left_states = _canonical_states(reference_left)
    reference_right_states = _canonical_states(reference_right)
    reference_graphs = (
        _weighted_graph(
            reference_left,
            reference_right,
            reference_left_states[1],
            reference_right_states[1],
            pair(reference_left[0], reference_right[0]),
        ),
        _weighted_graph(
            reference_left,
            reference_right,
            reference_left_states[4],
            reference_right_states[1],
            pair(reference_left[4], reference_right[0]),
        ),
    )
    labels_a = _labels(
        tuple((*reference_left, *reference_right)),
        reference_left_states[1],
        reference_right_states[1],
        {pair(reference_left[0], reference_right[0]): 3},
    )
    inertia_a = _exact_inertia(
        tuple((*reference_left, *reference_right)),
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
        tuple((*reference_left, *reference_right)),
        reference_left_states[4],
        reference_right_states[1],
        {pair(reference_left[4], reference_right[0]): 3},
    )
    inertia_b = _exact_inertia(
        tuple(reversed(tuple((*reference_left, *reference_right)))),
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
    if inertia_a["inertia"] != [8, 2, 0] or inertia_b["inertia"] != [8, 2, 0]:
        raise AssertionError("independent A5 exceptional inertia changed")

    multi_records = []
    for record_index, missing_edges, missing_fives, complements in multi_specs:
        vertices = _vertices_from_complements(complements)
        left = tuple(missing_fives[0])
        all_nonfaces = tuple((*missing_edges, *missing_fives))
        if not _permutation_preserves(all_nonfaces, vertices, left):
            raise AssertionError(f"record-{record_index} S5 quotient symmetry failed")
        left_states = _canonical_states(left)
        gates = []
        class_total = Counter()
        for right_raw in missing_fives[1:]:
            right = tuple(right_raw)
            right_states = _all_states(right)
            gate = _replay_gate(vertices, left, right, left_states, right_states)
            if (
                gate.single_test_count != 52_500
                or len(gate.passing_trials) != 300
                or gate.passing_per_cross_pair_histogram
                != ((0, 15), (24, 5), (36, 5))
                or gate.stronger_survivors != ((4, 0), (5, 0))
                or gate.raw_double_trials != 630_000
                or gate.double_survivors != 0
                or gate.baseline_failure_count != 0
            ):
                raise AssertionError(f"independent record-{record_index} gate changed")
            for left_index, right_index, edge in gate.passing_trials:
                graph = _weighted_graph(
                    left,
                    right,
                    left_states[left_index],
                    right_states[right_index],
                    edge,
                )
                matches = tuple(
                    index
                    for index, reference in enumerate(reference_graphs)
                    if _isomorphic(graph, reference)
                )
                if len(matches) != 1:
                    raise AssertionError("independent A5 block classification failed")
                class_total[("A5_A5_type_1_1", "A5_A5_type_4_1")[matches[0]]] += 1
            gates.append(gate)
        if class_total != {
            "A5_A5_type_1_1": 360,
            "A5_A5_type_4_1": 540,
        }:
            raise AssertionError(f"independent record-{record_index} class census changed")
        multi_records.append(
            {
                "record_index": record_index,
                "vertex_count": len(vertices),
                "left_S5_symmetry_verified_on_nonfaces_and_vertices": True,
                "normalizing_support_symmetry_order": 120,
                "gates": [_gate_record(gate) for gate in gates],
                "exceptional_class_histogram": dict(sorted(class_total.items())),
                "unrestricted_compact_coxeter_no_go": True,
            }
        )

    record32 = _record32_replay()
    record35 = json.loads(RECORD35_ARCHIVE.read_text(encoding="utf-8"))
    if (
        record35.get("certificate_sha256")
        != "fcb5b0374c1ee68b32193491dcb1d905a7ac39ef61e7c56a7bbfe6bf6b56d23c"
        or record35.get("theorem", {}).get("compact_coxeter_realization_exists")
        is not False
    ):
        raise AssertionError("the independently chained record-35 certificate changed")

    result = {
        "lanner_quotient_replay": quotient,
        "record30": {
            "vertex_count": len(vertices30),
            "internal_S4_symmetries_verified": True,
            "normalizing_support_symmetry_orders": [24, 24],
            "gate": _gate_record(gate30),
            "unrestricted_compact_coxeter_no_go": True,
        },
        "record32": record32,
        "record35_chained_digest": record35["certificate_sha256"],
        "record42": {
            "vertex_count": len(vertices42),
            "internal_S4_x_S5_symmetries_verified": True,
            "normalizing_support_symmetry_orders": [24, 120],
            "gate": _gate_record(gate42),
            "exact_index_two_witness": inertia42,
            "unrestricted_compact_coxeter_no_go": True,
        },
        "multiwedge_records": multi_records,
        "A5_exceptional_exact_inertias": {
            "A5_A5_type_1_1": inertia_a,
            "A5_A5_type_4_1": inertia_b,
        },
        "monotone_two_edge_reduction": (
            "Any locally elliptic labeling with at least two cross edges "
            "would remain elliptic after weakening all other cross edges to "
            "m=2 and the chosen two to m=3. Hence only pairs of individually "
            "passing m=3 edges need direct testing."
        ),
        "principal_block_scope_reaudit": (
            "Outside entries cannot alter a decisive principal block. For each "
            "multi-support nerve, fix one audited residual support: its entire "
            "zero/one/at-least-two cross-edge split is covered, and all entries "
            "with an endpoint outside that block are irrelevant to its inertia."
        ),
        "selected_nerve_indices": [30, 32, 35, 42, 43, 44, 47],
        "reaudited_support_indices": [30, 42, 43, 44, 47],
        "graph_shadow_or_uniqueness_claimed": False,
        "all_selected_nerves_compact_coxeter_excluded": True,
        "principal_block_scope_reaudit_passed": True,
        "runtime_seconds": round(perf_counter() - started, 6),
    }
    return result


def compare_archive(result: Mapping[str, object]) -> None:
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    if archive["theorem"]["selected_nerve_indices"] != result["selected_nerve_indices"]:
        raise AssertionError("aggregate selected-nerve scope changed")
    if "graph shadows" not in archive["theorem"]["scope"]:
        raise AssertionError("aggregate archive lost selected-nerve-only scope")
    if archive["scope_reaudit"]["reaudited_indices"] != result[
        "reaudited_support_indices"
    ]:
        raise AssertionError("aggregate scope re-audit changed")
    census = {
        record["record_index"]: record
        for record in archive["six_nerve_census"]
    }
    for record_index in (30, 42):
        independent = result[f"record{record_index}"]
        archived = census[record_index]
        gate = independent["gate"]
        archived_gate = archived["support_gates"][0]
        if (
            archived["rank_nine_vertex_count"] != independent["vertex_count"]
            or archived_gate["single_edge_test_count"] != gate["single_test_count"]
            or archived_gate["single_edge_passing_count"]
            != gate["single_passing_count"]
            or archived["raw_two_edge_trial_count"] != gate["raw_double_trials"]
            or archived["monotonicity_residual_two_edge_trial_count"]
            != gate["residual_double_trials"]
            or archived["two_edge_passing_count"] != gate["double_survivors"]
        ):
            raise AssertionError(f"aggregate record-{record_index} differs from hostile replay")
    for independent in result["multiwedge_records"]:
        archived = census[independent["record_index"]]
        if (
            archived["rank_nine_vertex_count"] != independent["vertex_count"]
            or archived["exceptional_class_histogram"]
            != independent["exceptional_class_histogram"]
            or archived["two_edge_passing_count"]
            != sum(gate["double_survivors"] for gate in independent["gates"])
        ):
            raise AssertionError("aggregate multiwedge census differs from hostile replay")
    detailed32 = archive["record32_detailed_census"]
    independent32 = result["record32"]
    if (
        detailed32["overlapping_order_four_gate"]["compatible_join_count"]
        != independent32["A4_compatible_join_count"]
        or detailed32["overlapping_order_four_gate"]["passing_state_count"]
        != independent32["A4_passing_count"]
        or detailed32["cross_support_gate"]["passing_count"]
        != independent32["gate"]["single_passing_count"]
    ):
        raise AssertionError("aggregate record-32 detail differs from hostile replay")
    witnesses42 = census[42]["exact_index_two_witnesses"]
    if witnesses42[0]["signs"] != result["record42"]["exact_index_two_witness"]["signs"]:
        raise AssertionError("aggregate record-42 exact signs differ")
    archived_a5 = {
        witness["name"]: witness
        for witness in census[44]["exact_index_two_witnesses"]
    }
    for name, witness in result["A5_exceptional_exact_inertias"].items():
        if (
            archived_a5[name]["signs"] != witness["signs"]
            or archived_a5[name]["inertia"] != witness["inertia"]
        ):
            raise AssertionError("aggregate A5 exact inertia differs")


if __name__ == "__main__":
    replay_result = replay()
    if ARCHIVE.exists():
        compare_archive(replay_result)
        replay_result["frozen_archive_compared"] = True
    else:
        replay_result["frozen_archive_compared"] = False
    print(json.dumps(replay_result, indent=2, sort_keys=True))
