#!/usr/bin/env python3
"""Standalone replay of the Record 30 orbit-10 whole-nerve no-go.

The Lanner and finite-diagram engines come from the independent, no-production
two-edge replay.  This program independently audits the residual internal
pairs, reconstructs the full cross-component CSP, quotients the 648 sparse
states by the exact nerve action, and uses a second finite-field Gram test.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from typing import Mapping

import networkx as nx

import scripts.independent_h9_14_record30_a2_orbit10_two_cross as hostile


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a2_orbit10_whole_nerve_no_go.json"
Assignment = tuple[tuple[int, int], ...]


def automorphisms() -> tuple[tuple[int, ...], ...]:
    graph = nx.Graph()
    for vertex in range(hostile.N):
        graph.add_node(("v", vertex), color="facet")
    for index, face in enumerate(hostile.DOTTED + hostile.SUPPORTS):
        graph.add_node(("f", index), color=f"missing_{len(face)}")
        for vertex in face:
            graph.add_edge(("v", vertex), ("f", index))
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=lambda left, right: left["color"] == right["color"],
    )
    result = {
        tuple(move[("v", vertex)][1] for vertex in range(hostile.N))
        for move in matcher.isomorphisms_iter()
    }
    if len(result) != 1152:
        raise AssertionError("independent nerve group order changed")
    return tuple(sorted(result))


def move_state(state: hostile.State, move: tuple[int, ...]) -> hostile.State:
    return tuple(
        sorted(
            (hostile.pair(move[edge[0]], move[edge[1]]), label)
            for edge, label in state
        )
    )


def state_orbits(
    states: tuple[hostile.State, ...], moves: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    lookup = {state: index for index, state in enumerate(states)}
    remaining = set(range(len(states)))
    result = []
    while remaining:
        representative = min(remaining)
        orbit = {
            lookup[move_state(states[representative], move)] for move in moves
        }
        result.append((representative, tuple(sorted(orbit))))
        remaining.difference_update(orbit)
    return tuple(result)


def cross_edges() -> tuple[hostile.Pair, ...]:
    return tuple(
        sorted(
            {
                hostile.pair(left, right)
                for first, second in ((0, 1), (0, 2), (1, 2))
                for left in hostile.COMPONENTS[first]
                for right in hostile.COMPONENTS[second]
            }
        )
    )


def residual_internal_pairs() -> tuple[hostile.Pair, ...]:
    main = hostile.COMPONENTS[0]
    fixed = {
        hostile.pair(left, right)
        for support in hostile.SUPPORTS
        for left, right in combinations(support, 2)
    }
    dotted = set(hostile.DOTTED)
    return tuple(
        edge
        for edge in combinations(main, 2)
        if edge not in fixed and edge not in dotted
    )


def labels_with_assignment(
    state: hostile.State,
    assignment: Assignment,
    edges: tuple[hostile.Pair, ...],
) -> dict[hostile.Pair, int]:
    labels = hostile.base_labels(state)
    for number, label in assignment:
        labels[edges[number]] = label
    return labels


def all_vertices_finite(
    state: hostile.State,
    assignment: Assignment,
    edges: tuple[hostile.Pair, ...],
    vertices: tuple[tuple[int, ...], ...],
) -> bool:
    labels = labels_with_assignment(state, assignment, edges)
    return all(hostile.restriction_finite(vertex, labels) for vertex in vertices)


def single_choices(
    state: hostile.State,
    edges: tuple[hostile.Pair, ...],
    vertices: tuple[tuple[int, ...], ...],
) -> frozenset[tuple[int, int]]:
    result = set()
    for number, _edge in enumerate(edges):
        for label in hostile.LABELS:
            assignment = ((number, label),)
            if all_vertices_finite(state, assignment, edges, vertices):
                result.add((number, label))
    return frozenset(result)


def component_category(edge: hostile.Pair) -> tuple[int, int]:
    memberships = {
        vertex: component
        for component, vertices in enumerate(hostile.COMPONENTS)
        for vertex in vertices
    }
    return tuple(sorted((memberships[edge[0]], memberships[edge[1]])))


def connected_two_seeds(
    state: hostile.State,
    choices: frozenset[tuple[int, int]],
    edges: tuple[hostile.Pair, ...],
    vertices: tuple[tuple[int, ...], ...],
) -> set[Assignment]:
    result = set()
    for first, second in combinations(sorted(choices), 2):
        if first[0] == second[0]:
            continue
        if component_category(edges[first[0]]) == component_category(edges[second[0]]):
            continue
        assignment = tuple(sorted((first, second)))
        if all_vertices_finite(state, assignment, edges, vertices):
            result.add(assignment)
    return result


def extend(
    state: hostile.State,
    current: set[Assignment],
    choices: frozenset[tuple[int, int]],
    edges: tuple[hostile.Pair, ...],
    vertices: tuple[tuple[int, ...], ...],
) -> set[Assignment]:
    result = set()
    tested = set()
    for assignment in current:
        occupied = {number for number, _label in assignment}
        for choice in choices:
            if choice[0] in occupied:
                continue
            candidate = tuple(sorted((*assignment, choice)))
            if candidate in tested:
                continue
            tested.add(candidate)
            if all_vertices_finite(state, candidate, edges, vertices):
                result.add(candidate)
    return result


def replay() -> dict[str, object]:
    vertices = hostile.actual_vertices()
    states = hostile.sparse_states(vertices)
    moves = automorphisms()
    orbits = state_orbits(states, moves)
    edges = cross_edges()

    internal = residual_internal_pairs()
    internal_passing = []
    for state_index, state in enumerate(states):
        labels = hostile.base_labels(state)
        for edge in internal:
            for label in hostile.LABELS:
                labels[edge] = label
                if all(hostile.restriction_finite(vertex, labels) for vertex in vertices):
                    internal_passing.append([state_index, list(edge), label])
                labels[edge] = 2

    choice_sets = tuple(single_choices(state, edges, vertices) for state in states)
    if len(set(choice_sets)) != 1:
        raise AssertionError("independent cross choices vary by sparse state")
    common_choices = choice_sets[0]

    common_shells: dict[int, set[Assignment]] | None = None
    quotient_rows = []
    zero_residue_count = 0
    residue_histogram: Counter[tuple[int, int]] = Counter()
    for representative, orbit in orbits:
        state = states[representative]
        current = connected_two_seeds(
            state, common_choices, edges, vertices
        )
        shells = {}
        size = 2
        row_counts = []
        while current:
            shells[size] = current
            row_counts.append(len(current))
            for assignment in current:
                labels = labels_with_assignment(state, assignment, edges)
                residue = hostile.fixed_residue(labels)
                residue_histogram[(size, residue)] += 1
                zero_residue_count += residue == 0
            current = extend(
                state, current, common_choices, edges, vertices
            )
            size += 1
        if common_shells is None:
            common_shells = shells
        elif shells != common_shells:
            raise AssertionError("independent quotient shell sets disagree")
        quotient_rows.append(
            {
                "representative_state_index": representative,
                "orbit_size": len(orbit),
                "shell_counts": row_counts,
                "first_empty_shell": size,
            }
        )
    if common_shells is None:
        raise AssertionError("independent quotient produced no representative")

    shell_payloads = {
        size: [
            [
                state_index,
                [[[edges[number][0], edges[number][1]], label] for number, label in assignment],
            ]
            for state_index in range(len(states))
            for assignment in sorted(shell)
        ]
        for size, shell in common_shells.items()
    }
    result = {
        "kind": "h9_14_record30_a2_orbit10_whole_nerve_independent_replay",
        "implementation": {
            "imports_production_enumerator": False,
            "imports_coxeter_label_module": False,
            "support_join_order": list(hostile.JOIN_ORDER),
            "finite_recognizer": "standalone ADEBFH forest classifier",
            "gram_modulus": hostile.MODULUS,
            "sqrt2_residue": hostile.SQRT2,
            "sqrt5_residue": hostile.SQRT5,
        },
        "actual_rank_nine_vertex_count": len(vertices),
        "sparse_state_count": len(states),
        "terminal_states_sha256": hostile.canonical_sha(states),
        "nerve_automorphism_order": len(moves),
        "sparse_state_orbit_count": len(orbits),
        "sparse_state_orbit_sizes": [len(orbit) for _representative, orbit in orbits],
        "residual_internal_pair_count": len(internal),
        "residual_internal_trial_count": len(states) * len(internal) * 3,
        "residual_internal_passing_count": len(internal_passing),
        "cross_meeting_pair_count": len(edges),
        "common_single_choice_count": len(common_choices),
        "common_single_choices_sha256": hostile.canonical_sha(
            [[list(edges[number]), label] for number, label in sorted(common_choices)]
        ),
        "connected_shell_counts_per_state": {
            str(size): len(shell) for size, shell in sorted(common_shells.items())
        },
        "literal_connected_shell_counts": {
            str(size): len(payload) for size, payload in sorted(shell_payloads.items())
        },
        "literal_connected_shell_sha256": {
            str(size): hostile.canonical_sha(payload)
            for size, payload in sorted(shell_payloads.items())
        },
        "first_empty_connected_shell": 5,
        "quotient_rows": quotient_rows,
        "second_prime_residue_histogram": {
            f"{size}:{residue}": count
            for (size, residue), count in sorted(residue_histogram.items())
        },
        "zero_second_prime_fixed_minor_count": zero_residue_count,
        "whole_orbit10_nerve_compact_coxeter_impossible": (
            not internal_passing and zero_residue_count == 0
        ),
    }
    if ARCHIVE.exists():
        source = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        result["agrees_with_frozen_production"] = (
            result["actual_rank_nine_vertex_count"]
            == source["candidate"]["actual_rank_nine_vertex_count"]
            and result["sparse_state_count"]
            == source["sparse_lanner_gate"]["terminal_state_count"]
            and result["terminal_states_sha256"]
            == source["sparse_lanner_gate"]["terminal_states_sha256"]
            and result["nerve_automorphism_order"]
            == source["symmetry_and_metric_gate"]["nerve_automorphism_order"]
            and result["residual_internal_trial_count"]
            == source["residual_internal_label_gate"]["literal_single_addition_trial_count"]
            and result["residual_internal_passing_count"] == 0
            and result["common_single_choices_sha256"]
            == source["cross_component_gate"]["common_single_edge_passing_choices_sha256"]
            and result["literal_connected_shell_counts"]
            == source["cross_component_gate"]["literal_connected_shell_counts"]
            and result["literal_connected_shell_sha256"]
            == source["cross_component_gate"]["literal_connected_shell_sha256"]
            and result["zero_second_prime_fixed_minor_count"] == 0
        )
    return result


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
