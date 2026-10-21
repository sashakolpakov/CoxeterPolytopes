#!/usr/bin/env python3
"""Grow all locally finite cross-component supersets from the two-edge shell.

This is a discovery program.  Its search is lossless for the cross-component
label model because every connected assignment has a two-edge spanning-tree
restriction and ellipticity is downward closed when non-right entries are
removed.  Internal meeting pairs not already present in a sparse Lanner state
remain right-angled; this script does not yet enumerate those extra labels.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Mapping

import networkx as nx

import scripts.explore_h9_14_record30_a2_orbit10_two_cross as base


Assignment = tuple[tuple[int, int], ...]  # (cross-edge index, label)

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a2_orbit10_whole_nerve_no_go.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "879c927f7bedb2d6fdb59a0eaccab7f936ad688e67444b8e74dddb750f8253b3"
)


def canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return canonical_sha(payload)


def cross_edges() -> tuple[base.Pair, ...]:
    return tuple(sorted({edge for category in base._component_cross_edges() for edge in category}))


def nerve_automorphisms() -> tuple[tuple[int, ...], ...]:
    """Enumerate the full colour-preserving minimal-nonface action."""

    graph = nx.Graph()
    for vertex in range(base.N):
        graph.add_node(("v", vertex), kind="facet")
    for index, face in enumerate(base.MINIMAL_NONFACES):
        graph.add_node(("f", index), kind=f"nonface_{len(face)}")
        for vertex in face:
            graph.add_edge(("v", vertex), ("f", index))
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=lambda left, right: left["kind"] == right["kind"],
    )
    result = {
        tuple(move[("v", vertex)][1] for vertex in range(base.N))
        for move in matcher.isomorphisms_iter()
    }
    if len(result) != 1152:
        raise AssertionError("the orbit-10 nerve automorphism order changed")
    return tuple(sorted(result))


def move_state(state: base.LabelState, move: tuple[int, ...]) -> base.LabelState:
    return tuple(
        sorted((base.pair(move[edge[0]], move[edge[1]]), label) for edge, label in state)
    )


def state_orbits(
    states: tuple[base.LabelState, ...], automorphisms: tuple[tuple[int, ...], ...]
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    state_index = {state: index for index, state in enumerate(states)}
    remaining = set(range(len(states)))
    result = []
    while remaining:
        representative = min(remaining)
        orbit = {
            state_index[move_state(states[representative], move)]
            for move in automorphisms
        }
        if not orbit <= remaining | (set(range(len(states))) - remaining):
            raise AssertionError("a sparse-state orbit left the terminal set")
        result.append((representative, tuple(sorted(orbit))))
        remaining.difference_update(orbit)
    if sorted(index for _representative, orbit in result for index in orbit) != list(range(len(states))):
        raise AssertionError("the sparse-state quotient is not a partition")
    return tuple(result)


def seed_assignments(
    state: base.LabelState,
    vertices: tuple[tuple[int, ...], ...],
    affected: dict[base.Pair, tuple[int, ...]],
    edge_index: dict[base.Pair, int],
) -> set[Assignment]:
    labels = base._base_labels(state)
    valid_choices = individually_valid_choices(state, vertices, affected)

    result = set()
    for first, second in base._two_edge_placements():
        common = tuple(sorted(set(affected[first]) & set(affected[second])))
        for first_label in base.FINITE_LABELS:
            if (first, first_label) not in valid_choices:
                continue
            for second_label in base.FINITE_LABELS:
                if (second, second_label) not in valid_choices:
                    continue
                labels[first] = first_label
                labels[second] = second_label
                if all(
                    base._restriction_finite(vertices[index], labels)
                    for index in common
                ):
                    result.add(
                        tuple(
                            sorted(
                                (
                                    (edge_index[first], first_label),
                                    (edge_index[second], second_label),
                                )
                            )
                        )
                    )
                labels[first] = 2
                labels[second] = 2
    return result


def individually_valid_choices(
    state: base.LabelState,
    vertices: tuple[tuple[int, ...], ...],
    affected: dict[base.Pair, tuple[int, ...]],
) -> frozenset[tuple[base.Pair, int]]:
    labels = base._base_labels(state)
    passing = set()
    for edge in cross_edges():
        for label in base.FINITE_LABELS:
            labels[edge] = label
            if all(
                base._restriction_finite(vertices[index], labels)
                for index in affected[edge]
            ):
                passing.add((edge, label))
            labels[edge] = 2
    return frozenset(passing)


def residual_internal_pairs() -> tuple[base.Pair, ...]:
    fixed = {
        base.pair(left, right)
        for support in base.HIGHER_NONFACES
        for left, right in combinations(support, 2)
    }
    dotted = {base.pair(*edge) for edge in base.DOTTED_PAIRS}
    return tuple(
        edge
        for edge in combinations(base.MAIN_COMPONENT, 2)
        if edge not in fixed and edge not in dotted
    )


def audit_residual_internal_edges(
    states: tuple[base.LabelState, ...],
    vertices: tuple[tuple[int, ...], ...],
) -> dict[str, object]:
    residual = residual_internal_pairs()
    vertex_sets = tuple(frozenset(vertex) for vertex in vertices)
    affected = {
        edge: tuple(
            index for index, vertex in enumerate(vertex_sets) if set(edge) <= vertex
        )
        for edge in residual
    }
    passing = []
    for state_index, state in enumerate(states):
        labels = base._base_labels(state)
        for edge in residual:
            for label in base.FINITE_LABELS:
                labels[edge] = label
                if all(
                    base._restriction_finite(vertices[index], labels)
                    for index in affected[edge]
                ):
                    passing.append([state_index, list(edge), label])
                labels[edge] = 2
    if len(residual) != 19:
        raise AssertionError("the residual internal meeting-pair count changed")
    return {
        "residual_internal_meeting_pairs": [list(edge) for edge in residual],
        "literal_single_addition_trial_count": len(states)
        * len(residual)
        * len(base.FINITE_LABELS),
        "passing_single_addition_count": len(passing),
        "passing_single_additions_sha256": canonical_sha(passing),
        "all_residual_internal_pairs_forced_label_two_by_monotonicity": not passing,
    }


def extend_assignments(
    state: base.LabelState,
    current: set[Assignment],
    vertices: tuple[tuple[int, ...], ...],
    affected: dict[base.Pair, tuple[int, ...]],
    edges: tuple[base.Pair, ...],
) -> set[Assignment]:
    labels = base._base_labels(state)
    passing = set()
    checked: set[Assignment] = set()
    for assignment in current:
        occupied = {index for index, _ in assignment}
        for edge_number, edge in enumerate(edges):
            if edge_number in occupied:
                continue
            for label in base.FINITE_LABELS:
                candidate = tuple(sorted((*assignment, (edge_number, label))))
                if candidate in checked:
                    continue
                checked.add(candidate)
                for old_number, old_label in candidate:
                    labels[edges[old_number]] = old_label
                if all(
                    base._restriction_finite(vertices[index], labels)
                    for index in affected[edge]
                ):
                    passing.add(candidate)
                for old_number, _old_label in candidate:
                    labels[edges[old_number]] = 2
    return passing


def gram_witness(
    state: base.LabelState,
    assignment: Assignment,
    edges: tuple[base.Pair, ...],
) -> dict[str, object] | None:
    labels = base._base_labels(state)
    for edge_number, label in assignment:
        labels[edges[edge_number]] = label
    return base._fixed_minor_witness(labels)


def build() -> dict[str, object]:
    vertices = base._actual_vertices()
    states, stages = base.sparse_states(vertices)
    automorphisms = nerve_automorphisms()
    orbits = state_orbits(states, automorphisms)
    edges = cross_edges()
    edge_index = {edge: index for index, edge in enumerate(edges)}
    vertex_sets = tuple(frozenset(vertex) for vertex in vertices)
    affected = {
        edge: tuple(index for index, vertex in enumerate(vertex_sets) if set(edge) <= vertex)
        for edge in edges
    }

    internal_gate = audit_residual_internal_edges(states, vertices)
    literal_single_choices = tuple(
        individually_valid_choices(state, vertices, affected) for state in states
    )
    distinct_single_choices = set(literal_single_choices)
    if len(distinct_single_choices) != 1:
        raise AssertionError("cross single-edge choices depend on the sparse state")
    common_single_choices = next(iter(distinct_single_choices))
    dotted_endpoints = {vertex for edge in base.DOTTED_PAIRS for vertex in edge}
    if (
        len(common_single_choices) != 20
        or any(set(edge) - dotted_endpoints for edge, _label in common_single_choices)
    ):
        raise AssertionError("the exact K2,2,2 cross-support reduction changed")

    quotient_rows: list[dict[str, object]] = []
    representative_shells: dict[int, dict[int, set[Assignment]]] = {}
    for representative, orbit in orbits:
        state = states[representative]
        current = seed_assignments(state, vertices, affected, edge_index)
        state_rows = []
        representative_shells[representative] = {}
        size = 2
        while current:
            representative_shells[representative][size] = current
            witness_histogram: Counter[tuple[tuple[int, ...], int]] = Counter()
            unwitnessed = 0
            for assignment in current:
                witness = gram_witness(state, assignment, edges)
                if witness is None:
                    unwitnessed += 1
                else:
                    witness_histogram[
                        (
                            tuple(witness["excluded_facets"]),
                            int(witness["scaled_determinant_mod_1009"]),
                        )
                    ] += 1
            witness_payload = [
                [list(excluded), residue, multiplicity]
                for (excluded, residue), multiplicity in sorted(
                    witness_histogram.items()
                )
            ]
            state_rows.append(
                {
                    "edge_count": size,
                    "locally_finite_assignment_count": len(current),
                    "no_mod_1009_fixed_minor_witness_count": unwitnessed,
                    "distinct_nonzero_fixed_minor_witness_count": len(
                        witness_histogram
                    ),
                    "fixed_minor_witness_histogram_sha256": canonical_sha(
                        witness_payload
                    ),
                }
            )
            current = extend_assignments(state, current, vertices, affected, edges)
            size += 1
        if size != 5:
            raise AssertionError("the first empty connected cross shell changed")
        quotient_rows.append(
            {
                "representative_state_index": representative,
                "orbit_size": len(orbit),
                "orbit_state_indices_sha256": canonical_sha(orbit),
                "shells": state_rows,
                "first_empty_shell_edge_count": size,
            }
        )

    first_representative = orbits[0][0]
    common_shells = representative_shells[first_representative]
    if any(shells != common_shells for shells in representative_shells.values()):
        raise AssertionError("connected cross shells depend on sparse-state orbit")
    if {size: len(shell) for size, shell in common_shells.items()} != {
        2: 56,
        3: 128,
        4: 38,
    }:
        raise AssertionError("the common connected cross-shell counts changed")

    # Transport every representative assignment through the complete nerve
    # group.  This retains a literal census, not merely orbit-weighted totals.
    transported: dict[int, set[tuple[int, Assignment]]] = {}
    state_index = {state: index for index, state in enumerate(states)}
    edge_moves = {
        move: tuple(
            edge_index[base.pair(move[edge[0]], move[edge[1]])] for edge in edges
        )
        for move in automorphisms
    }
    for representative, _orbit in orbits:
        for size, assignments in representative_shells[representative].items():
            target = transported.setdefault(size, set())
            for move in automorphisms:
                moved_state_index = state_index[move_state(states[representative], move)]
                edge_move = edge_moves[move]
                for assignment in assignments:
                    moved_assignment = tuple(
                        sorted((edge_move[number], label) for number, label in assignment)
                    )
                    target.add((moved_state_index, moved_assignment))

    expected_state_indices = {
        index for _representative, orbit in orbits for index in orbit
    }
    if any(
        state_index_value not in expected_state_indices
        for shell in transported.values()
        for state_index_value, _assignment in shell
    ):
        raise AssertionError("transport escaped the selected state orbits")

    shell_payloads = {
        size: [
            [
                state_index_value,
                [[[edges[number][0], edges[number][1]], label] for number, label in assignment],
            ]
            for state_index_value, assignment in sorted(shell)
        ]
        for size, shell in transported.items()
    }
    category_of = {
        edge_index[edge]: category_index
        for category_index, category in enumerate(base._component_cross_edges())
        for edge in category
    }
    shell_two_legacy_payload = []
    for state_index_value, assignment in transported.get(2, set()):
        ordered = tuple(sorted(assignment, key=lambda item: category_of[item[0]]))
        shell_two_legacy_payload.append(
            [
                state_index_value,
                list(edges[ordered[0][0]]),
                ordered[0][1],
                list(edges[ordered[1][0]]),
                ordered[1][1],
            ]
        )
    shell_two_legacy_payload.sort()

    if canonical_sha(shell_two_legacy_payload) != (
        "151a6305767cb86e227efb0b1be7559c519c1f1841897d181e42599c586a6806"
    ):
        raise AssertionError("the frozen two-cross survivor set was not reproduced")
    if any(
        len(transported[size]) != len(states) * len(common_shells[size])
        for size in common_shells
    ):
        raise AssertionError("literal symmetry transport lost a cross assignment")
    if any(
        row["no_mod_1009_fixed_minor_witness_count"]
        for quotient in quotient_rows
        for row in quotient["shells"]
    ):
        raise AssertionError("a quotient cross assignment escaped the fixed-minor gate")

    common_choice_payload = [
        [list(edge), label] for edge, label in sorted(common_single_choices)
    ]
    compact_quotient_rows = [
        {
            "representative_state_index": row["representative_state_index"],
            "orbit_size": row["orbit_size"],
            "first_empty_shell_edge_count": row[
                "first_empty_shell_edge_count"
            ],
            "shell_counts": [
                shell["locally_finite_assignment_count"]
                for shell in row["shells"]
            ],
            "unwitnessed_counts": [
                shell["no_mod_1009_fixed_minor_witness_count"]
                for shell in row["shells"]
            ],
            "full_row_sha256": canonical_sha(row),
        }
        for row in quotient_rows
    ]
    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record30_a2_orbit10_whole_nerve_compact_no_go",
        "claim_status": "frozen_exact_theorem_with_independent_replay",
        "claim": (
            "The Record 30 a=2 orbit-10 nerve admits no compact Coxeter "
            "realization in H^9."
        ),
        "candidate": {
            "record_index": 30,
            "addition_schedule": 2,
            "positive_A4_orbit_index": 10,
            "minimal_nonfaces": [list(face) for face in base.MINIMAL_NONFACES],
            "actual_rank_nine_vertex_count": len(vertices),
            "actual_vertices_sha256": canonical_sha(vertices),
        },
        "sparse_lanner_gate": {
            "join_stages": list(stages),
            "terminal_state_count": len(states),
            "terminal_states_sha256": canonical_sha(states),
            "complete_order_four_and_five_lanner_placements": True,
            "all_actual_vertex_restrictions_tested": True,
        },
        "residual_internal_label_gate": internal_gate,
        "cross_component_gate": {
            "component_vertex_sets": [list(component) for component in base.COMPONENTS],
            "all_cross_meeting_pair_count": len(edges),
            "literal_single_edge_trial_count": len(states)
            * len(edges)
            * len(base.FINITE_LABELS),
            "literal_single_edge_passing_count": sum(map(len, literal_single_choices)),
            "single_edge_passing_choices_are_identical_for_all_sparse_states": True,
            "common_single_edge_passing_choices": common_choice_payload,
            "common_single_edge_passing_choices_sha256": canonical_sha(common_choice_payload),
            "support_reduces_to_K2_2_2_on_dotted_pair_endpoints": True,
            "no_node_3_through_10_in_a_passing_cross_edge": True,
            "connected_component_graph_requires_at_least_two_cross_edges": True,
            "connected_shell_counts_per_sparse_state": {
                str(size): len(shell) for size, shell in sorted(common_shells.items())
            },
            "first_empty_connected_shell_edge_count": 5,
            "no_connected_assignment_with_five_or_more_edges_by_downward_closure": True,
            "literal_connected_shell_counts": {
                str(size): len(shell) for size, shell in sorted(transported.items())
            },
            "literal_connected_shell_sha256": {
                str(size): canonical_sha(payload)
                for size, payload in sorted(shell_payloads.items())
            },
            "frozen_two_edge_legacy_survivor_sha256": canonical_sha(
                shell_two_legacy_payload
            ),
        },
        "symmetry_and_metric_gate": {
            "nerve_automorphism_order": len(automorphisms),
            "sparse_state_orbit_count": len(orbits),
            "sparse_state_orbit_size_histogram": {
                str(size): count
                for size, count in sorted(
                    Counter(len(orbit) for _representative, orbit in orbits).items()
                )
            },
            "orbit_partition_sha256": canonical_sha(orbits),
            "quotient_transcript_sha256": canonical_sha(quotient_rows),
            "quotient_rows": compact_quotient_rows,
            "fixed_minor_normalization": (
                "diagonal 4, m3=-2, m4=-2*sqrt(2), m5=-(1+sqrt(5))"
            ),
            "coefficient_ring": "Z[sqrt(2),sqrt(5)]",
            "modulus": base.MODULUS,
            "sqrt2_residue": base.SQRT2_MOD,
            "sqrt5_residue": base.SQRT5_MOD,
            "all_quotient_shell_assignments_have_rank_at_least_eleven": True,
            "all_literal_shell_assignments_excluded_by_symmetry_transport": True,
            "H9_gram_rank_upper_bound": 10,
        },
        "conclusion": {
            "whole_orbit10_nerve_compact_coxeter_impossible": True,
            "new_hyperbolic_polytope": False,
        },
        "scope": {
            "whole_orbit10_nerve_compact_excluded": True,
            "whole_record30_a2_schedule_excluded": False,
            "whole_record30_shadow_excluded": False,
            "finite_volume_noncompact_statement": False,
            "Burcroff_reduced_finite_label_alphabet": [2, 3, 4, 5],
        },
        "independent_replay": {
            "script": "scripts/independent_h9_14_record30_a2_orbit10_whole_nerve.py",
            "imports_production_enumerator": False,
            "imports_coxeter_label_module": False,
            "opposite_support_join_order": [3, 1, 2, 0],
            "second_modulus": 1129,
            "second_sqrt2_residue": 406,
            "second_sqrt5_residue": 474,
        },
    }
    document["certificate_sha256"] = certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record30_a2_orbit10_whole_nerve_compact_no_go":
        raise ValueError("unexpected Record 30 orbit-10 whole-nerve kind")
    rebuilt = build()
    if document != rebuilt:
        raise ValueError("the Record 30 orbit-10 whole-nerve certificate did not replay")
    digest = certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("the Record 30 orbit-10 whole-nerve digest is invalid")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("the Record 30 orbit-10 whole-nerve certificate is not frozen")
    return {
        "sparse_state_count": rebuilt["sparse_lanner_gate"]["terminal_state_count"],
        "residual_internal_passing_count": rebuilt["residual_internal_label_gate"]["passing_single_addition_count"],
        "literal_connected_shell_counts": rebuilt["cross_component_gate"]["literal_connected_shell_counts"],
        "first_empty_connected_shell": rebuilt["cross_component_gate"]["first_empty_connected_shell_edge_count"],
        "whole_orbit10_nerve_compact_coxeter_impossible": rebuilt["conclusion"]["whole_orbit10_nerve_compact_coxeter_impossible"],
        "certificate_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        print(json.dumps(verify(json.loads(ARCHIVE.read_text(encoding="utf-8"))), indent=2, sort_keys=True))
    else:
        print(json.dumps(build(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
