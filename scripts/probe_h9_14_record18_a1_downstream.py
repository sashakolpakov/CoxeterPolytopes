#!/usr/bin/env python3
"""Symmetry, topology, and compact-label probe for record-18 a=1 targets."""

from __future__ import annotations

import argparse
from collections import Counter
from typing import Sequence

import networkx as nx

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record18 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
)
from scripts.explore_h9_14_record18_a1_ridge_z3 import mask, solve, vertices
from scripts.explore_h9_14_record18_a2_target_labels import (
    minimal_transversals,
    sparse_label_gate,
)


A4 = (0, 2, 4, 5)
EXPECTED_F_VECTOR = (14, 86, 310, 734, 1198, 1358, 1036, 486, 108)


def incidence_graph(additions: Sequence[tuple[int, ...]]) -> nx.Graph:
    graph = nx.Graph()
    for vertex in range(14):
        graph.add_node(("v", vertex), kind="vertex")
    for left, right in TARGET_MISSING_EDGES:
        graph.add_edge(("v", left), ("v", right))
    four_node = ("q", 0)
    graph.add_node(four_node, kind="four")
    for vertex in A4:
        graph.add_edge(four_node, ("v", vertex))
    for index, row in enumerate(additions):
        node = ("h", index)
        graph.add_node(node, kind="five")
        for vertex in row:
            graph.add_edge(node, ("v", vertex))
    return graph


def isomorphic(left: nx.Graph, right: nx.Graph) -> bool:
    return nx.is_isomorphic(
        left,
        right,
        node_match=nx.algorithms.isomorphism.categorical_node_match("kind", None),
    )


def quotient(
    families: Sequence[tuple[tuple[int, ...], ...]],
) -> tuple[tuple[tuple[tuple[int, ...], ...], tuple[int, ...]], ...]:
    representatives: list[tuple[tuple[int, ...], ...]] = []
    graphs: list[nx.Graph] = []
    hashes: list[str] = []
    members: list[list[int]] = []
    for index, family in enumerate(families):
        graph = incidence_graph(family)
        fingerprint = nx.weisfeiler_lehman_graph_hash(graph, node_attr="kind")
        found = None
        for orbit, (other, other_hash) in enumerate(zip(graphs, hashes, strict=True)):
            if fingerprint == other_hash and isomorphic(graph, other):
                found = orbit
                break
        if found is None:
            representatives.append(family)
            graphs.append(graph)
            hashes.append(fingerprint)
            members.append([index])
        else:
            members[found].append(index)
    return tuple(
        (representative, tuple(orbit_members))
        for representative, orbit_members in zip(representatives, members, strict=True)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit", type=int)
    args = parser.parse_args()
    source = solve(mask(A4))
    raw = source["scheduled_target_additions"]
    families = tuple(
        tuple(tuple(int(vertex) for vertex in row) for row in family)
        for family in raw
    )
    if (
        source["ridge_and_scheduled_model_count"] != 1080
        or source["scheduled_pure_target_count"] != 1080
        or len(families) != 1080
    ):
        raise AssertionError("record-18 a1 [53] target count changed")
    orbits = quotient(families)
    print("raw_targets", len(families), flush=True)
    print("full_graph_orbits", len(orbits), flush=True)
    print("orbit_sizes", sorted(len(members) for _, members in orbits), flush=True)
    topology_green = 0
    label_zero = 0
    for index, (additions, members) in enumerate(orbits):
        if args.orbit is not None and index != args.orbit:
            continue
        supports = (A4, *additions)
        blocker = minimal_transversals((*TARGET_MISSING_EDGES, *supports))
        histogram = Counter(value.bit_count() for value in blocker)
        if histogram != {5: 108}:
            raise AssertionError("a ridge-scheduled model lost literal purity")
        topology = audit_global_facet_selection(
            TARGET_GRAPH6,
            blocker,
            [{"missing_four_count": 1, "f_vector": EXPECTED_F_VECTOR}],
            blocker,
        ).to_record()
        green = bool(topology["passes_connected_homology_audit"])
        topology_green += green
        label_gate = sparse_label_gate(supports, 250_000) if green else None
        terminal = label_gate["terminal_sparse_state_count"] if label_gate else -1
        label_zero += bool(label_gate and not label_gate["state_cap_exceeded"] and terminal == 0)
        stages = (
            [stage["passing_state_count"] for stage in label_gate["stages"]]
            if label_gate
            else []
        )
        print(
            f"orbit {index} size {len(members)} topology {int(green)} "
            f"label_stages {stages} terminal {terminal} additions {additions}",
            flush=True,
        )
    print("topology_green_orbits", topology_green)
    print("topology_green_label_zero_orbits", label_zero)


if __name__ == "__main__":
    main()
