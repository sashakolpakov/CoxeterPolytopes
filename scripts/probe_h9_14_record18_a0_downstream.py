#!/usr/bin/env python3
"""Topology and full-graph orbit probe for corrected record-18 a=0 targets.

This consumes the literal ``SURVIVOR`` rows emitted by the standalone
balanced-cut C++ scout.  It intentionally remains a scout: theorem promotion
requires a separately implemented replay of both the generator and every
downstream gate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import re
import subprocess
from typing import Iterable, Sequence

import networkx as nx

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record18 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
)
from scripts.explore_h9_14_record18_a2_target_labels import sparse_label_gate


N = 14
FULL = (1 << N) - 1
EXPECTED_F_VECTOR = (14, 86, 310, 735, 1204, 1372, 1052, 495, 110)
SURVIVOR_RE = re.compile(r"^SURVIVOR key=([^ ]+) rows=([0-9,]+)$")


def vertices(value: int) -> tuple[int, ...]:
    return tuple(index for index in range(N) if value >> index & 1)


def minimal_transversals(rows: Sequence[int]) -> tuple[int, ...]:
    answer = []
    for value in range(1, FULL + 1):
        if not all(value & row for row in rows):
            continue
        if any(
            all((value ^ (1 << vertex)) & row for row in rows)
            for vertex in vertices(value)
        ):
            continue
        answer.append(value)
    return tuple(answer)


def incidence_graph(rows: Sequence[int]) -> nx.Graph:
    graph = nx.Graph()
    for vertex in range(N):
        graph.add_node(("v", vertex), kind="vertex")
    for left, right in TARGET_MISSING_EDGES:
        graph.add_edge(("v", left), ("v", right))
    for index, row in enumerate(rows):
        node = ("h", index)
        graph.add_node(node, kind="five")
        for vertex in vertices(row):
            graph.add_edge(node, ("v", vertex))
    return graph


def isomorphic(left: nx.Graph, right: nx.Graph) -> bool:
    return nx.is_isomorphic(
        left,
        right,
        node_match=nx.algorithms.isomorphism.categorical_node_match("kind", None),
    )


def parse(binary: str) -> tuple[tuple[int, ...], ...]:
    completed = subprocess.run(
        [binary], check=True, capture_output=True, text=True
    )
    families = []
    for line in completed.stdout.splitlines():
        match = SURVIVOR_RE.match(line)
        if match is None:
            continue
        rows = tuple(sorted(int(value) for value in match.group(2).split(",") if value))
        if len(rows) != 8 or len(set(rows)) != 8:
            raise AssertionError("malformed a0 survivor family")
        families.append(rows)
    if len(families) != 196 or len(set(families)) != 196:
        raise AssertionError("corrected a0 survivor count changed")
    return tuple(families)


def quotient(
    families: Sequence[tuple[int, ...]],
) -> tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]:
    representatives: list[tuple[int, ...]] = []
    graphs: list[nx.Graph] = []
    members: list[list[int]] = []
    buckets: list[str] = []
    for index, rows in enumerate(families):
        graph = incidence_graph(rows)
        fingerprint = nx.weisfeiler_lehman_graph_hash(graph, node_attr="kind")
        found = None
        for orbit, (other, other_fingerprint) in enumerate(zip(graphs, buckets, strict=True)):
            if fingerprint == other_fingerprint and isomorphic(graph, other):
                found = orbit
                break
        if found is None:
            representatives.append(rows)
            graphs.append(graph)
            buckets.append(fingerprint)
            members.append([index])
        else:
            members[found].append(index)
    return tuple(
        (representative, tuple(orbit_members))
        for representative, orbit_members in zip(representatives, members, strict=True)
    )


def audit(rows: Sequence[int]) -> dict[str, object]:
    edge_masks = tuple(sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES)
    blocker = minimal_transversals((*edge_masks, *rows))
    histogram = Counter(value.bit_count() for value in blocker)
    if histogram != {5: 110}:
        raise AssertionError("literal corrected target lost blocker purity")
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        blocker,
        [{"missing_four_count": 0, "f_vector": EXPECTED_F_VECTOR}],
        blocker,
    ).to_record()
    return {
        "rows": [list(vertices(row)) for row in rows],
        "blocker_size_histogram": dict(sorted(histogram.items())),
        "global_audit": global_audit,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary")
    args = parser.parse_args()
    families = parse(args.binary)
    orbits = quotient(families)
    print(f"raw_survivors {len(families)}")
    print(f"full_graph_orbits {len(orbits)}")
    print("orbit_sizes", sorted(len(members) for _, members in orbits))
    passing = 0
    for index, (representative, members) in enumerate(orbits):
        result = audit(representative)
        topology = result["global_audit"]
        green = bool(topology["passes_connected_homology_audit"])
        passing += green
        ridge = topology["ridge_degree_histogram"]
        label_gate = (
            sparse_label_gate(tuple(vertices(row) for row in representative), 250_000)
            if green
            else None
        )
        label_counts = (
            [stage["passing_state_count"] for stage in label_gate["stages"]]
            if label_gate is not None
            else []
        )
        print(
            f"orbit {index} size {len(members)} topology {int(green)} "
            f"ridge {ridge} label_stages {label_counts} "
            f"label_terminal "
            f"{label_gate['terminal_sparse_state_count'] if label_gate else -1} "
            f"rows {representative}"
        )
    print(f"topology_green_orbits {passing}")


if __name__ == "__main__":
    main()
