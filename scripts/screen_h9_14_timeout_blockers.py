#!/usr/bin/env python3
"""Exact structural screen for the pre-screen H9/14 facet timeouts.

This is deliberately a graph/link arithmetic screen, not a SAT model.  It
identifies from the frozen facet-completion archive the 23 records that were
theorem-undecided before this screen, enumerates every inclusion-minimal
vertex cover of their missing-edge graphs, and retains only
Dehn--Sommerville schedules compatible with strict
positivity of ``g4`` in all fourteen vertex links.

The output is intended to expose shadows for which blocker duality decomposes
the number of five-element minimal transversals into a few low-order hitting
counts, as in the record-23 certificate.
"""

from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    minimal_graph_vertex_covers,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


ROOT = Path(__file__).resolve().parents[1]
FACET_ARCHIVE = ROOT / "data/enumerations/h9_14_facet_completion.json"
SUPPORT_ARCHIVE = ROOT / "data/enumerations/h9_14_nerve_support.json"

# These four baseline timeouts have since been decided independently or by an
# explicit audited selection.  Record 23 and record 54 were baseline UNSAT,
# so they do not occur in the timeout list in the first place.
RESOLVED_TIMEOUTS = {30, 42, 92, 93}


def vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask & (1 << vertex))


def active_components(graph6: str) -> tuple[tuple[int, ...], ...]:
    order, edges = decode_graph6(graph6)
    adjacency = [set() for _ in range(order)]
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    unseen = {vertex for vertex in range(order) if adjacency[vertex]}
    components: list[tuple[int, ...]] = []
    while unseen:
        seed = min(unseen)
        stack = [seed]
        component: set[int] = set()
        while stack:
            vertex = stack.pop()
            if vertex in component:
                continue
            component.add(vertex)
            stack.extend(adjacency[vertex] - component)
        unseen.difference_update(component)
        components.append(tuple(sorted(component)))
    return tuple(components)


def component_signature(graph6: str) -> str:
    order, edges = decode_graph6(graph6)
    degree = Counter()
    for first, second in edges:
        degree[first] += 1
        degree[second] += 1
    pieces = []
    for component in active_components(graph6):
        edge_count = sum(
            first in component and second in component for first, second in edges
        )
        degrees = sorted((degree[vertex] for vertex in component), reverse=True)
        pieces.append(f"{len(component)}v/{edge_count}e/{','.join(map(str, degrees))}")
    isolate_count = order - sum(map(len, active_components(graph6)))
    return "+".join(pieces) + f"+{isolate_count}I"


def inclusion_exclusion_profile(covers: tuple[int, ...]) -> Counter[tuple[int, int]]:
    """Count blocker terms by coefficient sign and supplement order.

    A nonempty subfamily J contributes ``(-1)^(|J|+1) tau_k`` with
    ``k=5-|union(J)|``.  Unions larger than five cannot lie in a five-set and
    contribute zero.
    """

    profile: Counter[tuple[int, int]] = Counter()
    for size in range(1, len(covers) + 1):
        sign = 1 if size % 2 else -1
        for chosen in combinations(covers, size):
            union = 0
            for cover in chosen:
                union |= cover
            if union.bit_count() <= 5:
                profile[(sign, 5 - union.bit_count())] += 1
    return profile


def main() -> None:
    completion = json.loads(FACET_ARCHIVE.read_text())
    support = json.loads(SUPPORT_ARCHIVE.read_text())
    support_by_graph = {
        record["graph6"]: record
        for record in support["records"]
        if record.get("survives_second_layer") is True
    }
    timeouts = [
        record
        for record in completion["records"]
        if record["solver"]["status"] == "unknown"
        and record["index"] not in RESOLVED_TIMEOUTS
    ]
    if len(timeouts) != 23:
        raise ValueError(f"expected 23 theorem-undecided timeouts, got {len(timeouts)}")

    print(
        "record\tedges\tcomponent_signature\tminimal_covers(size:count)"
        "\tcover_count\tIE_profile\tstrict-g4 schedules a:b:f8:sum"
    )
    for record in timeouts:
        graph6 = record["graph6"]
        support_record = support_by_graph[graph6]
        covers = minimal_graph_vertex_covers(graph6)
        cover_histogram = Counter(mask.bit_count() for mask in covers)
        schedules = []
        for schedule in support_record["dehn_sommerville_schedule"]:
            summed_g4 = summed_vertex_link_g4(schedule["f_vector"])
            if summed_g4 >= 14:
                schedules.append(
                    (
                        schedule["missing_four_count"],
                        schedule["missing_five_count_offset"],
                        schedule["f_vector"][8],
                        summed_g4,
                    )
                )
        histogram = ",".join(
            f"{size}:{cover_histogram[size]}" for size in sorted(cover_histogram)
        )
        profile = inclusion_exclusion_profile(covers)
        profile_text = ",".join(
            f"{'+' if sign > 0 else '-'}{profile[(sign, degree)]}t{degree}"
            for sign, degree in sorted(profile, key=lambda item: (-item[1], -item[0]))
        )
        schedule_text = ",".join(":".join(map(str, row)) for row in schedules)
        print(
            f"{record['index']}\t{record['edges']}\t{component_signature(graph6)}"
            f"\t{histogram}\t{len(covers)}\t{profile_text}\t{schedule_text}"
        )
        print(
            "  covers="
            + " ".join("{" + ",".join(map(str, vertices(mask))) + "}" for mask in covers)
        )


if __name__ == "__main__":
    main()
