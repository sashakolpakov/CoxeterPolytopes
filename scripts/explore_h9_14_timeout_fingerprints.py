"""Exact structural fingerprints for the unresolved H9/14 graph shadows.

This is a theorem-search utility, not a certificate.  It reads only frozen
archives and uses exact integer graph/clutter calculations to rank the
remaining bounded-timeout cases for solver-independent blocker attacks.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import networkx as nx

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    independent_masks,
    minimal_graph_vertex_covers,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import (
    decode_graph6,
    macaulay_raise,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BASELINE = ROOT / "data/enumerations/h9_14_facet_completion.json"
REALIZED_OR_PROMOTED_TIMEOUTS = {30, 42, 92, 93}


def vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask & (1 << vertex))


def active_graph_data(graph6: str):
    _, edges = decode_graph6(graph6)
    active = tuple(sorted({vertex for edge in edges for vertex in edge}))
    isolates = tuple(vertex for vertex in range(14) if vertex not in active)
    graph = nx.Graph()
    graph.add_nodes_from(active)
    graph.add_edges_from(edges)
    automorphisms = tuple(
        tuple(mapping[vertex] for vertex in active)
        for mapping in nx.algorithms.isomorphism.GraphMatcher(
            graph, graph
        ).isomorphisms_iter()
    )
    return active, isolates, automorphisms


def active_image(mask: int, active: tuple[int, ...], image: tuple[int, ...]) -> int:
    return sum(
        1 << image[index]
        for index, vertex in enumerate(active)
        if mask & (1 << vertex)
    )


def a4_orbit_count(
    faces: tuple[int, ...],
    active: tuple[int, ...],
    isolates: tuple[int, ...],
    automorphisms: tuple[tuple[int, ...], ...],
) -> int:
    isolate_mask = sum(1 << vertex for vertex in isolates)
    keys = {
        (
            min(active_image(face, active, image) for image in automorphisms),
            (face & isolate_mask).bit_count(),
        )
        for face in faces
    }
    return len(keys)


def base_link_g3_g4(graph6: str) -> tuple[tuple[int, ...], tuple[int, ...]]:
    independent = {size: independent_masks(graph6, size) for size in range(2, 6)}
    g3_values = []
    g4_values = []
    for vertex in range(14):
        counts = {
            size: sum(bool(face & (1 << vertex)) for face in independent[size])
            for size in independent
        }
        f0, f1, f2, f3 = (counts[size] for size in range(2, 6))
        g3_values.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4_values.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3_values), tuple(g4_values)


def a4_interval_feasible(face: int, local) -> bool:
    for figure in local.vertex_figures:
        incidence = int(bool(face & (1 << figure.vertex)))
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.edge_figures:
        subset = sum(1 << vertex for vertex in figure.edge)
        incidence = int(face & subset == subset)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.triple_figures:
        subset = sum(1 << vertex for vertex in figure.triple)
        incidence = int(face & subset == subset)
        if incidence > figure.missing_four_incidence_max:
            return False
    return True


def a1_shadow_gate(
    face: int,
    fives: tuple[int, ...],
    base_g3: tuple[int, ...],
    base_g4: tuple[int, ...],
    total_missing_fives: int,
) -> tuple[bool, bool, int, int]:
    shadow = tuple(extension for extension in fives if extension & face == face)
    additional = total_missing_fives - len(shadow)
    if additional < 0:
        return False, False, len(shadow), additional
    lower_sum = 0
    upper_sum = 0
    strict = True
    for vertex in range(14):
        incidence = int(bool(face & (1 << vertex)))
        previous = base_g3[vertex] - incidence
        shadow_incidence = sum(bool(item & (1 << vertex)) for item in shadow)
        before = base_g4[vertex] + 6 * incidence - shadow_incidence
        if previous < 1 or before < 1:
            strict = False
        if previous >= 0:
            lower_sum += max(0, before - macaulay_raise(previous, 3))
        upper_sum += before - 1
    aggregate = strict and lower_sum <= 5 * additional <= upper_sum
    return strict, aggregate, len(shadow), additional


def component_signature(graph6: str) -> str:
    _, edges = decode_graph6(graph6)
    graph = nx.Graph()
    graph.add_nodes_from(range(14))
    graph.add_edges_from(edges)
    pieces = []
    for component in nx.connected_components(graph):
        subgraph = graph.subgraph(component)
        pieces.append(
            (len(component), subgraph.number_of_edges(), tuple(sorted(dict(subgraph.degree()).values())))
        )
    return "+".join(
        f"v{vertices_}e{edges_}d{','.join(map(str, degrees))}"
        for vertices_, edges_, degrees in sorted(pieces, reverse=True)
    )


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    timeout_indices = [
        record["index"]
        for record in baseline["records"]
        if record["solver"]["status"] == "unknown"
        and record["index"] not in REALIZED_OR_PROMOTED_TIMEOUTS
    ]
    output = []
    for index in timeout_indices:
        record = survivors[index]
        graph6 = record["graph6"]
        fours = independent_masks(graph6, 4)
        fives = independent_masks(graph6, 5)
        covers = minimal_graph_vertex_covers(graph6)
        schedules = []
        for schedule in record["dehn_sommerville_schedule"]:
            total = summed_vertex_link_g4(schedule["f_vector"])
            schedules.append(
                {
                    "a": schedule["missing_four_count"],
                    "missing_fives": schedule["missing_five_count_offset"],
                    "facets": schedule["f_vector"][8],
                    "summed_link_g4": total,
                    "strict_sum_feasible": total >= 14,
                }
            )
        local = analyze_local_link_graph(graph6)
        active, isolates, automorphisms = active_graph_data(graph6)
        interval_fours = tuple(face for face in fours if a4_interval_feasible(face, local))
        base_g3, base_g4 = base_link_g3_g4(graph6)
        a1_schedule = next(
            (item for item in schedules if item["a"] == 1), None
        )
        shadow_strict_fours: tuple[int, ...] = ()
        shadow_aggregate_fours: tuple[int, ...] = ()
        shadow_size_histogram: Counter[int] = Counter()
        additional_a5_histogram: Counter[int] = Counter()
        if a1_schedule is not None:
            strict_values = []
            aggregate_values = []
            for face in interval_fours:
                strict, aggregate, shadow_size, additional = a1_shadow_gate(
                    face,
                    fives,
                    base_g3,
                    base_g4,
                    a1_schedule["missing_fives"],
                )
                shadow_size_histogram[shadow_size] += 1
                additional_a5_histogram[additional] += 1
                if strict:
                    strict_values.append(face)
                if aggregate:
                    aggregate_values.append(face)
            shadow_strict_fours = tuple(strict_values)
            shadow_aggregate_fours = tuple(aggregate_values)
        cover_histogram = Counter(mask.bit_count() for mask in covers)
        pattern_counts = Counter(
            sum(bool(face & cover) << position for position, cover in enumerate(covers))
            for face in fives
        )
        output.append(
            {
                "index": index,
                "graph6": graph6,
                "edges": record["edges"],
                "component_signature": component_signature(graph6),
                "ridge_core": record["ridge_support_core_count"],
                "independent_four_count": len(fours),
                "independent_five_count": len(fives),
                "single_a4_orbits": a4_orbit_count(
                    fours, active, isolates, automorphisms
                ),
                "a1_interval_feasible_a4_count": len(interval_fours),
                "a1_interval_feasible_a4_orbits": a4_orbit_count(
                    interval_fours, active, isolates, automorphisms
                ),
                "a1_shadow_strict_a4_count": len(shadow_strict_fours),
                "a1_shadow_strict_a4_orbits": a4_orbit_count(
                    shadow_strict_fours, active, isolates, automorphisms
                ),
                "a1_vertex_aggregate_a5_feasible_a4_count": (
                    len(shadow_aggregate_fours)
                ),
                "a1_vertex_aggregate_a5_feasible_a4_orbits": a4_orbit_count(
                    shadow_aggregate_fours, active, isolates, automorphisms
                ),
                "a1_shadow_size_histogram": dict(sorted(shadow_size_histogram.items())),
                "a1_additional_minimal_a5_histogram": dict(
                    sorted(additional_a5_histogram.items())
                ),
                "active_graph_automorphism_order": len(automorphisms),
                "minimal_cover_count": len(covers),
                "minimal_cover_size_histogram": dict(sorted(cover_histogram.items())),
                "cover_pattern_alphabet": len(pattern_counts),
                "cover_pattern_largest_cell": max(pattern_counts.values()),
                "vertex_a4_intervals": [
                    list(figure.missing_four_incidence_interval)
                    for figure in local.vertex_figures
                ],
                "edge_a4_interval_count": len(local.edge_figures),
                "triple_a4_capacity_count": len(local.triple_figures),
                "schedules": schedules,
            }
        )
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
