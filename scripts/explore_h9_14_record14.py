#!/usr/bin/env python3
"""Exact exploratory census for the exceptional branch of H9/14 record 14."""

from __future__ import annotations

from collections import Counter
from itertools import combinations
from pathlib import Path
import subprocess

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph


ROOT = Path(__file__).resolve().parents[1]
GENG = Path("/opt/homebrew/bin/geng")
VERTEX_COUNT = 14
SPECIALS = (12, 13)
GRAPH_EDGES = ((0, 12), (1, 12), (0, 13), (1, 13), (2, 13))


def mask(vertices: tuple[int, ...] | list[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def candidate_facets() -> tuple[int, ...]:
    edge_masks = tuple(mask(list(edge)) for edge in GRAPH_EDGES)
    return tuple(
        mask(list(chosen))
        for chosen in combinations(range(VERTEX_COUNT), 5)
        if all(mask(list(chosen)) & edge for edge in edge_masks)
    )


FACETS = candidate_facets()
FACET_INDEX = {facet: index for index, facet in enumerate(FACETS)}
ALL_SELECTED_BITS = (1 << len(FACETS)) - 1
DISJOINT_BITS = []
for chosen in range(1 << VERTEX_COUNT):
    value = 0
    for index, facet in enumerate(FACETS):
        if not chosen & facet:
            value |= 1 << index
    DISJOINT_BITS.append(value)

MASKS_BY_SIZE = {
    size: tuple(mask(list(chosen)) for chosen in combinations(range(14), size))
    for size in range(2, 15)
}
GRAPH_EDGE_MASKS = tuple(mask(list(edge)) for edge in GRAPH_EDGES)
TARGET_F_VECTORS = {
    0: (14, 86, 310, 735, 1204, 1372, 1052, 495, 110),
    1: (14, 86, 310, 734, 1198, 1358, 1036, 486, 108),
    2: (14, 86, 310, 733, 1192, 1344, 1020, 477, 106),
    3: (14, 86, 310, 732, 1186, 1330, 1004, 468, 104),
    4: (14, 86, 310, 731, 1180, 1316, 988, 459, 102),
}
LOCAL_LINK = analyze_local_link_graph("M???????????o?w??")


def hits(selected: int, chosen: int) -> bool:
    return selected & DISJOINT_BITS[chosen] == 0


def has_graph_edge(chosen: int) -> bool:
    return any(chosen & edge == edge for edge in GRAPH_EDGE_MASKS)


def is_induced_2k2_free(adjacency: tuple[int, ...]) -> bool:
    for chosen in combinations(range(len(adjacency)), 4):
        chosen_mask = mask(list(chosen))
        if sorted(
            (adjacency[vertex] & chosen_mask).bit_count() for vertex in chosen
        ) == [1, 1, 1, 1]:
            return False
    return True


def y_count(
    adjacency: tuple[int, ...],
    h: int,
    side: int,
) -> int:
    edge_count = sum(value.bit_count() for value in adjacency) // 2
    odd_triangles = sum(
        (
            ((adjacency[first] >> second) & 1)
            + ((adjacency[first] >> third) & 1)
            + ((adjacency[second] >> third) & 1)
        )
        % 2
        for first, second, third in combinations(range(11), 3)
    )
    other_side = ((1 << 11) - 1) ^ side ^ (1 << h)
    toggled_ones = sum(
        (((adjacency[h] >> first) & 1) ^ ((adjacency[first] >> second) & 1))
        for first in range(11)
        if side & (1 << first)
        for second in range(11)
        if other_side & (1 << second)
    )
    return edge_count + odd_triangles + 25 - 2 * toggled_ones


def twisted_y_triples(
    adjacency: tuple[int, ...],
    h: int,
    side: int,
) -> tuple[tuple[int, int, int], ...]:
    other_side = ((1 << 11) - 1) ^ side ^ (1 << h)
    result = []
    # Vertex 11 denotes the external base vertex (actual vertex 0).
    for triple in combinations(range(12), 3):
        if 11 in triple:
            pair = [vertex for vertex in triple if vertex != 11]
            value = (adjacency[pair[0]] >> pair[1]) & 1
        else:
            value = sum(
                (adjacency[first] >> second) & 1
                for first, second in combinations(triple, 2)
            ) % 2
        if h in triple:
            remainder = [vertex for vertex in triple if vertex != h]
            if (
                len(remainder) == 2
                and bool(side & (1 << remainder[0]))
                != bool(side & (1 << remainder[1]))
                and 11 not in remainder
            ):
                value ^= 1
        if value:
            result.append(triple)
    return tuple(result)


FIXED_SELECTED_MASKS = tuple(
    [
        mask([0, 1, 13, first, second])
        for first in range(3, 8)
        for second in range(8, 12)
    ]
    + [
        mask([0, 1, 2, first, second])
        for first in range(3, 8)
        for second in (12, *range(8, 12))
    ]
)
FIXED_SELECTED_BITS = sum(1 << FACET_INDEX[value] for value in FIXED_SELECTED_MASKS)


def selected_facets(
    y_triples: tuple[tuple[int, int, int], ...],
    h: int,
    side: int,
    vertex_two: int,
) -> int:
    other_side = ((1 << 11) - 1) ^ side ^ (1 << h)
    side_vertices = [v for v in range(11) if side & (1 << v)]
    residual_vertices = [
        v for v in range(11)
        if other_side & (1 << v) and v != vertex_two
    ]
    abstract_to_actual = {h: 1, vertex_two: 2}
    abstract_to_actual.update(
        {vertex: 3 + index for index, vertex in enumerate(side_vertices)}
    )
    abstract_to_actual.update(
        {vertex: 8 + index for index, vertex in enumerate(residual_vertices)}
    )
    abstract_to_actual[11] = 0
    result = FIXED_SELECTED_BITS
    for triple in y_triples:
        facet = mask([12, 13, *(abstract_to_actual[value] for value in triple)])
        result |= 1 << FACET_INDEX[facet]
    if result.bit_count() != len(y_triples) + len(FIXED_SELECTED_MASKS):
        raise AssertionError("record-14 facet types unexpectedly overlap")
    return result


def small_blocker_filter(selected: int, expected_a4: int) -> tuple[bool, tuple[int, ...]]:
    observed_edges = tuple(chosen for chosen in MASKS_BY_SIZE[2] if hits(selected, chosen))
    if set(observed_edges) != set(GRAPH_EDGE_MASKS):
        return False, ()
    for chosen in MASKS_BY_SIZE[3]:
        if not has_graph_edge(chosen) and hits(selected, chosen):
            return False, ()
    a4 = tuple(
        chosen
        for chosen in MASKS_BY_SIZE[4]
        if not has_graph_edge(chosen) and hits(selected, chosen)
    )
    if len(a4) != expected_a4:
        return False, a4
    # Every A4 is on W and meets {0,1}; the exceptional A5 through 12 has
    # already spent the only remaining special-link g4 unit.
    if any(chosen & mask([12, 13]) or not chosen & mask([0, 1]) for chosen in a4):
        return False, a4
    return True, a4


def complete_blocker(selected: int) -> tuple[int, ...]:
    blockers: list[int] = []
    for size in range(2, 15):
        for chosen in MASKS_BY_SIZE[size]:
            if any(chosen & blocker == blocker for blocker in blockers):
                continue
            if hits(selected, chosen):
                blockers.append(chosen)
    return tuple(blockers)


def a4_local_constraints(a4: tuple[int, ...]) -> bool:
    for figure in LOCAL_LINK.vertex_figures:
        incidence = sum(bool(face & (1 << figure.vertex)) for face in a4)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in LOCAL_LINK.edge_figures:
        edge = mask(list(figure.edge))
        incidence = sum(face & edge == edge for face in a4)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in LOCAL_LINK.triple_figures:
        triple = mask(list(figure.triple))
        incidence = sum(face & triple == triple for face in a4)
        if incidence > figure.missing_four_incidence_max:
            return False
    return True


def strict_vertex_link_g4(a4: tuple[int, ...], a5: tuple[int, ...]) -> bool:
    baseline = (5, 5, 6, *([7] * 9), 2, 1)
    _, edges = decode_graph6("M???????????o?w??")
    edge_set = {tuple(sorted(edge)) for edge in edges}
    for vertex in range(14):
        link_mask = mask(
            [
                other
                for other in range(14)
                if other != vertex
                and tuple(sorted((vertex, other))) not in edge_set
            ]
        )
        inherited_a4 = sum(face & ~link_mask == 0 for face in a4)
        incident_a5 = sum(bool(face & (1 << vertex)) for face in a5)
        if baseline[vertex] - inherited_a4 - incident_a5 < 1:
            return False
    return True


def face_vector(blockers: tuple[int, ...]) -> tuple[int, ...]:
    counts = [0] * 9
    for chosen in range(1, 1 << 14):
        size = chosen.bit_count()
        if size <= 9 and not any(chosen & blocker == blocker for blocker in blockers):
            counts[size - 1] += 1
    return tuple(counts)


def main() -> None:
    counts: Counter[str] = Counter()
    a4_histogram: Counter[tuple[int, int]] = Counter()
    blocker_histogram: Counter[tuple[int, tuple[tuple[int, int], ...]]] = Counter()
    rejection_counts: Counter[str] = Counter()
    survivors: list[tuple[object, ...]] = []
    blocker_cache: dict[int, tuple[int, ...]] = {}
    process = subprocess.Popen(
        [str(GENG), "-q", "-t", "11"],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    for line in process.stdout:
        graph6 = line.strip()
        _, edges = decode_graph6(graph6)
        adjacency_list = [0] * 11
        for first, second in edges:
            adjacency_list[first] |= 1 << second
            adjacency_list[second] |= 1 << first
        adjacency = tuple(adjacency_list)
        counts["triangle_free_graphs"] += 1
        if not is_induced_2k2_free(adjacency):
            continue
        counts["admissible_graphs"] += 1
        all_vertices = (1 << 11) - 1
        for h in range(11):
            degree = adjacency[h].bit_count()
            if degree > 5:
                continue
            available = [
                vertex
                for vertex in range(11)
                if vertex != h and not adjacency[h] & (1 << vertex)
            ]
            for extra in combinations(available, 5 - degree):
                side = adjacency[h] | mask(list(extra))
                current_y_count = y_count(adjacency, h, side)
                if current_y_count not in (57, 59, 61, 63, 65):
                    continue
                a = (65 - current_y_count) // 2
                counts[f"target_roles_degree_{degree}_a_{a}"] += 1
                y_triples = twisted_y_triples(adjacency, h, side)
                if len(y_triples) != current_y_count:
                    raise AssertionError("twisted y count formula changed")
                other_side = all_vertices ^ side ^ (1 << h)
                for vertex_two in range(11):
                    if not other_side & (1 << vertex_two):
                        continue
                    counts["colored_trials"] += 1
                    selected = selected_facets(y_triples, h, side, vertex_two)
                    if selected.bit_count() != 110 - 2 * a:
                        raise AssertionError("facet count does not match schedule")
                    passed, a4 = small_blocker_filter(selected, a)
                    a4_histogram[(a, len(a4))] += 1
                    if not passed:
                        continue
                    counts["small_blocker_survivors"] += 1
                    blockers = blocker_cache.get(selected)
                    if blockers is None:
                        blockers = complete_blocker(selected)
                        blocker_cache[selected] = blockers
                    histogram = Counter(value.bit_count() for value in blockers)
                    blocker_histogram[(a, tuple(sorted(histogram.items())))] += 1
                    if any(size > 5 for size in histogram):
                        rejection_counts["lanner_order"] += 1
                        continue
                    a4_blockers = tuple(
                        value for value in blockers if value.bit_count() == 4
                    )
                    a5_blockers = tuple(
                        value for value in blockers if value.bit_count() == 5
                    )
                    if not a4_local_constraints(a4_blockers):
                        rejection_counts["local_a4"] += 1
                        continue
                    if not strict_vertex_link_g4(a4_blockers, a5_blockers):
                        rejection_counts["local_g4"] += 1
                        continue
                    if face_vector(blockers) != TARGET_F_VECTORS[a]:
                        rejection_counts["face_vector"] += 1
                        continue
                    counts["full_blocker_survivors"] += 1
                    survivors.append(
                        (
                            graph6,
                            h,
                            tuple(v for v in range(11) if side & (1 << v)),
                            vertex_two,
                            degree,
                            a,
                            current_y_count,
                            tuple(sorted(histogram.items())),
                            tuple(a4),
                            tuple(blockers),
                            selected,
                        )
                    )
    if process.wait() != 0:
        raise RuntimeError("geng failed")
    print("counts", dict(sorted(counts.items())))
    print("a4_histogram", dict(sorted(a4_histogram.items())))
    print("rejection_counts", dict(sorted(rejection_counts.items())))
    print("blocker_histogram", dict(sorted(blocker_histogram.items())))
    print("survivor_count", len(survivors))
    survivor_multiplicities = Counter(survivor[-1] for survivor in survivors)
    print("unique_survivor_count", len(survivor_multiplicities))
    print(
        "survivor_multiplicity_histogram",
        dict(sorted(Counter(survivor_multiplicities.values()).items())),
    )
    print(
        "survivor_branch_histogram",
        dict(sorted(Counter((survivor[4], survivor[5]) for survivor in survivors).items())),
    )
    for survivor in survivors[:20]:
        print("survivor", survivor[:-1])


if __name__ == "__main__":
    main()
