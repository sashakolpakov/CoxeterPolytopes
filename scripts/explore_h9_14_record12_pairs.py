#!/usr/bin/env python3
"""Exploratory paired-branch census for H9/14 graph-shadow record 12."""

from __future__ import annotations

from collections import Counter, deque
from itertools import combinations, permutations
import sys

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


N = 10
FULL = (1 << 14) - 1
PAIR_INDEX = {pair: index for index, pair in enumerate(combinations(range(N), 2))}
BALANCED_WORDS = (
    "I???@owM?",
    "I???@ow}?",
    "I??GOPw]?",
    "I???Bpw]?",
    "I??GORw}?",
    "I???opw]?",
    "I??GPow}?",
    "I???orw}?",
    "I??@pow}?",
    "I??Worw}?",
)
TRIVIAL_WORDS = (
    "I??WopfXo",
    "I??FF?^Fo",
    "I??xprFpo",
)


def adjacency(word: str) -> tuple[int, ...]:
    result = [0] * N
    for first, second in decode_graph6(word)[1]:
        result[first] |= 1 << second
        result[second] |= 1 << first
    return tuple(result)


def edge_mask_from_adjacency(values: tuple[int, ...]) -> int:
    return sum(
        1 << PAIR_INDEX[first, second]
        for first, second in combinations(range(N), 2)
        if values[first] & (1 << second)
    )


def transform_edges(value: int, permutation: tuple[int, ...]) -> int:
    result = 0
    for first, second in combinations(range(N), 2):
        if value & (1 << PAIR_INDEX[first, second]):
            image = tuple(sorted((permutation[first], permutation[second])))
            result |= 1 << PAIR_INDEX[image]
    return result


def edges_from_mask(value: int) -> set[tuple[int, int]]:
    return {
        pair for pair, index in PAIR_INDEX.items() if value & (1 << index)
    }


def branch_triples(x_mask: int, cut_mask: int, y_mask: int) -> set[int]:
    edges = edges_from_mask(x_mask)
    result: set[int] = set()
    for isolate in range(N):
        if y_mask & (1 << isolate):
            result.add((1 << 10) | (1 << 11) | (1 << isolate))
    for first, second in combinations(range(N), 2):
        x = (first, second) in edges
        delta = bool(cut_mask & (1 << first)) != bool(
            cut_mask & (1 << second)
        )
        if x:
            result.add((1 << 11) | (1 << first) | (1 << second))
        if x != delta:
            result.add((1 << 10) | (1 << first) | (1 << second))
    for triple in combinations(range(N), 3):
        if sum(tuple(sorted(pair)) in edges for pair in combinations(triple, 2)) % 2:
            result.add(sum(1 << vertex for vertex in triple))
    return result


def selected_complements(
    x_mask: int, z_mask: int, p_mask: int, q_mask: int, y_mask: int
) -> tuple[int, ...]:
    p_triples = branch_triples(x_mask, p_mask, y_mask)
    q_triples = branch_triples(z_mask, q_mask, y_mask)
    result = set()
    # P branch: local isolates 0..9 -> globals 2..11, local specials -> 12,13.
    for triple in p_triples:
        global_triple = 0
        for local in range(N):
            if triple & (1 << local):
                global_triple |= 1 << (local + 2)
        if triple & (1 << 10):
            global_triple |= 1 << 12
        if triple & (1 << 11):
            global_triple |= 1 << 13
        result.add((1 << 0) | (1 << 1) | global_triple)
    # Q branch: its specials correspond to globals 0,1 and its base is 12,13.
    for triple in q_triples:
        global_triple = 0
        for local in range(N):
            if triple & (1 << local):
                global_triple |= 1 << (local + 2)
        if triple & (1 << 10):
            global_triple |= 1 << 0
        if triple & (1 << 11):
            global_triple |= 1 << 1
        result.add((1 << 12) | (1 << 13) | global_triple)
    return tuple(sorted(result))


def face_table(complements: tuple[int, ...]) -> list[bool]:
    table = [False] * (1 << 14)
    for complement in complements:
        facet = FULL ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return table


def minimal_nonfaces(table: list[bool]) -> dict[int, tuple[int, ...]]:
    result: dict[int, list[int]] = {}
    for value in range(1, 1 << 14):
        if table[value]:
            continue
        if all(table[value ^ (bits & -bits)] for bits in _bits(value)):
            result.setdefault(value.bit_count(), []).append(value)
    return {size: tuple(rows) for size, rows in sorted(result.items())}


def _bits(value: int):
    while value:
        bit = value & -value
        yield bit
        value ^= bit


def ridge_histogram(complements: tuple[int, ...]) -> Counter[int]:
    groups: Counter[int] = Counter()
    for complement in complements:
        outside = FULL ^ complement
        for bit in _bits(outside):
            groups[complement | bit] += 1
    return Counter(groups.values())


def connected(complements: tuple[int, ...]) -> bool:
    if not complements:
        return False
    seen = {complements[0]}
    queue = deque(seen)
    values = set(complements)
    while queue:
        value = queue.popleft()
        for bit in _bits(FULL ^ value):
            for removed in _bits(value):
                neighbor = (value | bit) ^ removed
                if neighbor in values and neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
    return len(seen) == len(values)


def face_vector(table: list[bool]) -> tuple[int, ...]:
    histogram = Counter(value.bit_count() for value, yes in enumerate(table) if yes)
    return tuple(histogram[size] for size in range(1, 10))


def full_audit(
    x_mask: int, z_mask: int, p_mask: int, q_mask: int, y_mask: int
) -> tuple[bool, dict[str, object]]:
    complements = selected_complements(x_mask, z_mask, p_mask, q_mask, y_mask)
    table = face_table(complements)
    nonfaces = minimal_nonfaces(table)
    histogram = {size: len(rows) for size, rows in nonfaces.items()}
    graph_edges = {
        (0, 12),
        (1, 12),
        (0, 13),
        (1, 13),
    }
    decoded_edges = {
        tuple(vertex for vertex in range(14) if row & (1 << vertex))
        for row in nonfaces.get(2, ())
    }
    a = histogram.get(4, 0)
    vector = face_vector(table)
    allowed_sizes = set(histogram) <= {2, 4, 5}
    ok = (
        allowed_sizes
        and decoded_edges == graph_edges
        and a <= 7
        and len(complements) == 130 - 2 * a
        and vector[-1] == len(complements)
        and ridge_histogram(complements) == Counter({2: vector[-2]})
        and connected(complements)
    )
    return ok, {
        "facets": len(complements),
        "nonfaces": histogram,
        "f_vector": vector,
        "ridge_histogram": dict(ridge_histogram(complements)),
        "connected": connected(complements),
    }


def balanced_group() -> tuple[tuple[int, ...], ...]:
    return tuple(
        first + tuple(value + 5 for value in second)
        for first in permutations(range(5))
        for second in permutations(range(5))
    )


def balanced_census() -> None:
    group = balanced_group()
    graph_masks = tuple(edge_mask_from_adjacency(adjacency(word)) for word in BALANCED_WORDS)
    automorphisms = tuple(
        tuple(permutation for permutation in group if transform_edges(value, permutation) == value)
        for value in graph_masks
    )
    print("balanced automorphism orders", tuple(map(len, automorphisms)))
    total_double_cosets = 0
    passing = []
    for x_index, x_mask in enumerate(graph_masks):
        aut_x = automorphisms[x_index]
        for z_index, z_mask in enumerate(graph_masks):
            embeddings = {transform_edges(z_mask, permutation) for permutation in group}
            total_double_cosets += 0
            while embeddings:
                embedded = min(embeddings)
                orbit = {transform_edges(embedded, permutation) for permutation in aut_x}
                embeddings.difference_update(orbit)
                total_double_cosets += 1
                ok, audit = full_audit(x_mask, embedded, 0b1111100000, 0b1111100000, 0)
                if ok:
                    passing.append((x_index, z_index, x_mask, embedded, audit))
    print("balanced double cosets", total_double_cosets)
    print("balanced passing", len(passing))
    print(
        "balanced passing hist",
        Counter(
            (
                item[-1]["nonfaces"].get(4, 0),
                item[-1]["nonfaces"].get(5, 0),
                item[-1]["facets"],
            )
            for item in passing
        ),
    )
    print("balanced f-vectors", Counter(item[-1]["f_vector"] for item in passing))
    if "--emit" in sys.argv:
        for item in passing:
            print("BALANCED_PAIR", item[2], item[3])


def _adjacency_from_edge_mask(value: int) -> tuple[int, ...]:
    result = [0] * N
    for (first, second), index in PAIR_INDEX.items():
        if value & (1 << index):
            result[first] |= 1 << second
            result[second] |= 1 << first
    return tuple(result)


def _automorphisms(value: int) -> tuple[tuple[int, ...], ...]:
    adjacency_values = _adjacency_from_edge_mask(value)
    degree_cells: dict[int, list[int]] = {}
    for vertex, neighbors in enumerate(adjacency_values):
        degree_cells.setdefault(neighbors.bit_count(), []).append(vertex)
    cells = tuple(degree_cells[key] for key in sorted(degree_cells))
    result = []
    # The largest product in the three exceptional graphs is 6!4!=17280.
    choices = [tuple(permutations(cell)) for cell in cells]

    def recurse(index: int, image: list[int]) -> None:
        if index == len(cells):
            permutation = [0] * N
            for cell, targets in zip(cells, image):
                for source, target in zip(cell, targets):
                    permutation[source] = target
            value_permutation = tuple(permutation)
            if transform_edges(value, value_permutation) == value:
                result.append(value_permutation)
            return
        for target in choices[index]:
            recurse(index + 1, image + [target])

    recurse(0, [])
    return tuple(result)


def _embedded_subgraphs(source: int, allowed: int) -> set[int]:
    source_adjacency = _adjacency_from_edge_mask(source)
    allowed_adjacency = _adjacency_from_edge_mask(allowed)
    source_order = tuple(
        sorted(range(N), key=lambda vertex: (-source_adjacency[vertex].bit_count(), vertex))
    )
    mapping: dict[int, int] = {}
    used = 0
    result: set[int] = set()

    def recurse(index: int) -> None:
        nonlocal used
        if index == N:
            permutation = tuple(mapping[vertex] for vertex in range(N))
            result.add(transform_edges(source, permutation))
            return
        vertex = source_order[index]
        source_degree = source_adjacency[vertex].bit_count()
        for image in range(N):
            if used & (1 << image) or allowed_adjacency[image].bit_count() < source_degree:
                continue
            if any(
                source_adjacency[vertex] & (1 << prior)
                and not (allowed_adjacency[image] & (1 << mapped))
                for prior, mapped in mapping.items()
            ):
                continue
            mapping[vertex] = image
            used |= 1 << image
            recurse(index + 1)
            used ^= 1 << image
            del mapping[vertex]

    recurse(0)
    return result


def trivial_census() -> None:
    graph_masks = tuple(edge_mask_from_adjacency(adjacency(word)) for word in TRIVIAL_WORDS)
    automorphisms = tuple(_automorphisms(value) for value in graph_masks)
    print("trivial automorphism orders", tuple(map(len, automorphisms)))
    complete = (1 << len(PAIR_INDEX)) - 1
    total_double_cosets = 0
    passing = []
    failing = []
    embedding_counts = {}
    for x_index, x_mask in enumerate(graph_masks):
        aut_x = automorphisms[x_index]
        for z_index, z_mask in enumerate(graph_masks):
            embeddings = _embedded_subgraphs(z_mask, complete ^ x_mask)
            embedding_counts[x_index, z_index] = len(embeddings)
            while embeddings:
                embedded = min(embeddings)
                orbit = {transform_edges(embedded, permutation) for permutation in aut_x}
                embeddings.difference_update(orbit)
                total_double_cosets += 1
                ok, audit = full_audit(x_mask, embedded, 0, 0, 0)
                if ok:
                    passing.append((x_index, z_index, x_mask, embedded, audit))
                else:
                    failing.append((x_index, z_index, audit))
    print("trivial embedding counts", embedding_counts)
    print("trivial double cosets", total_double_cosets)
    print("trivial passing", len(passing))
    print(
        "trivial passing hist",
        Counter(
            (
                item[-1]["nonfaces"].get(4, 0),
                item[-1]["nonfaces"].get(5, 0),
                item[-1]["facets"],
            )
            for item in passing
        ),
    )
    print("trivial f-vectors", Counter(item[-1]["f_vector"] for item in passing))
    print(
        "trivial failure profiles",
        Counter(
            (
                tuple(sorted(item[-1]["nonfaces"].items())),
                item[-1]["facets"],
                tuple(sorted(item[-1]["ridge_histogram"].items())),
                item[-1]["connected"],
            )
            for item in failing
        ),
    )


if __name__ == "__main__":
    if "--trivial" not in sys.argv:
        balanced_census()
    if "--balanced" not in sys.argv:
        trivial_census()
