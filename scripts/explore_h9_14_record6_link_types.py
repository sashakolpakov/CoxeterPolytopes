#!/usr/bin/env python3
"""Classify one special-vertex link in the H9/14 record-6 shadow.

Input is a graph6 stream of ordinary graphs on ten vertices, normally
``geng -q -t 10``.  Fixing the eleventh isolate as switching base identifies
such a graph ``X`` with the two-graph ``delta X`` of possible link-facet
complements.  The script verifies the induced-2K2 restriction, constructs
the complete eight-dimensional link boundary, and applies exact incidence,
Lanner-order, strict-g, and mod-two sphere-link audits.

This is a structural census, not by itself a record-6 nonexistence proof.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
import sys

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6, macaulay_raise
from coxeter_polytopes.combinatorics.simplicial_homology import (
    first_mod2_homology_sphere_link_failure,
)


VERTEX_COUNT = 11
FULL_MASK = (1 << VERTEX_COUNT) - 1


def mask(vertices: tuple[int, ...]) -> int:
    return sum(1 << vertex for vertex in vertices)


def graph_edge_masks(graph6: str) -> tuple[int, ...]:
    order, edges = decode_graph6(graph6)
    if order != 10:
        raise ValueError("record-6 link census expects ten-vertex graphs")
    # Local vertex zero is the switching base; geng vertices become 1,...,10.
    return tuple((1 << (first + 1)) | (1 << (second + 1)) for first, second in edges)


def two_graph_triples(edges: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        mask(triple)
        for triple in combinations(range(VERTEX_COUNT), 3)
        if sum((mask(triple) & edge) == edge for edge in edges) % 2
    )


def generated_face_table(complements: tuple[int, ...]) -> tuple[bool, ...]:
    table = [False] * (1 << VERTEX_COUNT)
    table[0] = True
    for complement in complements:
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def minimal_nonfaces(table: tuple[bool, ...]) -> tuple[int, ...]:
    return tuple(
        subset
        for subset in range(1, 1 << VERTEX_COUNT)
        if not table[subset]
        and all(
            table[subset ^ (1 << vertex)]
            for vertex in range(VERTEX_COUNT)
            if subset & (1 << vertex)
        )
    )


def dual_connected(complements: tuple[int, ...]) -> bool:
    groups: dict[int, list[int]] = {}
    for complement in complements:
        outside = FULL_MASK ^ complement
        while outside:
            bit = outside & -outside
            groups.setdefault(complement | bit, []).append(complement)
            outside ^= bit
    if any(len(group) != 2 for group in groups.values()):
        return False
    adjacency = {facet: set() for facet in complements}
    for first, second in groups.values():
        adjacency[first].add(second)
        adjacency[second].add(first)
    if not complements:
        return False
    seen = {complements[0]}
    stack = [complements[0]]
    while stack:
        current = stack.pop()
        fresh = adjacency[current] - seen
        seen.update(fresh)
        stack.extend(fresh)
    return len(seen) == len(complements)


def audit(graph6: str) -> tuple[str, int, int, int] | None:
    edges = graph_edge_masks(graph6)
    # The 0/2 four-set ridge rule for delta X is equivalent to excluding
    # triangles and induced 2K2s in X.  Verify both rather than trusting geng.
    for triple in combinations(range(1, VERTEX_COUNT), 3):
        value = mask(triple)
        if sum((value & edge) == edge for edge in edges) == 3:
            return None
    for four in combinations(range(1, VERTEX_COUNT), 4):
        value = mask(four)
        induced = tuple(edge for edge in edges if value & edge == edge)
        if len(induced) == 2 and not (induced[0] & induced[1]):
            return None

    triples = two_graph_triples(edges)
    for four in combinations(range(VERTEX_COUNT), 4):
        value = mask(four)
        if sum((triple & value) == triple for triple in triples) not in (0, 2):
            raise AssertionError("two-graph conversion violated a ridge equation")
    table = generated_face_table(triples)
    if not all(table[mask(pair)] for pair in combinations(range(VERTEX_COUNT), 2)):
        return None
    nonfaces = minimal_nonfaces(table)
    if not nonfaces or max(value.bit_count() for value in nonfaces) > 5:
        return None
    if not dual_connected(triples):
        return None
    face_vector = tuple(
        sum(table[value] for value in range(1 << VERTEX_COUNT) if value.bit_count() == size)
        for size in range(1, 9)
    )
    missing_triangles = sum(value.bit_count() == 3 for value in nonfaces)
    g3 = 4 - missing_triangles
    g4 = face_vector[3] - 6 * face_vector[2] + 21 * face_vector[1] - 56 * face_vector[0] + 126
    if not 0 <= missing_triangles <= 3 or not 1 <= g4 <= macaulay_raise(g3, 3):
        return None
    if first_mod2_homology_sphere_link_failure(table, 7) is not None:
        return None
    return graph6, len(triples), missing_triangles, g4


def main() -> None:
    input_count = 0
    ridge_type_count = 0
    survivors: list[tuple[str, int, int, int]] = []
    for raw in sys.stdin:
        graph6 = raw.strip()
        if not graph6:
            continue
        input_count += 1
        edges = graph_edge_masks(graph6)
        forbidden = False
        for four in combinations(range(1, VERTEX_COUNT), 4):
            value = mask(four)
            induced = tuple(edge for edge in edges if value & edge == edge)
            if len(induced) == 2 and not (induced[0] & induced[1]):
                forbidden = True
                break
        if forbidden:
            continue
        ridge_type_count += 1
        result = audit(graph6)
        if result is not None:
            survivors.append(result)
    print(
        json.dumps(
            {
                "input_triangle_free_unlabelled_graphs": input_count,
                "triangle_and_induced_2k2_free_types": ridge_type_count,
                "audited_link_types": len(survivors),
                "invariant_histogram": [
                    {
                        "facet_count": key[0],
                        "missing_triangle_count": key[1],
                        "g4": key[2],
                        "multiplicity": multiplicity,
                    }
                    for key, multiplicity in sorted(
                        Counter((row[1], row[2], row[3]) for row in survivors).items()
                    )
                ],
                "survivor_graph6": [row[0] for row in survivors],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
