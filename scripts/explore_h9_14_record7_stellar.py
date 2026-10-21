#!/usr/bin/env python3
"""Search a rigid stellar-subdivision/multiwedge family for H9/14 record 7.

This is a witness finder, not a theorem backend.  It starts with the boundary
of ``C(8,5)``, performs two stellar subdivisions at edges or triangles, and
then applies four simplicial wedges.  A hit is accepted only after the common
exact facet/face/ridge/link audit.  A miss excludes only this declared family.
"""

from __future__ import annotations

import argparse
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


ROOT = Path(__file__).resolve().parents[1]
RECORD_INDEX = 7


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _cyclic_facets() -> tuple[frozenset[int], ...]:
    facets = []
    for positions in combinations(range(8), 5):
        omitted = tuple(
            position for position in range(8) if position not in positions
        )
        if all(
            sum(first < selected < second for selected in positions) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            facets.append(frozenset(positions))
    return tuple(facets)


def _stellar(
    facets: tuple[frozenset[int], ...],
    face: frozenset[int],
    new_vertex: int,
) -> tuple[frozenset[int], ...]:
    star = tuple(facet for facet in facets if face <= facet)
    if not star:
        raise ValueError("stellar face is not present")
    result = set(facets) - set(star)
    result.update(
        (facet - {removed}) | {new_vertex}
        for facet in star
        for removed in face
    )
    return tuple(sorted(result, key=lambda item: tuple(sorted(item))))


def _faces(
    facets: tuple[frozenset[int], ...], sizes: tuple[int, ...]
) -> tuple[frozenset[int], ...]:
    return tuple(
        sorted(
            {
                frozenset(value)
                for facet in facets
                for size in sizes
                for value in combinations(sorted(facet), size)
            },
            key=lambda item: (len(item), tuple(sorted(item))),
        )
    )


def _minimal_nonfaces(
    vertex_count: int, facets: tuple[frozenset[int], ...]
) -> tuple[frozenset[int], ...]:
    face_masks = {0}
    for facet in facets:
        ordered = tuple(sorted(facet))
        for size in range(1, len(facet) + 1):
            face_masks.update(_mask(value) for value in combinations(ordered, size))
    result = []
    for value in range(1, 1 << vertex_count):
        if value in face_masks:
            continue
        if all(value ^ bit in face_masks for bit in _bits(value)):
            result.append(
                frozenset(vertex for vertex in range(vertex_count) if value >> vertex & 1)
            )
    return tuple(result)


def _bits(value: int) -> tuple[int, ...]:
    result = []
    while value:
        bit = value & -value
        result.append(bit)
        value ^= bit
    return tuple(result)


def _expanded_clutter(
    clutter: tuple[frozenset[int], ...], cloned: tuple[int, ...]
) -> tuple[frozenset[int], ...]:
    clone_labels = {vertex: 10 + index for index, vertex in enumerate(cloned)}
    result = []
    for nonface in clutter:
        expanded = set(nonface)
        expanded.update(
            clone_labels[vertex] for vertex in nonface if vertex in clone_labels
        )
        result.append(frozenset(expanded))
    return tuple(result)


def _edge_isomorphism(
    candidate_edges: tuple[frozenset[int], ...],
    target_edges: tuple[tuple[int, int], ...],
) -> dict[int, int] | None:
    candidate_active = sorted(set().union(*candidate_edges))
    target_edge_set = {frozenset(edge) for edge in target_edges}
    target_active = sorted(set().union(*target_edge_set))
    if len(candidate_active) != len(target_active):
        return None
    candidate_edge_set = set(candidate_edges)
    for image in permutations(target_active):
        mapping = dict(zip(candidate_active, image))
        if {
            frozenset(mapping[vertex] for vertex in edge)
            for edge in candidate_edge_set
        } == target_edge_set:
            candidate_isolates = [
                vertex for vertex in range(14) if vertex not in candidate_active
            ]
            target_isolates = [
                vertex for vertex in range(14) if vertex not in target_active
            ]
            mapping.update(zip(candidate_isolates, target_isolates))
            return mapping
    return None


def _minimal_five_transversals(
    clutter: tuple[frozenset[int], ...]
) -> tuple[int, ...]:
    result = []
    for candidate in combinations(range(14), 5):
        value = frozenset(candidate)
        if not all(value & row for row in clutter):
            continue
        if not all(
            any((value - {vertex}).isdisjoint(row) for row in clutter)
            for vertex in value
        ):
            continue
        result.append(_mask(candidate))
    return tuple(result)


def _relabel_masks(masks: tuple[int, ...], mapping: dict[int, int]) -> tuple[int, ...]:
    return tuple(
        sorted(
            _mask(mapping[vertex] for vertex in range(14) if value >> vertex & 1)
            for value in masks
        )
    )


def search(first_sizes: tuple[int, ...], second_sizes: tuple[int, ...]) -> dict:
    source = json.loads(
        (ROOT / "data/enumerations/h9_14_nerve_support.json").read_text(
            encoding="utf-8"
        )
    )
    record = [item for item in source["records"] if item["survives_second_layer"]][
        RECORD_INDEX
    ]
    _, target_edges = decode_graph6(record["graph6"])
    initial = _cyclic_facets()
    counters = {
        "subdivision_pairs": 0,
        "clone_choices": 0,
        "matching_edge_graphs": 0,
        "audited_candidates": 0,
    }
    for first_face in _faces(initial, first_sizes):
        first = _stellar(initial, first_face, 8)
        for second_face in _faces(first, second_sizes):
            counters["subdivision_pairs"] += 1
            second = _stellar(first, second_face, 9)
            clutter = _minimal_nonfaces(10, second)
            for cloned in combinations(range(10), 4):
                counters["clone_choices"] += 1
                expanded = _expanded_clutter(clutter, cloned)
                if any(len(row) == 3 or len(row) > 5 for row in expanded):
                    continue
                edges = tuple(row for row in expanded if len(row) == 2)
                if len(edges) != len(target_edges):
                    continue
                mapping = _edge_isomorphism(edges, target_edges)
                if mapping is None:
                    continue
                counters["matching_edge_graphs"] += 1
                selected = _minimal_five_transversals(expanded)
                if not selected:
                    continue
                relabelled = _relabel_masks(selected, mapping)
                try:
                    audit = audit_global_facet_selection(
                        record["graph6"],
                        record["ridge_support_core_masks"],
                        record["dehn_sommerville_schedule"],
                        relabelled,
                    )
                except ValueError:
                    continue
                counters["audited_candidates"] += 1
                if audit.passes_connected_homology_audit:
                    return {
                        "status": "witness",
                        "record_index": RECORD_INDEX,
                        "graph6": record["graph6"],
                        "first_face": sorted(first_face),
                        "second_face": sorted(second_face),
                        "cloned": list(cloned),
                        "mapping": {str(key): value for key, value in mapping.items()},
                        "expanded_minimal_nonfaces": [
                            sorted(value) for value in expanded
                        ],
                        "selected_facet_complement_masks": list(relabelled),
                        "audit": audit.to_record(),
                        "counters": counters,
                    }
    return {
        "status": "no_witness_in_declared_family",
        "record_index": RECORD_INDEX,
        "graph6": record["graph6"],
        "first_face_sizes": list(first_sizes),
        "second_face_sizes": list(second_sizes),
        "counters": counters,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--first-sizes", default="2,3")
    parser.add_argument("--second-sizes", default="2,3")
    args = parser.parse_args()
    first_sizes = tuple(int(value) for value in args.first_sizes.split(","))
    second_sizes = tuple(int(value) for value in args.second_sizes.split(","))
    print(json.dumps(search(first_sizes, second_sizes), sort_keys=True))


if __name__ == "__main__":
    main()
