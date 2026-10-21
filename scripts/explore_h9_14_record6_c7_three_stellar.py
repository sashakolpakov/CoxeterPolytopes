#!/usr/bin/env python3
"""Search three-stellar multiwedges of ``C(7,5)`` for H9/14 record 6.

The search exhausts labelled sequences of three edge/triangle stellar
subdivisions, deduplicates the resulting ten-vertex base complexes, and then
tests every multiset of four wedge carriers.  Fast missing-edge filters run
before blocker construction.  A miss is scoped to this construction family.
"""

from __future__ import annotations

from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6
from scripts.explore_h9_14_record6_multiwedge import _blocks, _expand
from scripts.explore_h9_14_record7_stellar import (
    _edge_isomorphism,
    _faces,
    _minimal_five_transversals,
    _minimal_nonfaces,
    _relabel_masks,
    _stellar,
)


ROOT = Path(__file__).resolve().parents[1]
RECORD_INDEX = 6


def _cyclic_facets(vertex_count: int) -> tuple[frozenset[int], ...]:
    facets = []
    for chosen in combinations(range(vertex_count), 5):
        omitted = tuple(
            vertex for vertex in range(vertex_count) if vertex not in chosen
        )
        if all(
            sum(first < selected < second for selected in chosen) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            facets.append(frozenset(chosen))
    return tuple(facets)


def _key(facets: tuple[frozenset[int], ...]) -> tuple[int, ...]:
    return tuple(sorted(sum(1 << vertex for vertex in facet) for facet in facets))


def _missing_edges(
    facets: tuple[frozenset[int], ...]
) -> tuple[frozenset[int], ...]:
    present = {
        frozenset(edge)
        for facet in facets
        for edge in combinations(sorted(facet), 2)
    }
    return tuple(
        frozenset(edge)
        for edge in combinations(range(10), 2)
        if frozenset(edge) not in present
    )


def _possible_target_edges(
    base_edges: tuple[frozenset[int], ...], blocks: dict[int, tuple[int, ...]]
) -> tuple[frozenset[int], ...] | None:
    expanded = []
    for edge in base_edges:
        value = frozenset(label for vertex in edge for label in blocks[vertex])
        if len(value) == 3 or len(value) > 5:
            return None
        if len(value) == 2:
            expanded.append(value)
    if len(expanded) != 3:
        return None
    active = set().union(*expanded)
    if len(active) != 3:
        return None
    if any(sum(vertex in edge for edge in expanded) != 2 for vertex in active):
        return None
    return tuple(expanded)


def search() -> dict[str, object]:
    source = json.loads(
        (ROOT / "data/enumerations/h9_14_nerve_support.json").read_text(
            encoding="utf-8"
        )
    )
    record = [item for item in source["records"] if item["survives_second_layer"]][
        RECORD_INDEX
    ]
    _, target_edges = decode_graph6(record["graph6"])
    initial = _cyclic_facets(7)
    profiles = tuple(combinations_with_replacement(range(10), 4))
    seen_bases: set[tuple[int, ...]] = set()
    counters = {
        "labelled_subdivision_sequences": 0,
        "distinct_base_complexes": 0,
        "multiwedge_profiles": 0,
        "edge_filtered_profiles": 0,
        "matching_edge_graphs": 0,
        "audited_candidates": 0,
    }
    for first_face in _faces(initial, (2, 3)):
        first = _stellar(initial, first_face, 7)
        for second_face in _faces(first, (2, 3)):
            second = _stellar(first, second_face, 8)
            for third_face in _faces(second, (2, 3)):
                counters["labelled_subdivision_sequences"] += 1
                base = _stellar(second, third_face, 9)
                key = _key(base)
                if key in seen_bases:
                    continue
                seen_bases.add(key)
                counters["distinct_base_complexes"] += 1
                base_edges = _missing_edges(base)
                for carriers in profiles:
                    counters["multiwedge_profiles"] += 1
                    blocks = _blocks(carriers)
                    edges = _possible_target_edges(base_edges, blocks)
                    if edges is None:
                        continue
                    counters["edge_filtered_profiles"] += 1
                    mapping = _edge_isomorphism(edges, target_edges)
                    if mapping is None:
                        continue
                    counters["matching_edge_graphs"] += 1
                    clutter = _minimal_nonfaces(10, base)
                    expanded = _expand(clutter, blocks)
                    if any(len(row) == 3 or len(row) > 5 for row in expanded):
                        continue
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
                            "stellar_faces": [
                                sorted(first_face),
                                sorted(second_face),
                                sorted(third_face),
                            ],
                            "wedge_carriers": list(carriers),
                            "blocks": {
                                str(vertex): list(block)
                                for vertex, block in blocks.items()
                            },
                            "mapping": {
                                str(key): value for key, value in mapping.items()
                            },
                            "expanded_minimal_nonfaces": [
                                sorted(value) for value in expanded
                            ],
                            "selected_facet_complement_masks": list(relabelled),
                            "audit": audit.to_record(),
                            "counters": counters,
                        }
    return {
        "status": "no_witness_in_declared_C7_three_stellar_multiwedge_family",
        "record_index": RECORD_INDEX,
        "graph6": record["graph6"],
        "stellar_face_sizes": [2, 3],
        "wedge_multiplicity": 4,
        "counters": counters,
    }


def main() -> None:
    print(json.dumps(search(), sort_keys=True))


if __name__ == "__main__":
    main()
