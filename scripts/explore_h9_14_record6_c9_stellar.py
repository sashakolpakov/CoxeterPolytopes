#!/usr/bin/env python3
"""Search one-stellar multiwedges of ``C(9,5)`` for H9/14 record 6.

This exactly enumerates every edge/triangle stellar subdivision of the cyclic
5-polytope on nine vertices and every multiset of four simplicial wedges.
It is a bounded witness search, not a whole-shadow non-existence proof.
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


def _cyclic9_facets() -> tuple[frozenset[int], ...]:
    facets = []
    for chosen in combinations(range(9), 5):
        omitted = tuple(vertex for vertex in range(9) if vertex not in chosen)
        if all(
            sum(first < selected < second for selected in chosen) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            facets.append(frozenset(chosen))
    return tuple(facets)


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
    initial = _cyclic9_facets()
    profiles = tuple(combinations_with_replacement(range(10), 4))
    counters = {
        "stellar_faces": 0,
        "multiwedge_profiles": 0,
        "matching_edge_graphs": 0,
        "audited_candidates": 0,
    }
    for stellar_face in _faces(initial, (2, 3)):
        counters["stellar_faces"] += 1
        base = _stellar(initial, stellar_face, 9)
        clutter = _minimal_nonfaces(10, base)
        for carriers in profiles:
            counters["multiwedge_profiles"] += 1
            blocks = _blocks(carriers)
            expanded = _expand(clutter, blocks)
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
                    "stellar_face": sorted(stellar_face),
                    "wedge_carriers": list(carriers),
                    "blocks": {
                        str(vertex): list(block) for vertex, block in blocks.items()
                    },
                    "mapping": {str(key): value for key, value in mapping.items()},
                    "expanded_minimal_nonfaces": [
                        sorted(value) for value in expanded
                    ],
                    "selected_facet_complement_masks": list(relabelled),
                    "audit": audit.to_record(),
                    "counters": counters,
                }
    return {
        "status": "no_witness_in_declared_C9_stellar_multiwedge_family",
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
