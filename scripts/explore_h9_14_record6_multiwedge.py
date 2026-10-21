#!/usr/bin/env python3
"""Search repeated-wedge stellar constructions for H9/14 record 6.

The record-7 finder restricts the four wedges to distinct base vertices.
This companion exhausts all multisets of four wedge carriers after two
edge/triangle stellar subdivisions of ``C(8,5)``.  It is a scoped witness
finder: a negative result applies only to this declared construction family.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations_with_replacement
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6
from scripts.explore_h9_14_record7_stellar import (
    _cyclic_facets,
    _edge_isomorphism,
    _faces,
    _mask,
    _minimal_five_transversals,
    _minimal_nonfaces,
    _relabel_masks,
    _stellar,
)


ROOT = Path(__file__).resolve().parents[1]
RECORD_INDEX = 6


def _blocks(carriers: tuple[int, ...]) -> dict[int, tuple[int, ...]]:
    multiplicities = Counter(carriers)
    next_label = 10
    result: dict[int, tuple[int, ...]] = {}
    for vertex in range(10):
        block = [vertex]
        for _ in range(multiplicities[vertex]):
            block.append(next_label)
            next_label += 1
        result[vertex] = tuple(block)
    if next_label != 14:
        raise AssertionError("four wedges must produce fourteen vertices")
    return result


def _expand(
    clutter: tuple[frozenset[int], ...], blocks: dict[int, tuple[int, ...]]
) -> tuple[frozenset[int], ...]:
    return tuple(
        frozenset(label for vertex in nonface for label in blocks[vertex])
        for nonface in clutter
    )


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
    initial = _cyclic_facets()
    counters = {
        "subdivision_pairs": 0,
        "multiwedge_profiles": 0,
        "matching_edge_graphs": 0,
        "audited_candidates": 0,
    }
    profiles = tuple(combinations_with_replacement(range(10), 4))
    for first_face in _faces(initial, (2, 3)):
        first = _stellar(initial, first_face, 8)
        for second_face in _faces(first, (2, 3)):
            counters["subdivision_pairs"] += 1
            second = _stellar(first, second_face, 9)
            clutter = _minimal_nonfaces(10, second)
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
                        "first_face": sorted(first_face),
                        "second_face": sorted(second_face),
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
        "status": "no_witness_in_declared_repeated_wedge_family",
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
