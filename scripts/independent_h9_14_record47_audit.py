#!/usr/bin/env python3
"""Hostile replay of the H9/14 record-47 incidence and Gale certificates."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    generated_face_table,
)
from coxeter_polytopes.combinatorics.gale_chirotope import (
    audit_uniform_rank_four_gale,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record47_selection.json"
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
GRAPH6 = "M?????????O?o?g?_"
EDGES = ((0, 11), (0, 12), (1, 12), (0, 13), (2, 13), (12, 13))
FIVES = (
    (3, 4, 5, 6, 7),
    (1, 2, 8, 9, 10),
    (1, 8, 9, 10, 11),
    (2, 8, 9, 10, 11),
)


def _mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


def _minimal_transversals(hyperedges: tuple[int, ...]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, 1 << 14):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if all(
            not all((candidate ^ bit) & edge for edge in hyperedges)
            for bit in (
                1 << vertex
                for vertex in range(14)
                if candidate & (1 << vertex)
            )
        ):
            answers.append(candidate)
    return tuple(answers)


def replay() -> dict[str, object]:
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    order, decoded_edges = decode_graph6(GRAPH6)
    if order != 14 or decoded_edges != EDGES:
        raise AssertionError("independent record-47 graph decode failed")

    blocker = _minimal_transversals(
        tuple(_mask(edge) for edge in EDGES)
        + tuple(_mask(face) for face in FIVES)
    )
    if Counter(value.bit_count() for value in blocker) != {5: 90}:
        raise AssertionError("independent blocker replay is not pure with 90 members")
    if blocker != tuple(archive["incidence"]["facet_complement_masks"]):
        raise AssertionError("independent blocker differs from archived facets")

    source_records = [
        record for record in support["records"] if record["survives_second_layer"]
    ]
    source = source_records[47]
    if source["graph6"] != GRAPH6:
        raise AssertionError("independent source-index replay failed")
    global_audit = audit_global_facet_selection(
        GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        blocker,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("independent global record-47 audit failed")

    signs = tuple(archive["gale_chirotope_signs"])
    gale_audit = audit_uniform_rank_four_gale(blocker, signs)
    if not gale_audit.passes:
        raise AssertionError("independent record-47 Gale audit failed")
    face_table = generated_face_table(blocker)
    packed = bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )
    packed_digest = sha256(packed).hexdigest()
    if packed_digest != archive["incidence"]["packed_face_table_sha256"]:
        raise AssertionError("independent face-table digest mismatch")
    return {
        "minimal_transversal_size_histogram": {"5": 90},
        "face_vector": list(global_audit.face_vector),
        "ridge_degree_histogram": {
            str(key): value
            for key, value in global_audit.ridge_degree_histogram.items()
        },
        "dual_component_count": global_audit.dual_component_count,
        "all_links_mod2_homology_spheres": (
            global_audit.all_links_are_mod2_homology_spheres
        ),
        "packed_face_table_sha256": packed_digest,
        "gale_sign_count": gale_audit.sign_count,
        "positive_circuit_count": len(gale_audit.positive_circuit_masks),
        "gale_audit_passes": gale_audit.passes,
    }


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
