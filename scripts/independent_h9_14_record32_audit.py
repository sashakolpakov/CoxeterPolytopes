#!/usr/bin/env python3
"""Hostile replay of the H9/14 record-32 convex completion certificate."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from itertools import combinations
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
ARCHIVE = ROOT / "data/enumerations/h9_14_record32_polytopality.json"
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
GRAPH6 = "M?????????O?O?g??"
EDGES = ((0, 11), (1, 12), (0, 13), (2, 13))
FOURS = ((2, 3, 4, 5), (3, 4, 5, 11))
FIVES = ((6, 7, 8, 9, 10),)
BASE_VERTICES = (0, 1, 2, 11, 12, 13, 14)
BASE_MISSING_EDGES = (
    (0, 11), (1, 12), (0, 13), (2, 13), (2, 14), (11, 14)
)
BASE_CYCLE = (11, 13, 14, 0, 2)
BASE_APICES = (1, 12)
CLONE_BLOCK = (3, 4, 5)
JOIN_BLOCK = (6, 7, 8, 9, 10)


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


def _all_faces(facets: set[frozenset[int]]) -> set[frozenset[int]]:
    faces = {frozenset()}
    for facet in facets:
        for size in range(1, len(facet) + 1):
            faces.update(frozenset(face) for face in combinations(facet, size))
    return faces


def _minimal_nonfaces(
    vertices: tuple[int, ...], facets: set[frozenset[int]]
) -> set[frozenset[int]]:
    faces = _all_faces(facets)
    result = set()
    for size in range(1, len(vertices) + 1):
        for subset in combinations(vertices, size):
            candidate = frozenset(subset)
            if candidate not in faces and all(
                candidate - {vertex} in faces for vertex in candidate
            ):
                result.add(candidate)
    return result


def replay() -> dict[str, object]:
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    order, decoded_edges = decode_graph6(GRAPH6)
    if order != 14 or decoded_edges != EDGES:
        raise AssertionError("independent record-32 graph decode failed")

    hyperedges = (
        tuple(_mask(face) for face in EDGES)
        + tuple(_mask(face) for face in FOURS)
        + tuple(_mask(face) for face in FIVES)
    )
    blocker = _minimal_transversals(hyperedges)
    if Counter(value.bit_count() for value in blocker) != {5: 110}:
        raise AssertionError("independent blocker is not pure with 110 members")
    if blocker != tuple(archive["incidence"]["facet_complement_masks"]):
        raise AssertionError("independent blocker differs from archived facets")

    records = [record for record in support["records"] if record["survives_second_layer"]]
    source = records[32]
    if source["graph6"] != GRAPH6:
        raise AssertionError("independent source-index replay failed")
    global_audit = audit_global_facet_selection(
        GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        blocker,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("independent global record-32 audit failed")

    # Rebuild the pentagonal-bipyramid base without importing the production
    # record-32 module.  Its ten triangles are one cycle edge plus one apex.
    base_facets = {
        frozenset((BASE_CYCLE[index], BASE_CYCLE[(index + 1) % 5], apex))
        for index in range(5)
        for apex in BASE_APICES
    }
    base_nonfaces = _minimal_nonfaces(BASE_VERTICES, base_facets)
    if base_nonfaces != {frozenset(edge) for edge in BASE_MISSING_EDGES}:
        raise AssertionError("independent pentagonal-bipyramid base changed")
    expanded = set()
    for face in base_nonfaces:
        lifted = set()
        for vertex in face:
            lifted.update(CLONE_BLOCK if vertex == 14 else (vertex,))
        expanded.add(frozenset(lifted))
    expanded.add(frozenset(JOIN_BLOCK))
    target = (
        {frozenset(face) for face in EDGES}
        | {frozenset(face) for face in FOURS}
        | {frozenset(face) for face in FIVES}
    )
    if expanded != target:
        raise AssertionError("independent multiwedge expansion changed")

    signs = tuple(archive["gale_chirotope_signs"])
    gale_audit = audit_uniform_rank_four_gale(blocker, signs)
    if not gale_audit.passes:
        raise AssertionError("independent record-32 Gale audit failed")
    face_table = generated_face_table(blocker)
    packed = bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )
    packed_digest = sha256(packed).hexdigest()
    if packed_digest != archive["incidence"]["packed_face_table_sha256"]:
        raise AssertionError("independent face-table digest mismatch")
    return {
        "minimal_transversal_size_histogram": {"5": 110},
        "face_vector": list(global_audit.face_vector),
        "ridge_degree_histogram": {
            str(key): value for key, value in global_audit.ridge_degree_histogram.items()
        },
        "dual_component_count": global_audit.dual_component_count,
        "all_links_mod2_homology_spheres": (
            global_audit.all_links_are_mod2_homology_spheres
        ),
        "base_face_vector": [7, 15, 10],
        "ordinary_polytopal": True,
        "packed_face_table_sha256": packed_digest,
        "gale_sign_count": gale_audit.sign_count,
        "positive_circuit_count": len(gale_audit.positive_circuit_masks),
        "gale_audit_passes": gale_audit.passes,
    }


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
