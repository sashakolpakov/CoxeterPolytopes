#!/usr/bin/env python3
"""Hostile replay of the record-32 ``a=0`` selection and Coxeter no-go."""

from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    generated_face_table,
)
from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    is_finite_coxeter_labeling,
    lanner_placement_codes,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record32_a0_selection.json"
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
GRAPH6 = "M?????????O?O?g??"
EDGES = ((0, 11), (1, 12), (0, 13), (2, 13))
FIVES = (
    (4, 5, 7, 8, 9),
    (1, 2, 3, 6, 10),
    (2, 3, 4, 6, 10),
    (3, 4, 6, 10, 11),
    (5, 7, 8, 9, 12),
    (5, 7, 8, 9, 13),
)
CORE = (5, 7, 8, 9)
OUTERS = (4, 12, 13)
SUPPORTS = tuple(tuple(sorted((*CORE, outer))) for outer in OUTERS)
FULL = (1 << 14) - 1
F_VECTOR = (14, 87, 318, 762, 1254, 1428, 1092, 513, 114)


def _mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


def _minimal_transversals(hyperedges: tuple[int, ...]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, 1 << 14):
        if not all(candidate & edge for edge in hyperedges):
            continue
        removable = False
        for vertex in range(14):
            bit = 1 << vertex
            if candidate & bit and all((candidate ^ bit) & edge for edge in hyperedges):
                removable = True
                break
        if not removable:
            answers.append(candidate)
    return tuple(answers)


def _packed_face_table(face_table: tuple[bool, ...]) -> bytes:
    return bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )


def _placement_buckets(
    support: tuple[int, ...],
) -> dict[tuple[int, ...], tuple[dict[tuple[int, int], int], ...]]:
    support_pairs = tuple(combinations(support, 2))
    core_pairs = tuple(combinations(CORE, 2))
    buckets: defaultdict[tuple[int, ...], list[dict[tuple[int, int], int]]] = (
        defaultdict(list)
    )
    for code in lanner_placement_codes(5):
        labels = dict(zip(support_pairs, code, strict=True))
        buckets[tuple(labels[edge] for edge in core_pairs)].append(labels)
    return {key: tuple(value) for key, value in buckets.items()}


def _common_core_replay(
    vertices: tuple[tuple[int, ...], ...],
) -> dict[str, object]:
    buckets = tuple(_placement_buckets(support) for support in SUPPORTS)
    core_codes = set(buckets[0])
    if any(set(table) != core_codes for table in buckets[1:]):
        raise AssertionError("independent common-core extension keys disagree")
    extension_histogram = Counter(len(options) for options in buckets[0].values())
    join_counts = (
        420,
        sum(len(buckets[0][key]) * len(buckets[1][key]) for key in core_codes),
        sum(
            len(buckets[0][key])
            * len(buckets[1][key])
            * len(buckets[2][key])
            for key in core_codes
        ),
    )
    union = tuple(sorted((*CORE, *OUTERS)))
    union_set = set(union)
    restrictions = tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in union_set)
                for vertex in vertices
            },
            key=lambda item: (len(item), item),
        )
    )
    states = []
    template = {edge: 2 for edge in combinations(union, 2)}
    for key in sorted(core_codes):
        for extensions in product(*(table[key] for table in buckets)):
            labels = dict(template)
            for extension in extensions:
                for edge, label in extension.items():
                    old = labels[edge]
                    if edge[0] in CORE and edge[1] in CORE and old not in (2, label):
                        raise AssertionError("independent common-core labels disagree")
                    labels[edge] = label
            states.append(labels)
    if len(states) != join_counts[-1]:
        raise AssertionError("independent common-core join count changed")

    staged_counts = []
    survivors = tuple(states)
    for order in sorted(set(map(len, restrictions))):
        tests = tuple(face for face in restrictions if len(face) == order)
        survivors = tuple(
            labels
            for labels in survivors
            if all(is_finite_coxeter_labeling(face, labels) for face in tests)
        )
        staged_counts.append((order, len(tests), len(survivors)))
    return {
        "join_counts": list(join_counts),
        "core_extension_count_histogram": {
            str(size): count for size, count in sorted(extension_histogram.items())
        },
        "restriction_size_histogram": {
            str(size): count for size, count in sorted(Counter(map(len, restrictions)).items())
        },
        "staged_survivor_counts": [list(row) for row in staged_counts],
        "final_survivor_count": len(survivors),
        "unrestricted_external_label_no_go": len(survivors) == 0,
    }


def replay() -> dict[str, object]:
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    support_manifest = json.loads(SUPPORT.read_text(encoding="utf-8"))
    order, decoded_edges = decode_graph6(GRAPH6)
    if order != 14 or decoded_edges != EDGES:
        raise AssertionError("independent record-32 a=0 graph decode failed")

    hyperedges = tuple(_mask(face) for face in (*EDGES, *FIVES))
    blocker = _minimal_transversals(hyperedges)
    if Counter(value.bit_count() for value in blocker) != {5: 114}:
        raise AssertionError("independent a=0 blocker is not 114 pure five-sets")
    if blocker != tuple(archive["incidence"]["facet_complement_masks"]):
        raise AssertionError("independent a=0 blocker differs from archive")

    records = [
        record
        for record in support_manifest["records"]
        if record["survives_second_layer"]
    ]
    source = records[32]
    if source["graph6"] != GRAPH6:
        raise AssertionError("independent a=0 source-index replay failed")
    global_audit = audit_global_facet_selection(
        GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        blocker,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("independent a=0 global audit failed")
    if global_audit.face_vector != F_VECTOR:
        raise AssertionError("independent a=0 f-vector changed")

    face_table = generated_face_table(blocker)
    packed_digest = sha256(_packed_face_table(face_table)).hexdigest()
    if packed_digest != archive["incidence"]["packed_face_table_sha256"]:
        raise AssertionError("independent a=0 face-table digest mismatch")
    vertices = tuple(
        tuple(facet for facet in range(14) if (FULL ^ complement) & (1 << facet))
        for complement in blocker
    )
    common_core = _common_core_replay(vertices)
    if common_core["join_counts"] != [420, 1332, 4620]:
        raise AssertionError("independent a=0 join census changed")
    if common_core["staged_survivor_counts"] != [
        [4, 13, 4620],
        [5, 12, 96],
        [6, 4, 0],
    ]:
        raise AssertionError("independent a=0 staged census changed")
    archived_gate = archive["common_four_core_coxeter_gate"]
    if archived_gate["join_counts"] != common_core["join_counts"]:
        raise AssertionError("independent a=0 joins differ from archive")
    archived_stages = [
        [
            stage["order"],
            stage["restriction_count"],
            stage["cumulative_passing_state_count"],
        ]
        for stage in archived_gate["cumulative_restriction_stages"]
    ]
    if archived_stages != common_core["staged_survivor_counts"]:
        raise AssertionError("independent a=0 stages differ from archive")
    if not archived_gate["unrestricted_external_label_no_go"]:
        raise AssertionError("archived a=0 monotone no-go is absent")

    return {
        "minimal_transversal_size_histogram": {"5": 114},
        "face_vector": list(global_audit.face_vector),
        "ridge_degree_histogram": {
            str(key): value for key, value in global_audit.ridge_degree_histogram.items()
        },
        "dual_component_count": global_audit.dual_component_count,
        "all_links_mod2_homology_spheres": (
            global_audit.all_links_are_mod2_homology_spheres
        ),
        "packed_face_table_sha256": packed_digest,
        "common_four_core": common_core,
        "ordinary_polytopality_claimed": False,
    }


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
