#!/usr/bin/env python3
"""Hostile replay of the H9/14 record-18 abstract candidate certificate.

This file intentionally imports neither the record-18 construction module nor
its certificate driver.  It reconstructs the blocker from hard-coded sets,
replays the graph-cover/ridge/cut structure, and invokes only generic face and
Gale auditors on the resulting witness.
"""

from __future__ import annotations

from collections import Counter, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    generated_face_table,
)
from coxeter_polytopes.combinatorics.gale_chirotope import (
    audit_uniform_rank_four_gale,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record18_selection.json"
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
GRAPH6 = "M???????????o?K?_"
EDGES = ((0, 12), (1, 12), (2, 13), (3, 13), (12, 13))
FIVES = (
    (0, 2, 4, 5, 6),
    (0, 5, 6, 7, 8),
    (0, 3, 7, 8, 9),
    (3, 7, 8, 9, 10),
    (2, 4, 5, 6, 11),
    (1, 2, 4, 10, 11),
    (1, 3, 9, 10, 11),
    (1, 4, 9, 10, 11),
)
EXPECTED_F_VECTOR = (14, 86, 310, 735, 1204, 1372, 1052, 495, 110)
EXPECTED_FACET_SHA256 = "969da679ad87ca730e822ceb7c5a4eee7f14b461062cf71543293c50c9872afc"
EXPECTED_FACE_SHA256 = "620c83d073c4a141ef2ec5cf9333326cdc8d7c7e94d5e48e647b19b030943eaf"
EXPECTED_SUMMARY_SHA256 = "3eb0f00ca9c98144c5e3fea65b00ef6e149a493ede96c529282d3b9ca6c69bb4"


def msk(values) -> int:
    return sum(1 << value for value in values)


def verts(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if value & (1 << vertex))


def decode(word: str) -> tuple[int, tuple[tuple[int, int], ...]]:
    codes = [ord(character) - 63 for character in word]
    if not codes or codes[0] == 63 or any(code < 0 or code > 63 for code in codes):
        raise ValueError("unsupported graph6 word")
    order = codes[0]
    bit_count = order * (order - 1) // 2
    expected_payload = (bit_count + 5) // 6
    if len(codes) != expected_payload + 1:
        raise ValueError("wrong graph6 payload length")
    bits = []
    for code in codes[1:]:
        bits.extend((code >> shift) & 1 for shift in range(5, -1, -1))
    if any(bits[bit_count:]):
        raise ValueError("nonzero graph6 padding")
    edges = []
    position = 0
    for upper in range(1, order):
        for lower in range(upper):
            if bits[position]:
                edges.append((lower, upper))
            position += 1
    return order, tuple(edges)


def minimal_transversals(hyperedges: tuple[int, ...]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, 1 << 14):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if all(
            not all((candidate ^ (1 << vertex)) & edge for edge in hyperedges)
            for vertex in verts(candidate)
        ):
            answers.append(candidate)
    return tuple(answers)


def cut_parts(
    domain: tuple[int, ...], edges: tuple[tuple[int, int], ...]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    adjacency = {vertex: set() for vertex in domain}
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    colour = {domain[0]: 0}
    queue = deque([domain[0]])
    while queue:
        vertex = queue.popleft()
        for neighbor in adjacency[vertex]:
            expected = 1 - colour[vertex]
            if neighbor in colour and colour[neighbor] != expected:
                raise AssertionError("hostile cut replay found an odd cycle")
            if neighbor not in colour:
                colour[neighbor] = expected
                queue.append(neighbor)
    if len(colour) != 10:
        raise AssertionError("hostile cut replay found a disconnected graph")
    parts = tuple(
        tuple(sorted(vertex for vertex in domain if colour[vertex] == side))
        for side in (0, 1)
    )
    parts = tuple(sorted(parts))
    expected = {
        tuple(sorted((first, second))) for first in parts[0] for second in parts[1]
    }
    if tuple(map(len, parts)) != (5, 5) or set(edges) != expected:
        raise AssertionError("hostile cut replay is not K5,5")
    return parts  # type: ignore[return-value]


def replay_structure(blocker: tuple[int, ...]) -> dict[str, object]:
    edge_masks = tuple(msk(edge) for edge in EDGES)
    graph_blocker = minimal_transversals(edge_masks)
    if tuple(map(verts, graph_blocker)) != ((2, 3, 12), (0, 1, 13), (12, 13)):
        raise AssertionError("hostile graph blocker changed")

    potential = {"y": [], "z": [], "w": []}
    special = msk((12, 13))
    for chosen in combinations(range(14), 5):
        candidate = msk(chosen)
        if not all(candidate & edge for edge in edge_masks):
            continue
        if candidate & special == special:
            potential["y"].append(candidate)
        elif candidate & msk((0, 1, 13)) == msk((0, 1, 13)) and not candidate & (1 << 12):
            potential["z"].append(candidate)
        elif candidate & msk((2, 3, 12)) == msk((2, 3, 12)) and not candidate & (1 << 13):
            potential["w"].append(candidate)
        else:
            raise AssertionError("hostile graph-cover classification failed")
    if tuple(len(potential[kind]) for kind in ("y", "z", "w")) != (220, 45, 45):
        raise AssertionError("hostile potential-facet counts changed")

    terms = tuple(
        (candidate, kind)
        for kind in ("y", "z", "w")
        for candidate in potential[kind]
    )
    groups = Counter()
    incidences = Counter()
    for chosen in combinations(range(14), 6):
        six_set = msk(chosen)
        kinds = tuple(kind for candidate, kind in terms if candidate & six_set == candidate)
        if not kinds:
            continue
        incidences.update(kinds)
        if six_set & special == special:
            groups["two_special"] += 1
        elif six_set & msk((0, 1, 13)) == msk((0, 1, 13)):
            groups["z"] += 1
        elif six_set & msk((2, 3, 12)) == msk((2, 3, 12)):
            groups["w"] += 1
        else:
            raise AssertionError("hostile symbolic-ridge classification failed")
    if groups != {"two_special": 495, "z": 120, "w": 120}:
        raise AssertionError("hostile symbolic-ridge group counts changed")
    if incidences != {"y": 1980, "z": 405, "w": 405}:
        raise AssertionError("hostile symbolic-ridge incidences changed")

    selected = Counter()
    z_pairs = []
    w_pairs = []
    for candidate in blocker:
        if candidate in potential["y"]:
            selected["y"] += 1
        elif candidate in potential["z"]:
            selected["z"] += 1
            z_pairs.append(tuple(v for v in verts(candidate) if v not in (0, 1, 13)))
        elif candidate in potential["w"]:
            selected["w"] += 1
            w_pairs.append(tuple(v for v in verts(candidate) if v not in (2, 3, 12)))
        else:
            raise AssertionError("hostile selected-facet classification failed")
    if selected != {"y": 60, "z": 25, "w": 25}:
        raise AssertionError("hostile selected-facet counts changed")
    z_parts = cut_parts(tuple(range(2, 12)), tuple(z_pairs))
    w_parts = cut_parts((0, 1, *range(4, 12)), tuple(w_pairs))
    if z_parts != ((2, 4, 5, 6, 11), (3, 7, 8, 9, 10)):
        raise AssertionError("hostile z cut changed")
    if w_parts != ((0, 5, 6, 7, 8), (1, 4, 9, 10, 11)):
        raise AssertionError("hostile w cut changed")
    return {
        "minimal_graph_covers": [list(verts(value)) for value in graph_blocker],
        "potential_type_counts": {kind: len(potential[kind]) for kind in ("y", "z", "w")},
        "symbolic_ridge_group_counts": dict(groups),
        "symbolic_ridge_incidence_counts": dict(incidences),
        "selected_type_counts": dict(selected),
        "z_cut_parts": [list(part) for part in z_parts],
        "w_cut_parts": [list(part) for part in w_parts],
    }


def replay() -> dict[str, object]:
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    archive_payload = dict(archive)
    archived_certificate_digest = archive_payload.pop("certificate_sha256")
    observed_certificate_digest = sha256(
        json.dumps(archive_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if observed_certificate_digest != archived_certificate_digest:
        raise AssertionError("hostile archive digest replay failed")
    if archive.get("status") != "abstract_audited_nerve_candidate":
        raise AssertionError("hostile replay rejected an overclaimed status")

    if decode(GRAPH6) != (14, EDGES):
        raise AssertionError("hostile graph6 decode failed")
    clutter = tuple(msk(edge) for edge in EDGES) + tuple(msk(face) for face in FIVES)
    blocker = minimal_transversals(clutter)
    if Counter(value.bit_count() for value in blocker) != {5: 110}:
        raise AssertionError("hostile blocker replay is not pure with 110 members")
    archived_facets = tuple(archive["incidence"]["facet_complement_masks"])
    if blocker != archived_facets:
        raise AssertionError("hostile blocker differs from archived facets")
    facet_digest = sha256(
        b"".join(value.to_bytes(2, "little") for value in blocker)
    ).hexdigest()
    if facet_digest != EXPECTED_FACET_SHA256:
        raise AssertionError("hostile facet digest changed")

    structure = replay_structure(blocker)
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    survivors = [record for record in support["records"] if record["survives_second_layer"]]
    source = survivors[18]
    if source["graph6"] != GRAPH6:
        raise AssertionError("hostile source-index replay failed")
    global_audit = audit_global_facet_selection(
        GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        blocker,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("hostile global record-18 audit failed")
    if global_audit.face_vector != EXPECTED_F_VECTOR:
        raise AssertionError("hostile record-18 face vector changed")

    signs = tuple(archive["gale_chirotope_signs"])
    gale_audit = audit_uniform_rank_four_gale(blocker, signs)
    if not gale_audit.passes:
        raise AssertionError("hostile record-18 Gale audit failed")
    face_table = generated_face_table(blocker)
    packed = bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )
    packed_digest = sha256(packed).hexdigest()
    if packed_digest != EXPECTED_FACE_SHA256:
        raise AssertionError("hostile face-table digest changed")

    summary: dict[str, object] = {
        "status": "abstract_audited_nerve_candidate",
        "minimal_transversal_size_histogram": {"5": 110},
        "facet_complement_sha256": facet_digest,
        "structure": structure,
        "face_vector": list(global_audit.face_vector),
        "ridge_degree_histogram": {
            str(key): value for key, value in global_audit.ridge_degree_histogram.items()
        },
        "dual_component_count": global_audit.dual_component_count,
        "all_links_mod2_homology_spheres": global_audit.all_links_are_mod2_homology_spheres,
        "packed_face_table_sha256": packed_digest,
        "gale_sign_count": gale_audit.sign_count,
        "positive_circuit_count": len(gale_audit.positive_circuit_masks),
        "gale_audit_passes": gale_audit.passes,
        "certificate_sha256": archived_certificate_digest,
    }
    # Exclude the enclosing archive digest to avoid a source/archive hash
    # cycle: this hostile checker itself is pinned by the archive.
    summary_payload = dict(summary)
    summary_payload.pop("certificate_sha256")
    digest = sha256(
        json.dumps(summary_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("hostile record-18 summary digest is not frozen")
    summary["summary_sha256"] = digest
    return summary


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
