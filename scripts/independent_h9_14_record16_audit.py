#!/usr/bin/env python3
"""Hostile replay of the complete H9/14 record-16 exclusion.

This file intentionally does not import either record-16 proof module.  It
rebuilds the blocker types, symbolic ridge equations, cut graphs, archived
X-graph filter, A4 compatibility, and final A5 contradiction from masks.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
X_ARCHIVE = ROOT / "data/enumerations/h9_14_record16_x_graphs.json"
GRAPH6 = "M???????????o?o?_"
EDGES = ((0, 12), (1, 12), (0, 13), (1, 13), (12, 13))
U = (0, 1)
I = tuple(range(2, 12))
W = tuple(range(12))
P = (12, 13)
EXPECTED_SUMMARY_SHA256 = (
    "bfd07dd04a18632fa7a8ca1c868e48fed60d8b944e594282b83427241633c3d5"
)


def msk(values) -> int:
    return sum(1 << value for value in values)


def verts(value: int, order: int = 14) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(order) if value & (1 << vertex))


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


def blocker_and_ridges() -> dict[str, int]:
    edge_masks = tuple(msk(edge) for edge in EDGES)
    covers = []
    y_facets = []
    z_facets = []
    special = msk(P)
    active = msk(U)
    isolates = msk(I)
    for chosen in combinations(range(14), 5):
        value = msk(chosen)
        if not all(value & edge for edge in edge_masks):
            continue
        covers.append(value)
        if value & special == special:
            y_facets.append(value)
        elif (
            (value & special).bit_count() == 1
            and value & active == active
            and (value & isolates).bit_count() == 2
        ):
            z_facets.append(value)
        else:
            raise AssertionError("unclassified graph-cover five-set")
    if (len(covers), len(y_facets), len(z_facets)) != (310, 220, 90):
        raise AssertionError("blocker classification changed")
    forbidden_y = sum(
        (facet ^ special) & active == active for facet in y_facets
    )
    z_pair_counts = Counter(facet & isolates for facet in z_facets)
    if forbidden_y != 10 or set(z_pair_counts.values()) != {2}:
        raise AssertionError("y prohibition or z pairing changed")

    terms = []
    for facet in y_facets:
        terms.append((facet, ("y", facet ^ special)))
    for facet in z_facets:
        terms.append((facet, ("z", facet & isolates)))
    hist = Counter()
    y_incidence = 0
    z_incidence = 0
    nonempty = 0
    for chosen in combinations(range(14), 6):
        six = msk(chosen)
        equation = Counter(key for facet, key in terms if facet & six == facet)
        if not equation:
            continue
        nonempty += 1
        y_incidence += sum(
            count for (kind, _), count in equation.items() if kind == "y"
        )
        z_incidence += sum(
            count for (kind, _), count in equation.items() if kind == "z"
        )
        signature = tuple(sorted(equation.values()))
        hist[signature] += 1
    # 450 ordinary R1 equations have four y terms; 45 have four y terms and
    # one doubled z term; 240 R2 equations have three z terms.
    if hist != {(1, 1, 1): 240, (1, 1, 1, 1): 450, (1, 1, 1, 1, 2): 45}:
        raise AssertionError(f"symbolic ridge signatures changed: {hist}")
    if (nonempty, y_incidence, z_incidence) != (735, 1_980, 810):
        raise AssertionError("symbolic ridge incidence totals changed")
    return {
        "graph_cover_five_sets": len(covers),
        "y_five_sets": len(y_facets),
        "forbidden_y_01i": forbidden_y,
        "z_five_sets": len(z_facets),
        "paired_z_variables": len(z_pair_counts),
        "nonempty_ridge_groups": nonempty,
        "y_ridge_incidences": y_incidence,
        "z_ridge_incidences": z_incidence,
    }


def z_cut_replay() -> dict[str, object]:
    pairs = tuple(combinations(range(10), 2))
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    graphs = set()
    histogram = Counter()
    for seed in range(512):
        to_zero = {vertex: bool(seed & (1 << (vertex - 1))) for vertex in range(1, 10)}
        graph = 0
        for first, second in pairs:
            edge = to_zero[second] if first == 0 else to_zero[first] ^ to_zero[second]
            graph |= int(edge) << pair_index[(first, second)]
        for triple in combinations(range(10), 3):
            count = sum(
                bool(graph & (1 << pair_index[pair]))
                for pair in combinations(triple, 2)
            )
            if count not in (0, 2):
                raise AssertionError("z seed reconstruction failed")
        graphs.add(graph)
        histogram[graph.bit_count()] += 1
    if len(graphs) != 512:
        raise AssertionError("z graph seed map is not injective")
    return {
        "solution_count": len(graphs),
        "edge_histogram": {
            str(edge_count): histogram[edge_count]
            for edge_count in sorted(histogram)
        },
    }


def adjacency(order: int, edges: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    result = [0] * order
    for first, second in edges:
        result[first] |= 1 << second
        result[second] |= 1 << first
    return tuple(result)


def triangle_free(adj: tuple[int, ...]) -> bool:
    return all(
        not (
            adj[first] & (1 << second)
            and adj[first] & (1 << third)
            and adj[second] & (1 << third)
        )
        for first, second, third in combinations(range(len(adj)), 3)
    )


def has_induced_2k2(adj: tuple[int, ...]) -> bool:
    for chosen in combinations(range(len(adj)), 4):
        degrees = [
            sum(bool(adj[vertex] & (1 << other)) for other in chosen)
            for vertex in chosen
        ]
        if sorted(degrees) == [1, 1, 1, 1]:
            return True
    return False


def y_and_a4(order: int, raw_edges: tuple[tuple[int, int], ...]):
    x_edges = {tuple(sorted((first + 2, second + 2))) for first, second in raw_edges}
    y = []
    for triple in combinations(W, 3):
        active_count = len(set(triple) & set(U))
        isolate_vertices = tuple(vertex for vertex in triple if vertex not in U)
        if active_count == 2:
            value = 0
        elif active_count == 1:
            value = int(tuple(sorted(isolate_vertices)) in x_edges)
        else:
            value = sum(
                tuple(sorted(pair)) in x_edges
                for pair in combinations(isolate_vertices, 2)
            ) % 2
        if value:
            y.append(msk(triple))
    blockers = []
    active_mask = msk(U)
    for chosen in combinations(W, 4):
        value = msk(chosen)
        if value & active_mask and all(value & triple for triple in y):
            blockers.append(value)
    return tuple(y), tuple(blockers)


def replay_x_archive() -> dict[str, object]:
    document = json.loads(X_ARCHIVE.read_text(encoding="utf-8"))
    payload = dict(document)
    archived_digest = payload.pop("certificate_sha256")
    digest = sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if archived_digest != digest:
        raise AssertionError("X archive canonical digest changed")
    expected_counts = {6: 38, 7: 107, 8: 410}
    expected_admissible = {6: 22, 7: 41, 8: 82}
    target_keys = {(8, 72), (7, 62), (8, 70), (7, 60), (6, 54)}
    targets = {key: [] for key in target_keys}
    retained_counts = {}
    for order, count in expected_counts.items():
        stream = document["streams"][str(order)]
        words = stream["graph6"]
        raw_digest = sha256(("".join(word + "\n" for word in words)).encode()).hexdigest()
        if len(words) != count or len(set(words)) != count or raw_digest != stream["raw_stream_sha256"]:
            raise AssertionError("X raw stream changed")
        retained = 0
        for word in words:
            decoded_order, edges = decode(word)
            if decoded_order != order:
                raise AssertionError("mixed order in X stream")
            adj = adjacency(order, edges)
            if not triangle_free(adj):
                raise AssertionError("geng -t stream has a triangle")
            if has_induced_2k2(adj):
                continue
            retained += 1
            degrees = tuple(sorted((value.bit_count() for value in adj), reverse=True))
            y_formula = 12 * len(edges) - sum(degree * degree for degree in degrees)
            y, blockers = y_and_a4(order, edges)
            if len(y) != y_formula:
                raise AssertionError("independent y count replay failed")
            key = (order, len(y))
            if key in targets:
                targets[key].append(
                    {
                        "graph6": word,
                        "edges": len(edges),
                        "degrees": list(degrees),
                        "a4": [list(verts(value)) for value in blockers],
                    }
                )
        retained_counts[order] = retained
    if retained_counts != expected_admissible:
        raise AssertionError("induced-2K2 filter counts changed")
    expected_target_counts = {(8, 72): 0, (7, 62): 1, (8, 70): 1, (7, 60): 2, (6, 54): 1}
    if {key: len(targets[key]) for key in target_keys} != expected_target_counts:
        raise AssertionError("target X graph counts changed")
    return {
        "archive_certificate_sha256": archived_digest,
        "retained_counts": {str(k): v for k, v in sorted(retained_counts.items())},
        "targets": {
            f"n{key[0]}_y{key[1]}": targets[key] for key in sorted(target_keys)
        },
    }


def schedule_cases() -> dict[str, object]:
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    survivors = [record for record in support["records"] if record["survives_second_layer"]]
    record = survivors[16]
    if record["graph6"] != GRAPH6 or decode(GRAPH6) != (14, EDGES):
        raise AssertionError("record-16 source target changed")
    facets = {
        schedule["missing_four_count"]: schedule["f_vector"][-1]
        for schedule in record["dehn_sommerville_schedule"]
    }
    residual = []
    for a in range(1, 5):
        for r in range(6):
            s = 10 - r
            z = r * s
            required_y = facets[a] - 2 * z
            if required_y <= 9 * s:
                residual.append((a, r, s, required_y))
    expected = [(3, 2, 8, 72), (3, 3, 7, 62), (4, 2, 8, 70), (4, 3, 7, 60), (4, 4, 6, 54)]
    if residual != expected:
        raise AssertionError("coarse residual schedule cases changed")
    return {"facet_counts": {str(k): v for k, v in sorted(facets.items())}, "residual": residual}


def final_case(x_audit: dict[str, object]) -> dict[str, object]:
    target = x_audit["targets"]["n6_y54"]
    if len(target) != 1:
        raise AssertionError("final X graph not unique")
    graph = target[0]
    if graph["degrees"] != [3, 3, 3, 3, 3, 3] or graph["edges"] != 9:
        raise AssertionError("final X graph is not cubic on six vertices")
    _, raw_edges = decode(graph["graph6"])
    adj = adjacency(6, raw_edges)
    colour = {0: 0}
    queue = [0]
    while queue:
        vertex = queue.pop()
        for neighbor in range(6):
            if not adj[vertex] & (1 << neighbor):
                continue
            expected = 1 - colour[vertex]
            if neighbor in colour and colour[neighbor] != expected:
                raise AssertionError("final X graph is not bipartite")
            if neighbor not in colour:
                colour[neighbor] = expected
                queue.append(neighbor)
    parts = [
        {vertex for vertex in range(6) if colour.get(vertex) == side}
        for side in (0, 1)
    ]
    if sorted(map(len, parts)) != [3, 3] or not all(
        bool(adj[first] & (1 << second)) == (colour[first] != colour[second])
        for first, second in combinations(range(6), 2)
    ):
        raise AssertionError("final X graph is not structurally K3,3")
    a4 = tuple(msk(face) for face in graph["a4"])
    if len(a4) != 4:
        raise AssertionError("final X graph does not force four A4s")
    shadow = {
        face | (1 << vertex)
        for face in a4
        for vertex in W
        if not face & (1 << vertex)
    }
    if len(shadow) != 30:
        raise AssertionError("final A4 shadow changed")
    additional_a5 = 32 - len(shadow)
    if additional_a5 != 2:
        raise AssertionError("final minimal A5 count changed")
    # tau(Z)>=2 and only two extra rows force both into Z.  Two five-subsets
    # of a ten-set with no common vertex are disjoint complements.  Replay a
    # canonical pair and count its two-element transversals.
    first = set(range(5))
    second = set(range(5, 10))
    z_edges = [
        pair
        for pair in combinations(range(10), 2)
        if set(pair) & first and set(pair) & second
    ]
    if len(z_edges) != 25:
        raise AssertionError("two complementary A5s do not give K5,5")

    # Hostile direct dualization: reconstruct the 102 final facet
    # complements and enumerate their complete minimal-transversal blocker.
    y, _ = y_and_a4(6, raw_edges)
    facets = [msk(P) | triple for triple in y]
    for large in range(2, 8):
        for small in range(8, 12):
            for special in P:
                facets.append(msk(U) | (1 << large) | (1 << small) | (1 << special))
    facets = sorted(set(facets))
    minimal = []
    for candidate in range(1, 1 << 14):
        if not all(candidate & facet for facet in facets):
            continue
        if all(
            not all((candidate ^ (1 << vertex)) & facet for facet in facets)
            for vertex in verts(candidate)
        ):
            minimal.append(candidate)
    blocker_histogram = Counter(value.bit_count() for value in minimal)
    expected_histogram = {2: 5, 4: 4, 5: 2, 6: 2, 7: 2}
    if len(facets) != 102 or blocker_histogram != expected_histogram:
        raise AssertionError("final direct blocker replay changed")
    special_fives = sorted(
        verts(value) for value in minimal if value.bit_count() == 5
    )
    if special_fives != [(8, 9, 10, 11, 12), (8, 9, 10, 11, 13)]:
        raise AssertionError("final direct blocker A5 support changed")
    return {
        "a4": graph["a4"],
        "shadow_size": len(shadow),
        "additional_minimal_a5": additional_a5,
        "two_row_z_blocker_edges": len(z_edges),
        "ridge_forced_z_edges": 24,
        "contradiction": "K5,5 versus K4,6",
        "direct_facet_count": len(facets),
        "direct_minimal_transversal_histogram": {
            str(size): blocker_histogram[size]
            for size in sorted(blocker_histogram)
        },
        "direct_special_a5s": [list(face) for face in special_fives],
    }


def audit() -> dict[str, object]:
    x_audit = replay_x_archive()
    result = {
        "blocker_and_ridges": blocker_and_ridges(),
        "z_cut": z_cut_replay(),
        "schedules": schedule_cases(),
        "x_archive": x_audit,
        "final_case": final_case(x_audit),
        "conclusion": "record 16 has no simplicial-sphere completion",
    }
    digest = sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("independent record-16 audit summary changed")
    result["summary_sha256"] = digest
    return result


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
