#!/usr/bin/env python3
"""Hostile replay of the exact H9/14 record-17 whole-shadow residual.

This script intentionally does not import ``h9_14_record17``.  It decodes
the archived graph6 quotient, reconstructs all facet complements directly,
dualizes every residual blocker, quotients by the actual coloured graph
automorphisms, and audits every link of all eight orbit representatives over
``GF(2)``.
"""

from __future__ import annotations

from collections import Counter, deque
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.simplicial_homology import (
    first_mod2_homology_sphere_link_failure,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
X_ARCHIVE = ROOT / "data/enumerations/h9_14_record17_x_graphs.json"

GRAPH6 = "M???????????o?g?_"
EDGES = ((0, 12), (1, 12), (0, 13), (2, 13), (12, 13))
EDGE_MASKS = tuple(sum(1 << value for value in edge) for edge in EDGES)
I = tuple(range(3, 12))
L = tuple(range(9))
W = tuple(range(12))
FULL = (1 << 14) - 1
P_BASE = (0, 1, 13)
Q_BASE = (0, 2, 12)
P_GROUND = (2,) + I
Q_GROUND = (1,) + I

EXPECTED_STAGES = {
    "target_y_count": {0: 1434, 1: 1376, 2: 1455, 3: 1869, 4: 521},
    "independent_triples": {0: 654, 1: 584, 2: 159, 3: 136, 4: 21},
    "exact_a4": {0: 36, 1: 82, 2: 76, 3: 53, 4: 19},
    "exact_a5_and_support": {0: 36, 1: 82, 2: 76, 3: 53, 4: 19},
    "full_blocker": {0: 1, 1: 17, 2: 3, 3: 0, 4: 0},
    "full_face_vector": {0: 1, 1: 17, 2: 3, 3: 0, 4: 0},
}
EXPECTED_ORBITS = (
    (1, "H???FB_", 31, 6, "719e1b41557b6f6d4934cc86f9e199ecfefa9469fb133cedb6f5e852a7d35514"),
    (1, "H???FB_", 440, 4, "72b5451251f7649dd45bacb9aba4c16a35979f9cfb1071f684b0cc4b240a51e7"),
    (2, "H???F?^", 248, 1, "13e150f8162107a6272332090ac3d85109340010a9e83e0225a8164ec8dec9ae"),
    (2, "H???Fbo", 496, 1, "bb878935209777eee9693a6ccf09822015a6aee96699b2ca65897a7b36f1f282"),
    (0, "H??EEBo", 496, 1, "3580d27d6ea97728df071e51f92da84fe8379bb8255214ba59c82dd63174ef7f"),
    (1, "H??EEBw", 31, 1, "7f077b63f89695f8689ff139007d94893bae9634c8fccb16bbc2efc02dc80071"),
    (1, "H??FFB_", 31, 6, "267972a560ae95971013d65f09e9f1de64a15df996798c06ae2f384555567cb6"),
    (2, "H??FF?^", 248, 1, "d764cd7e03333180fff96961bafa7a10b74665e632a39bb03ab62be759848d49"),
)
EXPECTED_SUMMARY_SHA256 = (
    "d30dfc3d25b5e4bb0c5067f96fabfa863f843167d02774596b1efffd26adb06d"
)


def msk(values) -> int:
    result = 0
    for value in values:
        result |= 1 << value
    return result


def verts(value: int, order: int = 14) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(order) if value & (1 << vertex))


def decode(word: str) -> tuple[int, tuple[tuple[int, int], ...]]:
    codes = [ord(character) - 63 for character in word]
    if not codes or codes[0] == 63 or any(code < 0 or code > 63 for code in codes):
        raise ValueError("unsupported graph6 word")
    order = codes[0]
    bit_count = order * (order - 1) // 2
    if len(codes) != (bit_count + 5) // 6 + 1:
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


def adjacency(order: int, edges) -> tuple[int, ...]:
    result = [0] * order
    for first, second in edges:
        result[first] |= 1 << second
        result[second] |= 1 << first
    return tuple(result)


def triangle_free(adj: tuple[int, ...]) -> bool:
    return all(
        not all(bool(adj[a] & (1 << b)) for a, b in combinations(triple, 2))
        for triple in combinations(range(len(adj)), 3)
    )


def induced_2k2_free(adj: tuple[int, ...]) -> bool:
    for chosen in combinations(range(len(adj)), 4):
        chosen_mask = msk(chosen)
        if sorted((adj[v] & chosen_mask).bit_count() for v in chosen) == [1, 1, 1, 1]:
            return False
    return True


def colour_isomorphic(
    adj: tuple[int, ...], p1: int, q1: int, p2: int, q2: int
) -> bool:
    order = len(adj)
    degree = tuple(row.bit_count() for row in adj)

    def colours(p: int, q: int) -> tuple[int, ...]:
        return tuple(int(bool(p & 1 << v)) + 2 * int(bool(q & 1 << v)) for v in range(order))

    first_colours = colours(p1, q1)

    def attempt(second_colours: tuple[int, ...]) -> bool:
        def sig(c, v):
            return (
                c[v],
                degree[v],
                tuple(sorted((c[u], degree[u]) for u in range(order) if adj[v] & 1 << u)),
            )

        first_sig = tuple(sig(first_colours, v) for v in range(order))
        second_sig = tuple(sig(second_colours, v) for v in range(order))
        if sorted(first_sig) != sorted(second_sig):
            return False
        mapping = {}
        used = set()

        def search() -> bool:
            if len(mapping) == order:
                return True
            choices = []
            for vertex in range(order):
                if vertex in mapping:
                    continue
                images = [
                    image
                    for image in range(order)
                    if image not in used
                    and first_sig[vertex] == second_sig[image]
                    and all(
                        bool(adj[vertex] & 1 << other)
                        == bool(adj[image] & 1 << mapped)
                        for other, mapped in mapping.items()
                    )
                ]
                if not images:
                    return False
                choices.append((len(images), vertex, images))
            _, vertex, images = min(choices)
            for image in images:
                mapping[vertex] = image
                used.add(image)
                if search():
                    return True
                used.remove(image)
                del mapping[vertex]
            return False

        return search()

    return attempt(colours(p2, q2)) or attempt(colours(q2, p2))


def archive_replay() -> tuple[dict[str, tuple[int, ...]], dict[str, object]]:
    document = json.loads(X_ARCHIVE.read_text(encoding="utf-8"))
    payload = dict(document)
    archived_digest = payload.pop("certificate_sha256")
    digest = sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    words = document["graph6"]
    raw = sha256(("".join(word + "\n" for word in words)).encode()).hexdigest()
    if (
        digest != archived_digest
        or raw != document["raw_stream_sha256"]
        or len(words) != 1897
        or len(set(words)) != 1897
    ):
        raise AssertionError("record-17 X archive changed")
    admissible = {}
    all_triangle_free = True
    for word in words:
        order, raw_edges = decode(word)
        adj = adjacency(order, raw_edges)
        all_triangle_free &= order == 9 and triangle_free(adj)
        if induced_2k2_free(adj):
            admissible[word] = adj
    if not all_triangle_free or len(admissible) != 156:
        raise AssertionError("record-17 X property replay changed")
    return admissible, {
        "certificate_sha256": archived_digest,
        "raw_stream_sha256": raw,
        "stream_count": len(words),
        "induced_2k2_free_count": len(admissible),
    }


def blocker_and_ridges() -> dict[str, object]:
    covers = []
    terms = []
    type_count = Counter()
    special = msk((12, 13))
    p_base = msk(P_BASE)
    q_base = msk(Q_BASE)
    for chosen in combinations(range(14), 5):
        value = msk(chosen)
        if not all(value & edge for edge in EDGE_MASKS):
            continue
        covers.append(value)
        if value & special == special:
            kind, key = "y", value ^ special
        elif value & p_base == p_base:
            kind, key = "p", value ^ p_base
        elif value & q_base == q_base:
            kind, key = "q", value ^ q_base
        else:
            raise AssertionError("unclassified graph-cover five-set")
        type_count[kind] += 1
        terms.append((value, (kind, key)))
    if (len(covers), type_count["y"], type_count["p"], type_count["q"]) != (310, 220, 45, 45):
        raise AssertionError("record-17 facet types changed")
    signatures = Counter()
    incidence = Counter()
    for chosen in combinations(range(14), 6):
        six = msk(chosen)
        equation = Counter(key for facet, key in terms if facet & six == facet)
        if not equation:
            continue
        signatures[tuple(sorted(equation.values()))] += 1
        for (kind, _), coefficient in equation.items():
            incidence[kind] += coefficient
    expected_signatures = {
        (1, 1, 1): 240,
        (1, 1, 1, 1): 414,
        (1, 1, 1, 1, 1): 72,
        (1, 1, 1, 1, 1, 1): 9,
    }
    if signatures != expected_signatures or incidence != {"y": 1980, "p": 405, "q": 405}:
        raise AssertionError("record-17 symbolic ridge system changed")
    return {
        "graph_cover_count": len(covers),
        "type_counts": dict(sorted(type_count.items())),
        "nonempty_ridge_groups": sum(signatures.values()),
        "signature_histogram": {str(key): value for key, value in sorted(signatures.items())},
        "incidences": dict(sorted(incidence.items())),
    }


def cut_edge(first: int, second: int, cut: int) -> bool:
    return bool(cut & 1 << first) != bool(cut & 1 << second)


def switched_triangle_free(adj: tuple[int, ...], cut: int) -> bool:
    for triple in combinations(L, 3):
        if all(
            bool(adj[first] & 1 << second) != cut_edge(first, second, cut)
            for first, second in combinations(triple, 2)
        ):
            return False
    return True


def difference_rule(adj: tuple[int, ...], p: int, q: int) -> bool:
    return all(
        bool(adj[first] & 1 << second) == cut_edge(first, second, p)
        for first, second in combinations(verts(p ^ q, 9), 2)
    )


def complements(adj: tuple[int, ...], p: int, q: int) -> tuple[int, ...]:
    y = []
    for first, second in combinations(L, 2):
        pair = (first + 3, second + 3)
        edge = bool(adj[first] & 1 << second)
        if edge:
            y.append(msk((0,) + pair))
        if edge != cut_edge(first, second, p):
            y.append(msk((1,) + pair))
        if edge != cut_edge(first, second, q):
            y.append(msk((2,) + pair))
    for vertex in L:
        if bool(p & 1 << vertex) != bool(q & 1 << vertex):
            y.append(msk((1, 2, vertex + 3)))
    for triple in combinations(L, 3):
        if sum(bool(adj[a] & 1 << b) for a, b in combinations(triple, 2)) % 2:
            y.append(msk(value + 3 for value in triple))

    p_colour = {2: False} | {v + 3: bool(p & 1 << v) for v in L}
    q_colour = {1: False} | {v + 3: bool(q & 1 << v) for v in L}
    result = [msk((12, 13)) | value for value in y]
    result.extend(
        msk(P_BASE) | msk(pair)
        for pair in combinations(P_GROUND, 2)
        if p_colour[pair[0]] != p_colour[pair[1]]
    )
    result.extend(
        msk(Q_BASE) | msk(pair)
        for pair in combinations(Q_GROUND, 2)
        if q_colour[pair[0]] != q_colour[pair[1]]
    )
    return tuple(sorted(set(result)))


def independent(value: int) -> bool:
    return not any(value & edge == edge for edge in EDGE_MASKS)


def is_face(value: int, facets: tuple[int, ...]) -> bool:
    return any(not value & complement for complement in facets)


def face_table(facets: tuple[int, ...]) -> tuple[bool, ...]:
    table = [False] * (1 << 14)
    table[0] = True
    for complement in facets:
        facet = FULL ^ complement
        value = facet
        while True:
            table[value] = True
            if not value:
                break
            value = (value - 1) & facet
    return tuple(table)


def minimal_rows(table: tuple[bool, ...]) -> dict[int, tuple[int, ...]]:
    result = {}
    for value in range(1, 1 << 14):
        if table[value]:
            continue
        if all(table[value ^ (1 << vertex)] for vertex in verts(value)):
            result.setdefault(value.bit_count(), []).append(value)
    return {size: tuple(values) for size, values in sorted(result.items())}


def vector(table: tuple[bool, ...]) -> tuple[int, ...]:
    count = Counter(value.bit_count() for value, present in enumerate(table) if present)
    return tuple(count[size] for size in range(1, 10))


def ridges_two(facets: tuple[int, ...]) -> bool:
    groups = Counter()
    for complement in facets:
        outside = FULL ^ complement
        while outside:
            bit = outside & -outside
            groups[complement | bit] += 1
            outside ^= bit
    return set(groups.values()) == {2}


def dual_connected(facets: tuple[int, ...]) -> bool:
    facet_set = set(facets)
    adjacency_map = {facet: set() for facet in facets}
    groups = {}
    for facet in facets:
        outside = FULL ^ facet
        while outside:
            bit = outside & -outside
            groups.setdefault(facet | bit, []).append(facet)
            outside ^= bit
    for group in groups.values():
        if len(group) == 2:
            first, second = group
            adjacency_map[first].add(second)
            adjacency_map[second].add(first)
    seen = {facets[0]}
    queue = deque(seen)
    while queue:
        current = queue.popleft()
        for neighbor in adjacency_map[current] - seen:
            seen.add(neighbor)
            queue.append(neighbor)
    return seen == facet_set


def source_schedules() -> tuple[list[dict[str, object]], list[int]]:
    document = json.loads(SUPPORT.read_text(encoding="utf-8"))
    survivors = [record for record in document["records"] if record["survives_second_layer"]]
    record = survivors[17]
    if record["graph6"] != GRAPH6 or decode(GRAPH6) != (14, EDGES):
        raise AssertionError("record-17 source target changed")
    return record["dehn_sommerville_schedule"], record["ridge_support_core_masks"]


def census(admissible: dict[str, tuple[int, ...]]) -> tuple[dict[str, object], list[dict[str, object]]]:
    schedules, support_core = source_schedules()
    expected_vector = {row["missing_four_count"]: tuple(row["f_vector"]) for row in schedules}
    expected_five = {row["missing_four_count"]: row["missing_five_count_offset"] for row in schedules}
    cuts = tuple(msk(chosen) for chosen in combinations(L, 5))
    independent_triples = tuple(msk(chosen) for chosen in combinations(range(14), 3) if independent(msk(chosen)))
    independent_fours = tuple(msk(chosen) for chosen in combinations(range(14), 4) if independent(msk(chosen)))
    independent_fives = tuple(msk(chosen) for chosen in combinations(range(14), 5) if independent(msk(chosen)))
    work_mask = msk(W)
    stages = {name: Counter() for name in EXPECTED_STAGES}
    coupled = 0
    candidates = []
    support_set = set(support_core)
    for word, adj in admissible.items():
        good_cuts = tuple(cut for cut in cuts if switched_triangle_free(adj, cut))
        for p in good_cuts:
            for q in good_cuts:
                if not difference_rule(adj, p, q):
                    continue
                coupled += 1
                facets = complements(adj, p, q)
                if not ridges_two(facets):
                    raise AssertionError("coupled equations missed a ridge")
                y_count = len(facets) - 50
                if y_count not in (52, 54, 56, 58, 60):
                    continue
                a = (60 - y_count) // 2
                stages["target_y_count"][a] += 1
                if not all(is_face(value, facets) for value in independent_triples):
                    continue
                stages["independent_triples"][a] += 1
                a4 = tuple(value for value in independent_fours if not is_face(value, facets))
                if len(a4) != a:
                    continue
                stages["exact_a4"][a] += 1
                missing_fives = tuple(value for value in independent_fives if not is_face(value, facets))
                minimal_a5 = tuple(
                    value
                    for value in missing_fives
                    if not any(value & four == four for four in a4)
                )
                support_ok = all(
                    not four & ~work_mask
                    and (four & 1 or four & msk((1, 2)) == msk((1, 2)))
                    for four in a4
                ) and all(not five & ~work_mask for five in minimal_a5)
                if len(missing_fives) != expected_five[a] or not support_ok:
                    continue
                stages["exact_a5_and_support"][a] += 1
                table = face_table(facets)
                rows = minimal_rows(table)
                if (
                    set(rows) - {2, 4, 5}
                    or set(rows.get(2, ())) != set(EDGE_MASKS)
                    or set(rows.get(4, ())) != set(a4)
                    or set(rows.get(5, ())) != set(minimal_a5)
                ):
                    continue
                stages["full_blocker"][a] += 1
                face_vector = vector(table)
                if face_vector != expected_vector[a]:
                    continue
                stages["full_face_vector"][a] += 1
                if not set(facets) <= support_set:
                    raise AssertionError("residual facet left the archived support core")
                candidates.append(
                    {
                        "a": a,
                        "word": word,
                        "p": p,
                        "q": q,
                        "facets": facets,
                        "rows": rows,
                        "table": table,
                        "face_vector": face_vector,
                    }
                )
    normalized = {name: {a: counts[a] for a in range(5)} for name, counts in stages.items()}
    if coupled != 11545 or normalized != EXPECTED_STAGES or len(candidates) != 21:
        raise AssertionError("record-17 hostile census changed")
    return {"coupled_state_count": coupled, "stage_counts": normalized}, candidates


def orbit_and_topology(
    admissible: dict[str, tuple[int, ...]], candidates: list[dict[str, object]]
) -> dict[str, object]:
    groups = []
    for candidate in candidates:
        matching = None
        for group in groups:
            representative = group[0]
            if candidate["word"] != representative["word"]:
                continue
            if colour_isomorphic(
                admissible[candidate["word"]],
                candidate["p"],
                candidate["q"],
                representative["p"],
                representative["q"],
            ):
                matching = group
                break
        if matching is None:
            groups.append([candidate])
        else:
            matching.append(candidate)
    observed = []
    topology = []
    for group in groups:
        representative = group[0]
        facets = representative["facets"]
        digest = sha256((",".join(map(str, facets))).encode()).hexdigest()
        observed.append(
            (
                representative["a"],
                representative["word"],
                representative["p"],
                len(group),
                digest,
            )
        )
        failure = first_mod2_homology_sphere_link_failure(representative["table"], 8)
        topology.append(
            {
                "a": representative["a"],
                "x_graph6": representative["word"],
                "p_local_vertices": list(verts(representative["p"], 9)),
                "facet_count": len(facets),
                "face_vector": list(representative["face_vector"]),
                "ridge_degrees_two": ridges_two(facets),
                "dual_connected": dual_connected(facets),
                "all_links_mod2_homology_spheres": failure is None,
                "minimal_nonface_histogram": {
                    str(size): len(values)
                    for size, values in representative["rows"].items()
                },
                "facet_sha256": digest,
            }
        )
    if tuple(observed) != EXPECTED_ORBITS:
        raise AssertionError(f"record-17 orbit quotient changed: {observed}")
    if not all(
        row["ridge_degrees_two"]
        and row["dual_connected"]
        and row["all_links_mod2_homology_spheres"]
        for row in topology
    ):
        raise AssertionError("record-17 residual topology audit failed")
    return {
        "labelled_representation_count": len(candidates),
        "orbit_count": len(groups),
        "orbit_branch_histogram": {"0": 1, "1": 4, "2": 3},
        "orbits": topology,
    }


def audit() -> dict[str, object]:
    admissible, archive = archive_replay()
    census_summary, candidates = census(admissible)
    result = {
        "blocker_and_ridges": blocker_and_ridges(),
        "x_archive": archive,
        "census": census_summary,
        "residual": orbit_and_topology(admissible, candidates),
        "conclusion": (
            "record 17 has exactly eight incidence-complex orbits after the "
            "strict a=5,6 exclusion: one with a=0, four with a=1, three with a=2"
        ),
        "nonclaim": "mod-two homology links do not certify PL or convex realizability",
    }
    digest = sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("independent record-17 audit summary changed")
    result["summary_sha256"] = digest
    return result


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
