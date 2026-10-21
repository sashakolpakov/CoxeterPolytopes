#!/usr/bin/env python3
"""Hostile archive-only replay of all 344 record-12 residual orbits.

This script intentionally does not import the production record-12 residual
module.  It reconstructs both branch clutters, every facet complement, the
full blocker, ridge degrees, dual connectivity, face vector, and the link of
every face.  Only the generic GF(2) homology routine is shared.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from itertools import combinations
import json
from math import comb
import multiprocessing
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.simplicial_homology import (
    first_mod2_homology_sphere_link_failure,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record12_residual.json"
FULL14 = (1 << 14) - 1
FULL10 = (1 << 10) - 1
MATRIX_ALL = (1 << 25) - 1
PAIR_INDEX = {
    pair: index for index, pair in enumerate(combinations(range(10), 2))
}
GRAPH_EDGES = {
    (1 << 0) | (1 << 12),
    (1 << 1) | (1 << 12),
    (1 << 0) | (1 << 13),
    (1 << 1) | (1 << 13),
}
EXPECTED_SUMMARY_SHA256 = (
    "db5142ca271b7794a1e04cbfeea10b29da36f8a4e2b21ba8a3cc40e92c0ba531"
)


def msk(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def bits(value: int):
    while value:
        bit = value & -value
        yield bit
        value ^= bit


def decode(word: str) -> tuple[int, tuple[tuple[int, int], ...]]:
    codes = [ord(character) - 63 for character in word]
    if not codes or codes[0] == 63 or any(code < 0 or code > 63 for code in codes):
        raise ValueError("unsupported graph6 word in record-12 hostile replay")
    order = codes[0]
    bit_count = order * (order - 1) // 2
    if len(codes) != (bit_count + 5) // 6 + 1:
        raise ValueError("wrong graph6 payload length")
    payload = []
    for code in codes[1:]:
        payload.extend((code >> shift) & 1 for shift in range(5, -1, -1))
    if any(payload[bit_count:]):
        raise ValueError("nonzero graph6 padding")
    edges = []
    position = 0
    for upper in range(1, order):
        for lower in range(upper):
            if payload[position]:
                edges.append((lower, upper))
            position += 1
    return order, tuple(edges)


def adjacency(word: str) -> tuple[int, ...]:
    order, edges = decode(word)
    if order != 10:
        raise ValueError("hostile branch graph does not have order ten")
    result = [0] * order
    for first, second in edges:
        result[first] |= 1 << second
        result[second] |= 1 << first
    return tuple(result)


def pair_mask_from_adjacency(rows: Sequence[int]) -> int:
    return sum(
        1 << PAIR_INDEX[first, second]
        for first, second in combinations(range(10), 2)
        if rows[first] & (1 << second)
    )


def pair_mask_from_matrix(matrix: int) -> int:
    if matrix < 0 or matrix & ~MATRIX_ALL:
        raise ValueError("hostile branch matrix is out of range")
    return sum(
        1 << PAIR_INDEX[first, second + 5]
        for first in range(5)
        for second in range(5)
        if matrix & (1 << (5 * first + second))
    )


def edge(value: int, first: int, second: int) -> bool:
    if first > second:
        first, second = second, first
    return bool(value & (1 << PAIR_INDEX[first, second]))


def cut_edge(cut: int, first: int, second: int) -> bool:
    return bool(cut & (1 << first)) != bool(cut & (1 << second))


def branch_triples(x_mask: int, cut: int, y: int) -> tuple[int, ...]:
    result = set()
    for isolate in range(10):
        if y & (1 << isolate):
            result.add(msk((10, 11, isolate)))
    for first, second in combinations(range(10), 2):
        x = edge(x_mask, first, second)
        delta = cut_edge(cut, first, second)
        if x:
            result.add(msk((11, first, second)))
        if x != delta:
            result.add(msk((10, first, second)))
    for triple in combinations(range(10), 3):
        if sum(edge(x_mask, *pair) for pair in combinations(triple, 2)) % 2:
            result.add(msk(triple))
    return tuple(sorted(result))


def blocker(edges: Iterable[int]) -> tuple[int, ...]:
    rows = [0]
    for hyperedge in sorted(set(edges)):
        candidates = []
        for row in rows:
            if row & hyperedge:
                candidates.append(row)
            else:
                for bit in bits(hyperedge):
                    candidates.append(row | bit)
        candidates = sorted(set(candidates), key=lambda row: (row.bit_count(), row))
        rows = []
        for row in candidates:
            if not any(prior & row == prior for prior in rows):
                rows.append(row)
    return tuple(rows)


def allowed_z(x: int, dp: int, dq: int, yi: int, yj: int) -> tuple[int, ...]:
    return tuple(
        z
        for z in (0, 1)
        if yi + yj + x + (x ^ dp) + z + (z ^ dq) in (0, 2)
    )


def allowed_adjacency(word: str, p: int, y: int) -> tuple[int, ...]:
    source = adjacency(word)
    q = p ^ y
    result = [0] * 10
    for first, second in combinations(range(10), 2):
        if 1 in allowed_z(
            int(bool(source[first] & (1 << second))),
            int(cut_edge(p, first, second)),
            int(cut_edge(q, first, second)),
            int(bool(y & (1 << first))),
            int(bool(y & (1 << second))),
        ):
            result[first] |= 1 << second
            result[second] |= 1 << first
    return tuple(result)


def embedding_exists(
    source: Sequence[int],
    source_cut: int,
    source_y: int,
    target: Sequence[int],
    target_cut: int,
    target_y: int,
) -> bool:
    order = tuple(sorted(range(10), key=lambda v: (-source[v].bit_count(), v)))
    mapping: dict[int, int] = {}
    used = 0

    def search(index: int) -> bool:
        nonlocal used
        if index == 10:
            return True
        vertex = order[index]
        colour = (
            int(bool(source_cut & (1 << vertex))),
            int(bool(source_y & (1 << vertex))),
        )
        for image in range(10):
            if used & (1 << image):
                continue
            if colour != (
                int(bool(target_cut & (1 << image))),
                int(bool(target_y & (1 << image))),
            ):
                continue
            if any(
                source[vertex] & (1 << prior)
                and not target[image] & (1 << mapped)
                for prior, mapped in mapping.items()
            ):
                continue
            mapping[vertex] = image
            used |= 1 << image
            if search(index + 1):
                return True
            used ^= 1 << image
            del mapping[vertex]
        return False

    return search(0)


def branch_and_pair_replay(document: Mapping[str, object]) -> dict[str, object]:
    branches = document["coloured_branch_orbits"]
    if not isinstance(branches, list) or len(branches) != 28:
        raise AssertionError("hostile record-12 branch count changed")
    seen = set()
    for index, state in enumerate(branches):
        if state["index"] != index:
            raise AssertionError("hostile branch indices are not canonical")
        key = (state["colour_partition"], state["coloured_canonical_graph6"])
        if key in seen:
            raise AssertionError("hostile branch orbit repeats")
        seen.add(key)
        x_mask = pair_mask_from_adjacency(adjacency(state["x_graph6"]))
        rows = blocker(branch_triples(x_mask, int(state["p_mask"]), int(state["y_mask"])))
        histogram = Counter(row.bit_count() for row in rows)
        text = ",".join(f"{size}:{count}" for size, count in sorted(histogram.items()))
        if text != state["blocker_histogram"] or set(histogram) - {4, 5}:
            raise AssertionError("hostile coloured branch blocker changed")

    pairs = []
    complement_only = 0
    for p_index, p_state in enumerate(branches):
        p = int(p_state["p_mask"])
        y = int(p_state["y_mask"])
        q = p ^ y
        target = allowed_adjacency(p_state["x_graph6"], p, y)
        for q_index, q_state in enumerate(branches):
            if int(q_state["y_mask"]).bit_count() != y.bit_count():
                continue
            source = adjacency(q_state["x_graph6"])
            source_cut = int(q_state["p_mask"])
            source_y = int(q_state["y_mask"])
            direct = embedding_exists(source, source_cut, source_y, target, q, y)
            complemented = False
            if not direct:
                complemented = embedding_exists(
                    source, source_cut ^ FULL10, source_y, target, q, y
                )
            if direct or complemented:
                pairs.append([p_index, q_index, complemented])
                complement_only += int(complemented)
    if pairs != document["pairable_ordered_branch_pairs"]:
        raise AssertionError("hostile complementary-cut pair list changed")
    kinds = Counter()
    for first, second, _ in pairs:
        p_state = branches[first]
        q_state = branches[second]
        if int(p_state["y_mask"]) or int(q_state["y_mask"]):
            kinds["nonzero_y"] += 1
        elif int(p_state["p_mask"]).bit_count() == 5 and int(q_state["p_mask"]).bit_count() == 5:
            kinds["balanced"] += 1
        elif int(p_state["p_mask"]) == 0 and int(q_state["p_mask"]) == 0:
            kinds["trivial"] += 1
        else:
            kinds["other"] += 1
    if kinds != {"balanced": 100, "trivial": 8} or complement_only:
        raise AssertionError("hostile cut normalization conclusion changed")
    return {
        "coloured_branch_orbits": len(branches),
        "pairable_ordered_type_pairs": len(pairs),
        "balanced_ordered_type_pairs": kinds["balanced"],
        "trivial_ordered_type_pairs": kinds["trivial"],
        "nonzero_y_pair_count": kinds["nonzero_y"],
        "complement_only_pair_count": complement_only,
    }


def selected_complements(x_matrix: int, z_matrix: int) -> tuple[int, ...]:
    x_mask = pair_mask_from_matrix(x_matrix)
    z_mask = pair_mask_from_matrix(z_matrix)
    cut = sum(1 << vertex for vertex in range(5, 10))
    p_triples = branch_triples(x_mask, cut, 0)
    q_triples = branch_triples(z_mask, cut, 0)
    result = set()
    for triple in p_triples:
        global_triple = 0
        for local in range(10):
            if triple & (1 << local):
                global_triple |= 1 << (local + 2)
        if triple & (1 << 10):
            global_triple |= 1 << 12
        if triple & (1 << 11):
            global_triple |= 1 << 13
        result.add(msk((0, 1)) | global_triple)
    for triple in q_triples:
        global_triple = 0
        for local in range(10):
            if triple & (1 << local):
                global_triple |= 1 << (local + 2)
        if triple & (1 << 10):
            global_triple |= 1 << 0
        if triple & (1 << 11):
            global_triple |= 1 << 1
        result.add(msk((12, 13)) | global_triple)
    return tuple(sorted(result))


def face_table(complements: Sequence[int]) -> tuple[bool, ...]:
    table = [False] * (1 << 14)
    for complement in complements:
        facet = FULL14 ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def minimal_nonfaces(table: Sequence[bool]) -> dict[int, tuple[int, ...]]:
    result: dict[int, list[int]] = {}
    for value in range(1, 1 << 14):
        if table[value]:
            continue
        if all(table[value ^ bit] for bit in bits(value)):
            result.setdefault(value.bit_count(), []).append(value)
    return {size: tuple(rows) for size, rows in sorted(result.items())}


def ridge_histogram(complements: Sequence[int]) -> Counter[int]:
    groups: Counter[int] = Counter()
    for complement in complements:
        for bit in bits(FULL14 ^ complement):
            groups[complement | bit] += 1
    return Counter(groups.values())


def dual_connected(complements: Sequence[int]) -> bool:
    values = set(complements)
    root = min(values)
    seen = {root}
    queue = deque([root])
    while queue:
        current = queue.popleft()
        for added in bits(FULL14 ^ current):
            for removed in bits(current):
                neighbor = (current | added) ^ removed
                if neighbor in values and neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
    return len(seen) == len(values)


def macaulay_raise(value: int, degree: int) -> int:
    if value < 0 or degree < 1:
        raise ValueError("invalid hostile Macaulay input")
    remaining = value
    next_top = None
    raised = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top is None or top + 1 < next_top)
            and comb(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= comb(top, lower)
        raised += comb(top + 1, lower + 1)
        next_top = top
    if remaining:
        raise AssertionError("hostile Macaulay expansion was incomplete")
    return raised


def link_g_entries(table: Sequence[bool], face: int) -> tuple[int, tuple[int, ...]]:
    face_size = face.bit_count()
    polytope_dimension = 9 - face_size
    remaining = tuple(vertex for vertex in range(14) if not face & (1 << vertex))
    link_f = {-1: 1}
    for dimension in range(polytope_dimension):
        link_f[dimension] = sum(
            table[face | msk(extension)]
            for extension in combinations(remaining, dimension + 1)
        )
    h = []
    for index in range(polytope_dimension // 2 + 1):
        h.append(
            sum(
                (-1) ** (index - face_order)
                * comb(polytope_dimension - face_order, polytope_dimension - index)
                * link_f[face_order - 1]
                for face_order in range(index + 1)
            )
        )
    return link_f[0], tuple(h[index] - h[index - 1] for index in range(1, len(h)))


def local_g_and_macaulay_audit(table: Sequence[bool]) -> None:
    strict_orders = {1: (2, 3, 4), 2: (2, 3), 3: (2,)}
    for face_size in (1, 2, 3):
        polytope_dimension = 9 - face_size
        for chosen in combinations(range(14), face_size):
            face = msk(chosen)
            if any(face & graph_edge == graph_edge for graph_edge in GRAPH_EDGES):
                continue
            if not table[face]:
                raise AssertionError("hostile independent small face is absent")
            link_vertices, g_entries = link_g_entries(table, face)
            if link_vertices < polytope_dimension + 2:
                raise AssertionError("hostile face link is a forbidden simplex")
            for order, observed in enumerate(g_entries, start=1):
                minimum = 1 if order in strict_orders[face_size] else 0
                if observed < minimum:
                    raise AssertionError("hostile strict local g inequality failed")
            for previous_order in range(1, len(g_entries)):
                previous = g_entries[previous_order - 1]
                observed = g_entries[previous_order]
                if previous < 0 or observed > macaulay_raise(previous, previous_order):
                    raise AssertionError("hostile local Macaulay inequality failed")


def audit_orbit(task: tuple[int, Mapping[str, object]]) -> dict[str, object]:
    index, orbit = task
    complements = selected_complements(int(orbit["x_matrix"]), int(orbit["z_matrix"]))
    digest = sha256(
        b"".join(value.to_bytes(2, "little") for value in complements)
    ).hexdigest()
    if digest != orbit["selected_complement_sha256"]:
        raise AssertionError("hostile selected-complement digest changed")
    table = face_table(complements)
    nonfaces = minimal_nonfaces(table)
    histogram = {size: len(rows) for size, rows in nonfaces.items()}
    a = len(nonfaces.get(4, ()))
    counts = Counter(value.bit_count() for value, yes in enumerate(table) if yes)
    vector = tuple(counts[size] for size in range(1, 10))
    if (
        set(nonfaces) - {2, 4, 5}
        or set(nonfaces.get(2, ())) != GRAPH_EDGES
        or a != orbit["missing_four_count"]
        or len(complements) != orbit["facet_count"]
        or vector != tuple(orbit["face_vector"])
        or ridge_histogram(complements) != Counter({2: vector[-2]})
        or not dual_connected(complements)
    ):
        raise AssertionError("hostile full incidence audit failed")
    local_g_and_macaulay_audit(table)
    failure = first_mod2_homology_sphere_link_failure(table, 8)
    if failure is not None:
        raise AssertionError(
            f"hostile residual {index} has a bad link at {failure.face_mask}"
        )
    return {
        "index": index,
        "missing_four_count": a,
        "facet_count": len(complements),
        "minimal_nonface_histogram": {
            str(size): count for size, count in histogram.items()
        },
        "ridge_degrees_two": True,
        "dual_connected": True,
        "strict_local_g_and_macaulay": True,
        "all_links_mod2_homology_spheres": True,
        "selected_complement_sha256": digest,
    }


def archive_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def audit(workers: int = 1) -> dict[str, object]:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    digest = archive_digest(document)
    if digest != document["certificate_sha256"]:
        raise AssertionError("hostile record-12 archive digest changed")
    branch = branch_and_pair_replay(document)
    tasks = list(enumerate(document["residual_orbits"]))
    if workers > 1:
        context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
            orbits = list(executor.map(audit_orbit, tasks, chunksize=1))
    else:
        orbits = [audit_orbit(task) for task in tasks]
    schedule = Counter(orbit["missing_four_count"] for orbit in orbits)
    if dict(sorted(schedule.items())) != {0: 48, 1: 95, 2: 129, 3: 54, 4: 18}:
        raise AssertionError("hostile record-12 schedule histogram changed")
    summary: dict[str, object] = {
        "archive_certificate_sha256": digest,
        "branch_and_cut_normalization": branch,
        "residual": {
            "orbit_count": len(orbits),
            "schedule_orbit_counts": {
                str(key): value for key, value in sorted(schedule.items())
            },
            "all_full_blockers_verified": True,
            "all_ridge_degrees_two": True,
            "all_dual_graphs_connected": True,
            "all_strict_local_g_and_macaulay_audits": True,
            "all_face_links_mod2_homology_spheres": True,
            "orbits": orbits,
        },
    }
    payload = json.dumps(summary, sort_keys=True, separators=(",", ":")).encode()
    summary["summary_sha256"] = sha256(payload).hexdigest()
    if (
        EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN"
        and summary["summary_sha256"] != EXPECTED_SUMMARY_SHA256
    ):
        raise AssertionError("hostile record-12 residual summary changed")
    return summary


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args(argv)
    print(json.dumps(audit(args.workers), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
