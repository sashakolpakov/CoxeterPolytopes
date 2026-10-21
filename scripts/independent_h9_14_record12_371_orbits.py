#!/usr/bin/env python3
"""Hostile quotient and topology replay for the corrected 371 record-12 orbits.

No production record-12 census/residual module is imported.  The quotient is
recomputed by a direct canonical form for two-coloured 5x5 matrices under
independent layer complements, layer swap, row/column permutations, and
bipartition transpose.  Every quotient representative then receives a full
incidence and all-face-link audit.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb
import multiprocessing
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.simplicial_homology import (
    first_mod2_homology_sphere_link_failure,
)


PRODUCTION_PAIRS = Path("/tmp/record12_safe_pairs_all_v2.txt")
HOSTILE_PAIRS = Path("/tmp/record12_hostile_pairs.txt")
ORBIT_STREAM = Path("/tmp/record12_safe_full_orbits_v2.txt")
EXPECTED_PRODUCTION_PAIRS_SHA256 = (
    "6b86880f60359a6ae638acf55a730282f2b46e5cfceb74c45d8ae7f0035d635a"
)
EXPECTED_SORTED_PAIR_SHA256 = (
    "a189faa3d1faeb1f4db21630d2fbb816f8d3d4f1421d04eb251112d37945f336"
)
EXPECTED_ORBIT_STREAM_SHA256 = (
    "4d978594b6be052328b4056d6b05c6249e80e8b156a09d2eba8057f3b7daa23c"
)

N = 10
FULL10 = (1 << 10) - 1
FULL14 = (1 << 14) - 1
MATRIX_ALL = (1 << 25) - 1
BALANCED_CUT = sum(1 << vertex for vertex in range(5, 10))
PAIR_INDEX = {pair: index for index, pair in enumerate(combinations(range(10), 2))}
PERMUTATIONS5 = tuple(permutations(range(5)))
GRAPH_EDGES = {
    (1 << 0) | (1 << 12),
    (1 << 1) | (1 << 12),
    (1 << 0) | (1 << 13),
    (1 << 1) | (1 << 13),
}
SCHEDULES = {
    0: (14, 87, 320, 780, 1320, 1554, 1224, 585, 130),
    1: (14, 87, 320, 779, 1314, 1540, 1208, 576, 128),
    2: (14, 87, 320, 778, 1308, 1526, 1192, 567, 126),
    3: (14, 87, 320, 777, 1302, 1512, 1176, 558, 124),
    4: (14, 87, 320, 776, 1296, 1498, 1160, 549, 122),
}


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def bits(value: int):
    while value:
        bit = value & -value
        yield bit
        value ^= bit


def msk(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def matrix_from_pair_mask(value: int) -> int:
    result = 0
    for row in range(5):
        for column in range(5):
            if value & (1 << PAIR_INDEX[row, column + 5]):
                result |= 1 << (5 * row + column)
    if result.bit_count() != value.bit_count():
        raise AssertionError("hostile balanced graph has a non-cross edge")
    return result


def transpose(value: int) -> int:
    result = 0
    for row in range(5):
        for column in range(5):
            if value & (1 << (5 * row + column)):
                result |= 1 << (5 * column + row)
    return result


def canonical_bipartite_pair(x: int, z: int) -> int:
    """Canonicalize under S5xS5 for a fixed bipartition/layer orientation.

    For each of the 120 row permutations, column permutations are eliminated
    exactly by sorting the five combined layer-column vectors in decreasing
    order: assigning the largest vector to the least significant target
    column minimizes the row-major integer.
    """

    x_rows = tuple((x >> (5 * row)) & 31 for row in range(5))
    z_rows = tuple((z >> (5 * row)) & 31 for row in range(5))
    best = (1 << 50) - 1
    for row_order in PERMUTATIONS5:
        columns = []
        for source_column in range(5):
            value = 0
            for target_row, source_row in enumerate(row_order):
                value |= ((x_rows[source_row] >> source_column) & 1) << (5 * target_row)
                value |= ((z_rows[source_row] >> source_column) & 1) << (
                    25 + 5 * target_row
                )
            columns.append(value)
        encoded = 0
        for target_column, column in enumerate(sorted(columns, reverse=True)):
            encoded |= column << target_column
        best = min(best, encoded)
    return best


def canonical_pair(x: int, z: int) -> int:
    best = (1 << 50) - 1
    for swapped in (False, True):
        first, second = (z, x) if swapped else (x, z)
        for complement_first in (False, True):
            for complement_second in (False, True):
                xx = first ^ (MATRIX_ALL if complement_first else 0)
                zz = second ^ (MATRIX_ALL if complement_second else 0)
                best = min(
                    best,
                    canonical_bipartite_pair(xx, zz),
                    canonical_bipartite_pair(transpose(xx), transpose(zz)),
                )
    return best


def sorted_line_sha256(lines: Sequence[str]) -> str:
    return sha256("".join(line + "\n" for line in sorted(lines)).encode()).hexdigest()


def quotient_replay() -> tuple[dict[str, object], tuple[tuple[int, int, int], ...]]:
    if file_sha256(PRODUCTION_PAIRS) != EXPECTED_PRODUCTION_PAIRS_SHA256:
        raise AssertionError("hostile production pair stream changed")
    if file_sha256(ORBIT_STREAM) != EXPECTED_ORBIT_STREAM_SHA256:
        raise AssertionError("hostile orbit stream changed")
    production_lines = PRODUCTION_PAIRS.read_text(encoding="ascii").splitlines()
    hostile_lines = HOSTILE_PAIRS.read_text(encoding="ascii").splitlines()
    if (
        len(production_lines) != 42_000
        or len(set(production_lines)) != 42_000
        or set(production_lines) != set(hostile_lines)
        or sorted_line_sha256(hostile_lines) != EXPECTED_SORTED_PAIR_SHA256
    ):
        raise AssertionError("hostile and production pair sets differ")

    observed = Counter()
    key_to_a: dict[int, int] = {}
    for line in hostile_lines:
        fields = line.split()
        if fields[0] != "PAIR" or len(fields) != 10:
            raise AssertionError("malformed hostile pair row")
        _p_index, _q_index, raw_x, raw_z, cut, y, facets, a, b = map(
            int, fields[1:]
        )
        if (cut, y, facets, b) != (BALANCED_CUT, 0, 130 - 2 * a, 12 - 2 * a):
            raise AssertionError("hostile passing pair left the balanced schedule")
        x = matrix_from_pair_mask(raw_x)
        z = matrix_from_pair_mask(raw_z)
        key = canonical_pair(x, z)
        observed[key] += 1
        previous = key_to_a.setdefault(key, a)
        if previous != a:
            raise AssertionError("hostile quotient mixed Dehn--Sommerville branches")

    archived = {}
    representatives = []
    for line in ORBIT_STREAM.read_text(encoding="ascii").splitlines():
        marker, raw_x, raw_z, raw_multiplicity = line.split()
        if marker != "ORBIT":
            raise AssertionError("malformed hostile orbit row")
        x, z, multiplicity = int(raw_x), int(raw_z), int(raw_multiplicity)
        key = x | (z << 25)
        if canonical_pair(x, z) != key:
            raise AssertionError("archived hostile orbit representative is not canonical")
        if key in archived:
            raise AssertionError("archived hostile orbit representative repeats")
        archived[key] = multiplicity
        representatives.append((x, z, multiplicity))
    if observed != Counter(archived) or len(observed) != 371:
        raise AssertionError("hostile D8 x stabilizer quotient changed")
    schedule = Counter(key_to_a[key] for key in key_to_a)
    expected_schedule = {0: 48, 1: 95, 2: 144, 3: 63, 4: 21}
    if dict(sorted(schedule.items())) != expected_schedule:
        raise AssertionError("hostile 371-orbit schedule histogram changed")
    return (
        {
            "passing_pair_rows": len(hostile_lines),
            "semantic_pair_set_sha256": EXPECTED_SORTED_PAIR_SHA256,
            "full_symmetry_orbits": len(observed),
            "orbit_multiplicity_sum": sum(observed.values()),
            "schedule_orbit_histogram": {
                str(a): count for a, count in sorted(schedule.items())
            },
            "all_passing_pairs_balanced_y_zero": True,
            "canonical_action": (
                "independent layer complements and swap (D8), S5xS5 row/column "
                "permutations, and bipartition transpose"
            ),
        },
        tuple(representatives),
    )


def pair_mask_from_matrix(matrix: int) -> int:
    return sum(
        1 << PAIR_INDEX[row, column + 5]
        for row in range(5)
        for column in range(5)
        if matrix & (1 << (5 * row + column))
    )


def edge(mask: int, first: int, second: int) -> bool:
    if first > second:
        first, second = second, first
    return bool(mask & (1 << PAIR_INDEX[first, second]))


def cut_edge(cut: int, first: int, second: int) -> bool:
    return bool(cut & (1 << first)) != bool(cut & (1 << second))


def branch_triples(mask: int, cut: int) -> tuple[int, ...]:
    result = set()
    for first, second in combinations(range(10), 2):
        x = edge(mask, first, second)
        delta = cut_edge(cut, first, second)
        if x:
            result.add((1 << 11) | (1 << first) | (1 << second))
        if x != delta:
            result.add((1 << 10) | (1 << first) | (1 << second))
    for triple in combinations(range(10), 3):
        if sum(edge(mask, *pair) for pair in combinations(triple, 2)) % 2:
            result.add(msk(triple))
    return tuple(sorted(result))


def selected_complements(x: int, z: int) -> tuple[int, ...]:
    p_triples = branch_triples(pair_mask_from_matrix(x), BALANCED_CUT)
    q_triples = branch_triples(pair_mask_from_matrix(z), BALANCED_CUT)
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
            if subset == 0:
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
        raise AssertionError("hostile Macaulay expansion incomplete")
    return raised


def link_g_entries(table: Sequence[bool], face: int) -> tuple[int, tuple[int, ...]]:
    face_size = face.bit_count()
    dimension = 9 - face_size
    remaining = tuple(vertex for vertex in range(14) if not face & (1 << vertex))
    link_f = {-1: 1}
    for degree in range(dimension):
        link_f[degree] = sum(
            table[face | msk(extension)]
            for extension in combinations(remaining, degree + 1)
        )
    h = []
    for index in range(dimension // 2 + 1):
        h.append(
            sum(
                (-1) ** (index - order)
                * comb(dimension - order, dimension - index)
                * link_f[order - 1]
                for order in range(index + 1)
            )
        )
    return link_f[0], tuple(h[index] - h[index - 1] for index in range(1, len(h)))


def local_g_and_macaulay(table: Sequence[bool]) -> None:
    strict_orders = {1: (2, 3, 4), 2: (2, 3), 3: (2,)}
    for face_size in (1, 2, 3):
        dimension = 9 - face_size
        for chosen in combinations(range(14), face_size):
            face = msk(chosen)
            if any(face & graph_edge == graph_edge for graph_edge in GRAPH_EDGES):
                continue
            if not table[face]:
                raise AssertionError("hostile small face absent")
            link_vertices, g = link_g_entries(table, face)
            if link_vertices < dimension + 2:
                raise AssertionError("hostile forbidden simplex link")
            for order, observed in enumerate(g, start=1):
                if observed < (1 if order in strict_orders[face_size] else 0):
                    raise AssertionError("hostile strict local g failed")
            for previous_order in range(1, len(g)):
                if g[previous_order] > macaulay_raise(g[previous_order - 1], previous_order):
                    raise AssertionError("hostile local Macaulay failed")


def all_link_euler(table: Sequence[bool]) -> None:
    signed = [
        (-1 if value.bit_count() % 2 else 1) if is_face else 0
        for value, is_face in enumerate(table)
    ]
    supersets = signed[:]
    for vertex in range(14):
        bit = 1 << vertex
        for value in range(1 << 14):
            if not value & bit:
                supersets[value] += supersets[value | bit]
    for face, is_face in enumerate(table):
        if not is_face:
            continue
        reduced_euler = -((-1 if face.bit_count() % 2 else 1) * supersets[face])
        link_dimension = 8 - face.bit_count()
        expected = -1 if link_dimension % 2 else 1
        if reduced_euler != expected:
            raise AssertionError("hostile face-link Euler characteristic failed")


def audit_orbit(task: tuple[int, tuple[int, int, int]]) -> dict[str, object]:
    index, (x, z, multiplicity) = task
    complements = selected_complements(x, z)
    table = face_table(complements)
    nonfaces = minimal_nonfaces(table)
    a = len(nonfaces.get(4, ()))
    b = len(nonfaces.get(5, ()))
    counts = Counter(value.bit_count() for value, yes in enumerate(table) if yes)
    vector = tuple(counts[size] for size in range(1, 10))
    if (
        set(nonfaces) - {2, 4, 5}
        or set(nonfaces.get(2, ())) != GRAPH_EDGES
        or a not in SCHEDULES
        or b != 12 - 2 * a
        or len(complements) != 130 - 2 * a
        or vector != SCHEDULES[a]
        or ridge_histogram(complements) != Counter({2: vector[-2]})
        or not dual_connected(complements)
    ):
        raise AssertionError(f"hostile incidence audit failed at orbit {index}")
    local_g_and_macaulay(table)
    all_link_euler(table)
    failure = first_mod2_homology_sphere_link_failure(table, 8)
    if failure is not None:
        raise AssertionError(
            f"hostile homology-link failure at orbit {index}, face {failure.face_mask}"
        )
    digest = sha256(
        b"".join(value.to_bytes(2, "little") for value in complements)
    ).hexdigest()
    return {
        "index": index,
        "x_matrix": x,
        "z_matrix": z,
        "input_multiplicity": multiplicity,
        "missing_four_count": a,
        "minimal_nonface_histogram": {"2": 4, "4": a, "5": b},
        "facet_count": len(complements),
        "face_vector": list(vector),
        "selected_complement_sha256": digest,
        "ridge_degrees_two": True,
        "dual_connected": True,
        "strict_local_g_and_macaulay": True,
        "all_face_links_euler_spheres": True,
        "all_face_links_mod2_homology_spheres": True,
    }


def archived_representatives() -> tuple[tuple[int, int, int], ...]:
    """Load the frozen 371 rows for a topology-only hostile replay.

    Their quotient status is established independently by the companion C++
    replay.  This loader pins the exact byte stream and checks only its syntax
    and uniqueness, keeping the expensive topology pass separate from the
    42,000-row quotient calculation.
    """

    if file_sha256(ORBIT_STREAM) != EXPECTED_ORBIT_STREAM_SHA256:
        raise AssertionError("hostile orbit stream changed")
    representatives = []
    seen = set()
    for line in ORBIT_STREAM.read_text(encoding="ascii").splitlines():
        fields = line.split()
        if len(fields) != 4 or fields[0] != "ORBIT":
            raise AssertionError("malformed hostile orbit row")
        x, z, multiplicity = map(int, fields[1:])
        row = (x, z, multiplicity)
        if (x, z) in seen:
            raise AssertionError("hostile orbit representative repeats")
        seen.add((x, z))
        representatives.append(row)
    if len(representatives) != 371 or sum(row[2] for row in representatives) != 42_000:
        raise AssertionError("hostile orbit archive has wrong size or multiplicity")
    return tuple(representatives)


def audit_representatives(
    representatives: Sequence[tuple[int, int, int]], workers: int
) -> dict[str, object]:
    tasks = list(enumerate(representatives))
    if workers > 1:
        context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as executor:
            orbits = list(executor.map(audit_orbit, tasks, chunksize=1))
    else:
        orbits = [audit_orbit(task) for task in tasks]
    schedule = Counter(orbit["missing_four_count"] for orbit in orbits)
    if dict(sorted(schedule.items())) != {0: 48, 1: 95, 2: 144, 3: 63, 4: 21}:
        raise AssertionError("hostile topology schedule changed")
    return {
        "orbit_count": len(orbits),
        "schedule_orbit_histogram": {
            str(a): count for a, count in sorted(schedule.items())
        },
        "all_full_blockers_verified": True,
        "all_face_vectors_match_schedule": True,
        "all_ridge_degrees_two": True,
        "all_dual_graphs_connected": True,
        "all_strict_local_g_and_macaulay_audits": True,
        "all_face_links_euler_spheres": True,
        "all_face_links_mod2_homology_spheres": True,
        "orbits": orbits,
    }


def replay(workers: int = 1, topology_only: bool = False) -> dict[str, object]:
    if topology_only:
        representatives = archived_representatives()
        quotient: dict[str, object] = {
            "mode": "topology-only; quotient checked by independent C++ companion",
            "orbit_stream_sha256": EXPECTED_ORBIT_STREAM_SHA256,
            "full_symmetry_orbits": 371,
            "orbit_multiplicity_sum": 42_000,
        }
    else:
        quotient, representatives = quotient_replay()
    summary: dict[str, object] = {
        "quotient": quotient,
        "topology": audit_representatives(representatives, workers),
    }
    payload = json.dumps(summary, sort_keys=True, separators=(",", ":")).encode()
    summary["summary_sha256"] = sha256(payload).hexdigest()
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--topology-only",
        action="store_true",
        help="use the SHA-pinned 371 archive after the companion C++ quotient replay",
    )
    args = parser.parse_args()
    result = replay(args.workers, args.topology_only)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
