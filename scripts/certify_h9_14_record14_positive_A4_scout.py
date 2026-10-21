#!/usr/bin/env python3
"""Freeze and verify the compact-only record-14 positive-A4 scout."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb, factorial
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/explore_h9_14_record14_positive_orbits.cpp"
ARCHIVE = ROOT / "data/enumerations/h9_14_record14_positive_A4_scout.json"
GRAPH6 = "M???????????o?w??"
N = 14
MISSING_EDGES = (
    (0, 12),
    (1, 12),
    (0, 13),
    (1, 13),
    (2, 13),
)
EDGE_SET = frozenset(tuple(sorted(edge)) for edge in MISSING_EDGES)
EXPECTED_BRANCH_SUMMARY = (
    (1, 579, 579, 369, 5),
    (2, 167_331, 160_245, 60_810, 36),
    (3, 32_183_329, 28_293_741, 580_272, 59),
    (4, 4_634_399_376, 3_571_339_779, 28_854, 13),
)


def _is_missing_edge(left: int, right: int) -> bool:
    return tuple(sorted((left, right))) in EDGE_SET


def _independent(vertices: Iterable[int]) -> bool:
    chosen = tuple(vertices)
    return all(not _is_missing_edge(left, right) for left, right in combinations(chosen, 2))


def _independent_sets(size: int) -> tuple[tuple[int, ...], ...]:
    return tuple(chosen for chosen in combinations(range(N), size) if _independent(chosen))


def _macaulay_raise(value: int, degree: int) -> int:
    remaining = value
    next_top: int | None = None
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
        raise AssertionError("inexact Macaulay expansion")
    return raised


def _strict_interval(
    previous: int, beta: int, candidate_count: int, degree: int
) -> list[int]:
    lower = max(0, beta - _macaulay_raise(previous, degree))
    upper = min(candidate_count, beta - 1)
    if lower > upper:
        raise AssertionError("empty graph-derived local interval")
    return [lower, upper]


def derive_schedules() -> list[dict[str, object]]:
    independent = {size: len(_independent_sets(size)) for size in range(2, 10)}
    h1 = N - 9
    h2 = independent[2] - 8 * N + 36
    h3 = independent[3] - 7 * independent[2] + 28 * N - 84
    h4_without_fours = (
        independent[4]
        - 6 * independent[3]
        + 21 * independent[2]
        - 56 * N
        + 126
    )
    rows = []
    for missing_fours in range(independent[4] + 1):
        h4 = h4_without_fours - missing_fours
        global_g4 = h4 - h3
        if global_g4 < 1:
            break
        h = [1, h1, h2, h3, h4, h4, h3, h2, h1, 1]
        f = [
            sum(comb(9 - index, 9 - size) * h[index] for index in range(size + 1))
            for size in range(1, 10)
        ]
        summed_link_g4 = 5 * f[4] - 24 * f[3] + 63 * f[2] - 112 * f[1] + 126 * N
        rows.append(
            {
                "a": missing_fours,
                "missing_independent_fives": independent[5] - f[4],
                "global_g4": global_g4,
                "f_vector": f,
                "h_vector": h,
                "summed_vertex_link_g4": summed_link_g4,
                "strict_vertex_sum_feasible": summed_link_g4 >= N,
            }
        )
    return rows


def derive_local_link_constraints() -> dict[str, object]:
    vertices = []
    for vertex in range(N):
        link = [
            other
            for other in range(N)
            if other != vertex and not _is_missing_edge(vertex, other)
        ]
        inherited = sum(_is_missing_edge(*pair) for pair in combinations(link, 2))
        candidate_edges = comb(len(link), 2) - inherited
        candidate_triangles = sum(_independent(chosen) for chosen in combinations(link, 3))
        g1 = len(link) - 9
        g2 = candidate_edges - 8 * len(link) + 36
        beta = candidate_triangles - 7 * candidate_edges + 28 * len(link) - 84
        vertices.append(
            {
                "vertex": vertex,
                "link_vertex_count": len(link),
                "inherited_missing_edge_count": inherited,
                "candidate_link_triangle_count": candidate_triangles,
                "g1": g1,
                "g2": g2,
                "g3_without_missing_fours": beta,
                "missing_four_incidence_interval": _strict_interval(
                    g2, beta, candidate_triangles, 2
                ),
            }
        )

    edges = []
    for left, right in combinations(range(N), 2):
        if _is_missing_edge(left, right):
            continue
        link = [
            other
            for other in range(N)
            if other not in (left, right)
            and not _is_missing_edge(other, left)
            and not _is_missing_edge(other, right)
        ]
        candidate_edges = sum(
            not _is_missing_edge(first, second)
            for first, second in combinations(link, 2)
        )
        g1 = len(link) - 8
        beta = candidate_edges - 7 * len(link) + 28
        edges.append(
            {
                "edge": [left, right],
                "link_vertex_count": len(link),
                "candidate_link_edge_count": candidate_edges,
                "g1": g1,
                "g2_without_missing_fours": beta,
                "missing_four_incidence_interval": _strict_interval(
                    g1, beta, candidate_edges, 1
                ),
            }
        )

    triples = []
    for triple in combinations(range(N), 3):
        if not _independent(triple):
            continue
        candidate_vertices = sum(
            other not in triple
            and all(not _is_missing_edge(other, vertex) for vertex in triple)
            for other in range(N)
        )
        triples.append(
            {
                "triple": list(triple),
                "candidate_link_vertex_count": candidate_vertices,
                "missing_four_incidence_interval": [0, candidate_vertices - 8],
            }
        )

    histogram = Counter(
        row["missing_four_incidence_interval"][1]
        for row in (*vertices, *edges, *triples)
    )
    return {
        "vertices": vertices,
        "edges": edges,
        "triples": triples,
        "upper_endpoint_histogram": {
            str(upper): count for upper, count in sorted(histogram.items())
        },
    }


def _interval_a4_count(local: Mapping[str, object]) -> int:
    constraints = []
    for row in local["vertices"]:  # type: ignore[index]
        constraints.append((frozenset((row["vertex"],)), row["missing_four_incidence_interval"]))
    for row in local["edges"]:  # type: ignore[index]
        constraints.append((frozenset(row["edge"]), row["missing_four_incidence_interval"]))
    for row in local["triples"]:  # type: ignore[index]
        constraints.append((frozenset(row["triple"]), row["missing_four_incidence_interval"]))
    return sum(
        all(interval[0] <= int(required <= frozenset(face)) <= interval[1]
            for required, interval in constraints)
        for face in _independent_sets(4)
    )


def _active_automorphism_order() -> int:
    active = (0, 1, 2, 12, 13)
    target = frozenset(EDGE_SET)
    order = 0
    for images in permutations(active):
        move = dict(zip(active, images, strict=True))
        moved = frozenset(
            tuple(sorted((move[left], move[right]))) for left, right in EDGE_SET
        )
        order += moved == target
    return order


def _source_sha256() -> str:
    return sha256(SOURCE.read_bytes()).hexdigest()


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"record-14 positive scout {name} is malformed")
    return value


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record14_positive_A4_scout":
        raise ValueError("wrong record-14 positive scout kind")
    if document.get("graph6") != GRAPH6 or document.get("missing_edges") != [
        list(edge) for edge in MISSING_EDGES
    ]:
        raise ValueError("record-14 graph source changed")
    if (
        document.get("theorem_certificate") is not False
        or document.get("independent_replay_complete") is not False
        or document.get("whole_shadow_compact_coxeter_theorem") is not False
        or document.get("compact_only") is not True
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
        or document.get("collision_free_symmetry_quotient") is not True
        or document.get("all_positive_strict_schedules_enumerated") is not True
    ):
        raise ValueError("record-14 scout scope was over-promoted")

    expected_schedules = derive_schedules()
    if document.get("schedules") != expected_schedules:
        raise ValueError("record-14 graph-derived schedules changed")
    if document.get("strict_compact_surviving_a") != [0, 1, 2, 3, 4] or document.get(
        "strict_compact_excluded_a"
    ) != [5, 6]:
        raise ValueError("record-14 strict schedule split changed")

    local = derive_local_link_constraints()
    if document.get("local_link_constraints") != local:
        raise ValueError("record-14 exact graph-derived local table changed")
    if len(local["vertices"]) != 14 or len(local["edges"]) != 86 or len(local["triples"]) != 310:
        raise ValueError("record-14 local table sizes changed")
    if _interval_a4_count(local) != 579 or document.get("interval_A4_count") != 579:
        raise ValueError("record-14 interval A4 ground set changed")
    if document.get("independent_A4_count") != 735:
        raise ValueError("record-14 independent A4 count changed")
    if document.get("active_automorphism_order") != _active_automorphism_order() != 2:
        raise ValueError("record-14 active automorphism order changed")

    raw_branches = document.get("branches")
    if not isinstance(raw_branches, list):
        raise ValueError("record-14 branches are malformed")
    branches = tuple(_mapping(row, "branch") for row in raw_branches)
    summary = tuple(
        (
            int(row["a"]),
            int(row["raw"]),
            int(row["interval_feasible"]),
            int(row["strict_aggregate_viable"]),
            int(row["orbit_count"]),
        )
        for row in branches
    )
    if summary != EXPECTED_BRANCH_SUMMARY:
        raise ValueError("record-14 positive branch summary changed")

    group_order = 2 * factorial(9)
    for row in branches:
        a = int(row["a"])
        if int(row["raw"]) != comb(579, a):
            raise ValueError(f"record-14 a={a} raw count is inconsistent")
        raw_orbits = row.get("orbits")
        if not isinstance(raw_orbits, list) or len(raw_orbits) != int(row["orbit_count"]):
            raise ValueError(f"record-14 a={a} orbit list is malformed")
        keys = []
        sizes = []
        for raw_orbit in raw_orbits:
            orbit = _mapping(raw_orbit, f"a={a} orbit")
            key = tuple(map(int, orbit["key"]))
            representative = orbit.get("representative")
            if not isinstance(representative, list) or len(representative) != a:
                raise ValueError(f"record-14 a={a} representative is malformed")
            if any(
                not isinstance(face, list)
                or len(face) != 4
                or len(set(map(int, face))) != 4
                for face in representative
            ):
                raise ValueError(f"record-14 a={a} representative A4 changed")
            keys.append(key)
            sizes.append(int(orbit["size"]))
        if len(set(keys)) != len(keys):
            raise ValueError(f"record-14 a={a} canonical key collision")
        if any(size <= 0 or group_order % size for size in sizes):
            raise ValueError(f"record-14 a={a} orbit size is invalid")
        if sum(sizes) != int(row["strict_aggregate_viable"]):
            raise ValueError(f"record-14 a={a} orbit weights do not sum")
        gates = row.get("gate_histogram")
        if not isinstance(gates, list) or sum(int(gate["count"]) for gate in gates) != int(
            row["interval_feasible"]
        ):
            raise ValueError(f"record-14 a={a} gate weights do not sum")
        if sum(
            int(gate["count"])
            for gate in gates
            if gate["strict"] is True and gate["aggregate"] is True
        ) != int(row["strict_aggregate_viable"]):
            raise ValueError(f"record-14 a={a} viable gate weights do not sum")

    a4 = branches[-1]
    if int(a4["ground_labelled_family_count"]) != comb(579, 4):
        raise ValueError("record-14 a=4 matrix quotient did not recover the raw ground set")
    if int(a4["column_paths"]) - int(a4["equal_active_duplicate_paths"]) - int(
        a4["repeated_full_row_paths"]
    ) != int(a4["distinct_matrix_orbits"]):
        raise ValueError("record-14 a=4 matrix quotient identity failed")
    provenance = _mapping(document.get("source_provenance"), "source provenance")
    if provenance.get("path") != str(SOURCE.relative_to(ROOT)) or provenance.get(
        "sha256"
    ) != _source_sha256():
        raise ValueError("record-14 positive scout source provenance changed")

    return {
        "record_index": 14,
        "schedule_count": len(expected_schedules),
        "local_constraint_count": 410,
        "strict_schedule_split": {"surviving": [0, 1, 2, 3, 4], "excluded": [5, 6]},
        "positive_branch_orbit_count": sum(int(row["orbit_count"]) for row in branches),
        "positive_branch_closed_count": sum(
            int(row["strict_aggregate_viable"]) == 0 for row in branches
        ),
        "collision_free_a4_ground_recovery": True,
        "compact_only": True,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def fresh_record() -> dict[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("no C++20 compiler is available")
    with tempfile.TemporaryDirectory(prefix="record14-positive-") as directory:
        executable = Path(directory) / "positive"
        subprocess.run(
            [
                compiler,
                "-std=c++20",
                "-O3",
                "-Wall",
                "-Wextra",
                "-pedantic",
                str(SOURCE),
                "-o",
                str(executable),
            ],
            cwd=ROOT,
            check=True,
        )
        completed = subprocess.run(
            [str(executable), "4"],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    document = json.loads(completed.stdout)
    document["source_provenance"] = {
        "path": str(SOURCE.relative_to(ROOT)),
        "sha256": _source_sha256(),
    }
    return document


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", type=Path)
    args = parser.parse_args(argv)
    document = fresh_record() if args.fresh else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    result = verify(document)
    if args.write is not None:
        rendered = json.dumps(document, indent=2, sort_keys=True) + "\n"
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(rendered, encoding="utf-8")
        result["written"] = str(args.write)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
