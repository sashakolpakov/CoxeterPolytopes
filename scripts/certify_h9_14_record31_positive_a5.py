#!/usr/bin/env python3
"""Freeze and verify the compact-only record-31 positive-A4/a=5 closure."""

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
PRODUCTION_SOURCE = ROOT / "scripts/explore_h9_14_record31_positive_orbits.cpp"
INDEPENDENT_SOURCE = ROOT / "scripts/independent_h9_14_record31_a5_audit.cpp"
ARCHIVE = ROOT / "data/enumerations/h9_14_record31_positive_a5.json"
GRAPH6 = "M?????????O?O?o??"
N = 14
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EDGE_SET = frozenset(tuple(sorted(edge)) for edge in MISSING_EDGES)
EXPECTED_POSITIVE_SUMMARY = (
    (1, 663, 663, 495, 6),
    (2, 219_453, 210_177, 83_946, 49),
    (3, 48_352_811, 42_603_687, 1_415_316, 143),
    (4, 7_978_213_815, 6_218_797_221, 49_476, 20),
)
EXPECTED_A5 = {
    "raw": 1_051_528_580_817,
    "source_a4_full_orbits": 62_234,
    "extension_paths": 41_012_206,
    "interval_extension_paths": 34_615_149,
    "shadow_feasible_extension_paths": 158_736,
    "orbit_count": 0,
}
EXPECTED_INDEPENDENT = {
    "active_mask_count": 13,
    "active_multisets": 6_188,
    "active_multisets_passing_active_caps": 1_000,
    "terminal_column_paths": 15_889_463,
    "distinct_row_paths": 15_128_612,
    "ground_row_paths": 15_128_612,
    "interval_paths": 12_797_343,
    "shadow_feasible_paths": 31_078,
    "viable_paths": 0,
}


def _is_missing_edge(left: int, right: int) -> bool:
    return tuple(sorted((left, right))) in EDGE_SET


def _independent(vertices: Iterable[int]) -> bool:
    chosen = tuple(vertices)
    return all(
        not _is_missing_edge(left, right)
        for left, right in combinations(chosen, 2)
    )


def _independent_sets(size: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        chosen for chosen in combinations(range(N), size) if _independent(chosen)
    )


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
    independent = {size: len(_independent_sets(size)) for size in range(2, 6)}
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
    rows: list[dict[str, object]] = []
    for missing_fours in range(independent[4] + 1):
        h4 = h4_without_fours - missing_fours
        global_g4 = h4 - h3
        if global_g4 < 1:
            break
        h = [1, h1, h2, h3, h4, h4, h3, h2, h1, 1]
        f = [
            sum(
                comb(9 - index, 9 - size) * h[index]
                for index in range(size + 1)
            )
            for size in range(1, 10)
        ]
        summed_link_g4 = (
            5 * f[4] - 24 * f[3] + 63 * f[2] - 112 * f[1] + 126 * N
        )
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
        candidate_triangles = sum(
            _independent(chosen) for chosen in combinations(link, 3)
        )
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


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"record-31 {name} is malformed")
    return value


def _source_provenance(path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path.read_bytes()).hexdigest(),
    }


def _certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _active_automorphism_order() -> int:
    active = (0, 1, 11, 12, 13)
    order = 0
    for images in permutations(active):
        move = dict(zip(active, images, strict=True))
        moved = frozenset(
            tuple(sorted((move[left], move[right]))) for left, right in EDGE_SET
        )
        order += moved == EDGE_SET
    return order


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_positive_a5_compact_branch_no_go":
        raise ValueError("wrong record-31 positive-a5 certificate kind")
    if document.get("graph6") != GRAPH6 or document.get("missing_edges") != [
        list(edge) for edge in MISSING_EDGES
    ]:
        raise ValueError("record-31 graph source changed")
    scope = _mapping(document.get("scope"), "scope")
    expected_scope = {
        "compact_only": True,
        "a5_compact_schedule_excluded": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    if scope != expected_scope:
        raise ValueError("record-31 a5 scope changed or was over-promoted")

    production = _mapping(document.get("production"), "production replay")
    if production.get("kind") != "h9_14_record31_positive_A4_scout":
        raise ValueError("record-31 production kind changed")
    schedules = derive_schedules()
    if production.get("schedules") != schedules:
        raise ValueError("record-31 Dehn--Sommerville schedules changed")
    if production.get("strict_compact_surviving_a") != [0, 1, 2, 3, 4, 5] or production.get(
        "strict_compact_excluded_a"
    ) != [6, 7]:
        raise ValueError("record-31 strict compact schedule split changed")
    local = derive_local_link_constraints()
    if production.get("local_link_constraints") != local:
        raise ValueError("record-31 exact local link table changed")
    if len(local["vertices"]) != 14 or len(local["edges"]) != 87 or len(
        local["triples"]
    ) != 319:
        raise ValueError("record-31 local table sizes changed")
    if production.get("independent_A4_count") != 771 or production.get(
        "interval_A4_count"
    ) != 663:
        raise ValueError("record-31 A4 ground counts changed")
    if production.get("active_automorphism_order") != _active_automorphism_order() != 2:
        raise ValueError("record-31 active automorphism order changed")

    raw_branches = production.get("branches")
    if not isinstance(raw_branches, list) or len(raw_branches) != 5:
        raise ValueError("record-31 positive branch table is malformed")
    branches = tuple(_mapping(row, "positive branch") for row in raw_branches)
    summary = tuple(
        (
            int(row["a"]),
            int(row["raw"]),
            int(row["interval_feasible"]),
            int(row["strict_aggregate_viable"]),
            int(row["orbit_count"]),
        )
        for row in branches[:4]
    )
    if summary != EXPECTED_POSITIVE_SUMMARY:
        raise ValueError("record-31 a=1,...,4 branch census changed")
    group_order = 2 * factorial(9)
    for row in branches[:4]:
        a = int(row["a"])
        if int(row["raw"]) != comb(663, a):
            raise ValueError(f"record-31 a={a} raw count changed")
        raw_orbits = row.get("orbits")
        if not isinstance(raw_orbits, list) or len(raw_orbits) != int(row["orbit_count"]):
            raise ValueError(f"record-31 a={a} orbit list changed")
        sizes = [int(_mapping(orbit, "orbit")["size"]) for orbit in raw_orbits]
        if any(size <= 0 or group_order % size for size in sizes):
            raise ValueError(f"record-31 a={a} orbit-stabilizer check failed")
        if sum(sizes) != int(row["strict_aggregate_viable"]):
            raise ValueError(f"record-31 a={a} orbit weights do not sum")

    a4 = branches[3]
    if (
        int(a4["active_row_multisets"]) != 1820
        or int(a4["ground_labelled_family_count"]) != comb(663, 4)
        or int(a4["interval_full_graph_orbits"]) != 62_234
        or int(a4["column_paths"])
        - int(a4["equal_active_duplicate_paths"])
        - int(a4["repeated_full_row_paths"])
        != int(a4["distinct_matrix_orbits"])
    ):
        raise ValueError("record-31 a=4 lossless quotient identity changed")

    a5 = branches[4]
    if int(a5.get("a", -1)) != 5 or any(
        int(a5.get(key, -1)) != value for key, value in EXPECTED_A5.items()
    ):
        raise ValueError("record-31 a=5 production closure changed")
    if a5.get("orbits") != [] or int(a5["extension_paths"]) != 62_234 * (663 - 4):
        raise ValueError("record-31 a=5 extension coverage identity failed")
    shadow_histogram = _mapping(a5.get("shadow_count_histogram"), "a5 shadow histogram")
    if sum(map(int, shadow_histogram.values())) != int(a5["interval_extension_paths"]):
        raise ValueError("record-31 a=5 shadow histogram does not cover interval paths")
    gate_histogram = a5.get("gate_histogram")
    if not isinstance(gate_histogram, list) or sum(
        int(_mapping(row, "a5 gate row")["path_count"]) for row in gate_histogram
    ) != int(a5["shadow_feasible_extension_paths"]):
        raise ValueError("record-31 a=5 gate histogram is incomplete")
    if any(
        _mapping(row, "a5 gate row").get("strict") is True
        or _mapping(row, "a5 gate row").get("aggregate") is True
        for row in gate_histogram
    ):
        raise ValueError("record-31 a=5 branch unexpectedly survives")

    independent = _mapping(document.get("independent"), "independent replay")
    if independent.get("kind") != "h9_14_record31_a5_independent_direct_matrix_audit":
        raise ValueError("record-31 independent replay kind changed")
    if any(
        int(independent.get(key, -1)) != value
        for key, value in EXPECTED_INDEPENDENT.items()
    ):
        raise ValueError("record-31 independent direct-matrix census changed")
    if independent.get("a5_excluded_by_necessary_gate") is not True:
        raise ValueError("record-31 independent replay does not close a=5")
    independent_gates = independent.get("gate_histogram")
    if not isinstance(independent_gates, list) or sum(
        int(_mapping(row, "independent gate row")["path_count"])
        for row in independent_gates
    ) != int(independent["interval_paths"]):
        raise ValueError("record-31 independent gate coverage changed")
    if any(
        _mapping(row, "independent gate row").get("strict") is True
        or _mapping(row, "independent gate row").get("aggregate") is True
        for row in independent_gates
    ):
        raise ValueError("record-31 independent replay has a survivor")

    provenance = document.get("source_provenance")
    expected_provenance = [
        _source_provenance(PRODUCTION_SOURCE),
        _source_provenance(INDEPENDENT_SOURCE),
    ]
    if provenance != expected_provenance:
        raise ValueError("record-31 a5 source provenance changed")
    if document.get("certificate_sha256") != _certificate_digest(document):
        raise ValueError("record-31 a5 certificate digest changed")
    return {
        "record_index": 31,
        "schedule_count": 8,
        "strict_schedule_split": {"surviving": [0, 1, 2, 3, 4, 5], "excluded": [6, 7]},
        "newly_closed_schedule": 5,
        "remaining_open_schedules": [0, 1, 2, 3, 4],
        "production_a5_surviving_orbits": 0,
        "independent_a5_surviving_paths": 0,
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def _run_cpp(source: Path, argument: str | None = None) -> dict[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("no C++20 compiler is available")
    with tempfile.TemporaryDirectory(prefix="record31-a5-") as directory:
        executable = Path(directory) / "audit"
        subprocess.run(
            [
                compiler,
                "-std=c++20",
                "-O3",
                "-Wall",
                "-Wextra",
                "-pedantic",
                str(source),
                "-o",
                str(executable),
            ],
            cwd=ROOT,
            check=True,
        )
        command = [str(executable)]
        if argument is not None:
            command.append(argument)
        completed = subprocess.run(
            command,
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    parsed = json.loads(completed.stdout)
    if not isinstance(parsed, dict):
        raise ValueError(f"{source.name} did not emit a JSON object")
    return parsed


def fresh_record() -> dict[str, object]:
    production = _run_cpp(PRODUCTION_SOURCE, "5")
    independent = _run_cpp(INDEPENDENT_SOURCE)
    document: dict[str, object] = {
        "kind": "h9_14_record31_positive_a5_compact_branch_no_go",
        "claim": (
            "For graph-shadow record 31, the compact Dehn--Sommerville "
            "schedule with exactly five missing four-faces is impossible."
        ),
        "graph6": GRAPH6,
        "missing_edges": [list(edge) for edge in MISSING_EDGES],
        "scope": {
            "compact_only": True,
            "a5_compact_schedule_excluded": True,
            "whole_shadow_compact_coxeter_theorem": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "nonclaims": [
            "Schedules a=0,1,2,3,4 remain open for other completions of record 31.",
            "The selected record-31 nerve has a separate compact-label no-go, but that does not close the other schedules.",
            "No parabolic cusp or Vinberg finite-volume endpoint criterion is audited here.",
        ],
        "production": production,
        "independent": independent,
        "source_provenance": [
            _source_provenance(PRODUCTION_SOURCE),
            _source_provenance(INDEPENDENT_SOURCE),
        ],
    }
    document["certificate_sha256"] = _certificate_digest(document)
    verify(document)
    return document


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true", help="rerun both C++ audits")
    parser.add_argument("--emit", action="store_true", help="emit the full certificate JSON")
    args = parser.parse_args(argv)
    document = fresh_record() if args.fresh else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    summary = verify(document)
    print(json.dumps(document if args.emit else summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
