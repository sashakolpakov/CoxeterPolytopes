#!/usr/bin/env python3
"""Freeze the exact record-30 positive-A4 schedules through a=5."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Mapping


ROOT = Path(__file__).resolve().parents[1]
A1_A4_SOURCE = ROOT / "scripts/explore_h9_14_record18_positive_orbits.cpp"
A5_SOURCE = ROOT / "scripts/independent_h9_14_record30_a5.cpp"
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_positive_scout.json"
A1_A4_SHA256 = "1142408977f20f442d5178f132b02bfe5ede77434d0c519f251e5b5de2339dbd"
A5_SHA256 = "fe92bcdc6020469b2f5c346905551a3657fc2a285aa66ecf7e2f5b345e217922"
EXPECTED_CERTIFICATE_SHA256 = "02a2238457b33d0cba59bb9fb85579a5e27c4520ea86f17a4ce0614fcabb111e"
EXPECTED_BRANCHES = (
    (1, 742, 742, 742, 3),
    (2, 274_911, 270_375, 97_699, 21),
    (3, 67_811_380, 64_513_204, 524_552, 32),
    (4, 12_528_152_455, 11_245_592_071, 28_896, 7),
)


def canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _compile_run(source: Path, flags: tuple[str, ...], argument: str | None = None) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="record30_positive_") as temporary:
        executable = Path(temporary) / "audit"
        subprocess.run(
            ["clang++", "-std=c++20", "-O3", "-pthread", "-Wall", "-Wextra", "-Werror", *flags, str(source), "-o", str(executable)],
            check=True,
            cwd=ROOT,
        )
        command = [str(executable)] + ([] if argument is None else [argument])
        completed = subprocess.run(command, check=True, cwd=ROOT, text=True, capture_output=True)
    return json.loads(completed.stdout)


def _schedule_table() -> list[dict[str, object]]:
    edges = {(0, 11), (1, 12), (2, 13)}
    independent = {}
    for size in range(2, 10):
        independent[size] = sum(
            not any(set(edge) <= set(chosen) for edge in edges)
            for chosen in combinations(range(14), size)
        )
    result = []
    vectors = (
        (14, 88, 328, 806, 1364, 1596, 1248, 594, 132),
        (14, 88, 328, 805, 1358, 1582, 1232, 585, 130),
        (14, 88, 328, 804, 1352, 1568, 1216, 576, 128),
        (14, 88, 328, 803, 1346, 1554, 1200, 567, 126),
        (14, 88, 328, 802, 1340, 1540, 1184, 558, 124),
        (14, 88, 328, 801, 1334, 1526, 1168, 549, 122),
        (14, 88, 328, 800, 1328, 1512, 1152, 540, 120),
        (14, 88, 328, 799, 1322, 1498, 1136, 531, 118),
    )
    if (independent[4], independent[5]) != (806, 1372):
        raise AssertionError("record-30 independent-set census changed")
    for a, vector in enumerate(vectors):
        result.append({
            "a": a,
            "b": 8 + 6 * a,
            "facets": vector[-1],
            "global_g4": 8 - a,
            "summed_vertex_link_g4": 6 * (8 - a),
            "strict_vertex_sum_feasible": 6 * (8 - a) >= 14,
            "f_vector": list(vector),
        })
    return result


def verify(document: Mapping[str, object], *, allow_unfrozen: bool = False) -> dict[str, object]:
    if sha256(A1_A4_SOURCE.read_bytes()).hexdigest() != A1_A4_SHA256:
        raise ValueError("record-30 a1-a4 source digest changed")
    if sha256(A5_SOURCE.read_bytes()).hexdigest() != A5_SHA256:
        raise ValueError("record-30 a5 source digest changed")
    if document.get("kind") != "h9_14_record30_positive_scout":
        raise ValueError("wrong record-30 positive scout kind")
    if document.get("schedule_table") != _schedule_table():
        raise ValueError("record-30 schedule table changed")
    a1_a4 = document.get("a1_a4")
    if not isinstance(a1_a4, Mapping):
        raise ValueError("record-30 a1-a4 transcript missing")
    branches = a1_a4.get("branches")
    if not isinstance(branches, list):
        raise ValueError("record-30 positive branches malformed")
    summary = tuple(
        (row["a"], row["raw"], row["interval_feasible"], row["strict_aggregate_viable"], row["orbit_count"])
        for row in branches
    )
    if summary != EXPECTED_BRANCHES:
        raise ValueError("record-30 a1-a4 census changed")
    a5 = document.get("a5")
    if not isinstance(a5, Mapping) or (
        a5.get("raw") != 1_849_155_302_358
        or a5.get("interval_feasible") != 1_510_284_730_038
        or a5.get("strict_aggregate_viable") != 672
        or a5.get("full_orbit_count") != 1
        or a5.get("a5_compact_schedule_excluded") is not True
        or a5.get("orbits", [{}])[0].get("residual_link_witness")
        != {"subset": [0, 1], "h": [1, 3, 4, 4], "g": [2, 1, 0]}
    ):
        raise ValueError("record-30 a5 quotient or obstruction changed")
    if document.get("strict_sum_excluded_schedules") != [6, 7]:
        raise ValueError("record-30 strict-sum schedule cut changed")
    if (
        document.get("whole_shadow_compact_coxeter_excluded") is not False
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
        or document.get("compact_only") is not True
    ):
        raise ValueError("record-30 positive scout scope was over-promoted")
    critical = dict(document)
    digest = critical.pop("certificate_sha256", None)
    rebuilt = canonical_sha(critical)
    if digest != rebuilt:
        raise ValueError("record-30 positive scout internal digest invalid")
    if not allow_unfrozen and EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-30 positive scout certificate digest changed")
    return {
        "positive_orbit_counts": [3, 21, 32, 7, 1],
        "a5_compact_schedule_excluded": True,
        "strict_sum_excluded_schedules": [6, 7],
        "whole_shadow_compact_coxeter_excluded": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def fresh_record() -> dict[str, object]:
    document = {
        "kind": "h9_14_record30_positive_scout",
        "compact_only": True,
        "graph6": "M?????????O?O?G??",
        "missing_edges": [[0, 11], [1, 12], [2, 13]],
        "source_sha256": {"a1_a4": A1_A4_SHA256, "a5": A5_SHA256},
        "schedule_table": _schedule_table(),
        "a1_a4": _compile_run(A1_A4_SOURCE, ("-DH9_RECORD30",), "4"),
        "a5": _compile_run(A5_SOURCE, ()),
        "strict_sum_excluded_schedules": [6, 7],
        "whole_shadow_compact_coxeter_excluded": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "scope": (
            "Exact compact-only positive-A4 incidence quotient. Schedule a=5 "
            "is excluded by a forced residual edge-link g-vector; a=6,7 by "
            "the strict vertex-link sum. No whole-shadow or finite-volume "
            "noncompact conclusion is made."
        ),
    }
    document["certificate_sha256"] = canonical_sha(document)
    verify(document, allow_unfrozen=True)
    return document


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    if args.refresh:
        document = fresh_record()
        ARCHIVE.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    else:
        document = json.loads(ARCHIVE.read_text())
        verify(document)
        if args.recompute and fresh_record() != document:
            raise ValueError("fresh record-30 positive replay differs from archive")
    print(json.dumps(verify(document, allow_unfrozen=True), sort_keys=True))


if __name__ == "__main__":
    main()
