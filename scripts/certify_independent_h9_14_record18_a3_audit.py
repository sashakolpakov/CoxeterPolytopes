#!/usr/bin/env python3
"""Freeze the standalone all-orbit record-18 a=3 addition audit."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
POSITIVE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_positive_orbit_audit.cpp"
)
ADDITION_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_addition_cover_audit.cpp"
)
OUTPUT = ROOT / "data/enumerations/h9_14_record18_a3_independent_audit.json"

EXPECTED_ADDITION_ORBIT_HISTOGRAM = {"2": 29, "3": 10, "4": 1}
EXPECTED_AGGREGATE = {
    "raw_unordered_addition_count": 14_933_639_143,
    "vertex_interval_feasible_addition_count": 2_392_528_542,
    "no_small_transversal_addition_count": 6_011,
    "scheduled_pure_104_facet_target_count": 0,
}


def _compile(compiler: str, source: Path, target: Path) -> None:
    subprocess.run(
        [
            compiler,
            "-std=c++20",
            "-O3",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            str(source),
            "-o",
            str(target),
        ],
        cwd=ROOT,
        check=True,
    )


def _run_json(command: Sequence[str]) -> Mapping[str, object]:
    completed = subprocess.run(
        list(command),
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    value = json.loads(completed.stdout)
    if not isinstance(value, Mapping):
        raise ValueError("independent record-18 executable returned malformed JSON")
    return value


def _reconstruct(key: Sequence[int], a: int) -> tuple[int, ...]:
    if len(key) != a + 8:
        raise ValueError("positive orbit key has the wrong length")
    rows = []
    for row in range(a):
        value = int(key[row])
        for isolate, pattern in enumerate(key[a:], start=4):
            if (int(pattern) >> row) & 1:
                value |= 1 << isolate
        if value.bit_count() != 4:
            raise AssertionError("positive key did not reconstruct an A4")
        rows.append(value)
    if len(set(rows)) != a:
        raise AssertionError("positive key reconstructed duplicate A4 rows")
    return tuple(rows)


def build() -> dict[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("no C++20 compiler is available")
    with tempfile.TemporaryDirectory(prefix="record18-a3-independent-") as directory:
        positive_executable = Path(directory) / "positive"
        addition_executable = Path(directory) / "addition"
        _compile(compiler, POSITIVE_SOURCE, positive_executable)
        _compile(compiler, ADDITION_SOURCE, addition_executable)
        positive = _run_json([str(positive_executable)])
        branches = positive.get("branches")
        if not isinstance(branches, list) or len(branches) != 4:
            raise ValueError("independent positive branch output is malformed")
        a3 = branches[2]
        if not isinstance(a3, Mapping) or a3.get("a") != 3:
            raise ValueError("independent positive a=3 branch is malformed")
        raw_orbits = a3.get("viable_orbit_keys")
        if not isinstance(raw_orbits, list) or len(raw_orbits) != 40:
            raise ValueError("independent positive a=3 orbit list is malformed")

        rows = []
        for index, raw_orbit in enumerate(raw_orbits):
            if not isinstance(raw_orbit, Mapping):
                raise ValueError("independent positive a=3 orbit row is malformed")
            key = tuple(map(int, raw_orbit["key"]))
            fours = _reconstruct(key, 3)
            result = dict(
                _run_json(
                    [str(addition_executable), "3", *(str(face) for face in fours)]
                )
            )
            if result.get("representative_A4_masks") != list(fours):
                raise AssertionError("a=3 addition replay changed its input A4s")
            result["positive_orbit_index"] = index
            result["positive_orbit_key"] = list(key)
            result["positive_orbit_literal_size"] = int(raw_orbit["literal_size"])
            rows.append(result)
            print(
                f"a3 orbit {index + 1}/40: add={result['additional_A5_count']} "
                f"no-small={result['no_small_transversal_addition_count']} "
                f"pure104={result['scheduled_pure_target_count']}",
                file=sys.stderr,
                flush=True,
            )

    addition_histogram = {
        str(additions): sum(
            int(row["additional_A5_count"]) == additions for row in rows
        )
        for additions in (2, 3, 4)
    }
    aggregate = {
        "raw_unordered_addition_count": sum(
            int(row["raw_unordered_addition_count"]) for row in rows
        ),
        "vertex_interval_feasible_addition_count": sum(
            int(row["vertex_interval_feasible_addition_count"]) for row in rows
        ),
        "no_small_transversal_addition_count": sum(
            int(row["no_small_transversal_addition_count"]) for row in rows
        ),
        "scheduled_pure_104_facet_target_count": sum(
            int(row["scheduled_pure_target_count"]) for row in rows
        ),
    }
    if addition_histogram != EXPECTED_ADDITION_ORBIT_HISTOGRAM:
        raise AssertionError("independent a=3 addition-count histogram changed")
    if aggregate != EXPECTED_AGGREGATE:
        raise AssertionError(
            f"independent a=3 aggregate differs: {aggregate!r}"
        )
    return {
        "schema_version": "1.0.0",
        "kind": "independent_h9_14_record18_a3_addition_audit",
        "claim": (
            "A standalone coverage-descending earliest-pivot search and "
            "three-graph-cover blocker dualization independently exclude all "
            "forty record-18 a=3 positive-A4 orbits at a necessary incidence gate."
        ),
        "backend": (
            "warning-clean standalone C++20; no project headers, Python project "
            "modules, production representatives, or production addition tables"
        ),
        "input_orbit_source": (
            "fresh independent incidence-matrix positive-orbit executable"
        ),
        "input_positive_orbit_count": 40,
        "addition_orbit_histogram": addition_histogram,
        "aggregate": aggregate,
        "orbit_audits": rows,
        "a3_schedule_incidence_empty": True,
        "compact_only": True,
        "theorem_certificate": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "remaining_unfinished_schedules": [0, 1, 2],
        "scope": (
            "Necessary compact-incidence result only. The a=0,1,2 schedules "
            "remain open here; no finite-volume noncompact conclusion follows."
        ),
    }


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "independent_h9_14_record18_a3_addition_audit":
        raise ValueError("wrong independent record-18 a=3 audit kind")
    if (
        document.get("compact_only") is not True
        or document.get("theorem_certificate") is not False
        or document.get("whole_shadow_compact_coxeter_theorem") is not False
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
        or document.get("remaining_unfinished_schedules") != [0, 1, 2]
    ):
        raise ValueError("independent record-18 a=3 scope was over-promoted")
    if document.get("input_positive_orbit_count") != 40 or (
        document.get("addition_orbit_histogram") != EXPECTED_ADDITION_ORBIT_HISTOGRAM
    ):
        raise ValueError("independent record-18 a=3 input orbit census changed")
    if document.get("aggregate") != EXPECTED_AGGREGATE:
        raise ValueError("independent record-18 a=3 aggregate changed")
    rows = document.get("orbit_audits")
    if not isinstance(rows, list) or len(rows) != 40:
        raise ValueError("independent record-18 a=3 orbit transcript is malformed")
    if len({tuple(map(int, row["positive_orbit_key"])) for row in rows}) != 40:
        raise ValueError("independent record-18 a=3 orbit keys are not unique")
    if any(int(row["scheduled_pure_target_count"]) != 0 for row in rows):
        raise ValueError("independent record-18 a=3 pure target survived")
    if document.get("a3_schedule_incidence_empty") is not True:
        raise ValueError("independent record-18 a=3 conclusion changed")
    return {
        "input_positive_orbit_count": 40,
        "complete_addition_transcript_replayed": True,
        "no_small_transversal_addition_count": 6_011,
        "scheduled_pure_104_facet_target_count": 0,
        "a3_schedule_incidence_empty": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args(argv)
    if args.write:
        document = build()
        OUTPUT.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    else:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    result = verify(document)
    if args.fresh and not args.write:
        rebuilt = build()
        if rebuilt != document:
            raise AssertionError("fresh independent a=3 audit differs from archive")
        result["fresh_no_project_replay_matches_archive"] = True
    print(json.dumps({"output": str(OUTPUT), **result}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
