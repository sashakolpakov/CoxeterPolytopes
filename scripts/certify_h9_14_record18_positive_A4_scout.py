#!/usr/bin/env python3
"""Verify or freshly replay the record-18 positive-A4 production scout.

This is intentionally not a theorem certificate.  It verifies the lossless
local/forced-shadow census and its exact orbit weights; an independent
implementation remains mandatory before any whole-shadow promotion.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from math import comb
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record18_positive_A4_scout.json"
SOURCE = ROOT / "scripts/explore_h9_14_record18_positive_orbits.cpp"

EXPECTED_SCHEDULES = (
    (0, 8, 110, 42),
    (1, 14, 108, 36),
    (2, 20, 106, 30),
    (3, 26, 104, 24),
    (4, 32, 102, 18),
    (5, 38, 100, 12),
    (6, 44, 98, 6),
)
EXPECTED_BRANCH_SUMMARY = (
    (1, 495, 495, 145, 3),
    (2, 122_265, 120_969, 9_144, 14),
    (3, 20_092_215, 19_449_199, 188_944, 40),
    (4, 2_471_342_445, 2_310_730_275, 14_812, 13),
)
EXPECTED_A4_QUOTIENT = {
    "active_row_multisets": 17_550,
    "column_paths": 2_228_923,
    "equal_active_duplicate_paths": 448_772,
    "repeated_full_row_paths": 90_666,
    "s8_orbits": 1_689_485,
    "interval_s8_orbits": 280_869,
    "viable_s8_orbits": 54,
}


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"record-18 positive scout {name} is malformed")
    return value


def _branches(document: Mapping[str, object]) -> tuple[Mapping[str, object], ...]:
    value = document.get("branches")
    if not isinstance(value, list):
        raise ValueError("record-18 positive scout branches are malformed")
    return tuple(_mapping(row, "branch") for row in value)


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record18_positive_A4_scout":
        raise ValueError("wrong record-18 positive scout kind")
    if (
        document.get("theorem_certificate") is not False
        or document.get("independent_replay_complete") is not False
        or document.get("whole_shadow_compact_coxeter_theorem") is not False
        or document.get("compact_only") is not True
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("record-18 positive scout scope was over-promoted")
    if tuple(map(tuple, document.get("schedules", ()))) != EXPECTED_SCHEDULES:
        raise ValueError("record-18 schedule table changed")
    if document.get("strict_compact_surviving_a") != [0, 1, 2, 3, 4] or (
        document.get("strict_compact_excluded_a") != [5, 6]
    ):
        raise ValueError("record-18 strict schedule cut changed")
    if document.get("independent_A4_count") != 735 or (
        document.get("interval_A4_count") != 495
        or document.get("active_automorphism_order") != 8
    ):
        raise ValueError("record-18 A4 ground set or symmetry changed")

    branches = _branches(document)
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
        raise ValueError("record-18 positive branch summary changed")
    for row in branches:
        a = int(row["a"])
        if int(row["raw"]) != comb(495, a):
            raise ValueError(f"record-18 a={a} raw family count is inconsistent")
        orbits = row.get("orbits")
        if not isinstance(orbits, list) or len(orbits) != int(row["orbit_count"]):
            raise ValueError(f"record-18 a={a} orbit list is malformed")
        keys = []
        sizes = []
        for raw_orbit in orbits:
            orbit = _mapping(raw_orbit, f"a={a} orbit")
            key = tuple(map(int, orbit["key"]))
            representative = orbit["representative"]
            if not isinstance(representative, list) or len(representative) != a:
                raise ValueError(f"record-18 a={a} representative is malformed")
            if any(
                not isinstance(face, list)
                or len(face) != 4
                or len(set(map(int, face))) != 4
                for face in representative
            ):
                raise ValueError(f"record-18 a={a} representative A4 changed")
            keys.append(key)
            sizes.append(int(orbit["size"]))
        if len(set(keys)) != len(keys):
            raise ValueError(f"record-18 a={a} duplicate orbit key")
        if sum(sizes) != int(row["strict_aggregate_viable"]):
            raise ValueError(f"record-18 a={a} orbit weights do not sum")
        expected_histogram = Counter(map(int, sizes))
        archived_histogram = {
            int(size): int(count)
            for size, count in _mapping(
                row.get("orbit_size_histogram"), f"a={a} orbit histogram"
            ).items()
        }
        if expected_histogram != archived_histogram:
            raise ValueError(f"record-18 a={a} orbit histogram changed")
        gates = row.get("gate_histogram")
        if not isinstance(gates, list) or sum(int(gate["count"]) for gate in gates) != int(
            row["interval_feasible"]
        ):
            raise ValueError(f"record-18 a={a} shadow-gate weights do not sum")
        if sum(
            int(gate["count"])
            for gate in gates
            if gate["strict"] is True and gate["aggregate"] is True
        ) != int(row["strict_aggregate_viable"]):
            raise ValueError(f"record-18 a={a} viable gate weights do not sum")

    a4 = branches[3]
    if any(int(a4[key]) != value for key, value in EXPECTED_A4_QUOTIENT.items()):
        raise ValueError("record-18 a=4 quotient transcript changed")
    if int(a4["column_paths"]) - int(a4["equal_active_duplicate_paths"]) - int(
        a4["repeated_full_row_paths"]
    ) != int(a4["s8_orbits"]):
        raise ValueError("record-18 a=4 path quotient identity failed")
    full_group_order = 8 * 40_320
    if any(full_group_order % int(orbit["size"]) for orbit in a4["orbits"]):
        raise ValueError("record-18 a=4 orbit size does not divide the graph group")

    return {
        "record_index": 18,
        "production_positive_census_complete": True,
        "positive_orbit_count": sum(int(row["orbit_count"]) for row in branches),
        "a4_quotient_weighting_verified": True,
        "independent_replay_complete": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def fresh_replay() -> dict[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("no C++20 compiler is available")
    with tempfile.TemporaryDirectory(prefix="record18-positive-") as directory:
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
            check=True,
            cwd=ROOT,
        )
        completed = subprocess.run(
            [str(executable), "4"],
            check=True,
            cwd=ROOT,
            text=True,
            capture_output=True,
        )
    fresh = json.loads(completed.stdout)
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    for key in (
        "kind",
        "compact_only",
        "independent_A4_count",
        "interval_A4_count",
        "active_automorphism_order",
        "branches",
        "scope",
    ):
        if fresh.get(key) != archived.get(key):
            raise AssertionError(f"fresh record-18 replay differs at {key}")
    return {
        "fresh_source_replay_matches_archive": True,
        "compiler": Path(compiler).name,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args(argv)
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    result = verify(document)
    if args.fresh:
        result["fresh_replay"] = fresh_replay()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
