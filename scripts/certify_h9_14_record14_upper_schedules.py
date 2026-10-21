#!/usr/bin/env python3
"""Freeze Record 14's complete compact a=3,4,5,6 schedule closure."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Mapping, Sequence

from scripts.certify_h9_14_record14_positive_A4_scout import (
    ARCHIVE as POSITIVE_ARCHIVE,
    derive_schedules,
    verify as verify_positive,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    build_problem,
    enumerate_models,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record14_upper_schedules_no_go.json"
CPP_WRAPPER = ROOT / "scripts/explore_h9_14_record14_addition_cover.cpp"
CPP_ENGINE = ROOT / "scripts/explore_h9_14_record18_addition_cover.cpp"
INDEPENDENT_SOURCE = ROOT / "scripts/explore_h9_14_record14_fixed_orbit_z3.py"
EXPECTED_SCHEDULE_ROWS = (
    (0, 8, 110, 42),
    (1, 14, 108, 36),
    (2, 20, 106, 30),
    (3, 26, 104, 24),
    (4, 32, 102, 18),
    (5, 38, 100, 12),
    (6, 44, 98, 6),
)
EXPECTED_CERTIFICATE_DIGEST = (
    "9392e15c3d1c446150febb215573a7ecbe1a0843a31e6ed541b0d9896bd7fa8c"
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _mask(face: Sequence[int]) -> int:
    return sum(1 << vertex for vertex in face)


def _schedule_rows() -> tuple[tuple[int, int, int, int], ...]:
    return tuple(
        (
            int(row["a"]),
            int(row["missing_independent_fives"]),
            int(row["f_vector"][-1]),
            int(row["summed_vertex_link_g4"]),
        )
        for row in derive_schedules()
    )


def _compile_cpp(directory: Path) -> Path:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("a C++20 compiler is required")
    executable = directory / "record14-addition-cover"
    subprocess.run(
        [
            compiler,
            "-std=c++20",
            "-O3",
            "-Wall",
            "-Wextra",
            "-pedantic",
            str(CPP_WRAPPER),
            "-o",
            str(executable),
        ],
        cwd=ROOT,
        check=True,
    )
    return executable


def fresh_record() -> dict[str, object]:
    positive = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(positive)
    if _schedule_rows() != EXPECTED_SCHEDULE_ROWS:
        raise AssertionError("Record-14 Dehn--Sommerville schedules changed")
    rows = []
    with tempfile.TemporaryDirectory(prefix="record14-upper-") as raw_directory:
        executable = _compile_cpp(Path(raw_directory))
        for a in (3, 4):
            branch = next(row for row in positive["branches"] if row["a"] == a)
            for orbit_index, orbit in enumerate(branch["orbits"]):
                representative = orbit["representative"]
                masks = [_mask(face) for face in representative]
                completed = subprocess.run(
                    [
                        str(executable),
                        str(a),
                        *map(str, masks),
                        "--fast-scheduled-gate",
                    ],
                    cwd=ROOT,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                production = json.loads(completed.stdout)
                problem = build_problem(a, orbit_index)
                independent = enumerate_models(problem)
                if production["representative_A4_masks"] != masks:
                    raise AssertionError("C++ source representative changed")
                if production["scheduled_pure_target_count"] != 0:
                    raise AssertionError("upper-schedule C++ target survived")
                if not independent["enumeration_complete"]:
                    raise AssertionError("upper-schedule independent replay cut off")
                if independent["literal_pure_target_count"] != 0:
                    raise AssertionError("upper-schedule Z3 target survived")
                rows.append(
                    {
                        "a": a,
                        "source_positive_A4_orbit_index": orbit_index,
                        "source_orbit_key": orbit["key"],
                        "source_orbit_size": orbit["size"],
                        "representative_A4s": representative,
                        "forced_A5_shadow_count": production["forced_A5_shadow_count"],
                        "additional_A5_count": production["additional_A5_count"],
                        "scheduled_facet_count": production["scheduled_facet_count"],
                        "candidate_count": production["candidate_count"],
                        "base_small_transversal_count": production["base_small_transversal_count"],
                        "production_search_node_count": production["search_node_count"],
                        "production_no_small_transversal_family_count": production[
                            "no_small_transversal_addition_count"
                        ],
                        "production_scheduled_order_five_model_count": sum(
                            int(count)
                            for facets, count in production[
                                "order_five_transversal_count_histogram"
                            ].items()
                            if int(facets) == production["scheduled_facet_count"]
                        ),
                        "production_literal_pure_target_count": 0,
                        "independent_base_order_five_transversal_count": len(
                            problem["fives"]
                        ),
                        "independent_six_set_constraint_count": independent[
                            "six_set_constraint_count"
                        ],
                        "independent_z3_model_count": independent["z3_model_count"],
                        "independent_literal_pure_target_count": independent[
                            "literal_pure_target_count"
                        ],
                        "independent_enumeration_complete": independent[
                            "enumeration_complete"
                        ],
                    }
                )
    evidence = {
        "schedule_rows": [list(row) for row in EXPECTED_SCHEDULE_ROWS],
        "strict_vertex_sum_closed_schedules": [5, 6],
        "enumerated_schedules": [3, 4],
        "source_orbit_count_by_schedule": {"3": 59, "4": 13},
        "source_orbit_count": len(rows),
        "rows": rows,
        "all_rows_zero_production_targets": all(
            row["production_literal_pure_target_count"] == 0 for row in rows
        ),
        "all_rows_zero_independent_targets": all(
            row["independent_literal_pure_target_count"] == 0 for row in rows
        ),
        "production_search_node_total": sum(
            row["production_search_node_count"] for row in rows
        ),
        "production_no_small_transversal_family_total": sum(
            row["production_no_small_transversal_family_count"] for row in rows
        ),
        "independent_z3_model_total": sum(
            row["independent_z3_model_count"] for row in rows
        ),
    }
    document = {
        "kind": "h9_14_record14_upper_schedules_compact_incidence_no_go",
        "record_index": 14,
        "graph6": positive["graph6"],
        "source_positive_A4_archive": {
            "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
            "sha256": _sha(POSITIVE_ARCHIVE),
        },
        "source_files": {
            str(path.relative_to(ROOT)): _sha(path)
            for path in (CPP_WRAPPER, CPP_ENGINE, INDEPENDENT_SOURCE)
        },
        "evidence": evidence,
        "theorem": {
            "compact_record14_schedules_closed": [3, 4, 5, 6],
            "remaining_compact_schedules": [0, 1, 2],
            "whole_shadow_compact_coxeter_no_go": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "statement": (
                "Every compact Record-14 completion in schedules a=3,4,5,6 "
                "is impossible. Schedules a=5,6 fail strict summed vertex-link "
                "g4 positivity; all 72 source orbits in a=3,4 have zero exact "
                "incidence target in two lossless enumerators."
            ),
        },
    }
    document["certificate_digest"] = _digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record14_upper_schedules_compact_incidence_no_go":
        raise ValueError("wrong Record-14 upper-schedule certificate kind")
    if _schedule_rows() != EXPECTED_SCHEDULE_ROWS:
        raise ValueError("Record-14 schedule derivation changed")
    positive = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(positive)
    source = document.get("source_positive_A4_archive")
    if not isinstance(source, Mapping) or source.get("sha256") != _sha(POSITIVE_ARCHIVE):
        raise ValueError("Record-14 positive-A4 source changed")
    sources = document.get("source_files")
    if not isinstance(sources, Mapping):
        raise ValueError("Record-14 upper source map is malformed")
    for path in (CPP_WRAPPER, CPP_ENGINE, INDEPENDENT_SOURCE):
        if sources.get(str(path.relative_to(ROOT))) != _sha(path):
            raise ValueError(f"Record-14 source hash changed: {path.name}")
    evidence = document.get("evidence")
    if not isinstance(evidence, Mapping):
        raise ValueError("Record-14 evidence is malformed")
    rows = evidence.get("rows")
    if not isinstance(rows, list) or len(rows) != 72:
        raise ValueError("Record-14 upper orbit coverage changed")
    expected_sources = []
    for a in (3, 4):
        branch = next(row for row in positive["branches"] if row["a"] == a)
        expected_sources.extend(
            (a, index, orbit["key"], orbit["size"], orbit["representative"])
            for index, orbit in enumerate(branch["orbits"])
        )
    observed_sources = [
        (
            row["a"],
            row["source_positive_A4_orbit_index"],
            row["source_orbit_key"],
            row["source_orbit_size"],
            row["representative_A4s"],
        )
        for row in rows
    ]
    if observed_sources != expected_sources:
        raise ValueError("Record-14 upper source-orbit list is not literal")
    if any(
        not row["independent_enumeration_complete"]
        or row["production_literal_pure_target_count"] != 0
        or row["independent_literal_pure_target_count"] != 0
        for row in rows
    ):
        raise ValueError("Record-14 upper schedule has a residual target")
    theorem = document.get("theorem")
    if not isinstance(theorem, Mapping) or theorem.get(
        "compact_record14_schedules_closed"
    ) != [3, 4, 5, 6] or theorem.get("remaining_compact_schedules") != [0, 1, 2]:
        raise ValueError("Record-14 upper theorem scope changed")
    stored_digest = document.get("certificate_digest")
    body = dict(document)
    body.pop("certificate_digest", None)
    digest = _digest(body)
    if stored_digest != digest:
        raise ValueError("Record-14 upper certificate digest mismatch")
    if EXPECTED_CERTIFICATE_DIGEST not in ("TO_BE_FROZEN", digest):
        raise ValueError("Record-14 upper certificate is not frozen")
    return {
        "record_index": 14,
        "closed_schedules": [3, 4, 5, 6],
        "enumerated_source_orbit_count": 72,
        "production_target_count": 0,
        "independent_target_count": 0,
        "remaining_schedules": [0, 1, 2],
        "whole_shadow_compact_coxeter_no_go": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", type=Path)
    args = parser.parse_args(argv)
    document = fresh_record() if args.fresh else json.loads(OUTPUT.read_text(encoding="utf-8"))
    result = verify(document)
    if args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        result["written"] = str(args.write)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
