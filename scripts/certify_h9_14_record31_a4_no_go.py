#!/usr/bin/env python3
"""Freeze the compact-only exclusion of Record 31's entire ``a=4`` schedule."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_h9_14_record31_positive_a5 import (
    ARCHIVE as POSITIVE_ARCHIVE,
    verify as verify_positive,
)

LOW_SOURCE = ROOT / "scripts/explore_h9_14_record31_a4_additions.py"
HIGH_SOURCE = ROOT / "scripts/explore_h9_14_record31_a4_add3.cpp"
INDEPENDENT_SOURCE = ROOT / "scripts/independent_h9_14_record31_a4_additions_audit.cpp"
ARCHIVE = ROOT / "data/enumerations/h9_14_record31_a4_compact_no_go.json"
EXPECTED_CANDIDATES = (
    6, 21, 764, 56, 126, 246, 56, 56, 252, 449,
    126, 126, 246, 56, 126, 56, 126, 252, 462, 462,
)
EXPECTED_SMALL = (
    22, 21, 43, 42, 36, 30, 36, 36, 36, 30,
    27, 36, 30, 36, 30, 36, 42, 35, 43, 43,
)
EXPECTED_CAP = (
    6, 21, 3_951_530, 945, 2_310, 8_505, 540, 945, 8_541, 27_693,
    126, 2_310, 8_505, 540, 2_310, 540, 2_310, 8_541, 1_794_605, 2_819_841,
)
EXPECTED_BLOCKERS = {
    18: {"105": 150, "109": 120},
    19: {
        "94": 160,
        "95": 120,
        "96": 340,
        "99": 80,
        "100": 260,
        "102": 150,
        "104": 10,
        "115": 1,
    },
}


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Record 31 a=4 {name} is malformed")
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


def _positive_document() -> dict[str, object]:
    document = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a4_compact_schedule_no_go":
        raise ValueError("wrong Record 31 a=4 certificate kind")
    expected_scope = {
        "compact_only": True,
        "a4_compact_schedule_excluded": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    if _mapping(document.get("scope"), "scope") != expected_scope:
        raise ValueError("Record 31 a=4 scope changed or was over-promoted")

    positive = _positive_document()
    upstream = _mapping(document.get("upstream"), "upstream")
    if upstream != {
        "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
        "certificate_sha256": positive["certificate_sha256"],
        "a4_orbit_count": 20,
    }:
        raise ValueError("Record 31 a=4 upstream orbit certificate changed")
    upstream_orbits = positive["production"]["branches"][3]["orbits"]

    low = _mapping(document.get("production_low_additions"), "low-addition audit")
    if (
        low.get("kind") != "h9_14_record31_a4_low_addition_scout"
        or low.get("compact_only") is not True
        or int(low.get("orbit_representative_count", -1)) != 20
        or low.get("deferred_three_addition_indices") != [2, 18, 19]
        or int(low.get("literal_low_addition_target_count", -1)) != 0
    ):
        raise ValueError("Record 31 low-addition production summary changed")
    low_rows = low.get("representative_audits")
    if not isinstance(low_rows, list) or len(low_rows) != 20:
        raise ValueError("Record 31 low-addition rows changed")
    for index, raw in enumerate(low_rows):
        row = _mapping(raw, "low-addition row")
        expected_family = [
            sum(1 << vertex for vertex in face)
            for face in upstream_orbits[index]["representative"]
        ]
        if (
            int(row["orbit_representative_index"]) != index
            or row["missing_four_masks"] != expected_family
            or int(row["individual_candidate_count"]) != EXPECTED_CANDIDATES[index]
            or int(row["base_small_transversal_count"]) != EXPECTED_SMALL[index]
        ):
            raise ValueError(f"Record 31 low-addition row {index} changed")
        if index in (2, 18, 19):
            if row.get("deferred_to_pivoted_cpp") is not True:
                raise ValueError("three-addition row was not deferred")
        elif (
            row.get("deferred_to_pivoted_cpp") is not False
            or int(row["vertex_cap_feasible_count"]) != EXPECTED_CAP[index]
            or int(row["no_transversal_of_order_at_most_four_count"]) != 0
            or row["facet_count_histogram"] != {}
            or int(row["scheduled_target_count"]) != 0
        ):
            raise ValueError(f"Record 31 low-addition exclusion {index} changed")

    high = _mapping(document.get("production_three_additions"), "three-addition audit")
    if high.get("kind") != "h9_14_record31_a4_add3_pivot_scout" or high.get(
        "compact_only"
    ) is not True:
        raise ValueError("Record 31 pivoted production kind changed")
    high_rows_raw = high.get("representative_audits")
    if not isinstance(high_rows_raw, list) or len(high_rows_raw) != 3:
        raise ValueError("Record 31 pivoted production rows changed")
    high_rows = {
        int(_mapping(row, "three-addition row")["orbit_representative_index"]):
        _mapping(row, "three-addition row")
        for row in high_rows_raw
    }
    if set(high_rows) != {2, 18, 19}:
        raise ValueError("Record 31 pivoted production index coverage changed")
    expected_terminal = {
        2: (0, {}),
        18: (270, EXPECTED_BLOCKERS[18]),
        19: (1121, EXPECTED_BLOCKERS[19]),
    }
    for index, (no_small, histogram) in expected_terminal.items():
        row = high_rows[index]
        if (
            int(row["individual_candidate_count"]) != EXPECTED_CANDIDATES[index]
            or int(row["base_small_transversal_count"]) != EXPECTED_SMALL[index]
            or int(row["no_small_transversal_addition_count"]) != no_small
            or row["facet_count_histogram"] != histogram
            or int(row["scheduled_114_facet_target_count"]) != 0
        ):
            raise ValueError(f"Record 31 pivoted terminal row {index} changed")

    independent = _mapping(document.get("independent"), "independent audit")
    if (
        independent.get("kind") != "h9_14_record31_a4_independent_literal_audit"
        or independent.get("compact_only") is not True
        or int(independent.get("scheduled_target_count_across_representatives", -1)) != 0
    ):
        raise ValueError("Record 31 independent summary changed")
    independent_rows = independent.get("representative_audits")
    if not isinstance(independent_rows, list) or len(independent_rows) != 20:
        raise ValueError("Record 31 independent row coverage changed")
    for index, raw in enumerate(independent_rows):
        row = _mapping(raw, "independent row")
        expected_no_small = 270 if index == 18 else 1121 if index == 19 else 0
        expected_histogram = EXPECTED_BLOCKERS.get(index, {})
        if (
            int(row["orbit_representative_index"]) != index
            or int(row["individual_candidate_count"]) != EXPECTED_CANDIDATES[index]
            or int(row["base_small_transversal_count"]) != EXPECTED_SMALL[index]
            or int(row["vertex_cap_feasible_count"]) != EXPECTED_CAP[index]
            or int(row["no_transversal_of_order_at_most_four_count"]) != expected_no_small
            or row["blocker_count_histogram"] != expected_histogram
            or int(row["scheduled_114_facet_target_count"]) != 0
        ):
            raise ValueError(f"Record 31 independent row {index} changed")
        if index not in (2, 18, 19):
            production = _mapping(low_rows[index], "low-addition comparison row")
            if int(production["vertex_cap_feasible_count"]) != int(
                row["vertex_cap_feasible_count"]
            ):
                raise ValueError("low-addition production/independent mismatch")
        else:
            production = high_rows[index]
            if int(production["no_small_transversal_addition_count"]) != int(
                row["no_transversal_of_order_at_most_four_count"]
            ):
                raise ValueError("three-addition production/independent mismatch")

    downstream = _mapping(document.get("downstream"), "downstream")
    if downstream != {
        "scheduled_114_facet_target_count": 0,
        "ridge_topology_target_count": 0,
        "unrestricted_compact_label_target_count": 0,
        "reason": "Every completion fails before the ridge/topology and compact-label gates: it either retains a transversal of order at most four or has blocker count different from the scheduled 114 facets.",
    }:
        raise ValueError("Record 31 a=4 downstream-vacuity statement changed")

    expected_provenance = [
        _source_provenance(LOW_SOURCE),
        _source_provenance(HIGH_SOURCE),
        _source_provenance(INDEPENDENT_SOURCE),
    ]
    if document.get("source_provenance") != expected_provenance:
        raise ValueError("Record 31 a=4 source provenance changed")
    if document.get("certificate_sha256") != _certificate_digest(document):
        raise ValueError("Record 31 a=4 certificate digest changed")
    return {
        "record_index": 31,
        "newly_closed_schedule": 4,
        "a4_orbit_count": 20,
        "scheduled_target_count": 0,
        "remaining_open_schedules": [0, 1, 2, 3],
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def _run_cpp(source: Path) -> dict[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("no C++20 compiler is available")
    with tempfile.TemporaryDirectory(prefix="record31-a4-") as directory:
        executable = Path(directory) / "audit"
        subprocess.run(
            [compiler, "-std=c++20", "-O3", "-Wall", "-Wextra", "-pedantic",
             str(source), "-o", str(executable)],
            cwd=ROOT,
            check=True,
        )
        completed = subprocess.run(
            [str(executable)], cwd=ROOT, check=True, capture_output=True, text=True
        )
    parsed = json.loads(completed.stdout)
    if not isinstance(parsed, dict):
        raise ValueError(f"{source.name} did not emit a JSON object")
    return parsed


def _run_low() -> dict[str, object]:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        (str(ROOT), str(ROOT / "src"), environment.get("PYTHONPATH", ""))
    )
    completed = subprocess.run(
        [shutil.which("python3") or "python3", str(LOW_SOURCE)],
        cwd=ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    parsed = json.loads(completed.stdout)
    if not isinstance(parsed, dict):
        raise ValueError("low-addition source did not emit a JSON object")
    return parsed


def fresh_record() -> dict[str, object]:
    positive = _positive_document()
    document: dict[str, object] = {
        "kind": "h9_14_record31_a4_compact_schedule_no_go",
        "claim": (
            "Every completion of H9/14 graph-shadow Record 31 in the compact "
            "Dehn--Sommerville schedule a=4 is incidence-impossible."
        ),
        "scope": {
            "compact_only": True,
            "a4_compact_schedule_excluded": True,
            "whole_shadow_compact_coxeter_theorem": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "upstream": {
            "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
            "certificate_sha256": positive["certificate_sha256"],
            "a4_orbit_count": 20,
        },
        "production_low_additions": _run_low(),
        "production_three_additions": _run_cpp(HIGH_SOURCE),
        "independent": _run_cpp(INDEPENDENT_SOURCE),
        "downstream": {
            "scheduled_114_facet_target_count": 0,
            "ridge_topology_target_count": 0,
            "unrestricted_compact_label_target_count": 0,
            "reason": (
                "Every completion fails before the ridge/topology and compact-label "
                "gates: it either retains a transversal of order at most four or "
                "has blocker count different from the scheduled 114 facets."
            ),
        },
        "nonclaims": [
            "The compact schedules a=0,1,2,3 remain open for other completions of Record 31.",
            "No finite-volume noncompact or parabolic-cusp criterion is audited.",
            "The selected Record-31 nerve has a separate compact-label no-go; it is not used to infer this all-completion a=4 result.",
        ],
        "source_provenance": [
            _source_provenance(LOW_SOURCE),
            _source_provenance(HIGH_SOURCE),
            _source_provenance(INDEPENDENT_SOURCE),
        ],
    }
    document["certificate_sha256"] = _certificate_digest(document)
    verify(document)
    return document


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--emit", action="store_true")
    args = parser.parse_args(argv)
    document = fresh_record() if args.fresh else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    summary = verify(document)
    print(json.dumps(document if args.emit else summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
