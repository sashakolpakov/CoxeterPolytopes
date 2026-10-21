#!/usr/bin/env python3
"""Certify the compact-incidence exclusion of Record 30's a=3 schedule."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_h9_14_record30_positive_scout import (
    ARCHIVE as POSITIVE_ARCHIVE,
    verify as verify_positive,
)


PRODUCTION_SOURCE = ROOT / "scripts/explore_h9_14_record30_a3_additions.cpp"
INDEPENDENT_SOURCE = ROOT / "scripts/independent_h9_14_record30_a3_audit.cpp"
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a3_compact_no_go.json"

EXPECTED_SOURCE_SHA256 = {
    str(PRODUCTION_SOURCE.relative_to(ROOT)):
        "5b17fa1448d3322063f8a4c0d9dc49edfa982fb53685354caa273174bd6ee759",
    str(INDEPENDENT_SOURCE.relative_to(ROOT)):
        "c7c85107780728fd83857ca7e8eb30d73e14c617a5e7893cfd859d46d5220089",
}
EXPECTED_NO_SMALL = (
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    6, 1, 1, 0, 0, 1, 0, 0,
    735_086, 1_158_881,
    840, 840, 840,
    0, 0, 0, 0, 0, 0, 0, 0,
)
EXPECTED_PROJECTION_SHA256 = (
    "7afaebc75008dd16c6bf499df24de767c5f0a9cd9bc672b563e175c7de38c5a5"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "e3afaf33526c003307d797bd0789bd2b4275d3d74c83f0a8d9792bb1ad0b9d54"
)


def canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Record 30 a=3 {name} is malformed")
    return value


def _source_provenance(path: Path) -> dict[str, str]:
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": sha256(path.read_bytes()).hexdigest(),
    }


def _certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return canonical_sha(payload)


def _positive_document() -> dict[str, object]:
    document = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(document)
    return document


def _a3_input(
    positive: Mapping[str, object],
) -> tuple[str, list[list[int]]]:
    a1_a4 = _mapping(positive.get("a1_a4"), "positive transcript")
    branches = a1_a4.get("branches")
    if not isinstance(branches, list):
        raise ValueError("Record 30 positive branches are malformed")
    branch = next(
        (_mapping(row, "positive branch") for row in branches
         if _mapping(row, "positive branch").get("a") == 3),
        None,
    )
    if branch is None:
        raise ValueError("Record 30 positive a=3 branch is missing")
    orbits = branch.get("orbits")
    if not isinstance(orbits, list) or len(orbits) != 32:
        raise ValueError("Record 30 positive a=3 orbit count changed")
    masks: list[list[int]] = []
    rows = [str(len(orbits))]
    for index, raw_orbit in enumerate(orbits):
        orbit = _mapping(raw_orbit, "positive orbit")
        representative = orbit.get("representative")
        if not isinstance(representative, list) or len(representative) != 3:
            raise ValueError(f"Record 30 a=3 representative {index} changed")
        encoded = [sum(1 << int(vertex) for vertex in face)
                   for face in representative]
        if any(mask.bit_count() != 4 for mask in encoded):
            raise ValueError("Record 30 a=3 representative is not an A4 family")
        masks.append(encoded)
        rows.append(" ".join(map(str, (index, *encoded))))
    return "\n".join(rows) + "\n", masks


def _run_cpp(source: Path, input_text: str) -> dict[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("no C++20 compiler is available")
    with tempfile.TemporaryDirectory(prefix="record30-a3-") as directory:
        executable = Path(directory) / "audit"
        subprocess.run(
            [compiler, "-std=c++20", "-O3", "-Wall", "-Wextra",
             "-Werror", "-pedantic", str(source), "-o", str(executable)],
            cwd=ROOT,
            check=True,
        )
        completed = subprocess.run(
            [str(executable)],
            cwd=ROOT,
            input=input_text,
            text=True,
            capture_output=True,
            check=True,
        )
    parsed = json.loads(completed.stdout)
    if not isinstance(parsed, dict):
        raise ValueError(f"{source.name} did not emit a JSON object")
    return parsed


def _common_projection(row: Mapping[str, object]) -> dict[str, object]:
    return {
        "orbit_index": int(row["orbit_index"]),
        "forced_A5_shadow_count": int(row["forced_A5_shadow_count"]),
        "additional_A5_count": int(row["additional_A5_count"]),
        "candidate_count": int(row["candidate_count"]),
        "raw_unordered_addition_count": int(row["raw_unordered_addition_count"]),
        "base_small_transversal_count": int(row["base_small_transversal_count"]),
        "no_small_transversal_family_count": int(
            row["no_small_transversal_family_count"]
        ),
        "order_five_transversal_count_histogram": row[
            "order_five_transversal_count_histogram"
        ],
        "scheduled_pure_target_count": int(row["scheduled_pure_target_count"]),
    }


def _verified_rows(
    document: Mapping[str, object], masks: list[list[int]]
) -> list[dict[str, object]]:
    production = _mapping(document.get("production"), "production audit")
    independent = _mapping(document.get("independent"), "independent audit")
    if (
        production.get("kind") != "h9_14_record30_a3_addition_census"
        or production.get("compact_only") is not True
        or int(production.get("orbit_count", -1)) != 32
    ):
        raise ValueError("Record 30 a=3 production header changed")
    if (
        independent.get("kind")
        != "h9_14_record30_a3_independent_pivot_audit"
        or independent.get("compact_only") is not True
        or int(independent.get("orbit_count", -1)) != 32
        or int(independent.get("scheduled_pure_target_count_across_orbits", -1))
        != 0
    ):
        raise ValueError("Record 30 a=3 independent header changed")
    production_rows = production.get("orbits")
    independent_rows = independent.get("orbits")
    if (
        not isinstance(production_rows, list)
        or not isinstance(independent_rows, list)
        or len(production_rows) != 32
        or len(independent_rows) != 32
    ):
        raise ValueError("Record 30 a=3 audit row coverage changed")

    projection: list[dict[str, object]] = []
    for index, (raw_production, raw_independent) in enumerate(
        zip(production_rows, independent_rows, strict=True)
    ):
        production_row = _mapping(raw_production, "production row")
        independent_row = _mapping(raw_independent, "independent row")
        if production_row.get("A4_masks") != masks[index]:
            raise ValueError(f"Record 30 a=3 upstream orbit {index} changed")
        common = _common_projection(production_row)
        if common != _common_projection(independent_row):
            raise ValueError(
                f"Record 30 a=3 production/independent mismatch at orbit {index}"
            )
        histogram = common["order_five_transversal_count_histogram"]
        if not isinstance(histogram, Mapping):
            raise ValueError("Record 30 a=3 order-five histogram is malformed")
        if sum(int(value) for value in histogram.values()) != EXPECTED_NO_SMALL[index]:
            raise ValueError(f"Record 30 a=3 histogram total {index} changed")
        if (
            common["orbit_index"] != index
            or common["no_small_transversal_family_count"]
            != EXPECTED_NO_SMALL[index]
            or "126" in histogram
            or common["scheduled_pure_target_count"] != 0
            or production_row.get("scheduled_order_five_blocker_audit_histogram")
            != []
            or independent_row.get(
                "scheduled_order_five_blocker_count_histogram"
            ) != {}
        ):
            raise ValueError(f"Record 30 a=3 exclusion row {index} changed")
        projection.append(common)
    return projection


def verify(
    document: Mapping[str, object], *, allow_unfrozen: bool = False
) -> dict[str, object]:
    if document.get("kind") != "h9_14_record30_a3_compact_schedule_no_go":
        raise ValueError("wrong Record 30 a=3 certificate kind")
    expected_scope = {
        "compact_only": True,
        "a3_compact_schedule_excluded": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    if _mapping(document.get("scope"), "scope") != expected_scope:
        raise ValueError("Record 30 a=3 scope changed or was over-promoted")

    positive = _positive_document()
    input_text, masks = _a3_input(positive)
    del input_text
    upstream = _mapping(document.get("upstream"), "upstream")
    if upstream != {
        "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
        "certificate_sha256": positive["certificate_sha256"],
        "a3_orbit_count": 32,
    }:
        raise ValueError("Record 30 a=3 upstream certificate changed")

    projection = _verified_rows(document, masks)
    projection_sha = canonical_sha(projection)
    if (
        not allow_unfrozen
        and EXPECTED_PROJECTION_SHA256 != "TO_BE_FROZEN"
        and projection_sha != EXPECTED_PROJECTION_SHA256
    ):
        raise ValueError("Record 30 a=3 replay projection changed")

    expected_downstream = {
        "scheduled_126_facet_target_count": 0,
        "ridge_topology_target_count": 0,
        "unrestricted_compact_label_target_count": 0,
        "reason": (
            "Every completion either retains a transversal of order at most "
            "four or has a number of order-five minimal transversals different "
            "from the scheduled 126 facets."
        ),
    }
    if _mapping(document.get("downstream"), "downstream") != expected_downstream:
        raise ValueError("Record 30 a=3 downstream-vacuity statement changed")

    expected_provenance = [
        _source_provenance(PRODUCTION_SOURCE),
        _source_provenance(INDEPENDENT_SOURCE),
    ]
    if document.get("source_provenance") != expected_provenance:
        raise ValueError("Record 30 a=3 source provenance changed")
    for row in expected_provenance:
        if EXPECTED_SOURCE_SHA256[row["path"]] != row["sha256"]:
            raise ValueError("Record 30 a=3 frozen source digest changed")

    digest = _certificate_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("Record 30 a=3 certificate digest is invalid")
    if (
        not allow_unfrozen
        and EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("Record 30 a=3 frozen certificate changed")
    return {
        "record_index": 30,
        "newly_closed_schedule": 3,
        "a3_orbit_count": 32,
        "no_small_transversal_family_count": sum(EXPECTED_NO_SMALL),
        "scheduled_126_facet_target_count": 0,
        "remaining_open_schedules": [0, 1, 2, 4],
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "projection_sha256": projection_sha,
    }


def fresh_record() -> dict[str, object]:
    positive = _positive_document()
    input_text, _ = _a3_input(positive)
    document: dict[str, object] = {
        "kind": "h9_14_record30_a3_compact_schedule_no_go",
        "claim": (
            "Every completion of H9/14 graph-shadow Record 30 in the compact "
            "Dehn--Sommerville schedule a=3 is incidence-impossible."
        ),
        "scope": {
            "compact_only": True,
            "a3_compact_schedule_excluded": True,
            "whole_shadow_compact_coxeter_theorem": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "upstream": {
            "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
            "certificate_sha256": positive["certificate_sha256"],
            "a3_orbit_count": 32,
        },
        "production": _run_cpp(PRODUCTION_SOURCE, input_text),
        "independent": _run_cpp(INDEPENDENT_SOURCE, input_text),
        "downstream": {
            "scheduled_126_facet_target_count": 0,
            "ridge_topology_target_count": 0,
            "unrestricted_compact_label_target_count": 0,
            "reason": (
                "Every completion either retains a transversal of order at most "
                "four or has a number of order-five minimal transversals different "
                "from the scheduled 126 facets."
            ),
        },
        "nonclaims": [
            "The compact schedules a=0,1,2,4 remain open for other completions of Record 30.",
            "No finite-volume noncompact or parabolic-cusp criterion is audited.",
            "This schedule exclusion is combinatorial and does not construct a new hyperbolic Coxeter polytope.",
        ],
        "source_provenance": [
            _source_provenance(PRODUCTION_SOURCE),
            _source_provenance(INDEPENDENT_SOURCE),
        ],
    }
    document["certificate_sha256"] = _certificate_digest(document)
    verify(document, allow_unfrozen=True)
    return document


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--emit", action="store_true")
    args = parser.parse_args(argv)
    document = (
        fresh_record()
        if args.fresh
        else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    )
    summary = verify(document, allow_unfrozen=args.fresh)
    print(json.dumps(document if args.emit else summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
