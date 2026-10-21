#!/usr/bin/env python3
"""Freeze the compact-only exclusion of Record 31's entire ``a=3`` schedule."""

from __future__ import annotations

import argparse
from collections import Counter
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


PRODUCTION_SOURCE = ROOT / "scripts/explore_h9_14_record31_a3_additions.cpp"
DOWNSTREAM_SOURCE = ROOT / "scripts/explore_h9_14_record31_a3_downstream.py"
INDEPENDENT_SOURCE = (
    ROOT / "scripts/independent_h9_14_record31_a3_additions_audit.cpp"
)
INDEPENDENT_DOWNSTREAM_SOURCE = (
    ROOT / "scripts/independent_h9_14_record31_a3_downstream_audit.py"
)
HOSTILE_LANNER_DEPENDENCY = (
    ROOT / "scripts/independent_h9_14_record32_downstream_audit.py"
)
FACET_COMPLETION_DEPENDENCY = (
    ROOT / "src/coxeter_polytopes/combinatorics/facet_completion.py"
)
LABEL_DEPENDENCY = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py"
)
SUPPORT_ARCHIVE = ROOT / "data/enumerations/h9_14_nerve_support.json"
ARCHIVE = ROOT / "data/enumerations/h9_14_record31_a3_compact_no_go.json"

EXPECTED_TARGETS_BY_REPRESENTATIVE = {
    5: 1,
    10: 11,
    27: 1,
    49: 1,
    53: 1,
    59: 10,
    104: 2,
    108: 60,
    113: 14,
    116: 3,
    135: 6,
    136: 35,
    138: 68,
    139: 45,
    141: 85,
    142: 30,
}
EXPECTED_RIDGE_FAMILIES = {
    (60, 1984, 6212, 8248, 10112): 10,
    (60, 1984, 6336, 9984): 27,
    (60, 1984, 7712, 8220, 8640): 53,
}
EXPECTED_PRODUCTION_LABEL_TRANSCRIPTS = {
    (60, 1984, 6212, 8248, 10112): (
        (1, 420, 420), (4, 1332, 96), (2, 10368, 0),
    ),
    (60, 1984, 6336, 9984): (
        (1, 420, 420), (3, 1542, 84), (2, 2232, 0),
    ),
    (60, 1984, 7712, 8220, 8640): (
        (1, 420, 420), (4, 1542, 84), (2, 7992, 216), (0, 23328, 0),
    ),
}
EXPECTED_INDEPENDENT_LABEL_TRANSCRIPTS = {
    (60, 1984, 6212, 8248, 10112): (
        (4, 420, 420), (1, 1332, 96), (2, 10368, 0),
    ),
    (60, 1984, 6336, 9984): (
        (3, 108, 108), (1, 1542, 84), (2, 2232, 0),
    ),
    (60, 1984, 7712, 8220, 8640): (
        (4, 108, 108), (1, 1542, 84), (2, 7992, 216), (3, 23328, 0),
    ),
}


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"Record 31 a=3 {name} is malformed")
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
    branch = document["production"]["branches"][2]
    if branch["a"] != 3 or branch["orbit_count"] != 143:
        raise ValueError("Record 31 a=3 upstream orbit branch changed")
    return document


def _normalized_production_targets(
    production: Mapping[str, object],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    rows = production.get("representative_audits")
    if not isinstance(rows, list) or len(rows) != 143:
        raise ValueError("production row coverage changed")
    normalized = []
    for raw in rows:
        row = _mapping(raw, "production row")
        families = []
        for family in row["pure_target_minimal_nonfaces"]:
            families.append(
                tuple(
                    sorted(
                        sum(1 << int(vertex) for vertex in face)
                        for face in family
                    )
                )
            )
        normalized.append(tuple(sorted(families)))
    return tuple(normalized)


def _normalized_independent_targets(
    independent: Mapping[str, object],
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    rows = independent.get("representative_audits")
    if not isinstance(rows, list) or len(rows) != 143:
        raise ValueError("independent row coverage changed")
    return tuple(
        tuple(
            sorted(
                tuple(sorted(map(int, family)))
                for family in _mapping(raw, "independent row")[
                    "pure_target_minimal_nonface_masks"
                ]
            )
        )
        for raw in rows
    )


def _label_transcript(gate: Mapping[str, object]) -> tuple[tuple[int, int, int], ...]:
    stages = gate.get("stages")
    if not isinstance(stages, list):
        raise ValueError("label stages are malformed")
    return tuple(
        (
            int(_mapping(stage, "label stage")["support_index"]),
            int(_mapping(stage, "label stage")["compatible_extension_count"]),
            int(_mapping(stage, "label stage")["passing_state_count"]),
        )
        for stage in stages
    )


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a3_compact_schedule_no_go":
        raise ValueError("wrong Record 31 a=3 certificate kind")
    expected_scope = {
        "compact_only": True,
        "a3_compact_schedule_excluded": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    if _mapping(document.get("scope"), "scope") != expected_scope:
        raise ValueError("Record 31 a=3 scope changed or was over-promoted")

    positive = _positive_document()
    upstream = _mapping(document.get("upstream"), "upstream")
    if upstream != {
        "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
        "certificate_sha256": positive["certificate_sha256"],
        "a3_orbit_count": 143,
    }:
        raise ValueError("Record 31 a=3 upstream certificate changed")

    production = _mapping(document.get("production"), "production")
    if (
        production.get("kind") != "h9_14_record31_a3_pivoted_addition_scout"
        or production.get("compact_only") is not True
        or int(production.get("pure_116_facet_target_count_across_representatives", -1))
        != 373
    ):
        raise ValueError("Record 31 a=3 production summary changed")
    production_rows_raw = production.get("representative_audits")
    if not isinstance(production_rows_raw, list) or len(production_rows_raw) != 143:
        raise ValueError("Record 31 a=3 production row count changed")
    production_rows = [_mapping(row, "production row") for row in production_rows_raw]
    addition_orbits = Counter(
        int(row["additional_minimal_A5_count"]) for row in production_rows
    )
    addition_targets = Counter()
    observed_targets_by_index = {}
    total_scheduled = 0
    for index, row in enumerate(production_rows):
        if int(row["orbit_representative_index"]) != index:
            raise ValueError("production row indices changed")
        count = int(row["pure_116_facet_target_count"])
        if count:
            observed_targets_by_index[index] = count
        addition_targets[int(row["additional_minimal_A5_count"])] += count
        total_scheduled += sum(
            int(value)
            for value in _mapping(
                row["scheduled_order_five_count_blocker_size_histogram"],
                "production scheduled blocker histogram",
            ).values()
        )
    if (
        addition_orbits != Counter({1: 55, 2: 40, 3: 37, 4: 11})
        or addition_targets != Counter({1: 3, 2: 22, 3: 25, 4: 323})
        or observed_targets_by_index != EXPECTED_TARGETS_BY_REPRESENTATIVE
        or total_scheduled != 378
        or production_rows[141][
            "scheduled_order_five_count_blocker_size_histogram"
        ] != {"116": 85, "118": 5}
    ):
        raise ValueError("Record 31 a=3 production census changed")

    independent = _mapping(document.get("independent"), "independent audit")
    if (
        independent.get("kind")
        != "h9_14_record31_a3_independent_partition_audit"
        or independent.get("compact_only") is not True
        or int(independent.get(
            "scheduled_order_five_state_count_across_representatives", -1
        )) != 378
        or int(independent.get(
            "pure_116_facet_target_count_across_representatives", -1
        )) != 373
    ):
        raise ValueError("Record 31 a=3 independent summary changed")
    independent_rows_raw = independent.get("representative_audits")
    if not isinstance(independent_rows_raw, list) or len(independent_rows_raw) != 143:
        raise ValueError("Record 31 a=3 independent row count changed")
    independent_rows = [
        _mapping(row, "independent row") for row in independent_rows_raw
    ]
    for index, (production_row, independent_row) in enumerate(
        zip(production_rows, independent_rows, strict=True)
    ):
        if (
            int(independent_row["orbit_representative_index"]) != index
            or int(independent_row["additional_minimal_A5_count"])
            != int(production_row["additional_minimal_A5_count"])
            or int(independent_row["individual_candidate_count"])
            != int(production_row["individual_candidate_count"])
            or int(independent_row["base_small_transversal_count"])
            != int(production_row["base_small_transversal_count"])
            or int(independent_row["pure_116_facet_target_count"])
            != int(production_row["pure_116_facet_target_count"])
            or independent_row["scheduled_blocker_size_histogram"]
            != production_row["scheduled_order_five_count_blocker_size_histogram"]
        ):
            raise ValueError(f"production/independent metadata mismatch at {index}")
    if _normalized_production_targets(production) != _normalized_independent_targets(
        independent
    ):
        raise ValueError("production/independent literal target families differ")

    downstream = _mapping(document.get("downstream"), "production downstream")
    expected_downstream_counts = {
        "literal_target_count": 373,
        "target_orbit_count": 33,
        "ridge_pseudomanifold_orbit_count": 3,
        "connected_mod2_homology_sphere_orbit_count": 3,
        "unrestricted_compact_label_no_go_orbit_count": 3,
        "compact_label_survivor_orbit_count": 0,
    }
    if (
        downstream.get("kind")
        != "h9_14_record31_a3_downstream_compact_audit"
        or downstream.get("compact_only") is not True
        or any(int(downstream.get(key, -1)) != value
               for key, value in expected_downstream_counts.items())
    ):
        raise ValueError("Record 31 a=3 downstream summary changed")
    orbit_rows_raw = downstream.get("target_orbits")
    if not isinstance(orbit_rows_raw, list) or len(orbit_rows_raw) != 33:
        raise ValueError("Record 31 a=3 target orbit count changed")
    if sum(int(_mapping(row, "target orbit")["literal_multiplicity"])
           for row in orbit_rows_raw) != 373:
        raise ValueError("Record 31 a=3 target orbit multiplicities changed")
    observed_ridge = {}
    for raw in orbit_rows_raw:
        row = _mapping(raw, "target orbit")
        if not row["ridge_pseudomanifold"]:
            continue
        family = tuple(map(int, row["minimal_nonface_masks"]))
        topology = _mapping(row["global_topology_audit"], "global topology")
        gate = _mapping(row["sparse_compact_label_gate"], "sparse label gate")
        if (
            int(row["literal_multiplicity"]) != 1
            or topology.get("face_vector")
            != [14, 87, 319, 768, 1269, 1449, 1110, 522, 116]
            or topology.get("ridge_degree_histogram") != {"2": 522}
            or topology.get("passes_connected_homology_audit") is not True
            or int(gate["rank_nine_vertex_count"]) != 116
            or int(gate["terminal_sparse_state_count"]) != 0
            or gate.get("state_cap_used") is not False
            or gate.get("unrestricted_compact_label_no_go") is not True
            or _label_transcript(gate)
            != EXPECTED_PRODUCTION_LABEL_TRANSCRIPTS.get(family)
        ):
            raise ValueError("Record 31 a=3 ridge/topology/label row changed")
        observed_ridge[family] = int(row["representative_source_index"])
    if observed_ridge != EXPECTED_RIDGE_FAMILIES:
        raise ValueError("Record 31 a=3 ridge family identity changed")

    hostile_replay = _mapping(
        document.get("independent_downstream"), "independent downstream"
    )
    if (
        hostile_replay.get("kind")
        != "independent_h9_14_record31_a3_literal_downstream_audit"
        or hostile_replay.get("compact_only") is not True
        or hostile_replay.get("classical_lanner_placement_counts")
        != {"4": 108, "5": 420}
        or int(hostile_replay.get("literal_target_count", -1)) != 373
        or int(hostile_replay.get("ridge_pseudomanifold_literal_count", -1)) != 3
        or hostile_replay.get(
            "all_ridge_survivors_have_unrestricted_compact_label_no_go"
        ) is not True
        or int(hostile_replay.get("compact_label_survivor_count", -1)) != 0
    ):
        raise ValueError("Record 31 a=3 independent downstream summary changed")
    replay_rows_raw = hostile_replay.get("ridge_survivor_label_audits")
    if not isinstance(replay_rows_raw, list) or len(replay_rows_raw) != 3:
        raise ValueError("independent ridge survivor rows changed")
    replay_families = {}
    for raw in replay_rows_raw:
        row = _mapping(raw, "independent label row")
        family = tuple(map(int, row["minimal_nonface_masks"]))
        if (
            int(row["rank_nine_vertex_count"]) != 116
            or int(row["terminal_sparse_state_count"]) != 0
            or row.get("state_cap_used") is not False
            or row.get("unrestricted_compact_label_no_go") is not True
            or _label_transcript(row)
            != EXPECTED_INDEPENDENT_LABEL_TRANSCRIPTS.get(family)
        ):
            raise ValueError("independent Record 31 label transcript changed")
        replay_families[family] = int(row["source_representative_index"])
    if replay_families != EXPECTED_RIDGE_FAMILIES:
        raise ValueError("independent Record 31 ridge family set changed")

    expected_provenance = [
        _source_provenance(path)
        for path in (
            PRODUCTION_SOURCE,
            DOWNSTREAM_SOURCE,
            INDEPENDENT_SOURCE,
            INDEPENDENT_DOWNSTREAM_SOURCE,
            HOSTILE_LANNER_DEPENDENCY,
            FACET_COMPLETION_DEPENDENCY,
            LABEL_DEPENDENCY,
        )
    ]
    if document.get("source_provenance") != expected_provenance:
        raise ValueError("Record 31 a=3 source provenance changed")
    if document.get("support_archive") != _source_provenance(SUPPORT_ARCHIVE):
        raise ValueError("Record 31 a=3 support archive provenance changed")
    if document.get("certificate_sha256") != _certificate_digest(document):
        raise ValueError("Record 31 a=3 certificate digest changed")
    return {
        "record_index": 31,
        "newly_closed_schedule": 3,
        "a3_orbit_count": 143,
        "literal_pure_target_count": 373,
        "target_orbit_count": 33,
        "ridge_pseudomanifold_orbit_count": 3,
        "unrestricted_compact_label_survivor_count": 0,
        "remaining_open_schedules": [0, 1, 2],
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }


def _representative_input(positive: Mapping[str, object]) -> str:
    orbits = positive["production"]["branches"][2]["orbits"]
    lines = [str(len(orbits))]
    for index, orbit in enumerate(orbits):
        masks = [
            sum(1 << int(vertex) for vertex in face)
            for face in orbit["representative"]
        ]
        if len(masks) != 3:
            raise ValueError("Record 31 a=3 representative shape changed")
        lines.append(" ".join(map(str, (index, *masks))))
    return "\n".join(lines) + "\n"


def _run_cpp(source: Path, standard_input: str) -> dict[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("no C++20 compiler is available")
    with tempfile.TemporaryDirectory(prefix="record31-a3-") as directory:
        executable = Path(directory) / "audit"
        subprocess.run(
            [compiler, "-std=c++20", "-O3", "-Wall", "-Wextra", "-pedantic",
             str(source), "-o", str(executable)],
            cwd=ROOT,
            check=True,
        )
        completed = subprocess.run(
            [str(executable)],
            cwd=ROOT,
            input=standard_input,
            check=True,
            capture_output=True,
            text=True,
        )
    parsed = json.loads(completed.stdout)
    if not isinstance(parsed, dict):
        raise ValueError(f"{source.name} did not emit a JSON object")
    return parsed


def _run_downstream(source: Path, payload: Mapping[str, object]) -> dict[str, object]:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        (str(ROOT), str(ROOT / "src"), environment.get("PYTHONPATH", ""))
    )
    with tempfile.TemporaryDirectory(prefix="record31-a3-downstream-") as directory:
        input_path = Path(directory) / "input.json"
        input_path.write_text(json.dumps(payload), encoding="utf-8")
        completed = subprocess.run(
            [shutil.which("python3") or "python3", str(source), str(input_path)],
            cwd=ROOT,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )
    parsed = json.loads(completed.stdout)
    if not isinstance(parsed, dict):
        raise ValueError(f"{source.name} did not emit a JSON object")
    return parsed


def fresh_record() -> dict[str, object]:
    positive = _positive_document()
    standard_input = _representative_input(positive)
    production = _run_cpp(PRODUCTION_SOURCE, standard_input)
    independent = _run_cpp(INDEPENDENT_SOURCE, standard_input)
    document: dict[str, object] = {
        "kind": "h9_14_record31_a3_compact_schedule_no_go",
        "claim": (
            "Every completion of H9/14 graph-shadow Record 31 in the compact "
            "Dehn--Sommerville schedule a=3 fails an exact necessary condition: "
            "all but three target orbits fail the ridge-pseudomanifold gate, and "
            "each of the three ridge survivors has an unrestricted compact "
            "Coxeter-label no-go on its actual rank-nine vertices."
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
            "a3_orbit_count": 143,
        },
        "production": production,
        "independent": independent,
        "downstream": _run_downstream(DOWNSTREAM_SOURCE, production),
        "independent_downstream": _run_downstream(
            INDEPENDENT_DOWNSTREAM_SOURCE, independent
        ),
        "support_archive": _source_provenance(SUPPORT_ARCHIVE),
        "nonclaims": [
            "The compact schedules a=0,1,2 remain open for other completions of Record 31.",
            "A mod-two homology-sphere audit is not promoted to an abstract polytopality theorem.",
            "No finite-volume noncompact or parabolic-cusp criterion is audited.",
            "No new hyperbolic Coxeter polytope is produced.",
        ],
        "source_provenance": [
            _source_provenance(path)
            for path in (
                PRODUCTION_SOURCE,
                DOWNSTREAM_SOURCE,
                INDEPENDENT_SOURCE,
                INDEPENDENT_DOWNSTREAM_SOURCE,
                HOSTILE_LANNER_DEPENDENCY,
                FACET_COMPLETION_DEPENDENCY,
                LABEL_DEPENDENCY,
            )
        ],
    }
    document["certificate_sha256"] = _certificate_digest(document)
    verify(document)
    return document


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--emit", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    document = (
        fresh_record()
        if args.fresh
        else json.loads(ARCHIVE.read_text(encoding="utf-8"))
    )
    if args.write:
        ARCHIVE.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    summary = verify(document)
    print(json.dumps(document if args.emit else summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
