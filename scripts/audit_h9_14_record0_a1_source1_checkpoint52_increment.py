#!/usr/bin/env python3
"""Independently audit source-1 compact CEGAR cores 51 and 52.

This verifier does not import the exploratory CEGAR or its core extractor.
It replays the two appended local contradictions with the frozen independent
finite/Lanner classifier and rebuilds each complete ``S3 x S8`` source orbit
by generator BFS.  The result is a bounded incremental certificate, not a
terminal compact no-go theorem.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from scripts.audit_h9_14_record0_a1_source1_checkpoint48_core import (
    EXPECTED_SOURCE_GROUP_ORDER,
    FIXED_FOUR,
    FULL,
    _core_orbit,
    _independent_five,
    _project,
    _source_generators,
)
from scripts.audit_h9_14_record0_a1_source3_terminal_cegar import (
    _audit_lanner_tables,
    _checkpoint_stage_projection,
    _digest,
    _label_state_count,
    _mask,
)


ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT_50 = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_label_core_cegar_checkpoint_50_more.json"
)
CHECKPOINT_52 = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_label_core_cegar_checkpoint_52_more.json"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_checkpoint52_increment_audit.json"
)
BASE_AUDITOR = (
    ROOT / "scripts/audit_h9_14_record0_a1_source1_checkpoint48_core.py"
)

CHECKPOINT_50_SHA256 = (
    "5522deaca61bd37518ea59694263b359d3831db073351b7715366c8835265e9e"
)
CHECKPOINT_52_SHA256 = (
    "8577321be0061fca0aa9401069b8dad2fea7fabbb6cb583524234b8932a28f48"
)
BASE_AUDITOR_SHA256 = (
    "ad3d8cde9dc9bda1af69b3f748cb19b030890df43a2123a025396f4722d91ef5"
)
SCHEMA = "coxeter.h9_14.record0_a1_source1.checkpoint52_increment_audit.v1"
PREVIOUS_CUMULATIVE_ORBIT_SIZE_SUM = 367_836
EXPECTED_CUMULATIVE_ORBIT_SIZE_SUM = 387_996

EXPECTED = {
    50: {
        "digest": "e0b91dc79bf9f4bfd7fc8db11d9a618a92d07480070d0d08e78a156d0f2a89a4",
        "support_count": 4,
        "actual_vertex_count": 2,
        "stored_global_passing": [108, 84, 36, 6, 0],
        "independent_passing": [108, 110, 48, 0],
        "independent_compatible": [108, 1542, 686, 624],
        "orbit_size": 10080,
        "stabilizer": 24,
        "orbit_sha256": "c1fe586ce6ce9697c976a04f78bae6a6a6d774562633f48f9ae6386ceae48735",
    },
    51: {
        "digest": "e47a250ee632eec60f7f8d8690f2736888d7eb1abc6548d9f60125370d087bbc",
        "support_count": 6,
        "actual_vertex_count": 2,
        "stored_global_passing": [108, 84, 36, 4, 4, 0],
        "independent_passing": [108, 110, 62, 103, 10, 0],
        "independent_compatible": [108, 1542, 394, 266, 40, 34],
        "orbit_size": 10080,
        "stabilizer": 24,
        "orbit_sha256": "b89b8356f2372c2bb59300ad3292088c0ec80746204990636a6d15fec9694361",
    },
}


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _audit_row(row: dict[str, Any]) -> dict[str, Any]:
    iteration = int(row["iteration"])
    expected = EXPECTED[iteration]
    core = row["local_core"]
    supports = tuple(_mask(values) for values in core["core_supports"])
    additions = tuple(
        sorted(_mask(values) for values in core["core_additional_A5s"])
    )
    actual_vertices = tuple(
        _mask(values) for values in core["core_actual_vertices"]
    )
    restrictions = tuple(
        sorted(_mask(values) for values in core["core_face_restrictions"])
    )
    facet_complements = tuple(
        sorted(_mask(values) for values in core["core_facet_complements"])
    )
    support_union = 0
    for support in supports:
        support_union |= support

    actual_terminal, actual_stages = _label_state_count(
        supports, actual_vertices
    )
    cut_terminal, cut_stages = _label_state_count(supports, restrictions)
    expected_stages = _checkpoint_stage_projection(core["replay_stages"])
    actual_projection = _checkpoint_stage_projection(actual_stages)
    cut_projection = _checkpoint_stage_projection(cut_stages)

    orbit = _core_orbit((additions, restrictions), _source_generators())
    orbit_sha256 = _digest(sorted(orbit))
    passing = [stage["passing_state_count"] for stage in cut_stages]
    compatible = [
        stage["compatible_extension_count"] for stage in cut_stages
    ]
    checks = {
        "representative_digest_matches": row["representative_digest"]
        == expected["digest"],
        "fixed_A4_occurs_once": supports.count(FIXED_FOUR) == 1,
        "remaining_supports_are_exactly_stored_independent_A5s": (
            tuple(sorted(value for value in supports if value != FIXED_FOUR))
            == additions
            and len(supports) == expected["support_count"]
            and all(_independent_five(value) for value in additions)
        ),
        "all_actual_vertices_have_rank_nine": (
            len(actual_vertices) == expected["actual_vertex_count"]
            and all(value.bit_count() == 9 for value in actual_vertices)
        ),
        "facet_complements_match_actual_vertices": facet_complements
        == tuple(sorted(FULL ^ value for value in actual_vertices)),
        "face_cut_has_same_projected_semantics": _project(
            actual_vertices, support_union
        )
        == _project(restrictions, support_union),
        "actual_and_cut_replays_match_archived_stages_and_die": (
            actual_projection == expected_stages
            and cut_projection == expected_stages
            and actual_terminal == cut_terminal == 0
        ),
        "stored_global_label_counts_match": row["label_stage_state_counts"]
        == expected["stored_global_passing"],
        "independent_local_state_counts_match": (
            passing == expected["independent_passing"]
            and compatible == expected["independent_compatible"]
        ),
        "complete_generator_orbit_size_matches": len(orbit)
        == row["full_core_orbit_size"]
        == row["full_face_core_orbit_size"]
        == expected["orbit_size"],
        "complete_generator_orbit_hash_matches": orbit_sha256
        == expected["orbit_sha256"],
        "stabilizer_order_matches": EXPECTED_SOURCE_GROUP_ORDER // len(orbit)
        == row["core_stabilizer_order"]
        == expected["stabilizer"],
        "every_orbit_addition_remains_an_independent_A5": all(
            all(_independent_five(value) for value in image_additions)
            for image_additions, _ in orbit
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"independent source-1 core {iteration + 1} audit failed: {checks}"
        )
    return {
        "iteration": iteration,
        "human_core_number": iteration + 1,
        "representative_digest": row["representative_digest"],
        "support_profile": f"fixed A4 + {len(additions)} A5 supports",
        "local_actual_vertex_count": len(actual_vertices),
        "stored_global_label_stage_state_counts": row[
            "label_stage_state_counts"
        ],
        "independent_local_passing_state_counts": passing,
        "independent_local_compatible_extension_counts": compatible,
        "full_source_orbit_size": len(orbit),
        "core_stabilizer_order": EXPECTED_SOURCE_GROUP_ORDER // len(orbit),
        "full_source_orbit_sha256": orbit_sha256,
        "solver_branches": row["solver_branch_count"],
        "solver_conflicts": row["solver_conflict_count"],
        "solver_wall_time_seconds": row["solver_wall_time_seconds"],
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(CHECKPOINT_50) != CHECKPOINT_50_SHA256:
        raise ArithmeticError("source-1 checkpoint-50 hash changed")
    if _sha256(CHECKPOINT_52) != CHECKPOINT_52_SHA256:
        raise ArithmeticError("source-1 checkpoint-52 hash changed")
    if _sha256(BASE_AUDITOR) != BASE_AUDITOR_SHA256:
        raise ArithmeticError("frozen independent base auditor hash changed")

    checkpoint_50 = json.loads(CHECKPOINT_50.read_text(encoding="utf-8"))
    checkpoint_52 = json.loads(CHECKPOINT_52.read_text(encoding="utf-8"))
    rows = checkpoint_52["iterations"][-2:]
    if [int(row["iteration"]) for row in rows] != [50, 51]:
        raise ArithmeticError("checkpoint does not end in iterations 50 and 51")
    summaries = [_audit_row(row) for row in rows]
    lanner = _audit_lanner_tables()
    cumulative = sum(
        int(row["full_core_orbit_size"])
        for row in checkpoint_52["iterations"]
    )
    previous_cumulative = sum(
        int(row["full_core_orbit_size"])
        for row in checkpoint_50["iterations"]
    )
    checks = {
        "checkpoint_identity_is_record0_a1_source1": (
            checkpoint_52["record_index"] == 0
            and checkpoint_52["schedule_a"] == 1
            and checkpoint_52["source_orbit_index"] == 1
        ),
        "first_50_iterations_are_byte_for_byte_structurally_unchanged": (
            checkpoint_52["iterations"][:50] == checkpoint_50["iterations"]
        ),
        "exactly_52_bounded_core_iterations_are_present": (
            len(checkpoint_52["iterations"]) == 52
            and checkpoint_52["excluded_local_core_orbit_count"] == 52
        ),
        "checkpoint_remains_OPTIMAL_core_cutoff_not_a_theorem": (
            checkpoint_52["terminal_solver_status"] == "OPTIMAL"
            and checkpoint_52["stopped_reason"] == "core_cutoff"
            and checkpoint_52["fixed_source_compact_coxeter_no_go"] is False
            and checkpoint_52["theorem_certificate"] is False
        ),
        "both_incremental_core_audits_pass": all(
            all(summary["checks"].values()) for summary in summaries
        ),
        "independent_classical_Lanner_tables_pass": all(
            lanner["checks"].values()
        ),
        "cumulative_orbit_size_sum_matches": (
            previous_cumulative == PREVIOUS_CUMULATIVE_ORBIT_SIZE_SUM
            and cumulative == EXPECTED_CUMULATIVE_ORBIT_SIZE_SUM
            and cumulative - previous_cumulative == 20_160
        ),
        "finite_volume_noncompact_is_not_claimed": (
            checkpoint_52["finite_volume_noncompact_checked"] is False
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"checkpoint-52 incremental audit failed: {checks}")
    return {
        "schema": SCHEMA,
        "scope": (
            "Independent proof for newly appended iterations 50 and 51 and "
            "their complete source-1 orbits. The 52-core checkpoint remains "
            "bounded and is not a fixed-source compact no-go theorem."
        ),
        "inputs": {
            "checkpoint_50": str(CHECKPOINT_50.relative_to(ROOT)),
            "checkpoint_50_sha256": CHECKPOINT_50_SHA256,
            "checkpoint_52": str(CHECKPOINT_52.relative_to(ROOT)),
            "checkpoint_52_sha256": CHECKPOINT_52_SHA256,
            "base_independent_auditor": str(BASE_AUDITOR.relative_to(ROOT)),
            "base_independent_auditor_sha256": BASE_AUDITOR_SHA256,
        },
        "bounded_checkpoint_status": {
            "terminal_solver_status": checkpoint_52[
                "terminal_solver_status"
            ],
            "stopped_reason": checkpoint_52["stopped_reason"],
            "excluded_local_core_orbits": checkpoint_52[
                "excluded_local_core_orbit_count"
            ],
            "previous_cumulative_orbit_size_sum": previous_cumulative,
            "incremental_orbit_size_sum": cumulative - previous_cumulative,
            "cumulative_orbit_size_sum": cumulative,
            "orbit_size_sum_is_not_a_cross_core_distinct_clause_count": True,
            "fixed_source_compact_coxeter_no_go": False,
            "theorem_certificate": False,
            "finite_volume_noncompact_checked": False,
        },
        "new_cores": summaries,
        "independent_lanner_tables": lanner,
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        if archived != certificate:
            raise SystemExit("checkpoint-52 incremental archive mismatch")
        print("archive_check=PASS")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
