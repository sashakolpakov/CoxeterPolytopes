#!/usr/bin/env python3
"""Independently audit source-1 compact CEGAR cores 49 and 50.

The finite/Lanner replay and generator-orbit machinery are inherited from
the already frozen independent checkpoint-48 auditor.  No exploratory CEGAR
or core-extraction routine is imported.  This remains a bounded incremental
audit, not a terminal fixed-source no-go certificate.
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
CHECKPOINT = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_label_core_cegar_checkpoint_50_more.json"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_checkpoint50_increment_audit.json"
)
BASE_AUDITOR = (
    ROOT / "scripts/audit_h9_14_record0_a1_source1_checkpoint48_core.py"
)
CHECKPOINT_SHA256 = (
    "5522deaca61bd37518ea59694263b359d3831db073351b7715366c8835265e9e"
)
BASE_AUDITOR_SHA256 = (
    "ad3d8cde9dc9bda1af69b3f748cb19b030890df43a2123a025396f4722d91ef5"
)
SCHEMA = "coxeter.h9_14.record0_a1_source1.checkpoint50_increment_audit.v1"

EXPECTED = {
    48: {
        "digest": "030da2f5300465070a04f77c382eb641ac82e6dbf1bac2276c56fdbadd29bcae",
        "passing": [108, 110, 16, 0],
        "compatible": [108, 1542, 76, 56],
        "orbit_size": 1008,
        "stabilizer": 240,
        "orbit_sha256": "b8c5321b73ab5716d056369fad95d4b5d62e832c044fe8b637d8feb3b37fd001",
    },
    49: {
        "digest": "14161bd66d1218225e6c55aab23b50c03cd9fb61dd407857ec7706509def0659",
        "passing": [108, 350, 61, 17, 0],
        "compatible": [108, 1542, 1282, 32, 73],
        "orbit_size": 2016,
        "stabilizer": 120,
        "orbit_sha256": "5538a6b0026f050a418636cc4543b64c75122a11d11390f4c2e70c8fab9db126",
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

    generators = _source_generators()
    orbit = _core_orbit((additions, restrictions), generators)
    orbit_sha256 = _digest(sorted(orbit))
    passing = [stage["passing_state_count"] for stage in cut_stages]
    compatible = [stage["compatible_extension_count"] for stage in cut_stages]
    checks = {
        "representative_digest_matches": row["representative_digest"]
        == expected["digest"],
        "fixed_A4_occurs_once": supports.count(FIXED_FOUR) == 1,
        "remaining_supports_are_exactly_stored_independent_A5s": (
            tuple(sorted(value for value in supports if value != FIXED_FOUR))
            == additions
            and all(_independent_five(value) for value in additions)
        ),
        "all_actual_vertices_have_rank_nine": all(
            value.bit_count() == 9 for value in actual_vertices
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
        "independent_state_counts_match": passing == expected["passing"]
        and compatible == expected["compatible"],
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
        "support_profile": (
            f"fixed A4 + {len(additions)} A5 supports"
        ),
        "local_actual_vertex_count": len(actual_vertices),
        "independent_passing_state_counts": passing,
        "independent_compatible_extension_counts": compatible,
        "full_source_orbit_size": len(orbit),
        "core_stabilizer_order": EXPECTED_SOURCE_GROUP_ORDER // len(orbit),
        "full_source_orbit_sha256": orbit_sha256,
        "solver_branches": row["solver_branch_count"],
        "solver_conflicts": row["solver_conflict_count"],
        "solver_wall_time_seconds": row["solver_wall_time_seconds"],
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(CHECKPOINT) != CHECKPOINT_SHA256:
        raise ArithmeticError("source-1 checkpoint-50 hash changed")
    if _sha256(BASE_AUDITOR) != BASE_AUDITOR_SHA256:
        raise ArithmeticError("frozen checkpoint-48 auditor hash changed")
    checkpoint = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    rows = checkpoint["iterations"][-2:]
    if [int(row["iteration"]) for row in rows] != [48, 49]:
        raise ArithmeticError("checkpoint does not end in iterations 48 and 49")
    summaries = [_audit_row(row) for row in rows]
    lanner = _audit_lanner_tables()
    checks = {
        "checkpoint_identity_is_record0_a1_source1": (
            checkpoint["record_index"] == 0
            and checkpoint["schedule_a"] == 1
            and checkpoint["source_orbit_index"] == 1
        ),
        "exactly_50_bounded_core_iterations_are_present": (
            len(checkpoint["iterations"]) == 50
            and checkpoint["excluded_local_core_orbit_count"] == 50
        ),
        "checkpoint_remains_OPTIMAL_core_cutoff_not_a_theorem": (
            checkpoint["terminal_solver_status"] == "OPTIMAL"
            and checkpoint["stopped_reason"] == "core_cutoff"
            and checkpoint["fixed_source_compact_coxeter_no_go"] is False
            and checkpoint["theorem_certificate"] is False
        ),
        "both_incremental_core_audits_pass": all(
            all(summary["checks"].values()) for summary in summaries
        ),
        "independent_classical_Lanner_tables_pass": all(
            lanner["checks"].values()
        ),
        "finite_volume_noncompact_is_not_claimed": (
            checkpoint["finite_volume_noncompact_checked"] is False
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"checkpoint-50 incremental audit failed: {checks}")
    return {
        "schema": SCHEMA,
        "scope": (
            "Independent proof for newly appended iterations 48 and 49 and "
            "their complete source-1 orbits. The 50-core checkpoint remains "
            "bounded and is not a fixed-source no-go theorem."
        ),
        "inputs": {
            "checkpoint": str(CHECKPOINT.relative_to(ROOT)),
            "checkpoint_sha256": CHECKPOINT_SHA256,
            "base_independent_auditor": str(BASE_AUDITOR.relative_to(ROOT)),
            "base_independent_auditor_sha256": BASE_AUDITOR_SHA256,
        },
        "bounded_checkpoint_status": {
            "terminal_solver_status": checkpoint["terminal_solver_status"],
            "stopped_reason": checkpoint["stopped_reason"],
            "excluded_local_core_orbits": checkpoint[
                "excluded_local_core_orbit_count"
            ],
            "cumulative_orbit_image_clause_count": sum(
                int(row["full_core_orbit_size"])
                for row in checkpoint["iterations"]
            ),
            "fixed_source_compact_coxeter_no_go": False,
            "theorem_certificate": False,
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
            raise SystemExit("checkpoint-50 incremental archive mismatch")
        print("archive_check=PASS")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
