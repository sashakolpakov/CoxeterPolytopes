#!/usr/bin/env python3
"""Verify or recompute the record-32 whole-shadow production candidate.

The archived document is intentionally *not* a theorem certificate: its
independent-replay and theorem-promotion flags must remain false here.  The
optional expensive replay reconstructs every Python-side branch and label
gate.  The independent no-project-import replay is a separate promotion gate.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record32_coxeter_labels import (
    audit_record32_coxeter_labels,
)
from scripts.explore_h9_14_record32_a0_target_gate import audit_family
from scripts.explore_h9_14_record32_a1_completions import (
    A1_REPRESENTATIVES,
    audit_representative as audit_a1_representative,
)
from scripts.explore_h9_14_record32_a1_cross_gate import (
    TARGETS as A1_CROSS_TARGETS,
    audit_target as audit_a1_cross_target,
)
from scripts.explore_h9_14_record32_a2_additions import (
    A2_REPRESENTATIVES,
    audit_representative as audit_a2_representative,
)
from scripts.explore_h9_14_record32_a3_additions import audit_a3_representative
from scripts.explore_h9_14_record32_all_completions import (
    EXPECTED_SCHEDULES,
    positive_branch_scout,
    schedule_audit,
    state_space_audit,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record32_whole_shadow_candidate.json"
A0_TARGETS = ROOT / "data/enumerations/h9_14_record32_a0_all_targets_scout.json"
A0_CPP = ROOT / "scripts/explore_h9_14_record32_a0.cpp"
EXPECTED_CANDIDATE_FILE_SHA256 = (
    "54cd477eb4621c6744a8583c490acfadb2997a64c37feeba20c21b83d366c802"
)
EXPECTED_A0_TARGET_FILE_SHA256 = (
    "ae84f9c8dde0f66cd4d8fc8c0b9dee2afab985c7d484f2e052709ad3c3734938"
)
EXPECTED_A0_LABEL_TRANSCRIPT_SHA256 = (
    "02bae944190d58ce7a309c66ce893a0705ceb03721bfea040ae495373e3d43f3"
)
EXPECTED_A0_CPP_SHA256 = (
    "90ae323727f70aa1175dc903e801cfe8a3666a0a7a5f3b7fb226d69d014bcad7"
)


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _projection_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _require_mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"record-32 candidate {name} is malformed")
    return value


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("schema_version") != "0.1.0" or document.get("kind") != (
        "h9_14_record32_whole_shadow_compact_candidate"
    ):
        raise ValueError("wrong record-32 candidate kind")
    status = _require_mapping(document.get("status"), "status")
    scope = _require_mapping(document.get("scope"), "scope")
    conclusion = _require_mapping(
        document.get("candidate_conclusion"), "candidate conclusion"
    )
    if status != {
        "production_census_complete": True,
        "independent_replay_complete": False,
        "theorem_promoted": False,
    }:
        raise ValueError("record-32 candidate was promoted without replay")
    if (
        scope.get("record_index") != 32
        or scope.get("graph6") != "M?????????O?O?g??"
        or scope.get("compact_only") is not True
        or scope.get("abstract_targets_are_only_necessary_incidence_objects")
        is not True
        or scope.get("ordinary_polytopality_of_abstract_targets_claimed") is not False
        or scope.get("PL_sphere_status_of_abstract_targets_claimed") is not False
        or scope.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("record-32 candidate scope changed")
    if document.get("schedules") != [list(row) for row in EXPECTED_SCHEDULES]:
        raise ValueError("record-32 schedule archive changed")
    if document.get("strict_compact_surviving_a") != [0, 1, 2, 3] or document.get(
        "strict_compact_excluded_a"
    ) != [4, 5]:
        raise ValueError("record-32 strict schedule cut changed")

    a0 = _require_mapping(document.get("a0"), "a0")
    if (
        a0.get("active_row_multisets_after_caps") != 22_513
        or a0.get("active_row_multisets_after_cover_filter") != 14_731
        or a0.get("column_multiset_path_count") != 179_785_655
        or a0.get("literal_small_transversal_failure_count") != 0
        or a0.get("blocker_valid_path_count") != 176_272_067
        or a0.get("scheduled_target_path_count") != 1_263
        or a0.get("exact_target_family_count") != 1_025
        or a0.get("pure_blocker_count") != 941
        or a0.get("ridge_pseudomanifold_indices")
        != [6, 27, 43, 57, 132, 139, 436, 439, 481, 574, 626, 629]
        or a0.get("unrestricted_label_no_go_indices")
        != a0.get("ridge_pseudomanifold_indices")
        or _projection_sha256(
            a0.get("label_transcripts_support_compatible_passing")
        )
        != EXPECTED_A0_LABEL_TRANSCRIPT_SHA256
        or a0.get("target_file_sha256") != EXPECTED_A0_TARGET_FILE_SHA256
        or _file_sha256(A0_TARGETS) != EXPECTED_A0_TARGET_FILE_SHA256
    ):
        raise ValueError("record-32 a0 candidate archive changed")
    if _file_sha256(A0_CPP) != EXPECTED_A0_CPP_SHA256:
        raise ValueError("standalone record-32 a0 source changed")
    if "from coxeter_polytopes" in A0_CPP.read_text(encoding="utf-8"):
        raise ValueError("standalone record-32 a0 source gained a project import")

    a1 = _require_mapping(document.get("a1"), "a1")
    orbit_audits = a1.get("orbit_audits")
    cross_gates = a1.get("ridge_target_cross_gates")
    if not isinstance(orbit_audits, list) or not isinstance(cross_gates, list):
        raise ValueError("record-32 a1 candidate archive is malformed")
    expected_a1 = (
        (10, 2, 664, 1_161, 1_102, 9, 1, 1),
        (9, 3, 672, 43_229, 42_865, 271, 3, 1),
        (8, 4, 672, 2_323_448, 2_278_508, 12_873, 0, 0),
        (8, 4, 1_252, 6_068_459, 5_934_396, 19_836, 6, 0),
    )
    actual_a1 = tuple(
        (
            row.get("shadow"),
            row.get("additional_A5s"),
            row.get("candidates"),
            row.get("column_paths"),
            row.get("distinct_row_paths"),
            row.get("no_small_paths"),
            row.get("target_orbits"),
            row.get("ridge_valid_target_orbits"),
        )
        for row in map(lambda value: _require_mapping(value, "a1 orbit"), orbit_audits)
    )
    if actual_a1 != expected_a1:
        raise ValueError("record-32 a1 orbit census changed")
    if [
        (
            _require_mapping(row, "a1 cross gate").get("sparse_state_count"),
            _require_mapping(row, "a1 cross gate").get("single_m3_test_count"),
            _require_mapping(row, "a1 cross gate").get("single_m3_passing_count"),
            _require_mapping(row, "a1 cross gate").get(
                "unrestricted_compact_label_no_go"
            ),
        )
        for row in cross_gates
    ] != [(10_368, 207_360, 0, True), (864, 17_280, 0, True)]:
        raise ValueError("record-32 a1 cross gates changed")

    a2 = _require_mapping(document.get("a2"), "a2")
    a3 = _require_mapping(document.get("a3"), "a3")
    if (
        a2.get("strict_shadow_orbit_count") != 6
        or a2.get("scheduled_target_count_across_orbit_representatives") != 1
        or a2.get("unique_target_is_selected_record32_nerve") is not True
        or a2.get("selected_nerve_unrestricted_compact_label_no_go") is not True
        or a3.get("strict_shadow_orbit_count") != 1
        or a3.get("cap_feasible_count") != 540
        or a3.get("no_small_transversal_count") != 0
        or a3.get("scheduled_target_count") != 0
    ):
        raise ValueError("record-32 positive candidate closure changed")
    if conclusion != {
        "production_whole_graph_shadow_compact_coxeter_excluded": True,
        "independently_replayed": False,
        "whole_graph_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }:
        raise ValueError("record-32 candidate conclusion was over-promoted")
    return {
        "record_index": 32,
        "production_candidate_complete": True,
        "independent_replay_complete": False,
        "whole_graph_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "candidate_file_sha256": _file_sha256(OUTPUT),
    }


def recompute() -> dict[str, object]:
    """Expensively replay all Python-side production gates."""

    schedules = schedule_audit()
    state = state_space_audit()
    positive = positive_branch_scout(3)
    target_document = json.loads(A0_TARGETS.read_text(encoding="utf-8"))
    families = tuple(
        tuple(map(int, row)) for row in target_document["target_family_masks"]
    )
    incidence = tuple(
        audit_family(index, family, 250_000, True)
        for index, family in enumerate(families)
    )
    ridge_indices = tuple(
        row["family_index"] for row in incidence if row["all_ridge_degrees_two"]
    )
    label_audits = tuple(
        audit_family(index, families[index], 250_000, False)
        for index in ridge_indices
    )
    a1 = tuple(
        audit_a1_representative(
            index, representative, incidence_only=True
        )
        for index, representative in enumerate(A1_REPRESENTATIVES)
    )
    cross = tuple(audit_a1_cross_target(target) for target in A1_CROSS_TARGETS)
    a2 = tuple(
        audit_a2_representative(index, representative)
        for index, representative in enumerate(A2_REPRESENTATIVES)
    )
    a3 = audit_a3_representative()
    selected_label = audit_record32_coxeter_labels()
    replay = {
        "schedules": schedules,
        "state_space": state,
        "positive": positive,
        "a0": {
            "family_count": len(families),
            "pure_count": sum(row["pure_order_five_blocker"] for row in incidence),
            "ridge_indices": list(ridge_indices),
            "label_no_go_count": sum(
                row["sparse_label_gate"]["unrestricted_compact_label_no_go"]
                for row in label_audits
            ),
        },
        "a1": {
            "target_orbit_counts": [
                row["fixed_A4_stabilizer_target_orbit_count"] for row in a1
            ],
            "ridge_counts": [
                sum(target["pure_ridge_pseudomanifold"] for target in row["target_audits"])
                for row in a1
            ],
            "cross_single_m3_passing_counts": [
                row["single_m3_passing_count"] for row in cross
            ],
        },
        "a2_target_counts": [row["scheduled_target_count"] for row in a2],
        "a3_target_count": a3["scheduled_target_count"],
        "selected_a2_label_verified": selected_label.verified,
    }
    if replay["a0"] != {
        "family_count": 1_025,
        "pure_count": 941,
        "ridge_indices": [6, 27, 43, 57, 132, 139, 436, 439, 481, 574, 626, 629],
        "label_no_go_count": 12,
    }:
        raise AssertionError("record-32 production a0 replay changed")
    if replay["a1"] != {
        "target_orbit_counts": [1, 3, 0, 6],
        "ridge_counts": [1, 1, 0, 0],
        "cross_single_m3_passing_counts": [0, 0],
    }:
        raise AssertionError("record-32 production a1 replay changed")
    if replay["a2_target_counts"] != [0, 1, 0, 0, 0, 0] or (
        replay["a3_target_count"] != 0
        or replay["selected_a2_label_verified"] is not True
    ):
        raise AssertionError("record-32 production positive replay changed")
    return replay


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args(argv)
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    result = verify(document)
    archive_digest = _file_sha256(OUTPUT)
    if EXPECTED_CANDIDATE_FILE_SHA256 not in ("TO_BE_FROZEN", archive_digest):
        raise ValueError("record-32 candidate archive file hash changed")
    if args.recompute:
        result["production_recompute"] = recompute()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
