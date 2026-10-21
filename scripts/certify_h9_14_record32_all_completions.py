#!/usr/bin/env python3
"""Certify the independently replayed record-32 compact whole-shadow no-go."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Mapping, Sequence

from scripts.certify_h9_14_record32_all_completions_candidate import (
    OUTPUT as CANDIDATE_OUTPUT,
    verify as verify_candidate,
)
from scripts.certify_h9_14_record32_coxeter_labels import (
    OUTPUT as SELECTED_LABEL_OUTPUT,
    verify as verify_selected_label,
)
from scripts.certify_h9_14_record32_whole_shadow_hostile import (
    DEFAULT_ARCHIVE as HOSTILE_OUTPUT,
    validate_archive as validate_hostile,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record32_all_completions.json"
CERTIFIER_SOURCE = Path(__file__).resolve()
EXPECTED_HOSTILE_FILE_SHA256 = "381c97ac36f85d146f96601e19de0200efe68cd48bb7ad98728035f4cc0a0417"
EXPECTED_CERTIFICATE_SHA256 = "486b063921841de9aa3c4611ea171801c3368c1c1c7c5cdcef938fabedfb69d9"
EXPECTED_ARCHIVE_FILE_SHA256 = "600d107e9739085897349cf905d411df1af260afb75facb589c0edadf74a52b0"


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _normalized(value: object) -> object:
    return json.loads(json.dumps(value, sort_keys=True))


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    payload = json.loads(
        json.dumps(payload, sort_keys=True, separators=(",", ":"))
    )
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _normalized_certifier_source(source: str) -> str:
    """Remove only the three unavoidable outer-freeze literal values."""

    normalized = source
    for name in (
        "EXPECTED_HOSTILE_FILE_SHA256",
        "EXPECTED_CERTIFICATE_SHA256",
        "EXPECTED_ARCHIVE_FILE_SHA256",
    ):
        pattern = re.compile(
            rf'^(?P<prefix>{name}[ \t]*=[ \t]*)"[^"\n]*"(?P<suffix>[ \t]*)$',
            re.MULTILINE,
        )
        normalized, count = pattern.subn(
            lambda match: (
                f'{match.group("prefix")}"<OUTER_FREEZE>"'
                f'{match.group("suffix")}'
            ),
            normalized,
        )
        if count != 1:
            raise ValueError(
                f"expected one top-level {name} assignment, got {count}"
            )
    return normalized


def _certifier_source_record() -> dict[str, str]:
    normalized = _normalized_certifier_source(
        CERTIFIER_SOURCE.read_text(encoding="utf-8")
    )
    return {
        "path": str(CERTIFIER_SOURCE.relative_to(ROOT)),
        "sha256_excluding_outer_freeze_assignments": sha256(
            normalized.encode("utf-8")
        ).hexdigest(),
    }


def build() -> dict[str, object]:
    candidate = json.loads(CANDIDATE_OUTPUT.read_text(encoding="utf-8"))
    hostile = json.loads(HOSTILE_OUTPUT.read_text(encoding="utf-8"))
    selected = json.loads(SELECTED_LABEL_OUTPUT.read_text(encoding="utf-8"))
    candidate_summary = verify_candidate(candidate)
    validate_hostile(hostile)
    selected_summary = verify_selected_label(selected)
    sources = hostile["independent_sources"]
    hostile_file_sha256 = _file_sha256(HOSTILE_OUTPUT)
    if EXPECTED_HOSTILE_FILE_SHA256 not in (
        "TO_BE_FROZEN",
        hostile_file_sha256,
    ):
        raise ValueError("record-32 hostile archive file hash changed")
    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record32_whole_shadow_compact_coxeter_no_go",
        "scope": {
            "dimension": 9,
            "facet_count": 14,
            "record_index": 32,
            "graph6": "M?????????O?O?g??",
            "all_six_dehn_sommerville_schedules_covered": True,
            "compact_only": True,
            "abstract_targets_are_only_necessary_incidence_objects": True,
            "ordinary_polytopality_of_abstract_targets_claimed": False,
            "PL_sphere_status_of_abstract_targets_claimed": False,
            "ideal_or_parabolic_faces_enumerated": False,
            "finite_volume_noncompact_chamber_excluded": False,
        },
        "historical_production_candidate": {
            "path": str(CANDIDATE_OUTPUT.relative_to(ROOT)),
            "file_sha256": _file_sha256(CANDIDATE_OUTPUT),
            "production_census_complete": candidate_summary[
                "production_candidate_complete"
            ],
            "independent_replay_complete_at_that_stage": candidate_summary[
                "independent_replay_complete"
            ],
            "theorem_promoted_at_that_stage": candidate_summary[
                "whole_graph_shadow_compact_coxeter_theorem"
            ],
            "superseded_by_this_independently_replayed_certificate": True,
        },
        "independent_hostile_replay": {
            "path": str(HOSTILE_OUTPUT.relative_to(ROOT)),
            "file_sha256": hostile_file_sha256,
            "source_sha256": {
                name: row["sha256"] for name, row in sorted(sources.items())
            },
            "independent_replay_complete": hostile["status"][
                "independent_replay_complete"
            ],
            "whole_graph_shadow_compact_coxeter_no_go": hostile["status"][
                "whole_graph_shadow_compact_coxeter_no_go"
            ],
        },
        "selected_a2_label_theorem_chain": {
            "path": str(SELECTED_LABEL_OUTPUT.relative_to(ROOT)),
            "file_sha256": _file_sha256(SELECTED_LABEL_OUTPUT),
            "certificate_sha256": selected_summary["certificate_sha256"],
            "unrestricted_compact_coxeter_no_go": selected_summary[
                "unrestricted_compact_coxeter_no_go"
            ],
        },
        "certifier_source": _certifier_source_record(),
        "branch_closure": {
            "dehn_sommerville_schedules_a_b_facets_link_g4_sum": candidate[
                "schedules"
            ],
            "strict_compact_excluded_a": [4, 5],
            "a3": {
                "strict_shadow_orbits": hostile["a3"]["strict_shadow_orbit_count"],
                "cap_feasible_addition_pairs": hostile["a3"]["cap_feasible"],
                "no_small_transversal_pairs": hostile["a3"]["no_small"],
                "closed_before_labels": hostile["a3"]["no_small"] == 0,
            },
            "a2": {
                "strict_shadow_orbits": len(
                    hostile["a2"]["representative_transcripts"]
                ),
                "scheduled_target_orbits": sum(
                    row["targets"]
                    for row in hostile["a2"]["representative_transcripts"]
                ),
                "sole_target_is_selected_a2_nerve": hostile["a2"][
                    "sole_target_is_older_selected_record32_nerve"
                ],
                "closed_by_selected_label_theorem": selected_summary[
                    "unrestricted_compact_coxeter_no_go"
                ],
            },
            "a1": {
                "strict_shadow_A4_orbits": len(
                    hostile["a1"]["fixed_A4_transcripts"]
                ),
                "scheduled_target_orbits": sum(
                    row["target_orbits"]
                    for row in hostile["a1"]["fixed_A4_transcripts"]
                ),
                "ridge_valid_target_orbits": hostile["a1"][
                    "ridge_valid_target_count"
                ],
                "sparse_state_counts": [
                    row["sparse_states"] for row in hostile["a1"]["cross_gates"]
                ],
                "single_m3_cross_test_counts": [
                    row["single_m3_tests"]
                    for row in hostile["a1"]["cross_gates"]
                ],
                "single_m3_cross_passing_counts": [
                    row["passing"] for row in hostile["a1"]["cross_gates"]
                ],
                "closed_by_unrestricted_cross_gates": all(
                    row["passing"] == 0 for row in hostile["a1"]["cross_gates"]
                ),
            },
            "a0": {
                "production_literal_representatives": hostile["a0"][
                    "production_literal_families"
                ],
                "independent_literal_representatives": hostile["a0"][
                    "independent_replay"
                ]["target_family_count"],
                "full_action_target_orbits": hostile["a0"][
                    "full_action_orbit_comparison"
                ]["production_orbit_count"],
                "target_orbit_key_sha256": hostile["a0"][
                    "full_action_orbit_comparison"
                ]["orbit_key_sha256"],
                "pure_blocker_families": hostile["a0"]["pure_blocker_families"],
                "ridge_valid_families": len(hostile["a0"]["ridge_valid_indices"]),
                "unrestricted_zero_terminal_families": sum(
                    value == 0
                    for value in hostile["a0"][
                        "ridge_valid_unrestricted_terminal_state_counts"
                    ]
                ),
                "closed_by_blocker_ridge_and_label_gates": all(
                    value == 0
                    for value in hostile["a0"][
                        "ridge_valid_unrestricted_terminal_state_counts"
                    ]
                ),
            },
        },
        "theorem": {
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    document["certificate_sha256"] = _canonical_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("schema_version") != "1.0.0" or document.get("kind") != (
        "h9_14_record32_whole_shadow_compact_coxeter_no_go"
    ):
        raise ValueError("wrong record-32 whole-shadow certificate kind")
    digest = _canonical_digest(document)
    if document.get("certificate_sha256") != digest:
        raise ValueError("record-32 whole-shadow certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 not in ("TO_BE_FROZEN", digest):
        raise ValueError("record-32 whole-shadow certificate is not frozen")
    rebuilt = build()
    if _normalized(document) != _normalized(rebuilt):
        raise ValueError("record-32 whole-shadow certificate does not replay")
    scope = document["scope"]
    theorem = document["theorem"]
    branches = document["branch_closure"]
    if (
        not scope["compact_only"]
        or scope["ordinary_polytopality_of_abstract_targets_claimed"]
        or scope["PL_sphere_status_of_abstract_targets_claimed"]
        or scope["ideal_or_parabolic_faces_enumerated"]
        or scope["finite_volume_noncompact_chamber_excluded"]
        or theorem
        != {
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
        }
    ):
        raise ValueError("record-32 whole-shadow theorem scope changed")
    if branches["dehn_sommerville_schedules_a_b_facets_link_g4_sum"] != [
        [0, 6, 114, 36],
        [1, 12, 112, 30],
        [2, 18, 110, 24],
        [3, 24, 108, 18],
        [4, 30, 106, 12],
        [5, 36, 104, 6],
    ] or branches["strict_compact_excluded_a"] != [4, 5]:
        raise ValueError("record-32 exact schedule partition changed")
    return {
        "record_index": 32,
        "a0_full_action_target_orbits": document["branch_closure"]["a0"][
            "full_action_target_orbits"
        ],
        "a1_ridge_valid_target_orbits": document["branch_closure"]["a1"][
            "ridge_valid_target_orbits"
        ],
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    args = parser.parse_args(argv)
    if args.build:
        document = build()
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    archive_hash = _file_sha256(OUTPUT)
    if EXPECTED_ARCHIVE_FILE_SHA256 not in ("TO_BE_FROZEN", archive_hash):
        raise ValueError("record-32 whole-shadow archive file hash changed")
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
