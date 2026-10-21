#!/usr/bin/env python3
"""Freeze Record 14's compact ``a=0`` and whole-shadow no-go."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from scripts.certify_h9_14_record14_a1_schedule_no_go import (
    EXPECTED_CERTIFICATE_DIGEST as A1_DIGEST,
    OUTPUT as A1_ARCHIVE,
    verify as verify_a1,
)
from scripts.explore_h9_14_record14_a0_downstream import full_source_orbit_key
from scripts.explore_h9_14_record14_a0_symmetry_z3 import enumerate_a0
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    audit_compact_labels,
    finite_code,
    mask,
    placement_options,
)
from scripts.independent_h9_14_record14_a0_audit import replay_a0


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record14_a0_whole_shadow_no_go.json"
SOURCE_PATHS = (
    ROOT / "scripts/explore_h9_14_record14_a0_symmetry_z3.py",
    ROOT / "scripts/explore_h9_14_record14_a1_symmetry_z3.py",
    ROOT / "scripts/explore_h9_14_record14_fixed_orbit_z3.py",
    ROOT / "scripts/independent_h9_14_record14_a0_audit.py",
    ROOT / "scripts/independent_h9_14_record14_a1_audit.py",
    ROOT / "scripts/independent_h9_14_record14_a2_audit.py",
    ROOT / "scripts/explore_h9_14_record14_a0_downstream.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
)
EXPECTED_TARGET_GROUPS = (
    (0, 1, 7, 16, 19, 23, 26, 30, 32, 34, 38, 44, 45, 51),
    (3, 4, 8, 9, 10, 13, 14, 17, 20, 22, 24, 27, 28, 31, 37, 39, 43, 46, 48, 50),
    (2, 5, 6, 11, 12, 15, 18, 21, 25, 29, 33, 35, 36, 40, 41, 42, 47, 49),
)
EXPECTED_LABEL_STATE_SEQUENCES = (
    (420, 96, 144, 16, 0),
    (420, 96, 8, 0),
    (420, 96, 8, 0),
)
EXPECTED_CERTIFICATE_DIGEST = "75fa29552f12864ca668e82550a0d63c9b643d5b114f85a15fd641c4b523cd0e"


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _target_masks(
    raw_targets: Sequence[Sequence[Sequence[int]]],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(sorted(mask(face) for face in family)) for family in raw_targets
    )


def _target_groups(
    targets: Sequence[Sequence[int]],
) -> tuple[tuple[tuple[int, ...], ...], str]:
    grouped: dict[tuple[tuple[int, ...], tuple[int, ...]], list[int]] = {}
    for index, family in enumerate(targets):
        grouped.setdefault(full_source_orbit_key(family), []).append(index)
    groups = tuple(tuple(indices) for _key, indices in sorted(grouped.items()))
    return groups, _digest(sorted(grouped))


def fresh_record() -> dict[str, object]:
    a1 = json.loads(A1_ARCHIVE.read_text(encoding="utf-8"))
    verify_a1(a1)
    production = enumerate_a0()
    independent = replay_a0()
    production_targets = _target_masks(
        production["generator_lex_slice_pure_targets"]
    )
    independent_targets = _target_masks(
        independent["generator_lex_slice_pure_targets"]
    )
    if production_targets != independent_targets:
        raise AssertionError("Record-14 a=0 literal target lists disagree")
    if (
        not production["enumeration_complete"]
        or not independent["enumeration_complete"]
    ):
        raise AssertionError("Record-14 a=0 enumeration cut off")
    if len(production_targets) != 52:
        raise AssertionError("Record-14 a=0 target count changed")

    groups, key_digest = _target_groups(production_targets)
    if groups != EXPECTED_TARGET_GROUPS:
        raise AssertionError("Record-14 a=0 full target quotient changed")
    label_audits = []
    for orbit_index, group in enumerate(groups):
        target_index = group[0]
        audit = audit_compact_labels((), production_targets[target_index])
        sequence = tuple(stage["passing_state_count"] for stage in audit["stages"])
        if sequence != EXPECTED_LABEL_STATE_SEQUENCES[orbit_index]:
            raise AssertionError("Record-14 a=0 compact-label transcript changed")
        label_audits.append(
            {
                "target_orbit_index": orbit_index,
                "representative_generator_lex_target_index": target_index,
                "generator_lex_target_indices_covered": list(group),
                "representative_additional_A5s": production[
                    "generator_lex_slice_pure_targets"
                ][target_index],
                **audit,
            }
        )
        finite_code.cache_clear()
        placement_options.cache_clear()
    if not all(
        audit["unrestricted_compact_label_no_go"] for audit in label_audits
    ):
        raise AssertionError("Record-14 a=0 compact-label survivor appeared")

    evidence = {
        "source_branch_count": 1,
        "source_symmetry_group": "S9 x S2",
        "source_symmetry_order": 725760,
        "additional_A5_count": 8,
        "scheduled_facet_count": 110,
        "candidate_count": production["candidate_count"],
        "base_small_transversal_count": production[
            "base_small_transversal_count"
        ],
        "base_order_five_transversal_count": production[
            "base_order_five_transversal_count"
        ],
        "six_set_constraint_count": production["six_set_constraint_count"],
        "stabilizer_generator_count": production["stabilizer_generator_count"],
        "generator_moved_candidate_pair_counts": production[
            "generator_moved_candidate_pair_counts"
        ],
        "production_z3_generator_lex_model_count": production[
            "generator_lex_slice_model_count"
        ],
        "independent_cp_sat_generator_lex_model_count": independent[
            "generator_lex_slice_model_count"
        ],
        "generator_lex_slice_pure_target_count": len(production_targets),
        "generator_lex_slice_pure_targets": production[
            "generator_lex_slice_pure_targets"
        ],
        "production_independent_literal_equality": True,
        "terminal_z3_status": production["terminal_solver_status"],
        "terminal_cp_sat_status": independent["solver_status"],
        "full_target_orbit_count": len(groups),
        "full_target_orbit_groups": [list(group) for group in groups],
        "target_orbit_size_histogram": dict(
            (str(size), count)
            for size, count in sorted(Counter(map(len, groups)).items())
        ),
        "target_orbit_key_digest": key_digest,
        "unrestricted_compact_label_audits": label_audits,
        "all_target_orbits_unrestricted_compact_label_excluded": True,
    }
    document: dict[str, object] = {
        "kind": "h9_14_record14_a0_whole_shadow_compact_no_go",
        "record_index": 14,
        "graph6": "M???????????o?w??",
        "schedule": {
            "a": 0,
            "missing_independent_fives": 8,
            "scheduled_facet_count": 110,
            "summed_vertex_link_g4": 42,
        },
        "dependencies": {
            "a1_schedule_certificate": {
                "path": str(A1_ARCHIVE.relative_to(ROOT)),
                "sha256": _sha(A1_ARCHIVE),
                "certificate_digest": A1_DIGEST,
            }
        },
        "source_files": {
            str(path.relative_to(ROOT)): _sha(path) for path in SOURCE_PATHS
        },
        "proof_scope": (
            "The a=0 source is unique. Two exact generator-lex Boolean "
            "censuses agree literally on all 52 necessary incidence targets; "
            "every S9 x S2 orbit meets the lossless slice. An exact unordered "
            "incidence-matrix quotient gives three full target orbits, and "
            "the unrestricted compact Lanner/actual-vertex label gate kills "
            "all three. The a=1 dependency closes every other schedule."
        ),
        "evidence": evidence,
        "theorem": {
            "compact_record14_schedule_a0_excluded": True,
            "compact_record14_schedules_closed": [0, 1, 2, 3, 4, 5, 6],
            "remaining_compact_schedules": [],
            "whole_shadow_compact_coxeter_no_go": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    document["certificate_digest"] = _digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record14_a0_whole_shadow_compact_no_go":
        raise ValueError("wrong Record-14 a=0 certificate kind")
    a1 = json.loads(A1_ARCHIVE.read_text(encoding="utf-8"))
    verify_a1(a1)
    dependencies = document.get("dependencies")
    if not isinstance(dependencies, Mapping):
        raise ValueError("Record-14 a=0 dependency map is malformed")
    a1_dep = dependencies.get("a1_schedule_certificate")
    if (
        not isinstance(a1_dep, Mapping)
        or a1_dep.get("sha256") != _sha(A1_ARCHIVE)
        or a1_dep.get("certificate_digest") != A1_DIGEST
    ):
        raise ValueError("Record-14 a=1 dependency changed")
    sources = document.get("source_files")
    if not isinstance(sources, Mapping):
        raise ValueError("Record-14 a=0 source map is malformed")
    for path in SOURCE_PATHS:
        if sources.get(str(path.relative_to(ROOT))) != _sha(path):
            raise ValueError(f"Record-14 a=0 source hash changed: {path.name}")
    evidence = document.get("evidence")
    if not isinstance(evidence, Mapping):
        raise ValueError("Record-14 a=0 evidence is malformed")
    if (
        evidence.get("production_z3_generator_lex_model_count") != 52
        or evidence.get("independent_cp_sat_generator_lex_model_count") != 52
        or evidence.get("generator_lex_slice_pure_target_count") != 52
        or evidence.get("production_independent_literal_equality") is not True
    ):
        raise ValueError("Record-14 a=0 literal equality changed")
    groups = evidence.get("full_target_orbit_groups")
    if not isinstance(groups, list) or tuple(tuple(group) for group in groups) != EXPECTED_TARGET_GROUPS:
        raise ValueError("Record-14 a=0 target quotient changed")
    labels = evidence.get("unrestricted_compact_label_audits")
    if not isinstance(labels, list) or len(labels) != 3:
        raise ValueError("Record-14 a=0 label coverage is malformed")
    if tuple(
        tuple(stage["passing_state_count"] for stage in audit["stages"])
        for audit in labels
    ) != EXPECTED_LABEL_STATE_SEQUENCES:
        raise ValueError("Record-14 a=0 label transcripts changed")
    if not all(
        audit["unrestricted_compact_label_no_go"]
        and audit["terminal_sparse_label_state_count"] == 0
        for audit in labels
    ):
        raise ValueError("Record-14 a=0 compact-label survivor appeared")
    theorem = document.get("theorem")
    if (
        not isinstance(theorem, Mapping)
        or theorem.get("remaining_compact_schedules") != []
        or theorem.get("whole_shadow_compact_coxeter_no_go") is not True
        or theorem.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("Record-14 whole-shadow theorem scope changed")
    stored_digest = document.get("certificate_digest")
    body = dict(document)
    body.pop("certificate_digest", None)
    digest = _digest(body)
    if stored_digest != digest:
        raise ValueError("Record-14 a=0 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_DIGEST not in ("TO_BE_FROZEN", digest):
        raise ValueError("Record-14 a=0 certificate is not frozen")
    return {
        "record_index": 14,
        "schedule_a": 0,
        "generator_lex_target_count": 52,
        "full_target_orbit_count": 3,
        "terminal_compact_label_state_count": 0,
        "closed_schedules": [0, 1, 2, 3, 4, 5, 6],
        "remaining_schedules": [],
        "whole_shadow_compact_coxeter_no_go": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", type=Path)
    args = parser.parse_args(argv)
    document = (
        fresh_record()
        if args.fresh
        else json.loads(OUTPUT.read_text(encoding="utf-8"))
    )
    result = verify(document)
    if args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        result["written"] = str(args.write)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
