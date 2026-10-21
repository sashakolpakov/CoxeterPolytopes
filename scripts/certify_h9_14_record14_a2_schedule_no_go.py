#!/usr/bin/env python3
"""Freeze the complete Record-14 compact ``a=2`` schedule no-go."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from scripts.certify_h9_14_record14_positive_A4_scout import (
    ARCHIVE as POSITIVE_ARCHIVE,
    verify as verify_positive,
)
from scripts.certify_h9_14_record14_upper_schedules import (
    EXPECTED_CERTIFICATE_DIGEST as UPPER_DIGEST,
    OUTPUT as UPPER_ARCHIVE,
    verify as verify_upper,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    audit_compact_labels,
    build_problem,
    enumerate_models,
    finite_code,
    mask,
    placement_options,
    vertices,
)
from scripts.independent_h9_14_record14_a2_audit import replay_orbit


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record14_a2_schedule_no_go.json"
PRODUCTION_SOURCE = ROOT / "scripts/explore_h9_14_record14_fixed_orbit_z3.py"
INDEPENDENT_SOURCE = ROOT / "scripts/independent_h9_14_record14_a2_audit.py"
LABEL_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py"
TARGET_SOURCE_ORBITS = (1, 5, 13, 16)
EXPECTED_TARGET_GROUPS = {
    1: ((0, 1, 2, 3, 4, 5),),
    5: ((0, 1, 2, 3, 4, 5),),
    13: ((0, 1, 2, 3, 4, 5),),
    16: ((0, 1, 2), (3, 4), (5,)),
}
EXPECTED_CERTIFICATE_DIGEST = (
    "d47d7c062ae1b40a3dca8bb1f3ffe6c3f87b5444e5c0d78a8515e2caa4aab7a1"
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _target_masks(raw_targets: Sequence[Sequence[Sequence[int]]]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(tuple(sorted(mask(face) for face in family)) for family in raw_targets)
    )


def _move(value: int, action: Sequence[int]) -> int:
    return sum(1 << action[vertex] for vertex in vertices(value))


def _fixed_a4_actions(fours: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    isolates = tuple(range(3, 12))
    target = set(fours)
    actions = []
    for images in permutations(isolates):
        for swap in (False, True):
            action = list(range(14))
            action[0], action[1] = ((1, 0) if swap else (0, 1))
            for source, image in zip(isolates, images, strict=True):
                action[source] = image
            moved = {_move(face, action) for face in fours}
            if moved == target:
                actions.append(tuple(action))
    return tuple(actions)


def _target_groups(
    fours: Sequence[int], targets: Sequence[Sequence[int]]
) -> tuple[int, tuple[tuple[int, ...], ...], str]:
    actions = _fixed_a4_actions(fours)
    keys = []
    for family in targets:
        keys.append(
            min(
                tuple(sorted(_move(face, action) for face in family))
                for action in actions
            )
        )
    grouped: dict[tuple[int, ...], list[int]] = {}
    for index, key in enumerate(keys):
        grouped.setdefault(key, []).append(index)
    groups = tuple(tuple(values) for _key, values in sorted(grouped.items()))
    return len(actions), groups, _digest(sorted(grouped))


def fresh_record() -> dict[str, object]:
    positive = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(positive)
    upper = json.loads(UPPER_ARCHIVE.read_text(encoding="utf-8"))
    verify_upper(upper)
    branch = next(row for row in positive["branches"] if row["a"] == 2)
    if len(branch["orbits"]) != 36:
        raise AssertionError("Record-14 a=2 source orbit count changed")
    rows = []
    label_audits = []
    total_targets = 0
    for orbit_index, source in enumerate(branch["orbits"]):
        problem = build_problem(2, orbit_index)
        production = enumerate_models(problem)
        independent = replay_orbit(orbit_index)
        production_targets = _target_masks(production["literal_pure_targets"])
        independent_targets = _target_masks(independent["literal_pure_targets"])
        if production_targets != independent_targets:
            raise AssertionError(f"Record-14 a=2 literal mismatch at orbit {orbit_index}")
        if not production["enumeration_complete"] or not independent["enumeration_complete"]:
            raise AssertionError("Record-14 a=2 enumeration cut off")
        target_count = len(production_targets)
        total_targets += target_count
        row = {
            "source_positive_A4_orbit_index": orbit_index,
            "source_orbit_key": source["key"],
            "source_orbit_size": source["size"],
            "representative_A4s": source["representative"],
            "forced_A5_shadow_count": len(problem["shadow"]),
            "additional_A5_count": problem["addition_count"],
            "candidate_count": len(problem["candidates"]),
            "base_small_transversal_count": len(problem["small"]),
            "base_order_five_transversal_count": len(problem["fives"]),
            "six_set_constraint_count": production["six_set_constraint_count"],
            "production_z3_model_count": production["z3_model_count"],
            "independent_cp_sat_model_count": independent["literal_model_count"],
            "literal_pure_target_count": target_count,
            "literal_pure_targets": production["literal_pure_targets"],
            "production_independent_literal_equality": True,
        }
        if target_count:
            action_count, groups, key_digest = _target_groups(
                problem["fours"], production_targets
            )
            if groups != EXPECTED_TARGET_GROUPS[orbit_index]:
                raise AssertionError("Record-14 a=2 target orbit partition changed")
            row.update(
                {
                    "fixed_A4_stabilizer_action_count": action_count,
                    "literal_target_orbit_groups": [list(group) for group in groups],
                    "literal_target_orbit_count": len(groups),
                    "target_orbit_key_digest": key_digest,
                }
            )
            for group_index, group in enumerate(groups):
                target_index = group[0]
                audit = audit_compact_labels(
                    problem["fours"], production_targets[target_index]
                )
                label_audits.append(
                    {
                        "source_positive_A4_orbit_index": orbit_index,
                        "target_orbit_index": group_index,
                        "representative_literal_target_index": target_index,
                        "literal_target_indices_covered": list(group),
                        "representative_additional_A5s": production[
                            "literal_pure_targets"
                        ][target_index],
                        **audit,
                    }
                )
                finite_code.cache_clear()
                placement_options.cache_clear()
        rows.append(row)
    if total_targets != 24:
        raise AssertionError("Record-14 a=2 literal target count changed")
    if len(label_audits) != 6 or not all(
        row["unrestricted_compact_label_no_go"] for row in label_audits
    ):
        raise AssertionError("Record-14 a=2 label closure changed")
    evidence = {
        "source_positive_A4_orbit_count": 36,
        "source_target_count_histogram": dict(
            sorted(Counter(row["literal_pure_target_count"] for row in rows).items())
        ),
        "target_bearing_source_orbits": list(TARGET_SOURCE_ORBITS),
        "literal_pure_target_count": total_targets,
        "literal_target_orbit_count": len(label_audits),
        "rows": rows,
        "unrestricted_compact_label_audits": label_audits,
        "all_target_orbits_unrestricted_compact_label_excluded": True,
    }
    document = {
        "kind": "h9_14_record14_a2_compact_schedule_no_go",
        "record_index": 14,
        "graph6": positive["graph6"],
        "schedule": {
            "a": 2,
            "missing_independent_fives": 20,
            "scheduled_facet_count": 106,
            "summed_vertex_link_g4": 30,
        },
        "dependencies": {
            "positive_A4_archive": {
                "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
                "sha256": _sha(POSITIVE_ARCHIVE),
            },
            "upper_schedule_certificate": {
                "path": str(UPPER_ARCHIVE.relative_to(ROOT)),
                "sha256": _sha(UPPER_ARCHIVE),
                "certificate_digest": UPPER_DIGEST,
            },
        },
        "source_files": {
            str(path.relative_to(ROOT)): _sha(path)
            for path in (PRODUCTION_SOURCE, INDEPENDENT_SOURCE, LABEL_SOURCE)
        },
        "evidence": evidence,
        "theorem": {
            "compact_record14_schedule_a2_excluded": True,
            "compact_record14_schedules_closed_with_upper_certificate": [2, 3, 4, 5, 6],
            "remaining_compact_schedules": [0, 1],
            "whole_shadow_compact_coxeter_no_go": False,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    document["certificate_digest"] = _digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record14_a2_compact_schedule_no_go":
        raise ValueError("wrong Record-14 a=2 certificate kind")
    positive = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(positive)
    upper = json.loads(UPPER_ARCHIVE.read_text(encoding="utf-8"))
    verify_upper(upper)
    dependencies = document.get("dependencies")
    if not isinstance(dependencies, Mapping):
        raise ValueError("Record-14 a=2 dependency map is malformed")
    positive_dep = dependencies.get("positive_A4_archive")
    upper_dep = dependencies.get("upper_schedule_certificate")
    if not isinstance(positive_dep, Mapping) or positive_dep.get("sha256") != _sha(POSITIVE_ARCHIVE):
        raise ValueError("Record-14 a=2 positive source changed")
    if (
        not isinstance(upper_dep, Mapping)
        or upper_dep.get("sha256") != _sha(UPPER_ARCHIVE)
        or upper_dep.get("certificate_digest") != UPPER_DIGEST
    ):
        raise ValueError("Record-14 upper dependency changed")
    sources = document.get("source_files")
    if not isinstance(sources, Mapping):
        raise ValueError("Record-14 a=2 source map is malformed")
    for path in (PRODUCTION_SOURCE, INDEPENDENT_SOURCE, LABEL_SOURCE):
        if sources.get(str(path.relative_to(ROOT))) != _sha(path):
            raise ValueError(f"Record-14 a=2 source hash changed: {path.name}")
    evidence = document.get("evidence")
    if not isinstance(evidence, Mapping):
        raise ValueError("Record-14 a=2 evidence is malformed")
    rows = evidence.get("rows")
    labels = evidence.get("unrestricted_compact_label_audits")
    if not isinstance(rows, list) or len(rows) != 36 or not isinstance(labels, list):
        raise ValueError("Record-14 a=2 coverage is malformed")
    branch = next(row for row in positive["branches"] if row["a"] == 2)
    if [
        (
            row["source_positive_A4_orbit_index"],
            row["source_orbit_key"],
            row["source_orbit_size"],
            row["representative_A4s"],
        )
        for row in rows
    ] != [
        (index, source["key"], source["size"], source["representative"])
        for index, source in enumerate(branch["orbits"])
    ]:
        raise ValueError("Record-14 a=2 literal source coverage changed")
    if [row["source_positive_A4_orbit_index"] for row in rows if row["literal_pure_target_count"]] != list(TARGET_SOURCE_ORBITS):
        raise ValueError("Record-14 a=2 target-bearing sources changed")
    if sum(row["literal_pure_target_count"] for row in rows) != 24:
        raise ValueError("Record-14 a=2 literal target total changed")
    if any(
        not row["production_independent_literal_equality"]
        or row["production_z3_model_count"] != row["independent_cp_sat_model_count"]
        or row["production_z3_model_count"] != row["literal_pure_target_count"]
        for row in rows
    ):
        raise ValueError("Record-14 a=2 independent equality changed")
    if len(labels) != 6 or not all(row["unrestricted_compact_label_no_go"] for row in labels):
        raise ValueError("Record-14 a=2 label closure changed")
    theorem = document.get("theorem")
    if not isinstance(theorem, Mapping) or theorem.get("remaining_compact_schedules") != [0, 1]:
        raise ValueError("Record-14 a=2 theorem scope changed")
    stored_digest = document.get("certificate_digest")
    body = dict(document)
    body.pop("certificate_digest", None)
    digest = _digest(body)
    if stored_digest != digest:
        raise ValueError("Record-14 a=2 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_DIGEST not in ("TO_BE_FROZEN", digest):
        raise ValueError("Record-14 a=2 certificate is not frozen")
    return {
        "record_index": 14,
        "schedule_a": 2,
        "source_orbit_count": 36,
        "literal_target_count": 24,
        "literal_target_orbit_count": 6,
        "terminal_compact_label_state_count": 0,
        "closed_schedules": [2, 3, 4, 5, 6],
        "remaining_schedules": [0, 1],
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
