#!/usr/bin/env python3
"""Freeze the complete Record-14 compact ``a=1`` schedule no-go."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
from typing import Mapping, Sequence

from scripts.certify_h9_14_record14_a2_schedule_no_go import (
    EXPECTED_CERTIFICATE_DIGEST as A2_DIGEST,
    OUTPUT as A2_ARCHIVE,
    verify as verify_a2,
)
from scripts.certify_h9_14_record14_positive_A4_scout import (
    ARCHIVE as POSITIVE_ARCHIVE,
    verify as verify_positive,
)
from scripts.explore_h9_14_record14_a1_symmetry_z3 import (
    _move,
    enumerate_orbit,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    audit_compact_labels,
    build_problem,
    finite_code,
    mask,
    placement_options,
)
from scripts.independent_h9_14_record14_a1_audit import replay_orbit


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record14_a1_schedule_no_go.json"
PRODUCTION_SOURCE = ROOT / "scripts/explore_h9_14_record14_a1_symmetry_z3.py"
INDEPENDENT_SOURCE = ROOT / "scripts/independent_h9_14_record14_a1_audit.py"
LABEL_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py"
TARGET_SOURCE_ORBITS = (0, 3)
EXPECTED_TARGET_GROUPS = {
    0: ((0, 10), (1, 3), (2, 8), (4, 7), (5, 6), (9,)),
    3: ((0, 5), (1, 3), (2, 7), (4, 6)),
}
EXPECTED_CERTIFICATE_DIGEST = (
    "7c579a3ac2260275a504d7a5a2fa2cf851cb2eb379c040069e130a07cc1da4f3"
)


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


def _fixed_source_actions(four: int) -> tuple[tuple[int, ...], ...]:
    """Enumerate the full Record-14 graph stabilizer of one A4 row."""

    inside = tuple(vertex for vertex in range(3, 12) if four & (1 << vertex))
    outside = tuple(vertex for vertex in range(3, 12) if not four & (1 << vertex))
    active_actions = ((0, 1), (1, 0)) if (
        bool(four & 1) == bool(four & 2)
    ) else ((0, 1),)
    actions = []
    for inside_images, outside_images, active in product(
        permutations(inside), permutations(outside), active_actions
    ):
        action = list(range(14))
        action[0], action[1] = active
        for source, image in zip(inside, inside_images, strict=True):
            action[source] = image
        for source, image in zip(outside, outside_images, strict=True):
            action[source] = image
        moved = _move(four, action)
        if moved != four:
            raise AssertionError("Record-14 fixed-source action escaped the A4 row")
        actions.append(tuple(action))
    return tuple(actions)


def _target_groups(
    four: int, targets: Sequence[Sequence[int]]
) -> tuple[int, tuple[tuple[int, ...], ...], str]:
    actions = _fixed_source_actions(four)
    keys = [
        min(
            tuple(sorted(_move(face, action) for face in family))
            for action in actions
        )
        for family in targets
    ]
    grouped: dict[tuple[int, ...], list[int]] = {}
    for index, key in enumerate(keys):
        grouped.setdefault(key, []).append(index)
    groups = tuple(tuple(values) for _key, values in sorted(grouped.items()))
    return len(actions), groups, _digest(sorted(grouped))


def fresh_record() -> dict[str, object]:
    positive = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(positive)
    a2 = json.loads(A2_ARCHIVE.read_text(encoding="utf-8"))
    verify_a2(a2)
    branch = next(row for row in positive["branches"] if row["a"] == 1)
    if len(branch["orbits"]) != 5:
        raise AssertionError("Record-14 a=1 source orbit count changed")

    rows = []
    label_audits = []
    total_slice_targets = 0
    for orbit_index, source in enumerate(branch["orbits"]):
        problem = build_problem(1, orbit_index)
        production = enumerate_orbit(orbit_index)
        independent = replay_orbit(orbit_index)
        production_targets = _target_masks(
            production["generator_lex_slice_pure_targets"]
        )
        independent_targets = _target_masks(
            independent["generator_lex_slice_pure_targets"]
        )
        if production_targets != independent_targets:
            raise AssertionError(
                f"Record-14 a=1 generator-lex target mismatch at orbit {orbit_index}"
            )
        if (
            not production["enumeration_complete"]
            or not independent["enumeration_complete"]
        ):
            raise AssertionError("Record-14 a=1 enumeration cut off")
        target_count = len(production_targets)
        total_slice_targets += target_count
        row: dict[str, object] = {
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
            "source_stabilizer_order": production["source_stabilizer_order"],
            "stabilizer_generator_count": production["stabilizer_generator_count"],
            "production_z3_generator_lex_model_count": production[
                "generator_lex_slice_model_count"
            ],
            "independent_cp_sat_generator_lex_model_count": independent[
                "generator_lex_slice_model_count"
            ],
            "generator_lex_slice_pure_target_count": target_count,
            "generator_lex_slice_pure_targets": production[
                "generator_lex_slice_pure_targets"
            ],
            "production_independent_literal_equality": True,
            "terminal_z3_status": production["terminal_solver_status"],
            "terminal_cp_sat_status": independent["solver_status"],
        }
        if target_count:
            action_count, groups, key_digest = _target_groups(
                problem["fours"][0], production_targets
            )
            if tuple(sorted(groups)) != tuple(
                sorted(EXPECTED_TARGET_GROUPS[orbit_index])
            ):
                raise AssertionError("Record-14 a=1 target orbit partition changed")
            row.update(
                {
                    "full_fixed_A4_stabilizer_action_count": action_count,
                    "generator_lex_target_orbit_groups": [
                        list(group) for group in groups
                    ],
                    "target_orbit_count": len(groups),
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
        rows.append(row)

    if total_slice_targets != 19:
        raise AssertionError("Record-14 a=1 generator-lex target total changed")
    if len(label_audits) != 10 or not all(
        row["unrestricted_compact_label_no_go"] for row in label_audits
    ):
        raise AssertionError("Record-14 a=1 label closure changed")

    evidence = {
        "source_positive_A4_orbit_count": 5,
        "source_target_count_histogram": dict(
            (str(count), multiplicity)
            for count, multiplicity in sorted(
                Counter(
                    row["generator_lex_slice_pure_target_count"] for row in rows
                ).items()
            )
        ),
        "target_bearing_source_orbits": list(TARGET_SOURCE_ORBITS),
        "generator_lex_slice_pure_target_count": total_slice_targets,
        "full_target_orbit_count": len(label_audits),
        "rows": rows,
        "unrestricted_compact_label_audits": label_audits,
        "all_target_orbits_unrestricted_compact_label_excluded": True,
    }
    document: dict[str, object] = {
        "kind": "h9_14_record14_a1_compact_schedule_no_go",
        "record_index": 14,
        "graph6": positive["graph6"],
        "schedule": {
            "a": 1,
            "missing_independent_fives": 14,
            "scheduled_facet_count": 108,
            "summed_vertex_link_g4": 36,
        },
        "dependencies": {
            "positive_A4_archive": {
                "path": str(POSITIVE_ARCHIVE.relative_to(ROOT)),
                "sha256": _sha(POSITIVE_ARCHIVE),
            },
            "a2_schedule_certificate": {
                "path": str(A2_ARCHIVE.relative_to(ROOT)),
                "sha256": _sha(A2_ARCHIVE),
                "certificate_digest": A2_DIGEST,
            },
        },
        "source_files": {
            str(path.relative_to(ROOT)): _sha(path)
            for path in (PRODUCTION_SOURCE, INDEPENDENT_SOURCE, LABEL_SOURCE)
        },
        "proof_scope": (
            "Two exact Boolean models impose the complete compact a=1 "
            "incidence, blocker, facet-count, and ridge necessities. Their "
            "generator-lex slices agree literally; every finite group orbit "
            "meets each slice. Full fixed-source stabilizer quotienting then "
            "gives ten target orbits, all killed by the unrestricted compact "
            "Lanner/actual-vertex label gate."
        ),
        "evidence": evidence,
        "theorem": {
            "compact_record14_schedule_a1_excluded": True,
            "compact_record14_schedules_closed_with_dependencies": [1, 2, 3, 4, 5, 6],
            "remaining_compact_schedules": [0],
            "whole_shadow_compact_coxeter_no_go": False,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    document["certificate_digest"] = _digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record14_a1_compact_schedule_no_go":
        raise ValueError("wrong Record-14 a=1 certificate kind")
    positive = json.loads(POSITIVE_ARCHIVE.read_text(encoding="utf-8"))
    verify_positive(positive)
    a2 = json.loads(A2_ARCHIVE.read_text(encoding="utf-8"))
    verify_a2(a2)
    dependencies = document.get("dependencies")
    if not isinstance(dependencies, Mapping):
        raise ValueError("Record-14 a=1 dependency map is malformed")
    positive_dep = dependencies.get("positive_A4_archive")
    a2_dep = dependencies.get("a2_schedule_certificate")
    if (
        not isinstance(positive_dep, Mapping)
        or positive_dep.get("sha256") != _sha(POSITIVE_ARCHIVE)
    ):
        raise ValueError("Record-14 a=1 positive source changed")
    if (
        not isinstance(a2_dep, Mapping)
        or a2_dep.get("sha256") != _sha(A2_ARCHIVE)
        or a2_dep.get("certificate_digest") != A2_DIGEST
    ):
        raise ValueError("Record-14 a=2 dependency changed")
    sources = document.get("source_files")
    if not isinstance(sources, Mapping):
        raise ValueError("Record-14 a=1 source map is malformed")
    for path in (PRODUCTION_SOURCE, INDEPENDENT_SOURCE, LABEL_SOURCE):
        if sources.get(str(path.relative_to(ROOT))) != _sha(path):
            raise ValueError(f"Record-14 a=1 source hash changed: {path.name}")
    evidence = document.get("evidence")
    if not isinstance(evidence, Mapping):
        raise ValueError("Record-14 a=1 evidence is malformed")
    rows = evidence.get("rows")
    labels = evidence.get("unrestricted_compact_label_audits")
    if not isinstance(rows, list) or len(rows) != 5 or not isinstance(labels, list):
        raise ValueError("Record-14 a=1 coverage is malformed")
    branch = next(row for row in positive["branches"] if row["a"] == 1)
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
        raise ValueError("Record-14 a=1 source coverage changed")
    if [
        row["source_positive_A4_orbit_index"]
        for row in rows
        if row["generator_lex_slice_pure_target_count"]
    ] != list(TARGET_SOURCE_ORBITS):
        raise ValueError("Record-14 a=1 target-bearing sources changed")
    if sum(row["generator_lex_slice_pure_target_count"] for row in rows) != 19:
        raise ValueError("Record-14 a=1 target total changed")
    if any(
        not row["production_independent_literal_equality"]
        or row["production_z3_generator_lex_model_count"]
        != row["independent_cp_sat_generator_lex_model_count"]
        or row["production_z3_generator_lex_model_count"]
        != row["generator_lex_slice_pure_target_count"]
        for row in rows
    ):
        raise ValueError("Record-14 a=1 independent equality changed")
    if len(labels) != 10 or not all(
        row["unrestricted_compact_label_no_go"] for row in labels
    ):
        raise ValueError("Record-14 a=1 label closure changed")
    theorem = document.get("theorem")
    if (
        not isinstance(theorem, Mapping)
        or theorem.get("remaining_compact_schedules") != [0]
        or theorem.get("whole_shadow_compact_coxeter_no_go") is not False
    ):
        raise ValueError("Record-14 a=1 theorem scope changed")
    stored_digest = document.get("certificate_digest")
    body = dict(document)
    body.pop("certificate_digest", None)
    digest = _digest(body)
    if stored_digest != digest:
        raise ValueError("Record-14 a=1 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_DIGEST not in ("TO_BE_FROZEN", digest):
        raise ValueError("Record-14 a=1 certificate is not frozen")
    return {
        "record_index": 14,
        "schedule_a": 1,
        "source_orbit_count": 5,
        "generator_lex_target_count": 19,
        "full_target_orbit_count": 10,
        "terminal_compact_label_state_count": 0,
        "closed_schedules": [1, 2, 3, 4, 5, 6],
        "remaining_schedules": [0],
        "whole_shadow_compact_coxeter_no_go": False,
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
