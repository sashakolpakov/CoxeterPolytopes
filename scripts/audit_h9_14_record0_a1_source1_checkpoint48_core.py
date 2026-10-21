#!/usr/bin/env python3
"""Independent local audit of the 48th source-1 compact CEGAR core.

This is deliberately not a terminal no-go verifier.  It checks only that the
newly appended iteration 47 is a lossless compact-label contradiction and
that its complete fixed-source orbit has the archived size.  The finite and
Lanner classifiers are imported from the independently frozen source-3
auditor, not from the exploratory source-1 CEGAR or core extractor.
"""

from __future__ import annotations

import argparse
from collections import deque
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
from typing import Any, Iterable, Sequence

from scripts.audit_h9_14_record0_a1_source3_terminal_cegar import (
    _audit_lanner_tables,
    _checkpoint_stage_projection,
    _digest,
    _label_state_count,
    _mask,
    _move,
)


ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_label_core_cegar_checkpoint_48_more.json"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_checkpoint48_core_audit.json"
)
INDEPENDENT_CLASSIFIER = (
    ROOT / "scripts/audit_h9_14_record0_a1_source3_terminal_cegar.py"
)

CHECKPOINT_SHA256 = (
    "1d1f2ae760364730db8fc3aaaed78b2a3514237182174697186270d85c332931"
)
INDEPENDENT_CLASSIFIER_SHA256 = (
    "6540c3c39afb6de6a710f50348c2ff0bc2d4a42cf2bda72d4ce96c37aea219cf"
)
SCHEMA = "coxeter.h9_14.record0_a1_source1.checkpoint48_core_audit.v1"

N = 14
FULL = (1 << N) - 1
FIXED_FOUR_VERTICES = (0, 2, 3, 4)
FIXED_FOUR = _mask(FIXED_FOUR_VERTICES)
MISSING_EDGES = (_mask((0, 13)), _mask((1, 13)))
EXPECTED_SOURCE_GROUP_ORDER = factorial(3) * factorial(8)
EXPECTED_ORBIT_SIZE = 168
EXPECTED_ORBIT_SHA256 = (
    "81cd3f7b7c749f7ded551b2565ec6987122054d74411288346bc07afcd85e3a9"
)
EXPECTED_REPRESENTATIVE_DIGEST = (
    "889d7c9ac49418f055d0bd856100b3957b7da1f326d7d5f0475508632023961b"
)

CoreSignature = tuple[tuple[int, ...], tuple[int, ...]]


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value & (1 << vertex))


def _swap(left: int, right: int) -> tuple[int, ...]:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def _source_generators() -> tuple[tuple[int, ...], ...]:
    inside = (2, 3, 4)
    outside = tuple(range(5, 13))
    return tuple(
        _swap(block[index], block[index + 1])
        for block in (inside, outside)
        for index in range(len(block) - 1)
    )


def _move_signature(
    signature: CoreSignature, action: Sequence[int]
) -> CoreSignature:
    additions, restrictions = signature
    return (
        tuple(sorted(_move(value, action) for value in additions)),
        tuple(sorted(_move(value, action) for value in restrictions)),
    )


def _core_orbit(
    signature: CoreSignature, generators: Sequence[Sequence[int]]
) -> set[CoreSignature]:
    orbit = {signature}
    queue = deque([signature])
    while queue:
        current = queue.popleft()
        for generator in generators:
            image = _move_signature(current, generator)
            if image not in orbit:
                orbit.add(image)
                queue.append(image)
    return orbit


def _independent_five(value: int) -> bool:
    return value.bit_count() == 5 and all(value & edge != edge for edge in MISSING_EDGES)


def _project(values: Iterable[int], support_union: int) -> tuple[int, ...]:
    return tuple(sorted({value & support_union for value in values}))


def build_certificate() -> dict[str, Any]:
    if _sha256(CHECKPOINT) != CHECKPOINT_SHA256:
        raise ArithmeticError("source-1 checkpoint hash changed")
    if _sha256(INDEPENDENT_CLASSIFIER) != INDEPENDENT_CLASSIFIER_SHA256:
        raise ArithmeticError("independent finite/Lanner classifier hash changed")

    checkpoint = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    row = checkpoint["iterations"][-1]
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
    signature = (additions, restrictions)
    orbit = _core_orbit(signature, generators)
    sorted_orbit = sorted(orbit)

    local_checks = {
        "checkpoint_identity_is_record0_a1_source1": (
            checkpoint["record_index"] == 0
            and checkpoint["schedule_a"] == 1
            and checkpoint["source_orbit_index"] == 1
            and checkpoint["representative_A4s"] == [list(FIXED_FOUR_VERTICES)]
        ),
        "checkpoint_is_bounded_not_terminal": (
            checkpoint["terminal_solver_status"] == "OPTIMAL"
            and checkpoint["stopped_reason"] == "core_cutoff"
            and checkpoint["fixed_source_compact_coxeter_no_go"] is False
            and checkpoint["theorem_certificate"] is False
        ),
        "exactly_48_core_iterations_are_present": (
            len(checkpoint["iterations"]) == 48
            and checkpoint["excluded_local_core_orbit_count"] == 48
            and row["iteration"] == 47
        ),
        "new_representative_digest_matches": (
            row["representative_digest"] == EXPECTED_REPRESENTATIVE_DIGEST
        ),
        "fixed_A4_occurs_once": supports.count(FIXED_FOUR) == 1,
        "remaining_supports_are_the_three_stored_A5s": (
            tuple(sorted(value for value in supports if value != FIXED_FOUR))
            == additions
            and len(additions) == 3
            and all(_independent_five(value) for value in additions)
        ),
        "four_actual_vertices_have_rank_nine": (
            len(actual_vertices) == 4
            and all(value.bit_count() == 9 for value in actual_vertices)
        ),
        "facet_complements_match_actual_vertices": (
            facet_complements
            == tuple(sorted(FULL ^ value for value in actual_vertices))
        ),
        "face_cut_has_the_same_projected_semantics": (
            _project(actual_vertices, support_union)
            == _project(restrictions, support_union)
        ),
        "independent_actual_vertex_replay_matches": (
            actual_projection == expected_stages and actual_terminal == 0
        ),
        "independent_face_cut_replay_matches": (
            cut_projection == expected_stages and cut_terminal == 0
        ),
        "stored_terminal_state_is_zero": (
            core["terminal_sparse_label_state_count"] == 0
        ),
    }

    generator_checks = [
        _move(FIXED_FOUR, generator) == FIXED_FOUR
        and {_move(edge, generator) for edge in MISSING_EDGES}
        == set(MISSING_EDGES)
        for generator in generators
    ]
    orbit_checks = {
        "nine_adjacent_swaps_generate_S3_times_S8": (
            len(generators) == 9
            and all(generator_checks)
            and checkpoint["source_group_order"] == EXPECTED_SOURCE_GROUP_ORDER
        ),
        "complete_generator_BFS_orbit_has_size_168": (
            len(orbit) == row["full_core_orbit_size"]
            == row["full_face_core_orbit_size"]
            == EXPECTED_ORBIT_SIZE
        ),
        "orbit_hash_matches": _digest(sorted_orbit) == EXPECTED_ORBIT_SHA256,
        "stabilizer_order_is_1440": (
            EXPECTED_SOURCE_GROUP_ORDER // len(orbit)
            == row["core_stabilizer_order"]
            == 1440
        ),
        "every_orbit_addition_remains_an_independent_A5": all(
            all(_independent_five(value) for value in image_additions)
            for image_additions, _ in orbit
        ),
        "every_orbit_clause_has_a_literal": all(
            bool(image_additions or image_restrictions)
            for image_additions, image_restrictions in orbit
        ),
    }

    lanner = _audit_lanner_tables()
    checks = local_checks | orbit_checks | {
        "independent_classical_Lanner_tables_pass": all(
            lanner["checks"].values()
        ),
        "finite_volume_noncompact_is_not_claimed": (
            checkpoint["finite_volume_noncompact_checked"] is False
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"checkpoint-48 local audit failed: {checks}")

    return {
        "schema": SCHEMA,
        "scope": (
            "Independent local proof for iteration 47 and its complete source-1 "
            "orbit only. This is not a terminal fixed-source compact no-go and "
            "makes no finite-volume noncompact claim."
        ),
        "inputs": {
            "checkpoint": str(CHECKPOINT.relative_to(ROOT)),
            "checkpoint_sha256": CHECKPOINT_SHA256,
            "independent_classifier": str(
                INDEPENDENT_CLASSIFIER.relative_to(ROOT)
            ),
            "independent_classifier_sha256": INDEPENDENT_CLASSIFIER_SHA256,
        },
        "bounded_checkpoint_status": {
            "terminal_solver_status": checkpoint["terminal_solver_status"],
            "stopped_reason": checkpoint["stopped_reason"],
            "excluded_local_core_orbits": checkpoint[
                "excluded_local_core_orbit_count"
            ],
            "cumulative_orbit_image_clause_count": sum(
                int(item["full_core_orbit_size"])
                for item in checkpoint["iterations"]
            ),
            "fixed_source_compact_coxeter_no_go": False,
            "theorem_certificate": False,
        },
        "new_core": {
            "iteration": 47,
            "representative_digest": row["representative_digest"],
            "selected_facet_complements": row[
                "selected_facet_complement_count"
            ],
            "actual_vertices": row["rank_nine_actual_vertex_count"],
            "local_support_profile": "fixed A4 + 3 A5 supports",
            "local_actual_vertex_count": len(actual_vertices),
            "independent_passing_state_counts": [
                stage["passing_state_count"] for stage in cut_stages
            ],
            "independent_compatible_extension_counts": [
                stage["compatible_extension_count"] for stage in cut_stages
            ],
            "full_source_orbit_size": len(orbit),
            "core_stabilizer_order": EXPECTED_SOURCE_GROUP_ORDER // len(orbit),
            "full_source_orbit_sha256": _digest(sorted_orbit),
            "solver_branches": row["solver_branch_count"],
            "solver_conflicts": row["solver_conflict_count"],
            "solver_wall_time_seconds": row["solver_wall_time_seconds"],
        },
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
            raise SystemExit("checkpoint-48 local audit archive mismatch")
        print("archive_check=PASS")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
