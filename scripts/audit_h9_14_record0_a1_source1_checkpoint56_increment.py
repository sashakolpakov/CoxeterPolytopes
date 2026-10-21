#!/usr/bin/env python3
"""Independently audit source-1 compact CEGAR cores 53--56.

This verifier does not import the exploratory CEGAR or its core extractor.
It replays the four appended local contradictions with the frozen independent
finite/Lanner classifier, rebuilds every complete ``S3 x S8`` source orbit by
generator BFS, and proves that all 401,100 installed clause signatures in the
56-core checkpoint are distinct.  The result is a bounded incremental
certificate, not a terminal compact no-go theorem.
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
CHECKPOINT_52 = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_label_core_cegar_checkpoint_52_more.json"
)
CHECKPOINT_53 = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_label_core_cegar_checkpoint_53_more.json"
)
CHECKPOINT_56 = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_label_core_cegar_checkpoint_56_more.json"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source1_checkpoint56_increment_audit.json"
)
BASE_AUDITOR = (
    ROOT / "scripts/audit_h9_14_record0_a1_source1_checkpoint48_core.py"
)
INDEPENDENT_CLASSIFIER = (
    ROOT / "scripts/audit_h9_14_record0_a1_source3_terminal_cegar.py"
)

CHECKPOINT_52_SHA256 = (
    "8577321be0061fca0aa9401069b8dad2fea7fabbb6cb583524234b8932a28f48"
)
CHECKPOINT_53_SHA256 = (
    "0567d3998fbabf0ae3a2e6bb3329ea6040de83474a177011a6815b149083807e"
)
CHECKPOINT_56_SHA256 = (
    "a611447e70b7964c3643c5d837917427dd92c63b66d52ba50650457b40b93a46"
)
BASE_AUDITOR_SHA256 = (
    "ad3d8cde9dc9bda1af69b3f748cb19b030890df43a2123a025396f4722d91ef5"
)
INDEPENDENT_CLASSIFIER_SHA256 = (
    "6540c3c39afb6de6a710f50348c2ff0bc2d4a42cf2bda72d4ce96c37aea219cf"
)
SCHEMA = "coxeter.h9_14.record0_a1_source1.checkpoint56_increment_audit.v1"

PREVIOUS_DISTINCT_CLAUSE_COUNT = 387_996
INCREMENTAL_DISTINCT_CLAUSE_COUNT = 13_104
EXPECTED_DISTINCT_CLAUSE_COUNT = 401_100
EXPECTED_FAMILY_CANONICAL_SHA256 = (
    "c0c0bb19fb3860ae4779c0b4459f1aa4031f471ca20535102d716296a25acd6a"
)
EXPECTED_DISTINCT_CLAUSES_SHA256 = (
    "09556bd47f128544e250eaa10761e94b048ac84715e825cbacbb5c8887407fbb"
)

EXPECTED = {
    52: {
        "digest": "2b6671d5510370e353e50aa19e229e75097d968a0796aeddb6dbcb0f2a6790cf",
        "support_count": 4,
        "actual_vertex_count": 6,
        "stored_global_passing": [108, 98, 16, 0],
        "independent_passing": [108, 98, 16, 0],
        "independent_compatible": [108, 1542, 76, 64],
        "orbit_size": 504,
        "stabilizer": 480,
        "orbit_sha256": "5b4133e716b26897816a2297c3978996db85397e520eeb0bde435d5492c99f61",
        "solver_branches": 793_961,
        "solver_conflicts": 96_087,
        "solver_wall_time_seconds": 146.091561,
    },
    53: {
        "digest": "10b118cba05ea78a4cfc67b1587950cc82863bea2e983dd8d09478379cb4e722",
        "support_count": 3,
        "actual_vertex_count": 6,
        "stored_global_passing": [108, 98, 16, 0],
        "independent_passing": [108, 110, 0],
        "independent_compatible": [108, 1542, 922],
        "orbit_size": 504,
        "stabilizer": 480,
        "orbit_sha256": "24b4b54878487513bdeeb43554989909f932c69f04cdd5c92b4ebb84f5d67de5",
        "solver_branches": 1_930_858,
        "solver_conflicts": 315_027,
        "solver_wall_time_seconds": 625.765564,
    },
    54: {
        "digest": "f50a6e1bd82af47b227fc4eb70186e0bae30ddc779084adaa41fb9ef0626a6da",
        "support_count": 4,
        "actual_vertex_count": 2,
        "stored_global_passing": [108, 84, 8, 0],
        "independent_passing": [108, 256, 110, 0],
        "independent_compatible": [108, 1542, 812, 395],
        "orbit_size": 2016,
        "stabilizer": 120,
        "orbit_sha256": "7adef54b708e265d991d69896b0f9a51be9d144a59a136dd63168695cb560149",
        "solver_branches": 3_842_878,
        "solver_conflicts": 724_416,
        "solver_wall_time_seconds": 1401.682469,
    },
    55: {
        "digest": "09a56ee03c235d984358f4f368d013e805fe15acb4cb4a5892d72fc33ca54257",
        "support_count": 5,
        "actual_vertex_count": 2,
        "stored_global_passing": [108, 84, 10, 2, 0],
        "independent_passing": [108, 350, 122, 21, 0],
        "independent_compatible": [108, 1542, 1113, 75, 64],
        "orbit_size": 10080,
        "stabilizer": 24,
        "orbit_sha256": "e7cc7a1de3dff93fdedb77dd54785bfd7ad35f4a79d6bb46959f40b78d8449eb",
        "solver_branches": 1_252_517,
        "solver_conflicts": 170_032,
        "solver_wall_time_seconds": 253.580417,
    },
}


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _core_signature(
    row: dict[str, Any],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    core = row["local_core"]
    return (
        tuple(sorted(_mask(values) for values in core["core_additional_A5s"])),
        tuple(sorted(_mask(values) for values in core["core_face_restrictions"])),
    )


def _eligible_addition(value: int) -> bool:
    """Whether an A5 belongs to the source-1 missing-five variable ground."""

    return _independent_five(value) and value & FIXED_FOUR != FIXED_FOUR


def _audit_row(row: dict[str, Any]) -> dict[str, Any]:
    iteration = int(row["iteration"])
    expected = EXPECTED[iteration]
    core = row["local_core"]
    supports = tuple(_mask(values) for values in core["core_supports"])
    additions, restrictions = _core_signature(row)
    actual_vertices = tuple(
        _mask(values) for values in core["core_actual_vertices"]
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
            and core["core_support_count"] == expected["support_count"]
            and all(_eligible_addition(value) for value in additions)
        ),
        "local_actual_vertices_have_rank_nine": (
            len(actual_vertices) == expected["actual_vertex_count"]
            and core["core_actual_vertex_count"]
            == expected["actual_vertex_count"]
            and all(value.bit_count() == 9 for value in actual_vertices)
        ),
        "complete_target_counts_remain_178": (
            row["selected_facet_complement_count"] == 178
            and row["rank_nine_actual_vertex_count"] == 178
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
            and core["terminal_sparse_label_state_count"] == 0
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
        "solver_metrics_match_frozen_checkpoint": (
            row["solver_branch_count"] == expected["solver_branches"]
            and row["solver_conflict_count"] == expected["solver_conflicts"]
            and row["solver_wall_time_seconds"]
            == expected["solver_wall_time_seconds"]
        ),
        "every_orbit_addition_remains_in_missing_five_ground": all(
            all(_eligible_addition(value) for value in image_additions)
            for image_additions, _ in orbit
        ),
        "every_orbit_clause_has_at_least_one_literal": all(
            bool(image_additions or image_restrictions)
            for image_additions, image_restrictions in orbit
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


def _audit_all_orbits(rows: list[dict[str, Any]]) -> dict[str, Any]:
    generators = _source_generators()
    distinct_clauses: set[
        tuple[tuple[int, ...], tuple[int, ...]]
    ] = set()
    canonical_representatives = []
    independently_rebuilt_sizes = []
    all_additions_remain_independent = True
    all_clauses_have_a_literal = True
    for row in rows:
        orbit = _core_orbit(_core_signature(row), generators)
        canonical_representatives.append(min(orbit))
        independently_rebuilt_sizes.append(len(orbit))
        distinct_clauses.update(orbit)
        all_additions_remain_independent &= all(
            all(_eligible_addition(value) for value in additions)
            for additions, _ in orbit
        )
        all_clauses_have_a_literal &= all(
            bool(additions or restrictions)
            for additions, restrictions in orbit
        )

    stored_sizes = [int(row["full_core_orbit_size"]) for row in rows]
    checks = {
        "all_56_stored_orbit_sizes_rebuild_independently": (
            independently_rebuilt_sizes == stored_sizes
        ),
        "all_56_orbit_families_are_pairwise_distinct": (
            len(set(canonical_representatives)) == len(rows) == 56
        ),
        "all_401100_clause_signatures_are_pairwise_distinct": (
            sum(independently_rebuilt_sizes)
            == len(distinct_clauses)
            == EXPECTED_DISTINCT_CLAUSE_COUNT
        ),
        "family_canonical_representative_hash_matches": (
            _digest(canonical_representatives)
            == EXPECTED_FAMILY_CANONICAL_SHA256
        ),
        "distinct_clause_signature_hash_matches": (
            _digest(sorted(distinct_clauses))
            == EXPECTED_DISTINCT_CLAUSES_SHA256
        ),
        "all_orbit_additions_remain_in_missing_five_ground": (
            all_additions_remain_independent
        ),
        "all_orbit_clauses_have_at_least_one_literal": (
            all_clauses_have_a_literal
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global 56-orbit census failed: {checks}")
    return {
        "orbit_family_count": len(rows),
        "independently_rebuilt_orbit_size_sum": sum(
            independently_rebuilt_sizes
        ),
        "distinct_clause_signature_count": len(distinct_clauses),
        "family_canonical_representatives_sha256": _digest(
            canonical_representatives
        ),
        "distinct_clause_signatures_sha256": _digest(
            sorted(distinct_clauses)
        ),
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    frozen_inputs = {
        CHECKPOINT_52: CHECKPOINT_52_SHA256,
        CHECKPOINT_53: CHECKPOINT_53_SHA256,
        CHECKPOINT_56: CHECKPOINT_56_SHA256,
        BASE_AUDITOR: BASE_AUDITOR_SHA256,
        INDEPENDENT_CLASSIFIER: INDEPENDENT_CLASSIFIER_SHA256,
    }
    for path, expected_hash in frozen_inputs.items():
        if _sha256(path) != expected_hash:
            raise ArithmeticError(f"frozen input hash changed: {path}")

    checkpoint_52 = json.loads(CHECKPOINT_52.read_text(encoding="utf-8"))
    checkpoint_53 = json.loads(CHECKPOINT_53.read_text(encoding="utf-8"))
    checkpoint_56 = json.loads(CHECKPOINT_56.read_text(encoding="utf-8"))
    rows = checkpoint_56["iterations"][-4:]
    if [int(row["iteration"]) for row in rows] != [52, 53, 54, 55]:
        raise ArithmeticError("checkpoint does not end in iterations 52--55")

    summaries = [_audit_row(row) for row in rows]
    global_orbits = _audit_all_orbits(checkpoint_56["iterations"])
    lanner = _audit_lanner_tables()
    previous_cumulative = sum(
        int(row["full_core_orbit_size"])
        for row in checkpoint_52["iterations"]
    )
    cumulative = sum(
        int(row["full_core_orbit_size"])
        for row in checkpoint_56["iterations"]
    )
    metadata_keys = (
        "record_index",
        "graph6",
        "schedule_a",
        "source_orbit_index",
        "representative_A4s",
        "source_group_order",
        "dual_model_counts",
    )
    checks = {
        "checkpoint_identity_is_record0_a1_source1": (
            checkpoint_56["record_index"] == 0
            and checkpoint_56["schedule_a"] == 1
            and checkpoint_56["source_orbit_index"] == 1
        ),
        "fixed_model_metadata_is_unchanged": all(
            checkpoint_52[key] == checkpoint_53[key] == checkpoint_56[key]
            for key in metadata_keys
        ),
        "first_52_iterations_are_structurally_unchanged": (
            checkpoint_53["iterations"][:52] == checkpoint_52["iterations"]
            and checkpoint_56["iterations"][:52]
            == checkpoint_52["iterations"]
        ),
        "checkpoint53_is_exact_prefix_of_checkpoint56": (
            checkpoint_56["iterations"][:53]
            == checkpoint_53["iterations"]
        ),
        "exactly_56_bounded_core_iterations_are_present": (
            len(checkpoint_56["iterations"]) == 56
            and checkpoint_56["excluded_local_core_orbit_count"] == 56
            and all(
                row["iteration"] == index and row["local_core"] is not None
                for index, row in enumerate(checkpoint_56["iterations"])
            )
        ),
        "checkpoint_remains_OPTIMAL_core_cutoff_not_a_theorem": (
            checkpoint_56["terminal_solver_status"] == "OPTIMAL"
            and checkpoint_56["stopped_reason"] == "core_cutoff"
            and checkpoint_56["fixed_source_compact_coxeter_no_go"] is False
            and checkpoint_56["theorem_certificate"] is False
        ),
        "four_incremental_core_audits_pass": all(
            all(summary["checks"].values()) for summary in summaries
        ),
        "independent_global_orbit_census_passes": all(
            global_orbits["checks"].values()
        ),
        "independent_classical_Lanner_tables_pass": all(
            lanner["checks"].values()
        ),
        "distinct_clause_counts_match": (
            previous_cumulative == PREVIOUS_DISTINCT_CLAUSE_COUNT
            and cumulative - previous_cumulative
            == INCREMENTAL_DISTINCT_CLAUSE_COUNT
            and cumulative == EXPECTED_DISTINCT_CLAUSE_COUNT
            and global_orbits["distinct_clause_signature_count"]
            == EXPECTED_DISTINCT_CLAUSE_COUNT
        ),
        "finite_volume_noncompact_is_not_claimed": (
            checkpoint_56["finite_volume_noncompact_checked"] is False
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"checkpoint-56 incremental audit failed: {checks}")

    return {
        "schema": SCHEMA,
        "scope": (
            "Independent proof for newly appended iterations 52--55, their "
            "complete source-1 orbits, and pairwise distinctness of all "
            "401100 installed clause signatures. The 56-core checkpoint "
            "remains bounded and is not a fixed-source compact no-go theorem."
        ),
        "inputs": {
            "checkpoint_52": str(CHECKPOINT_52.relative_to(ROOT)),
            "checkpoint_52_sha256": CHECKPOINT_52_SHA256,
            "checkpoint_53": str(CHECKPOINT_53.relative_to(ROOT)),
            "checkpoint_53_sha256": CHECKPOINT_53_SHA256,
            "checkpoint_56": str(CHECKPOINT_56.relative_to(ROOT)),
            "checkpoint_56_sha256": CHECKPOINT_56_SHA256,
            "base_independent_auditor": str(BASE_AUDITOR.relative_to(ROOT)),
            "base_independent_auditor_sha256": BASE_AUDITOR_SHA256,
            "independent_classifier": str(
                INDEPENDENT_CLASSIFIER.relative_to(ROOT)
            ),
            "independent_classifier_sha256": INDEPENDENT_CLASSIFIER_SHA256,
        },
        "bounded_checkpoint_status": {
            "terminal_solver_status": checkpoint_56[
                "terminal_solver_status"
            ],
            "stopped_reason": checkpoint_56["stopped_reason"],
            "excluded_local_core_orbits": checkpoint_56[
                "excluded_local_core_orbit_count"
            ],
            "previous_distinct_clause_count": previous_cumulative,
            "incremental_distinct_clause_count": (
                cumulative - previous_cumulative
            ),
            "distinct_clause_signature_count": cumulative,
            "cross_core_distinctness_independently_checked": True,
            "fixed_source_compact_coxeter_no_go": False,
            "theorem_certificate": False,
            "finite_volume_noncompact_checked": False,
        },
        "new_cores": summaries,
        "global_orbit_census": global_orbits,
        "independent_lanner_tables": lanner,
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        if archived != certificate:
            raise SystemExit("checkpoint-56 incremental archive mismatch")
        print("archive_check=PASS")
    elif args.write_archive:
        ARCHIVE.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"wrote={ARCHIVE.relative_to(ROOT)}")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
