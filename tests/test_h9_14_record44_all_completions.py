from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.h9_14_record44 import (
    TARGET_MISSING_FIVES,
)
from coxeter_polytopes.combinatorics.h9_14_record44_all_completions import (
    EXPECTED_SCHEDULES,
    audit_record44_all_completions,
)
from scripts.certify_h9_14_record44_all_completions import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    SUPPORT,
    verify,
)
from scripts.independent_h9_14_record44_all_completions_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay,
)


ROOT = Path(__file__).resolve().parents[1]
HOSTILE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record44_all_completions_audit.py"
)
EXPECTED_ARCHIVE_FILE_SHA256 = (
    "ad01bf7d5431def08403385c824371ca6be19257a00c3f9a9623191a7d8db019"
)


@pytest.fixture(scope="session")
def production_audit():
    source = json.loads(SUPPORT.read_text(encoding="utf-8"))
    record = next(
        row for row in source["records"] if row["graph6"] == "M?????????O?o?o@?"
    )
    schedules_json = json.dumps(
        record["dehn_sommerville_schedule"],
        sort_keys=True,
        separators=(",", ":"),
    )
    return audit_record44_all_completions(schedules_json)


@pytest.fixture(scope="session")
def hostile_replay():
    return replay()


def test_all_four_schedules_reduce_to_one_unique_completion(
    production_audit,
) -> None:
    assert production_audit.verified
    assert production_audit.exact_residual_orbit_count == 1
    assert production_audit.schedule.rows == EXPECTED_SCHEDULES
    assert production_audit.schedule.surviving_missing_four_counts == (0, 1)
    assert production_audit.schedule.excluded_missing_four_counts == (2, 3)
    assert production_audit.a1.independent_a4_count == 690
    assert production_audit.a1.interval_feasible_a4_count == 378
    assert production_audit.a1.strict_shadow_feasible_count == 0
    assert production_audit.a0.enumeration_path_count == 2782
    assert production_audit.a0.exact_mask_state_count == 2214
    assert production_audit.a0.canonical_quotient_orbit_count == 502
    assert 2782 > 2214 > 502
    assert production_audit.a0.column_distribution_problem_count == 7
    assert sum(
        multiplicity * orbit_count
        for multiplicity, orbit_count in (
            production_audit.a0.path_multiplicity_per_quotient_orbit_histogram
        )
    ) == 2782
    assert production_audit.a0.small_transversal_path_count == 0
    assert production_audit.a0.path_facet_count_histogram[-1] == (90, 1)
    assert production_audit.a0.exact_mask_facet_count_histogram[-1] == (90, 1)
    assert production_audit.a0.quotient_facet_count_histogram[-1] == (90, 1)
    assert production_audit.a0.target_path_count == 1
    assert production_audit.a0.target_exact_mask_state_count == 1
    assert production_audit.a0.target_quotient_orbit_count == 1
    assert production_audit.a0.quotient_transcript_sha256 == (
        "d4d139e4ae208509eaf5d692a8d8cdd1f8c58a9b1e27c36d26696e53301977cb"
    )
    assert tuple(
        tuple(vertex for vertex in range(14) if value & (1 << vertex))
        for value in production_audit.a0.unique_target_missing_fives
    ) == TARGET_MISSING_FIVES


def test_no_import_hostile_replay_closes_the_same_compact_theorem(
    hostile_replay,
) -> None:
    source = HOSTILE_SOURCE.read_text(encoding="utf-8")
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    assert hostile_replay["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert hostile_replay["strict_schedule_survivors"] == [0, 1]
    assert hostile_replay["a1_census"]["aggregate_macaulay_passing_count"] == 0
    assert hostile_replay["a0_census"]["enumeration_path_count"] == 2782
    assert hostile_replay["a0_census"]["exact_mask_clutter_count"] == 2214
    assert hostile_replay["a0_census"]["s4_s9_canonical_state_count"] == 502
    assert hostile_replay["a0_census"][
        "active_graph_automorphism_quotient_applied"
    ] is False
    assert hostile_replay["a0_census"][
        "scheduled_target_s4_s9_canonical_count"
    ] == 1
    assert hostile_replay["unique_target_audit"][
        "minimal_transversal_size_histogram"
    ] == {"5": 90}
    assert hostile_replay["compact_coxeter_archive"][
        "unrestricted_compact_coxeter_no_go"
    ] is True
    assert hostile_replay["compact_coxeter_archive"][
        "normalizing_support_symmetry_orders"
    ] == [120]
    assert "principal submatrix" in hostile_replay["compact_coxeter_archive"][
        "principal_block_lemma"
    ]
    assert hostile_replay["scope"]["finite_volume_noncompact_checked"] is False


def test_frozen_whole_shadow_certificate_rebuilds_with_compact_only_scope() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert EXPECTED_ARCHIVE_FILE_SHA256 != "TO_BE_FROZEN"
    assert sha256(OUTPUT.read_bytes()).hexdigest() == EXPECTED_ARCHIVE_FILE_SHA256
    summary = verify(document)
    assert summary == {
        "record_index": 44,
        "exact_residual_incidence_orbits": 1,
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert document["unrestricted_compact_label_chain"][
        "record44_projection_sha256"
    ] == "99ba161f503317df8b96ce9adec033e4dd29102f576915e9b797a21c3bb498db"
    chain = document["unrestricted_compact_label_chain"]
    assert chain["normalizing_support_symmetry_orders"] == [120]
    assert chain["unrestricted_compact_coxeter_no_go"] is True
    assert "principal submatrix" in chain["principal_block_lemma"]
    assert "inertia" in chain["principal_block_lemma"]
    assert document["scope"]["compact_only"] is True
    assert document["scope"]["finite_volume_noncompact_chamber_excluded"] is False
