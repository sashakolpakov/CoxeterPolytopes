from __future__ import annotations

import json

import pytest

from scripts.certify_h9_14_record0_a1_source1_joint_branches import (
    EXPECTED_BRANCH_INDICES,
    EXPECTED_EXPANDED_CUT_COUNT,
    EXPECTED_MODEL_COUNTS,
    EXPECTED_PARTITION_DIGEST,
    EXPECTED_PATTERN_DIGEST,
    EXPECTED_TRANSCRIPT,
    EXPECTED_TRANSCRIPT_DIGEST,
    OUTPUT,
    ROOT,
    aggregate_coverage,
    audit_checkpoint,
    independent_model_counts,
    partition_record,
    placement_orbits,
    replay_seed_cuts,
    stabilizer_record,
    verify_record,
)


def _load_certificate() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_independent_108_placement_partition_and_exact_stabilizers() -> None:
    partition = partition_record()
    assert partition["literal_placement_count"] == 108
    assert partition["orbit_count"] == 19
    assert partition["orbit_sizes"] == [6, 6, 6, 3, *([6] * 14), 3]
    assert sum(partition["orbit_sizes"]) == 108
    assert partition["partition_digest"] == EXPECTED_PARTITION_DIGEST
    assert partition["coverage_is_disjoint_and_complete"] is True
    assert len({code for orbit in placement_orbits() for code in orbit}) == 108

    for branch, orbit in enumerate(placement_orbits()):
        audit = stabilizer_record(branch)
        internal_order = audit["exact_internal_code_stabilizer_order"]
        assert len(orbit) * internal_order == 6
        assert audit["exact_labelled_A4_stabilizer_order"] == 40320 * internal_order
        assert audit["generators_generate_exact_stabilizer"] is True
        assert audit[
            "all_generators_preserve_fixed_source_A4_code_and_facet_ground"
        ] is True
        assert audit["lex_generator_step_counts"] == (
            [332] * 7 + ([402] if len(orbit) == 3 else [])
        )


def test_independent_model_shape_and_full_seed_cut_replay() -> None:
    assert independent_model_counts() == EXPECTED_MODEL_COUNTS
    replay = replay_seed_cuts()
    assert replay["legacy_concrete_seed_cut_count"] == 73
    assert replay["abstract_nonfinite_pattern_count"] == 23
    assert replay["abstract_nonfinite_pattern_order_histogram"] == {
        "3": 5,
        "4": 9,
        "5": 8,
        "6": 1,
    }
    assert replay["all_abstract_patterns_deletion_minimal_nonfinite"] is True
    assert replay["abstract_nonfinite_pattern_digest"] == EXPECTED_PATTERN_DIGEST
    assert replay["expansion_transcript"] == EXPECTED_TRANSCRIPT
    assert replay["expansion_transcript_digest"] == EXPECTED_TRANSCRIPT_DIGEST
    assert replay["expanded_global_monotone_face_cut_count"] == (
        EXPECTED_EXPANDED_CUT_COUNT
    )
    assert replay["aggregate_expanded_cut_stream_digest"] == (
        "a4c01774988940813a748b744e585fc22086fdbadc8a5796e90078f0d76262eb"
    )


def test_frozen_aggregate_replays_all_terminal_branch_metadata() -> None:
    document = _load_certificate()
    summary = verify_record(document)
    assert summary["closed_branch_count"] == 19
    assert summary["placement_orbit_size_sum"] == 108
    assert summary["compact_label_survivor_count"] == 0
    assert summary["fixed_source_compact_coxeter_no_go"] is True
    assert [row["branch_index"] for row in document["branch_evidence"]] == list(
        EXPECTED_BRANCH_INDICES
    )
    assert all(
        row["terminal_solver_status"] == "INFEASIBLE"
        and row["compact_label_survivor"] is None
        and row["branch_compact_coxeter_no_go"] is True
        and row["checkpoint_fixed_source_compact_coxeter_no_go"] is False
        for row in document["branch_evidence"]
    )


def test_fixed_source_closure_requires_every_branch() -> None:
    document = _load_certificate()
    partial = aggregate_coverage(document["branch_evidence"][:-1])
    assert partial["individual_branch_closure_count"] == 18
    assert partial["missing_branch_indices"] == [18]
    assert partial["present_placement_orbit_size_sum"] == 105
    assert partial["fixed_source_compact_coxeter_no_go"] is False


def test_branch0_legacy_missing_terminal_statistics_are_tolerated(
    tmp_path,
) -> None:
    source = (
        ROOT
        / "data/enumerations"
        / "h9_14_record0_a1_source1_branch00_joint_checkpoint.json"
    )
    document = json.loads(source.read_text(encoding="utf-8"))
    for key in (
        "terminal_solver_wall_time_seconds",
        "terminal_solver_branch_count",
        "terminal_solver_conflict_count",
    ):
        document.pop(key, None)
    fixture = tmp_path / source.name
    fixture.write_text(json.dumps(document), encoding="utf-8")
    audit = audit_checkpoint(fixture, replay_seed_cuts())
    assert audit["terminal_solver_status"] == "INFEASIBLE"
    assert audit["terminal_statistics_available"] is False
    assert audit["terminal_statistics"] is None


def test_missing_terminal_statistics_are_not_tolerated_on_other_branches(
    tmp_path,
) -> None:
    source = (
        ROOT
        / "data/enumerations"
        / "h9_14_record0_a1_source1_branch01_joint_checkpoint.json"
    )
    document = json.loads(source.read_text(encoding="utf-8"))
    for key in (
        "terminal_solver_wall_time_seconds",
        "terminal_solver_branch_count",
        "terminal_solver_conflict_count",
    ):
        document.pop(key, None)
    fixture = tmp_path / source.name
    fixture.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(ValueError, match="terminal statistics are missing"):
        audit_checkpoint(fixture, replay_seed_cuts())


def test_certificate_does_not_claim_other_sources_shadow_or_noncompact() -> None:
    theorem = _load_certificate()["theorem"]
    assert theorem["fixed_source_compact_coxeter_no_go"] is True
    assert theorem["whole_record0_a1_schedule_excluded"] is False
    assert theorem["whole_record0_shadow_excluded"] is False
    assert theorem["finite_volume_noncompact_checked"] is False
    assert theorem["finite_volume_noncompact_chamber_excluded"] is False
