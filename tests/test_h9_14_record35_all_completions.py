from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.h9_14_record35 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
)
from coxeter_polytopes.combinatorics.h9_14_record35_all_completions import (
    EXPECTED_INDEPENDENT_ACTIVE_MASKS,
    EXPECTED_SCHEDULES,
    audit_record35_all_completions,
)
from scripts.certify_h9_14_record35_all_completions import (
    EXPECTED_ARCHIVE_FILE_SHA256,
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    SUPPORT,
    _canonical_digest,
    _normalized_certifier_source,
    verify,
)
from scripts.independent_h9_14_record35_all_completions_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay,
)


ROOT = Path(__file__).resolve().parents[1]
HOSTILE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record35_all_completions_audit.py"
)


@pytest.fixture(scope="session")
def production_audit():
    source = json.loads(SUPPORT.read_text(encoding="utf-8"))
    survivors = [
        row for row in source["records"] if row["survives_second_layer"]
    ]
    record = survivors[TARGET_RECORD_INDEX]
    assert record["graph6"] == TARGET_GRAPH6
    schedules_json = json.dumps(
        record["dehn_sommerville_schedule"],
        sort_keys=True,
        separators=(",", ":"),
    )
    support_core_json = json.dumps(
        record["ridge_support_core_masks"],
        sort_keys=True,
        separators=(",", ":"),
    )
    return audit_record35_all_completions(schedules_json, support_core_json)


@pytest.fixture(scope="session")
def hostile_replay():
    return replay()


def test_certificate_digest_is_json_roundtrip_invariant() -> None:
    document = {
        "histogram": {2: 7, 10: 3},
        "nested": [{5: 1, 12: 4}],
    }
    document["certificate_sha256"] = _canonical_digest(document)
    parsed = json.loads(json.dumps(document))
    assert _canonical_digest(document) == _canonical_digest(parsed)
    assert parsed["certificate_sha256"] == _canonical_digest(parsed)


def test_certifier_self_hash_normalizes_only_two_quoted_rhs_literals() -> None:
    source = (
        'EXPECTED_CERTIFICATE_SHA256 = "first"\n'
        'preserved = "EXPECTED_CERTIFICATE_SHA256 = trailing"\n'
        'EXPECTED_ARCHIVE_FILE_SHA256 = "second"\n'
        'def preserved_function():\n'
        '    return "all other bytes"\n'
    )
    assert _normalized_certifier_source(source) == (
        'EXPECTED_CERTIFICATE_SHA256 = "<OUTER_FREEZE>"\n'
        'preserved = "EXPECTED_CERTIFICATE_SHA256 = trailing"\n'
        'EXPECTED_ARCHIVE_FILE_SHA256 = "<OUTER_FREEZE>"\n'
        'def preserved_function():\n'
        '    return "all other bytes"\n'
    )


def test_lossless_all_active_a0_enumeration(production_audit) -> None:
    assert production_audit.verified
    assert production_audit.exact_residual_incidence_orbit_count == 4
    assert production_audit.schedule.rows == EXPECTED_SCHEDULES
    assert production_audit.schedule.surviving_a == (0, 1, 2)
    assert production_audit.schedule.excluded_a == (3, 4)

    a0 = production_audit.a0
    assert a0.independent_active_masks == EXPECTED_INDEPENDENT_ACTIVE_MASKS
    assert tuple(mask.bit_count() for mask in a0.independent_active_masks).count(3) == 4
    assert a0.independent_active_mask_count == 18
    assert a0.active_row_configuration_count == 1_221
    assert a0.column_distribution_problem_count == 235
    assert a0.total_distribution_path_count == 375_582
    assert a0.repeated_row_path_count == 7_239
    assert a0.small_transversal_path_count == 0
    assert a0.blocker_valid_path_count == 368_343
    assert sum(count for _facets, count in a0.facet_count_histogram) == 368_343
    assert dict(a0.facet_count_histogram)[102] == 20
    assert a0.target_path_count == 20
    assert a0.exact_target_family_count == 18
    assert a0.topology_valid_family_indices == (0, 1, 2, 3, 4)
    assert a0.all_target_full_graph_orbit_count == 11
    assert a0.topology_valid_full_graph_orbit_count == 3
    assert all(
        row.minimal_transversal_histogram == ((5, 102),)
        for row in a0.family_audits
    )
    assert len(a0.label_audits) == 18
    assert all(row.compact_no_go for row in a0.label_audits)
    assert all(row.stages[-1].passing_state_count == 0 for row in a0.label_audits)
    statement = a0.to_record()["lossless_enumeration_scope"][
        "active_row_statement"
    ]
    assert "all independent active subsets" in statement
    assert "No singleton-active-row normalization" in statement


def test_positive_schedules_are_exhausted(production_audit) -> None:
    positive = production_audit.positive_branches
    assert positive.verified
    assert positive.individually_interval_feasible_a4_count == 490
    assert positive.a1_viable_a4_count == 112
    assert positive.a1_orbit_size_histogram == ((56, 2),)
    first, second = positive.a1_completion_censuses
    assert (
        first.forced_shadow_count,
        first.additional_five_count,
        first.individually_cap_feasible_candidate_count,
        first.raw_addition_family_count,
        first.vertex_cap_feasible_count,
        first.no_small_transversal_count,
        first.scheduled_target_count,
    ) == (9, 2, 252, 31_626, 5_626, 1, 1)
    assert (
        second.forced_shadow_count,
        second.additional_five_count,
        second.individually_cap_feasible_candidate_count,
        second.raw_addition_family_count,
        second.vertex_cap_feasible_count,
        second.no_small_transversal_count,
        second.scheduled_target_count,
    ) == (8, 3, 462, 16_328_620, 1_794_605, 840, 0)
    assert positive.a2_unordered_pair_count == 119_805
    assert positive.a2_family_interval_feasible_count == 111_699
    assert positive.a2_viable_pair_count == 28
    assert positive.a2_orbit_size_histogram == ((28, 1),)
    assert positive.a2_completion_census.raw_addition_family_count == 1_540
    assert positive.a2_completion_census.vertex_cap_feasible_count == 540
    assert positive.a2_completion_census.no_small_transversal_count == 0


def test_a1_unrestricted_label_gate_has_exact_index_two_witnesses(
    production_audit,
) -> None:
    incidence = dict(production_audit.a1_incidence_audit)
    assert incidence["minimal_transversal_size_histogram"] == {"5": 100}
    assert incidence["ridge_degree_histogram"] == {"2": 450}
    assert incidence["passes_connected_homology_audit"] is True

    gate = production_audit.a1_label
    assert gate.unrestricted_compact_no_go
    assert gate.rank_nine_vertex_count == 100
    assert (
        gate.order_four_placement_count,
        gate.overlapping_order_five_placement_count,
        gate.disjoint_order_five_placement_count,
    ) == (108, 420, 420)
    assert gate.compatible_overlapping_join_count == 1_542
    assert gate.left_cluster_passing_state_count == 84
    assert gate.block_restriction_size_histogram == ((7, 20), (8, 55), (9, 15))
    assert gate.sparse_block_state_count == 35_280
    assert gate.sparse_block_failure_count == 0
    assert gate.label_three_test_count == 1_058_400
    assert gate.label_three_passing_count == 5_040
    assert gate.passing_edges_per_sparse_state_histogram == ((1, 5_040),)
    assert gate.stronger_single_edge_survivors == ((4, 0), (5, 0))
    assert gate.exceptional_class_histogram == (
        ("A5_A5_type_1_1", 4_320),
        ("A5_A5_type_4_1", 720),
    )
    assert tuple(row.inertia for row in gate.exact_index_two_witnesses) == (
        (8, 2, 0),
        (8, 2, 0),
    )
    assert tuple(
        row.multiplicity_per_support_gate
        for row in gate.exact_index_two_witnesses
    ) == (4_320, 720)


def test_no_project_import_hostile_replay_closes_the_compact_theorem(
    hostile_replay,
) -> None:
    source = HOSTILE_SOURCE.read_text(encoding="utf-8")
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    assert EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN"
    assert hostile_replay["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert hostile_replay["strict_schedule_survivors"] == [0, 1, 2]
    assert hostile_replay["positive_branches"]["a1_orbit_sizes"] == [56, 56]
    assert hostile_replay["positive_branches"]["a2_orbit_sizes"] == [28]
    a0 = hostile_replay["a0_census"]
    assert a0["complete_active_row_scope"] == (
        "all independent active subsets allowed by exact incidence caps; "
        "no singleton-row normalization"
    )
    assert a0["total_distribution_path_count"] == 375_582
    assert a0["exact_target_family_count"] == 18
    assert a0["independent_active_masks"] == list(EXPECTED_INDEPENDENT_ACTIVE_MASKS)
    assert a0["independent_active_mask_size_histogram"] == {
        0: 1,
        1: 5,
        2: 8,
        3: 4,
    }
    assert a0["topology_valid_family_indices"] == [0, 1, 2, 3, 4]
    assert a0["topology_valid_orbit_count"] == 3
    assert len(a0["full_topology_audits"]) == 5
    assert all(
        row["topology_valid"]
        and row["all_links_have_spherical_euler_characteristic"]
        and row["all_links_are_mod2_homology_spheres"]
        for row in a0["full_topology_audits"]
    )
    assert len(hostile_replay["a0_all_target_label_audits"]) == 18
    assert all(
        row["compact_no_go"]
        for row in hostile_replay["a0_all_target_label_audits"]
    )
    gate = hostile_replay["a1_label_audit"]
    assert gate["placement_counts"] == [108, 420, 420]
    assert gate["label_three_test_count"] == 1_058_400
    assert gate["label_three_passing_count"] == 5_040
    assert [row["inertia"] for row in gate["exact_witnesses"]] == [
        [8, 2, 0],
        [8, 2, 0],
    ]
    assert gate["compact_no_go"] is True
    assert hostile_replay["a1_incidence_audit"]["full_topology_audit"][
        "topology_valid"
    ] is True
    assert hostile_replay["theorem"] == {
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "exact_topology_valid_residual_incidence_orbit_count": 4,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def test_frozen_whole_shadow_certificate_rebuilds_compact_only() -> None:
    assert EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
    assert EXPECTED_ARCHIVE_FILE_SHA256 != "TO_BE_FROZEN"
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert sha256(OUTPUT.read_bytes()).hexdigest() == EXPECTED_ARCHIVE_FILE_SHA256
    summary = verify(document)
    assert summary == {
        "record_index": 35,
        "exact_residual_incidence_orbits": 4,
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert document["scope"]["compact_only"] is True
    assert document["scope"][
        "topology_audit_is_necessary_not_a_PL_sphere_certificate"
    ] is True
    assert document["scope"]["finite_volume_noncompact_chamber_excluded"] is False
    assert document["scope"]["ideal_or_parabolic_faces_enumerated"] is False
