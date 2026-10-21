from __future__ import annotations

from copy import deepcopy
import json

import pytest

from scripts.certify_h9_14_record31_a2_no_go import (
    ARCHIVE,
    EXPECTED_FAILURES,
    EXPECTED_RESIDUALS,
    EXPECTED_TARGET_DISTRIBUTION,
    canonical_family_key,
    verify,
)


def _document() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def _target_keys(census: dict[str, object]) -> set[tuple[int, ...]]:
    return {
        canonical_family_key(family)
        for row in census["orbits"]
        for family in row["ridge_valid_target_minimal_nonfaces"]
    }


def _target_distribution(census: dict[str, object]) -> dict[int, int]:
    return {
        row["orbit_index"]: row["ridge_valid_target_orbit_count"]
        for row in census["orbits"]
        if row["ridge_valid_target_orbit_count"]
    }


def test_record31_entire_compact_a2_schedule_is_excluded() -> None:
    document = _document()
    assert verify(document) == {
        "strict_face_link_excluded_orbit_count": 31,
        "enumerated_residual_orbit_count": 18,
        "ridge_target_orbit_count": 12,
        "compact_label_survivor_count": 0,
        "certificate_sha256": (
            "b8d8d6e7bd3ea4f16425a8f13c2e5449957717c2e78fe5172db734b38e6608b5"
        ),
    }
    assert document["scope"] == {
        "compact_only": True,
        "entire_a2_schedule_excluded": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    assert document["schedule"] == {
        "a": 2,
        "missing_five_count": 21,
        "facet_count": 118,
        "f_vector": [14, 87, 319, 769, 1275, 1463, 1126, 531, 118],
    }


def test_record31_a2_strict_face_link_partition_replays_exactly() -> None:
    document = _document()
    production = document["production"]
    independent_incidence = document["independent_incidence"]
    independent_link = document["independent_link_audit"]

    for census in (production, independent_incidence):
        failed = [
            row["orbit_index"]
            for row in census["orbits"]
            if row["fixed_local_failure"]
        ]
        assert failed == EXPECTED_FAILURES
        assert [
            index for index in range(49) if index not in failed
        ] == EXPECTED_RESIDUALS
        witnesses = [
            witness
            for row in census["orbits"]
            for witness in row["fixed_local_failure_witnesses"]
        ]
        assert len(witnesses) == 158
        assert {len(witness["subset"]) for witness in witnesses} == {2, 3}
        assert all(witness["g_prefix"][-1] == 0 for witness in witnesses)
        assert all(
            witness["additional_A5_incidence_interval"] == [0, -1]
            for witness in witnesses
        )

    assert (
        independent_link["strict_face_link_excluded_orbit_indices"]
        == EXPECTED_FAILURES
    )
    assert independent_link["residual_orbit_indices"] == EXPECTED_RESIDUALS
    assert [
        row["orbit_index"]
        for row in independent_link["orbits"]
        if row["strict_face_link_failure"]
    ] == EXPECTED_FAILURES


def test_record31_a2_independent_incidence_target_orbits_are_equal() -> None:
    document = _document()
    production = document["production"]
    independent = document["independent_incidence"]

    assert production["enumeration_mode"] == "nauty_partial_quotient"
    assert (
        independent["enumeration_mode"]
        == "direct_binary_matrix_partial_quotient"
    )
    assert _target_distribution(production) == EXPECTED_TARGET_DISTRIBUTION
    assert _target_distribution(independent) == EXPECTED_TARGET_DISTRIBUTION
    assert len(_target_keys(production)) == 12
    assert _target_keys(production) == _target_keys(independent)

    # Reversed tie/candidate traversal must not merely replay the production
    # transcript while claiming independence.
    production_nodes = [row["search_node_count"] for row in production["orbits"]]
    independent_nodes = [row["search_node_count"] for row in independent["orbits"]]
    assert production_nodes != independent_nodes
    assert sum(production_nodes) == 17_015_439
    assert sum(independent_nodes) == 13_253_521


def test_record31_a2_all_twelve_topology_targets_have_independent_label_no_go() -> None:
    document = _document()
    downstream = document["downstream"]
    independent = document["independent_downstream"]

    assert downstream["input_target_orbit_count"] == 12
    assert downstream["connected_mod2_homology_sphere_orbit_count"] == 12
    assert downstream["unrestricted_compact_label_no_go_orbit_count"] == 12
    assert downstream["compact_label_survivor_orbit_count"] == 0
    for row in downstream["target_orbits"]:
        assert row["minimal_transversal_size_histogram"] == {"5": 118}
        assert row["ridge_degree_histogram"] == {"2": 531}
        topology = row["global_topology_audit"]
        assert topology["matches_dehn_sommerville_schedule"] is True
        assert topology["dual_graph_connected"] is True
        assert topology["all_links_have_spherical_euler_characteristic"] is True
        assert topology["all_links_are_mod2_homology_spheres"] is True
        label = row["sparse_compact_label_gate"]
        assert label["rank_nine_vertex_count"] == 118
        assert label["state_cap_used"] is False
        assert label["terminal_sparse_state_count"] == 0
        assert label["unrestricted_compact_label_no_go"] is True

    assert independent["classical_lanner_placement_counts"] == {
        "4": 108,
        "5": 420,
    }
    assert independent["literal_target_orbit_count"] == 12
    assert independent["compact_label_survivor_count"] == 0
    assert all(
        row["rank_nine_vertex_count"] == 118
        and row["state_cap_used"] is False
        and row["terminal_sparse_state_count"] == 0
        and row["unrestricted_compact_label_no_go"] is True
        for row in independent["ridge_survivor_label_audits"]
    )


def test_record31_a2_scope_overpromotion_is_rejected() -> None:
    document = deepcopy(_document())
    document["scope"]["finite_volume_noncompact_chamber_excluded"] = True
    with pytest.raises(ValueError, match="scope changed or was over-promoted"):
        verify(document)
