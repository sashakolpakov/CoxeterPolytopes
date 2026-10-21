from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_jacobian_kummer_h12_small_subtree_census import (
    OUTPUT,
    build,
)


@pytest.fixture(scope="module")
def payload() -> dict:
    rebuilt = build()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_complete_one_edge_census(payload: dict) -> None:
    census = payload["complete_one_edge_census"]
    assert census["raw_selected_wall_count"] == 112
    assert census["exact_diagram_type_count"] == 2
    assert [row["seed_wall_orbit_size"] for row in census["type_representatives"]] == [
        16,
        96,
    ]
    assert [row["facet_count"] for row in census["type_representatives"]] == [
        317,
        326,
    ]
    assert len(census["all_wall_records"]) == 112


def test_complete_three_vertex_path_census(payload: dict) -> None:
    census = payload["complete_three_vertex_coloured_path_census"]
    assert census["raw_unordered_generator_pair_count"] == 6_216
    assert census["exact_seed_action_orbit_count"] == 35
    assert census["diagram_type_count"] == 35
    assert census["all_35_pairwise_nonisomorphic"]
    assert census["diagram_type_facet_histogram"] == {
        "462": 4,
        "471": 7,
        "480": 24,
    }
    assert census["raw_pair_facet_histogram"] == {
        "462": 120,
        "471": 1_536,
        "480": 4_560,
    }
    hashes = {
        row["exact_vertex_edge_profile_multiset_sha256"]
        for row in census["orbit_representatives"]
    }
    assert len(hashes) == 35
    distance_forgotten_hashes = {
        row["distance_forgotten_vertex_relation_profile_multiset_sha256"]
        for row in census["orbit_representatives"]
    }
    ordinary_matrix_hashes = {
        row["ordinary_Coxeter_matrix_vertex_profile_multiset_sha256"]
        for row in census["orbit_representatives"]
    }
    assert len(distance_forgotten_hashes) == 35
    assert len(ordinary_matrix_hashes) == 35


def test_allcock_and_picard_gates(payload: dict) -> None:
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())
    assert payload["selected_stock"]["wall_orbit_sizes"] == [16, 96]
    assert payload["allcock_and_direct_certification"][
        "finite_subtree_consequence"
    ].startswith("each two- or three-vertex union")
    for census_name, record_name in (
        ("complete_one_edge_census", "all_wall_records"),
        ("complete_three_vertex_coloured_path_census", "orbit_representatives"),
    ):
        for row in payload[census_name][record_name]:
            assert row["full_linear_rank"] == 13
            assert row["fixed_boundary_core_Smith_diagonal"] == [1] * 13
            assert row["primitive_integral_reflection_root_count"] == row["facet_count"]


def test_bounded_lineage_claim_is_scoped(payload: dict) -> None:
    lineage = payload["bounded_lineage"]
    assert lineage["headers_462_471_480_absent"]
    assert lineage["evidence_grade"] == "B1: bounded same-header absence only"
    assert lineage["worldwide_priority"] == "open"
    assert payload["carrier_scope"]["weak_divisor_reflection_carrier_only"]
    assert not payload["carrier_scope"]["full_nef_or_effectivity_claim"]
