from __future__ import annotations

import json

from scripts.certify_h9_14_record31_a3_no_go import (
    ARCHIVE,
    EXPECTED_RIDGE_FAMILIES,
    _normalized_independent_targets,
    _normalized_production_targets,
    verify,
)


def test_record31_entire_compact_a3_schedule_is_excluded() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 31,
        "newly_closed_schedule": 3,
        "a3_orbit_count": 143,
        "literal_pure_target_count": 373,
        "target_orbit_count": 33,
        "ridge_pseudomanifold_orbit_count": 3,
        "unrestricted_compact_label_survivor_count": 0,
        "remaining_open_schedules": [0, 1, 2],
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
    assert document["production"][
        "pure_116_facet_target_count_across_representatives"
    ] == 373
    assert document["independent"][
        "scheduled_order_five_state_count_across_representatives"
    ] == 378
    assert document["independent"][
        "pure_116_facet_target_count_across_representatives"
    ] == 373
    assert document["downstream"]["target_orbit_count"] == 33
    assert document["downstream"]["compact_label_survivor_orbit_count"] == 0


def test_record31_a3_independent_literal_and_label_replays_match() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert _normalized_production_targets(
        document["production"]
    ) == _normalized_independent_targets(document["independent"])

    production_ridge = {
        tuple(row["minimal_nonface_masks"]): row["representative_source_index"]
        for row in document["downstream"]["target_orbits"]
        if row["ridge_pseudomanifold"]
    }
    independent_ridge = {
        tuple(row["minimal_nonface_masks"]): row["source_representative_index"]
        for row in document["independent_downstream"][
            "ridge_survivor_label_audits"
        ]
    }
    assert production_ridge == independent_ridge == EXPECTED_RIDGE_FAMILIES
    assert all(
        row["global_topology_audit"]["passes_connected_homology_audit"]
        and row["sparse_compact_label_gate"][
            "unrestricted_compact_label_no_go"
        ]
        for row in document["downstream"]["target_orbits"]
        if row["ridge_pseudomanifold"]
    )
    assert all(
        row["unrestricted_compact_label_no_go"]
        and row["terminal_sparse_state_count"] == 0
        for row in document["independent_downstream"][
            "ridge_survivor_label_audits"
        ]
    )
