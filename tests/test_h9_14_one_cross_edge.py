import json

from scripts.certify_h9_14_one_cross_edge import OUTPUT, verify


EXPECTED = {
    30: (2_659_392, 1_218, 575_424, 180),
    42: (9_057_484_800, 603_840, 78_520_320, 6_120),
    43: (1_701_000, 1_440, 19_440, 16),
    44: (8_164_800, 2_220, 95_040, 20),
    47: (1_020_600, 360, 16_200, 6),
}


def test_frozen_one_cross_edge_certificate_replays_exactly() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary == {
        "record_indices": [30, 42, 43, 44, 47],
        "extension_orbits": 609_078,
        "local_survivor_orbits": 6_342,
        "metric_survivor_orbits": 0,
        "certificate_sha256": (
            "b6b440c8649725bbc49a43662036dc8765423b252e7c8092d85f8f7436e3e4f7"
        ),
    }
    assert manifest["totals"] == {
        "full_labelled_trials": 9_071_030_592,
        "full_nerve_symmetry_orbits": 609_078,
        "locally_admissible_full_labelled_extensions": 79_226_424,
        "locally_admissible_full_nerve_symmetry_orbits": 6_342,
        "metric_survivor_orbits": 0,
    }

    for record in manifest["records"]:
        shell = record["one_edge_shell"]
        observed = (
            shell["full_labelled_trial_count"],
            shell["full_nerve_symmetry_orbit_count"],
            shell["locally_admissible_full_labelled_count"],
            shell["locally_admissible_full_nerve_symmetry_orbit_count"],
        )
        assert observed == EXPECTED[record["record_index"]]
        assert shell["all_unaffected_vertices_reused_from_frozen_exact_audit"]
        assert shell["all_affected_vertices_checked"]
        assert shell["all_higher_lanner_nonfaces_unchanged_and_checked"]
        assert len(record["survivor_orbit_representatives"]) == observed[-1]
        assert record["metric_principal_gate"]["not_excluded_orbit_count"] == 0
        assert all(
            survivor["metric_principal_test"]["status"]
            == "negative_index_at_least_two"
            for survivor in record["survivor_orbit_representatives"]
        )


def test_nontrivial_multiwedge_schur_witnesses_are_frozen() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    expected_expression = "(50 - 22*sqrt(5))/(16*(2 - sqrt(5))**2)"
    schur = {
        record["record_index"]: [
            survivor
            for survivor in record["survivor_orbit_representatives"]
            if survivor["metric_principal_test"]["method"]
            == "exact_two_lanner_schur_factor"
        ]
        for record in manifest["records"]
    }
    assert {index: len(rows) for index, rows in schur.items()} == {
        30: 0,
        42: 0,
        43: 1,
        44: 1,
        47: 0,
    }
    for index in (43, 44):
        metric = schur[index][0]["metric_principal_test"]
        assert metric["schur_factor_sign"] == 1
        assert metric["schur_factor"] == expected_expression


def test_scope_is_explicitly_one_edge_only() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert "not an unrestricted label census" in manifest["theorem"][
        "scope_warning"
    ]
    assert manifest["theorem"]["not_a_new_polytope"] is True
    assert all(
        record["scope"]["changed_cross_edge_count"] == 1
        and record["scope"]["unrestricted_label_census"] is False
        for record in manifest["records"]
    )
