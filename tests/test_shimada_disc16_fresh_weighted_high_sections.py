from __future__ import annotations

from collections import Counter
import json

from scripts.screen_shimada_disc16_fresh_weighted_high_sections import (
    CASE_DIRECTORY,
    EXPECTED_CASE_COUNT,
    INVENTORY,
    MANIFEST,
    _inventory_cases,
    _sha256,
    audit_archive,
)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_fresh_weighted_denominator_is_exactly_one_h13_plus_574_h18() -> None:
    inventory = _load(INVENTORY)
    cases = _inventory_cases(inventory)
    assert len(cases) == EXPECTED_CASE_COUNT == 575
    assert Counter(record["hyperbolic_dimension"] for _, record in cases) == {
        13: 1,
        18: 574,
    }
    assert all(
        chamber["class_source"] == "fresh_full_weighted_wall_group_enumeration"
        and chamber["groups_equal"] is False
        and chamber["archived_stabilizer_order"] == 1
        and chamber["weighted_wall_group_order"] == 2
        and record["order"] == 2
        for chamber, record in cases
    )


def test_manifest_closes_every_fresh_class_by_a_forbidden_facet_pair() -> None:
    manifest = _load(MANIFEST)
    summary = manifest["summary"]
    assert summary["inventory_case_count"] == 575
    assert summary["completed_case_count"] == 575
    assert summary["pending_case_count"] == 0
    assert summary["classification_histogram"] == {"not_coxeter": 575}
    assert summary["coxeter_candidate_count_requiring_volume_audit"] == 0
    assert "0" not in summary["bad_pair_count_histogram"]
    assert sum(summary["bad_pair_count_histogram"].values()) == 575
    assert summary["bad_pair_count_histogram"]["1"] == 29
    assert manifest["survivors_requiring_volume_audit"] == []
    assert all(
        row["status"] == "complete"
        and row["classification"] == "not_coxeter"
        and row["non_coxeter_pair_count"] >= 1
        for row in manifest["cases"]
    )


def test_h13_lattice_orientation_and_exact_five_pair_obstruction() -> None:
    case = _load(CASE_DIRECTORY / "W0600_C001.json")
    assert case["case"]["hyperbolic_dimension"] == 13
    assert case["case"]["published_wall_count"] == 20
    assert case["case"]["published_wall_orientation_sign_histogram"] == {
        "-1": 3,
        "1": 17,
    }
    lattice = case["saturated_fixed_divisor_lattice"]
    assert lattice["inertia"] == [1, 13, 0]
    assert lattice["determinant"] == -7168
    assert lattice["rational_constraint_clearing_denominator"] == 7
    mirror = case["restricted_mirror_geometry"]
    assert mirror["wall_orbit_count"] == 14
    assert (
        mirror[
            "saturated_fixed_divisor_lattice_is_stable_under_all_retained_mirror_reflections"
        ]
        is False
    )
    assert case["exact_reduction"]["hyperbolic_facet_count"] == 14
    assert case["exact_reduction"]["unresolved_conic_facet_count"] == 0
    assert case["coxeter_gate"]["non_coxeter_ratio_histogram"] == {
        "1/117": 1,
        "1/45": 1,
        "2/5": 1,
        "4/13": 1,
        "9/65": 1,
    }


def test_all_case_no_go_grams_and_selected_full_arithmetic_replay() -> None:
    inventory = _load(INVENTORY)
    manifest = audit_archive(
        _sha256(INVENTORY), _inventory_cases(inventory), CASE_DIRECTORY, MANIFEST
    )
    assert manifest["summary"]["classification_histogram"] == {
        "not_coxeter": 575
    }
