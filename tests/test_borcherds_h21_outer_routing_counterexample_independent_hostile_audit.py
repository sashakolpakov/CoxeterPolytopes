from __future__ import annotations

from pathlib import Path

from scripts.audit_borcherds_h21_outer_routing_counterexample_independent_hostile import (
    audit,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = (
    ROOT
    / "results"
    / "borcherds_h21_outer_routing_counterexample_independent_hostile_audit.md"
)


def test_independent_audit_passes() -> None:
    payload = audit()
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())


def test_special_absolute_maximality_gates() -> None:
    maximality = audit()["maximality_audit"]
    assert maximality["short_coordinate_smith"] == [1] * 22
    assert maximality["short_odd_graph_component_size"] == 42
    assert maximality["discriminant_action_group_order"] == 6
    assert maximality["stable_discriminant_subgroup_orders"] == [1, 4]
    assert "trace 20" in maximality["rationality_gap_closed_by_trace"]
    assert "absolute reflection maximality" in maximality["conclusion"]


def test_fano_contraction_obstruction() -> None:
    obstruction = audit()["physical_and_contraction_audit"]
    assert obstruction["common_square_classes"] == [1, 2]
    assert obstruction["point_count"] == obstruction["line_count"] == 21
    assert obstruction["degree"] == 5
    assert obstruction["fano_determinant"] == -2
    assert obstruction["contraction"]["picard_rank_after"] == 1
    assert obstruction["contraction"]["surface_after"] == "P2"
    assert obstruction["contraction"]["image_square"] == 1
    assert obstruction["contraction"]["image_degree"] == 1


def test_scope_preserves_characteristic_two_escape() -> None:
    scope = audit()["scope"]
    assert "absolute reflection maximality" in scope["certified"]
    assert "characteristic-two" in scope["preserved_escape"]
    assert "not" not in scope["preserved_escape"].lower()
    assert "nonexistence of every" in scope["not_claimed"]


def test_report_needles() -> None:
    text = " ".join(REPORT.read_text().split())
    for needle in (
        "Hostile verdict: PASS",
        "projective-rationality step",
        "trace 20",
        "M=L\\text{ or }L^\\vee",
        "determinant is \\(-2\\)",
        "characteristic-two",
        "did not edit any owner byte",
    ):
        assert needle in text
