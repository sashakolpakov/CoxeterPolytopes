from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_compact_246_property4_counterexample import (
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/compact_246_raw_maximal_property4_counterexample.md"
EXPECTED_CERTIFICATE_SHA256 = (
    "6b5d46c2eb09e07cacde56e9580fb77e23ff326ed5d19de29039be9642b87563"
)


def test_certificate_passes_and_is_frozen() -> None:
    payload = build_certificate()
    assert payload["all_checks_pass"]
    assert payload["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256


def test_archived_payload_is_exact_replay() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived == build_certificate()


def test_integral_arithmetic_triangle() -> None:
    payload = build_certificate()
    triangle = payload["coxeter_triangle"]
    assert triangle["wall_gram"] == [[-12, 6, 0], [6, -6, 3], [0, 3, -2]]
    assert triangle["signature"] == [1, 2]
    assert triangle["determinant"] == 36
    assert triangle["edge_labels"] == [4, 6, 2]
    assert triangle["principal_pair_determinants"] == [36, 3, 24]
    assert triangle["area_over_pi"] == "1/12"
    assert payload["arithmeticity"]["minimal_vinberg_field"] == "Q"


def test_all_reflections_are_integral_isometries() -> None:
    payload = build_certificate()
    for key in (
        "reflections_integral",
        "reflections_preserve_gram",
        "reflections_have_order_two",
    ):
        assert payload["checks"][key]


def test_rational_form_is_anisotropic_by_mod_three_descent() -> None:
    anisotropy = build_certificate()["arithmeticity"]["rational_anisotropy"]
    assert anisotropy["prime"] == 3
    assert anisotropy["residual_binary_zeros_mod_3"] == [[0, 0]]
    assert "Infinite descent" in anisotropy["argument"]


def test_absolute_maximality_has_only_one_area_candidate() -> None:
    maximality = build_certificate()["absolute_reflection_maximality"]
    assert maximality["proper_overgroup_area_upper_bound_over_pi"] == "1/24"
    assert [item["labels"] for item in maximality["all_triangle_types_in_area_window"]] == [
        [2, 3, 7],
        [2, 3, 8],
    ]
    assert maximality["all_triangle_types_in_area_window"][0]["source_area_ratio"] == "7/2"
    assert maximality["integer_index_triangle_candidates"] == [
        {"labels": [2, 3, 8], "area_over_pi": "1/24", "index": 2}
    ]


def test_index_two_candidate_does_not_double_to_target() -> None:
    doublings = build_certificate()["absolute_reflection_maximality"][
        "index_two_doubling_analysis"
    ]
    assert {item["union"] for item in doublings} == {
        "triangle_(3,3,4)",
        "non_Coxeter_triangle_with_angle_2pi/3",
        "quadrilateral",
    }
    assert all(item["union"] != "triangle_(2,4,6)" for item in doublings)


def test_rational_and_minimal_trace_obstructions_are_distinct() -> None:
    obstruction = build_certificate()["carrier_obstruction"]
    assert "anisotropic over Q" in obstruction["rational_host_gate"]
    assert "isotropic odd unimodular form I_(1,2)" in obstruction["rational_host_gate"]
    assert "nef isotropic ruling" in obstruction["rational_branches"]
    assert "rank-three trace block" in obstruction["minimalizable_full_trace_branch"]
    assert "labels 2 or 3" in obstruction["angle_contradiction"]
    assert "all empty" in obstruction["conclusion"]


def test_report_states_the_scope_and_escape_clause() -> None:
    text = " ".join(REPORT.read_text(encoding="utf-8").lower().split())
    required = (
        "raw maximal arithmetic counterexample",
        "delta(2,4,6)",
        "absolute reflection maximality",
        "proper invariant/eigen-nef slice",
        "direct rational",
        "certified weak-log rational",
        "minimalizable canonical-orthogonal full trace",
        "presently advertised exact-return outer three-branch compiler cannot be exhaustive",
        "bare carrier-relative",
        "exact-return minimalizable full-trace",
    )
    for token in required:
        assert token in text


def test_report_does_not_claim_abstract_polytope_nonexistence() -> None:
    text = REPORT.read_text(encoding="utf-8")
    forbidden = (
        "Delta(2,4,6) does not exist",
        "ABSW is false",
        "all algebraic carriers are impossible without qualification",
    )
    for token in forbidden:
        assert token not in text
