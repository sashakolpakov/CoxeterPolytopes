from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_rational_pfirst_bianticanonical_section_algebra_gate import (
    ARCHIVE,
    UPSTREAM_PINS,
    build,
    dot,
    ordinary_jet_count,
    plane_monomial_count,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_pfirst_bianticanonical_section_algebra_gate.md"


def frozen() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text())


def test_live_certificate_matches_frozen_archive() -> None:
    assert build() == frozen()


def test_rank_eleven_anticanonical_arithmetic() -> None:
    row = frozen()["lattice"]
    assert row["K_square"] == row["F_square"] == -1
    assert row["B_square"] == -4
    assert (row["K_dot_F"], row["K_dot_B"]) == (1, 2)
    assert (row["F_arithmetic_genus"], row["B_arithmetic_genus"]) == (1, 0)
    assert row["B_minus_2K"] == [2 * value for value in row["F_minus_K"]]


def test_common_proper_point_and_contraction_frame() -> None:
    frame = frozen()["common_formal_frame"]
    assert frame["all_centers_proper_and_distinct"] is True
    assert frame["proximity_parent_lists"] == [[] for _ in range(10)]
    assert frame["proximity_matrix"] == [
        [int(i == j) for j in range(10)] for i in range(10)
    ]
    assert frame["labelled_exceptional_intersection_matrix"] == [
        [-int(i == j) for j in range(10)] for i in range(10)
    ]
    assert frame["labelled_exceptional_contraction_subsets"] == 1024
    assert frame["labelled_exceptional_complete_contraction_orders"] == 3628800


def test_formal_jet_descriptors_are_identical_but_answers_flip() -> None:
    payload = frozen()
    frame = payload["common_formal_frame"]
    assert frame["minus_K"] == {
        "plane_degree": 3,
        "multiplicities": [1] * 10,
        "matrix_shape": [10, 10],
    }
    assert frame["minus_2K"] == {
        "plane_degree": 6,
        "multiplicities": [2] * 10,
        "matrix_shape": [30, 28],
    }
    elliptic = payload["geometric_branches"]["smooth_cubic_ten_points"]
    coble = payload["geometric_branches"][
        "irreducible_ten_node_rational_sextic"
    ]
    assert (elliptic["minus_K_matrix_rank"], coble["minus_K_matrix_rank"]) == (
        9,
        10,
    )
    assert elliptic["minus_2K_matrix_rank"] == coble["minus_2K_matrix_rank"] == 27
    assert elliptic["physical_boundary_gate"] is False
    assert coble["physical_boundary_gate"] is True
    assert elliptic["degree_two_generator"] == "square_of_degree_one_generator"
    assert coble["degree_two_generator"] == "divisor_prime_generator"


def test_all_four_bezout_budgets_are_strict() -> None:
    for branch in frozen()["bezout_replay"].values():
        for name, row in branch.items():
            if name == "conclusion":
                continue
            assert row["forced_local_lower_bound"] > row["global_intersection"]


def test_winger_positive_control_is_byte_pinned_and_geometric() -> None:
    dependency = frozen()["explicit_ten_node_existence_dependency"]
    assert dependency["pins"] == UPSTREAM_PINS
    assert {
        relative: sha256_file(ROOT / relative) for relative in UPSTREAM_PINS
    } == UPSTREAM_PINS
    projected = dependency["projected_certificate"]
    assert projected["center_count"] == 10
    assert projected["all_ten_nodes"] is True
    assert projected["absolutely_irreducible"] is True
    assert projected["normalization_genus"] == 0
    assert projected["strong_coble"] is True
    assert (projected["cubic_rank"], projected["sextic_rank"]) == (10, 27)


def test_public_integer_boundaries_reject_coercion() -> None:
    assert dot((1,) + (0,) * 10, (1,) + (0,) * 10) == 1
    with pytest.raises(ValueError):
        dot((1,), (1,))
    for bad in (True, -1, 3.5):
        with pytest.raises(ValueError):
            plane_monomial_count(bad)  # type: ignore[arg-type]
        with pytest.raises(ValueError):
            ordinary_jet_count(bad)  # type: ignore[arg-type]


def test_report_keeps_the_no_go_and_exact_return_scope_separate() -> None:
    text = REPORT.read_text()
    assert "same-frame bi-anticanonical dichotomy" in text
    assert "enriched numerical/proximity no-go" in text
    assert "fat-point Rees slice" in text
    assert "classification-free Cox-prime criterion" in text
    assert "finite-sided finite-volume Coxeter" in text
    assert "no proper discrete reflection overgroup" in text
    assert "same full nef cone or the same raw Coxeter" in text
    assert "universal nonemptiness" in text
    scope = frozen()["scope"]
    assert scope["same_full_nef_cone_claimed"] is False
    assert scope["same_raw_coxeter_chamber_claimed"] is False
    assert scope["global_property4_claimed"] is False
    assert scope["universal_incidence_nonemptiness_claimed"] is False
    hypotheses = frozen()["criterion"]["exact_return_handoff_hypotheses"]
    assert "finite_sided_finite_volume_coxeter_chamber" in hypotheses
    assert "all_side_reflections_preserve_adapted_host" in hypotheses
    assert "no_proper_discrete_reflection_overgroup_in_Oplus_adapted_host" in hypotheses
