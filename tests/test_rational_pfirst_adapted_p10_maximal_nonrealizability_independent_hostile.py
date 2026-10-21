"""Independent hostile tests for the adapted decorated-P10 no-go."""

from __future__ import annotations

from scripts.audit_rational_pfirst_adapted_p10_maximal_nonrealizability_independent_hostile import (
    OWNER_PINS,
    ROOT,
    ag_terminal_certificate,
    build_audit,
    cone_certificate,
    descent_regression,
    physical_and_scope_certificate,
    root_facet_certificate,
    sha256_file,
)


def test_owner_bytes_are_locked() -> None:
    for item in OWNER_PINS.values():
        assert sha256_file(ROOT / item["path"]) == item["sha256"]


def test_all_sixty_six_active_sets_recover_exact_cone() -> None:
    cone = cone_certificate()
    assert cone["active_ten_subsets_tested"] == 66
    assert (cone["finite_rays"], cone["ideal_rays"], cone["hyperideal_rays"]) == (
        17,
        2,
        0,
    )
    assert cone["strict_interior_square"] == 3940
    assert cone["all_facets_rank_ten"] is True
    assert cone["finite_volume_by_closed_ball_polytope_triangulation"] is True
    assert cone["wall_cone_contains_no_negative_vector"] is True


def test_coble_mukai_root_chamber_itself_has_finite_volume() -> None:
    root = root_facet_certificate()
    assert root["dual_squares"] == [4, 18, 42, 30, 20, 12, 6, 2, 0, 10]
    assert (root["finite_vertices"], root["ideal_vertices"], root["hyperideal_vertices"]) == (
        9,
        1,
        0,
    )
    assert root["finite_covolume_root_reflection_group"] is True


def test_unbounded_descent_has_no_bounded_search_gap() -> None:
    descent = descent_regression(18)
    assert descent["large_d_Cauchy_bound"] == "a^2<=9q<=18"
    assert descent["large_d_exhaustive_survivor"] == "F=(3,-1^10)"
    assert descent["bounded_regression_roots"] > 10000
    assert descent["every_nonsimple_step_strictly_decreases_measure"] is True
    assert descent["bounded_search_used_as_proof"] is False


def test_physical_decorated_and_uncoloured_meanings_are_distinct() -> None:
    physical = physical_and_scope_certificate()
    assert physical["physical_square_split"] == [-1] + [-2] * 10 + [-4]
    assert physical["canonical_pairings"] == [-1] + [0] * 10 + [2]
    assert physical["uncoloured_F_arithmetic_genus"] == 1
    assert physical["decorated_B_arithmetic_genus"] == 0
    assert physical["same_mirror"] is True
    assert physical["different_geometric_effectivity_demand"] is True


def test_terminal_no_go_is_ag_and_complex_scoped() -> None:
    no_go = ag_terminal_certificate()
    assert no_go["base_field"] == "C"
    assert no_go["displayed_effective_root_count"] == 10
    assert no_go["effective_irreducible_roots_in_each_allowed_type"] == 12
    assert no_go["decorated_full_nef_carrier_exists"] is False
    assert no_go["raw_maximality_uses_this_classification"] is False
    assert no_go["terminal_no_go_uses_this_classification"] is True


def test_final_verdict_preserves_all_nonclaims() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_SCOPED"
    separation = audit["logical_separation"]
    assert separation["raw_finite_volume_and_maximality_classification_free"]
    assert separation["decorated_complex_rational_full_nef_nonrealizability"]
    assert separation["uncoloured_P10_nonexistence_claimed"] is False
    assert separation["uncoloured_rational_full_nef_carrier_exists"] is True
    assert separation["positive_characteristic_nonrealizability_claimed"] is False
    assert separation["universal_adapted_dimension_bound_claimed"] is False
