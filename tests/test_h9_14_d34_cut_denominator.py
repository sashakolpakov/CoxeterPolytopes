import json
from pathlib import Path
import shutil

import pytest

import scripts.certify_h9_14_t1_t2_d34_unfolding_census as old_d34
import scripts.certify_h9_14_t2_surface_unfolding_pilot as old_pilot
from scripts.certify_h9_14_d34_cut_denominator import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def result() -> dict:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        pytest.fail("cddexec_gmp is required for the exact D34-cut replay")
    originals = (
        old_d34._corpus,
        old_pilot._load_factory_input,
        old_pilot._load_comparison_inputs,
    )

    def forbidden_archive_loader(*_args, **_kwargs):
        raise AssertionError("the variable compiler attempted to load a pinned archive")

    old_d34._corpus = forbidden_archive_loader
    old_pilot._load_factory_input = forbidden_archive_loader
    old_pilot._load_comparison_inputs = forbidden_archive_loader
    try:
        return build_certificate(executable)
    finally:
        (
            old_d34._corpus,
            old_pilot._load_factory_input,
            old_pilot._load_comparison_inputs,
        ) = originals


def test_archive_is_exact_replay(result: dict) -> None:
    assert json.loads(Path(OUTPUT).read_text()) == result
    assert result["input_artifacts"] == []
    assert result["all_checks_pass"] is True
    assert all(result["checks"].values())


def test_full_surface_free_denominator_counts(result: dict) -> None:
    assert result["global_counts"] == {
        "raw_labelled_m4_assignments": 37_180,
        "raw_diagram_orbits": 3_610,
        "norm_2_4_compatible_labelled_cuts": 4_789,
        "norm_2_4_compatible_diagram_orbits": 374,
        "scaled_dual_integral_Gram_lanes": 748,
        "rank10_Lorentzian_corank4_diagrams": 2,
        "stable_even_fixed_overlattice_candidates": 28,
        "minimal_involutive_weak_K3_packages": 2,
        "finite_volume_noncompact_polytopes": 2,
    }
    supports = {row["support_id"]: row for row in result["support_censuses"]}
    assert supports["S16"]["all_eight_m4_labellings"]["automorphism_orbit_count"] == 376
    assert supports["S17"]["all_eight_m4_labellings"]["automorphism_orbit_count"] == 3_234
    assert supports["S16"]["primitive_norm_2_4_compatible_labellings"]["automorphism_orbit_count"] == 136
    assert supports["S17"]["primitive_norm_2_4_compatible_labellings"]["automorphism_orbit_count"] == 238


def test_only_t2_and_t1_survive_and_all_lifts_are_recompiled(result: dict) -> None:
    supports = {row["support_id"]: row for row in result["support_censuses"]}
    s16 = supports["S16"]["passing_diagram"]
    s17 = supports["S17"]["passing_diagram"]
    assert s16["relation_histogram"] == {"2": 75, "3": 8, "4": 8}
    assert s17["relation_histogram"] == {"2": 74, "3": 9, "4": 8}
    assert s16["m4_edges_one_based"] == [
        [4, 5], [8, 9], [1, 10], [1, 11],
        [5, 12], [9, 13], [2, 14], [6, 14],
    ]
    assert s17["m4_edges_one_based"] == [
        [3, 4], [4, 5], [6, 7], [7, 8],
        [1, 12], [11, 12], [2, 14], [9, 14],
    ]
    assert [
        lane["stable_even_overlattice_candidate_count"]
        for lane in s16["scaled_dual_lanes"]
    ] == [16, 5]
    assert [
        lane["stable_even_overlattice_candidate_count"]
        for lane in s17["scaled_dual_lanes"]
    ] == [5, 2]
    assert s16["unique_minimal_involutive_weak_K3_package"]["rank18_unfolding"]["two_elementary_invariants"] == [18, 4, 0]
    assert s17["unique_minimal_involutive_weak_K3_package"]["rank18_unfolding"]["two_elementary_invariants"] == [18, 2, 0]


def test_fresh_polytope_and_literal_vinberg_certificates(result: dict) -> None:
    supports = {row["support_id"]: row for row in result["support_censuses"]}
    expected = {
        "S16": (110, 90, 20, 683, {"0,2": 97, "1,1": 362, "2,0": 224}),
        "S17": (127, 112, 15, 707, {"0,2": 53, "1,1": 300, "2,0": 354}),
    }
    for support_id, values in expected.items():
        cone = supports[support_id]["passing_diagram"]["fresh_cone_and_volume_certificate"]
        rays, finite, ideal, edges, histogram = values
        assert cone["extreme_ray_count"] == rays
        assert cone["finite_ray_count"] == finite
        assert cone["ideal_ray_count"] == ideal
        assert cone["all_14_facets_irredundant"] is True
        assert len(cone["facet_witnesses"]) == 14
        assert cone["strict_interior"]["square"] > 0
        assert cone["strict_interior"]["minimum_slack"] > 0
        vinberg = cone["literal_vinberg"]
        assert vinberg["decision"] == "finite_volume_noncompact"
        assert vinberg["elliptic_rank8_edge_count"] == edges
        assert vinberg["edge_endpoint_type_histogram"] == histogram
        assert vinberg["bad_edge_count"] == 0
        assert vinberg["bad_parabolic_neighborhood_count"] == 0


def test_urcv_scope_is_not_overstated(result: dict) -> None:
    boundary = result["URCV_boundary"]
    assert boundary["U_arithmetic"].startswith("closed")
    assert boundary["U_weak_surface"].startswith("closed")
    assert boundary["U_strong"].startswith("not proved")
    assert "negative" in boundary["C"]
    assert boundary["V"].startswith("closed afresh")
