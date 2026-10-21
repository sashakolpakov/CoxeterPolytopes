from __future__ import annotations

import json

from scripts.audit_coble_four_point_constructible_incidence_hostile import (
    OUTPUT,
    build_certificate,
)


def test_symbolic_replay_strengthens_the_core_geometry() -> None:
    algebra = build_certificate()["independent_four_point_algebra"]
    assert algebra["grid"]["rank_histogram_at_char3_lambda2"] == {"3": 24}
    assert algebra["grid"]["distinct_projective_kernel_sections"] == 24
    assert algebra["grid"]["all_coefficient_matrices_invertible"]
    assert algebra["polar_chord"]["rank_histogram_at_char3_lambda2"] == {"2": 6}
    assert algebra["polar_chord"]["distinct_residual_lines"] == 6
    assert algebra["polar_chord"]["residual_boundary_collisions"] == 0


def test_raw_dependencies_are_reconstructed_without_the_bridge() -> None:
    raw = build_certificate()["independent_raw_dependency_extraction"]
    mii = raw["rank18_48wall"]
    vi = raw["rank15_60wall"]
    assert mii["bridge_archive_read"] is False
    assert mii["frame_census"] == {
        "short_mirrors": 24,
        "all_eight_subsets": 735_471,
        "pairwise_orthogonal": 21_211,
        "odd_edges": 32,
        "odd_edge_vertex_covers": 1,
        "integral_frames": 1,
        "unique_frame_zero_based": list(range(8)),
    }
    assert mii["reconstructed_NS"]["inertia"] == [1, 17, 0]
    assert mii["dependency_extraction"]["all_S4_permutations_once"]
    assert vi["bridge_archive_read"] is False
    assert vi["reconstructed_NS"]["inertia"] == [1, 14, 0]
    assert vi["dependency_extraction"]["plane_branches"] == 5
    assert vi["dependency_extraction"]["every_branch_has_complete_pair_complement_design"]


def test_two_explicit_cores_really_are_strong_coble_root_seeds() -> None:
    certificate = build_certificate()
    assert certificate["strong_Coble"]["verdict_for_both_explicit_cores"] == "PASS"
    assert certificate["D_CobVin_exact_return"]["verdict"] == "CONDITIONAL_PASS"
    assert certificate["D_CobVin_exact_return"]["not_merely_an_abstract_isomorphism"]


def test_standard_free_extension_implication_has_an_exact_counterexample() -> None:
    audit = build_certificate()["standard_free_extension_audit"]
    counterexample = audit["counterexample"]
    assert audit["verdict"] == "FALSE_AS_STATED"
    assert counterexample["base_point_free"]
    assert counterexample["generically_separable"]
    assert counterexample["smooth"]
    assert counterexample["geometrically_irreducible"] is False
    assert counterexample["reason"] == "the system is composed with a pencil"


def test_strict_five_gate_rating_is_not_L3_plus_or_L4() -> None:
    certificate = build_certificate()
    gates = certificate["five_gate_rating"]
    assert gates["C_independent_P_only_class"]["verdict"] == "PARTIAL_PASS"
    assert gates["Car_nonempty_for_every_member"]["verdict"].startswith("FAIL_")
    assert gates["Return_exact_same_P"]["verdict"] == "CONDITIONAL_PASS"
    assert gates["AG_independent_classification"]["verdict"] == "CONDITIONAL_PASS"
    assert gates["Pushback_new_theorem_on_every_member"]["verdict"].startswith("FAIL_")
    assert certificate["level_recommendation"]["recommended_rating"] == (
        "L3_TWO_MEMBER_INVERSE_CALIBRATION; NOT_AN_EXTENSION_CLASS_THEOREM; NOT_L4"
    )
    assert certificate["level_recommendation"]["corrected_production_rating_accepted"] == (
        "L3_TWO_MEMBER_CORE_CALIBRATION;_NOT_L3_PLUS_OR_L4"
    )


def test_enlarged_mobius_core_is_independently_reconstructed_but_separate() -> None:
    enlarged = build_certificate()["enlarged_finite_geometry_class"]
    mi = enlarged["third_raw_member"]
    assert mi["surface_labels_or_stored_blocks_read"] is False
    assert mi["physical_taxonomy"] == {
        "boundary_minus4": 2,
        "actual_minus2": 30,
        "crossing_minus1": 10,
        "single_boundary_cusp_lifts": 164,
    }
    assert mi["raw_design"]["identification"] == "Steiner system S(3,4,10)"
    assert mi["raw_design"]["triple_multiplicity_histogram"] == {"1": 120}
    assert mi["independent_F9_realization"]["PGL2_F9_order"] == 720
    assert mi["independent_F9_realization"][
        "unique_projectivity_graphs_interpolating_Frobenius"
    ] == 30
    assert mi["general_mobius_factory"]["exact_Coble_root_locus"] == "q=3"
    assert enlarged["Möbius_factory_theorem"]["q_not_3_nonexistence_claim"] is False
    scope = enlarged["production_artifact_scope"]
    assert scope["production_core_matches_independent_replay"]
    assert scope["arbitrary_input_option"] is False
    assert scope["actual_frozen_rating"] == "L3_STRUCTURAL_P_FIRST_CALIBRATION"


def test_frozen_hostile_archive_matches_live_replay() -> None:
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == build_certificate()
