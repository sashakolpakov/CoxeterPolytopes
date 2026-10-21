import os
from pathlib import Path

import pytest
from sympy import Matrix

from scripts.analyze_shimada_disc11_enhanced_quotient import (
    _exact_matrix_order,
    analyze_enhanced_quotient,
)


DISC11_PATH = Path(
    os.environ.get(
        "SHIMADA_DISC11",
        "/tmp/shimada_threediscs/threediscs/disc11.txt",
    )
)


def test_exact_GL20_torsion_gate() -> None:
    assert _exact_matrix_order(Matrix.eye(20))[:2] == (True, 1)
    assert _exact_matrix_order(Matrix.diag(-1, *([1] * 19)))[:2] == (True, 2)
    assert _exact_matrix_order(Matrix([[2, 0], [0, 1]]))[0] is False
    assert _exact_matrix_order(Matrix([[1, 1], [0, 1]]))[0] is False


@pytest.mark.skipif(
    not DISC11_PATH.is_file(),
    reason="Shimada's 45 MB disc11.txt supplement is not vendored",
)
def test_enhanced_returns_and_bounded_C2_quotient() -> None:
    result = analyze_enhanced_quotient(DISC11_PATH, traversal_states=50)

    assert result["verified"] is True
    enhanced = result["enhanced_return"]
    assert enhanced["D68_gallery_representatives"] == [68, 10, 10, 69, 68]
    assert enhanced["D27_gallery_representatives"] == [27, 11, 11, 74, 27]
    assert enhanced["NS_isometry"] is True
    assert enhanced["involution"] is True
    assert enhanced["centralizes_tau"] is True
    assert enhanced["discriminant_action"] == "minus_identity"
    assert enhanced["D27_global_relative_equals_D68_return"] is True
    assert enhanced["D27_local_left_matrix_centralizes_tau"] is False
    assert enhanced["D68_and_D27_return_wall_classes_are_same_mirror"] is True
    assert enhanced["generated_Coxeter_subgroup"] == "W(A1)=C2"
    assert enhanced["generated_subgroup_order"] == 2

    quotient = result["bounded_state_quotient"]
    assert quotient["processed_states"] == 50
    assert quotient["queued_states"] == 26
    assert quotient["internal_crossings_from_processed_states"] == 208
    assert quotient["A_images_missing_at_hard_cutoff"] == 6
    assert quotient["A_closed_state_set"] == 82
    assert quotient["A_closed_orbits"] == 41
    assert quotient["processed_orbits"] == 27
    assert quotient["frontier_only_orbits"] == 14
    assert quotient["fixed_physical_states"] == 0
    assert quotient["same_representative_pair_comparisons"] == 32
    assert quotient["centralizing_same_representative_pairs"] == 23
    assert quotient["distinct_centralizing_relative_transports"] == 1
    assert quotient["only_centralizing_relative_is_D68_involution"] is True
    assert quotient["quotient_closed"] is False

    decision_ball = result["centralizer_decision_ball"]
    assert decision_ball["facet_star_complete"] is False
    assert decision_ball["complete_restricted_facet_radius"] is None
    assert decision_ball["one_witness_prefix_radius"] == 5
    assert decision_ball["one_witness_prefix_layer_counts"] == {
        "0": 1,
        "1": 3,
        "2": 5,
        "3": 7,
        "4": 10,
        "5": 14,
    }
    assert decision_ball["states_in_one_witness_prefix_layers"] == 40
    assert decision_ball["D68_type_states_in_one_witness_prefix"] == 2
    assert decision_ball["D68_representative_stabilizer_order"] == 2
    assert decision_ball["centralizing_automorphisms_in_sample"] == 2
    assert decision_ball["only_identity_and_A_in_one_witness_sample"] is True
    assert decision_ball["centralizer_in_complete_ball_is_exactly_C2"] is False
    assert sorted(
        item["exact_order"] for item in decision_ball["centralizing_records"]
    ) == [1, 2]
    assert decision_ball["new_centralizer_requires_distance_at_least"] is None
    assert decision_ball["one_witness_queue_closed"] is False
    assert decision_ball["traversal_closed"] is False
    assert decision_ball["infinite_order_return_found"] is False
    assert decision_ball["decision_status"] == (
        "pending_exact_ambient_carrier_enumeration"
    )

    pell = result["Pell_seed_factorization"]
    assert pell["verified"] is True
    assert [item["residual"] for item in pell["records"]] == ["E", "identity"]

    roots = result["root_wall_quotient"]
    assert roots["visible_root_rays"] == 11
    assert roots["fixed_root_rays"] == 9
    assert roots["transposed_root_ray_pairs"] == [[8, 9]]
    assert roots["scaled_wall_gram_invariant"] is True
    assert roots["wall_gram_inertia"] == {"positive": 9, "negative": 1, "zero": 1}
    assert roots["all_pairs_coxeter"] is True
    assert roots["outer_cone"]["finite_vertices"]["quotient_orbits"] == 3
    assert roots["outer_cone"]["hyperideal_rays"]["quotient_orbits"] == 6
    assert roots["outer_cone"]["quotient_still_has_hyperideal_rays"] is True
