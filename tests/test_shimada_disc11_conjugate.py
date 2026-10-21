from coxeter_polytopes.exact.golden import Golden
from scripts.shimada_disc11_conjugate import (
    ARCHIVED_AMPLE_VECTOR,
    CONJUGATING_ROW_ISOMETRY,
    ample_witness,
    analyze_d68_minus_two_roots,
    conjugating_row_isometry_inverse,
    d68_wall_orbits,
    d68_minus_two_roots,
    verify_d68_restricted_facet_certificate,
    verify_d68_root_facet_certificate,
    verify_conjugate_marking_certificate,
)


def test_conjugate_marking_has_exact_ample_witness_in_D68() -> None:
    certificate = verify_conjugate_marking_certificate()

    assert certificate["verified"] is True
    assert certificate["scope"] == (
        "nef_ample_nonemptiness_only_not_a_Coxeter_or_compactness_claim"
    )
    assert certificate["row_isometry_determinant"] == "1"
    assert certificate["row_isometry_verified"] is True
    assert certificate["discriminant_action_trivial_verified"] is True
    assert certificate["phi_polynomial_verified"] is True
    assert certificate["phi_self_adjoint_verified"] is True
    assert certificate["phi_eigenvector_verified"] is True
    assert certificate["witness_square"] == {"a": "-2474", "b": "2688"}
    assert certificate["defining_wall_count"] == 28
    assert certificate["all_wall_pairings_strictly_positive"] is True
    assert certificate["minimum_wall_pairing"] == {"a": "1", "b": "0"}


def test_D68_certificate_archive_fingerprints() -> None:
    assert len(CONJUGATING_ROW_ISOMETRY) == 20
    assert all(len(row) == 20 for row in CONJUGATING_ROW_ISOMETRY)
    assert len(conjugating_row_isometry_inverse()) == 20
    assert ample_witness() == list(ARCHIVED_AMPLE_VECTOR)
    assert len(d68_wall_orbits()) == 17
    assert sum(len(orbit) for _, orbit in d68_wall_orbits()) == 28
    assert ARCHIVED_AMPLE_VECTOR[0] == Golden(1005, 1470)
    assert ARCHIVED_AMPLE_VECTOR[-1] == Golden(-495, -727)


def test_D68_actual_curve_walls_pass_conjugated_field_filter() -> None:
    root_data = d68_minus_two_roots()
    result = analyze_d68_minus_two_roots()

    assert len(root_data) == 17
    assert result["input_root_count"] == 17
    assert result["trace_minus_two_verified"] == 17
    assert result["positive_at_both_places"] == 17
    assert result["positive_shell_distribution"] == {
        "2": 9,
        "4 - 2*phi": 8,
    }
    assert result["field_integral_reflection_roots"] == 17
    assert result["field_span_rank"] == 9
    assert result["pair_sign_distribution"] == {
        "negative": 30,
        "positive": 8,
        "zero": 98,
    }
    assert result["all_ample_pairings_positive"] is True
    assert result["minimum_ample_pairing"] == {"a": "1", "b": "0"}
    assert result["all_eigen_pairings_positive"] is True


def test_D68_restricted_root_facets_have_exact_face_witnesses() -> None:
    result = verify_d68_root_facet_certificate()

    assert result["verified"] is True
    assert result["ambient_minus_two_wall_count"] == 17
    assert result["ambient_roots_on_restricted_facets"] == 16
    assert result["distinct_restricted_root_facets"] == 8
    assert result["redundant_restricted_root_walls"] == [[0, 0]]
    assert result["orbit_zero_positive_dependence_verified"] is True
    assert len(result["facets"]) == 8
    assert all(facet["target_pairings_zero"] for facet in result["facets"])
    assert all(
        facet["all_other_D68_pairings_positive"] for facet in result["facets"]
    )


def test_D68_complete_restricted_facet_list_is_exact() -> None:
    result = verify_d68_restricted_facet_certificate()

    assert result["verified"] is True
    assert result["ambient_wall_count"] == 28
    assert result["distinct_restricted_halfspaces"] == 19
    assert result["irredundant_restricted_facets"] == 11
    assert result["root_facets"] == 8
    assert result["internal_Borcherds_facets"] == 3
    assert result["internal_facet_wall_classes"] == [
        [[9, 0], [16, 0]],
        [[11, 1]],
        [[14, 0]],
    ]
    assert len(result["positive_face_witnesses"]) == 11
    assert len(result["positive_redundancy_dependences"]) == 8
    assert all(
        item["verified"] for item in result["positive_redundancy_dependences"]
    )
