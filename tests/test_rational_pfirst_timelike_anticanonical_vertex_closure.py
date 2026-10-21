from scripts.certify_rational_pfirst_timelike_anticanonical_vertex_closure import (
    payload,
)


def test_exact_timelike_anticanonical_vertex_closure() -> None:
    data = payload()
    assert data["ranks_covered"] == [4, 5, 6, 7, 8, 9]
    assert data["exceptional_counts_covered"] == [3, 4, 5, 6, 7, 8]
    for record, n in zip(data["records"], range(3, 9)):
        assert record["exceptional_count"] == n
        assert record["rank"] == n + 1
        assert record["canonical_square"] == 9 - n
        assert record["anticanonical_square"] == 9 - n
        assert record["anticanonical_is_timelike_vertex"] is True
        assert record["facet_squares"] == [-1] + [-2] * n
        assert record["canonical_pairings"] == [-1] + [0] * n
        assert record["facet_basis_determinant"] == -1
        assert record["q2_gram_determinant_absolute"] == 9 - n
        assert record["vertex_squares"][0] == 9 - n
        assert record["ideal_vertex_count"] == 1
        assert record["reflection_descent_gap_coefficient"] == 9 - n
        assert record["characteristic_rigidity_gap_coefficient"] == 9 - n
        assert record["flex_cluster_carrier"]["residual_base_ideal"] == (
            f"(x^{9 - n},w)"
        )
        assert record["flex_cluster_carrier"]["anticanonical_is_nef"] is True
        assert record["flex_cluster_carrier"]["full_nef_return"] is True
    assert data["uniform_reflection_maximality"][
        "all_integral_reflections_generated"
    ] is True
    assert data["uniform_characteristic_rigidity"][
        "ordinary_physical_maximal_chamber_is_Stab_K_conjugate"
    ] is True
    assert data["theorem_scope"]["classification_free_exact_carrier"] == "PROVED"
    assert data["theorem_scope"]["twice_minus_four_or_minus_K_outside_P"] == "OPEN"
