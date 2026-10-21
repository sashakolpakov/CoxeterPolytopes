from scripts.certify_golden_lanner_hodge_coherent_carriers import analyze


def test_all_twelve_complete_lanner_returns_are_hodge_coherent() -> None:
    result = analyze()
    assert result["positive_count"] == 12
    assert result["all_certified_rows_pass"] is True
    assert len(result["rows"]) == 12
    for row in result["rows"]:
        assert row["same_polytope_input_certified"] is True
        assert row["rationally_isometric_to_K3_lattice"] is True
        assert row["period_plane_available"] is True
        assert row["golden_action_minimal_polynomial"] is True
        assert row["golden_action_self_adjoint"] is True
        assert row["sqrt5_Hodge_operator_square_is_5I"] is True
        assert row["sqrt5_Hodge_similarity_multiplier"] == 5
        assert row["kuga_satake_isogeny_follows"] is True
        assert row["X_times_X_algebraic_correspondence_certified"] is False
        assert row["kuga_satake_isogeny_used_as_carrier_gate"] is False
        assert row["picard_action_matches_upstream"] is True
        assert row["picard_action_minimal_polynomial"] is True
        assert row["picard_action_self_adjoint"] is True
        assert row["full_H2_action_unital_and_faithful"] is True
        assert row["full_H2_action_self_adjoint"] is True
        assert row["distinguished_eigenspace_signature"][0] == 2
        assert row["conjugate_eigenspace_signature"][0] == 0
        assert row["witt_cancellation_for_every_fixed_primitive_embedding"] is True
        assert row["single_real_eigenspace_has_no_nonzero_Q_vectors"] is True
        assert row["exact_NS_Hodge_period_locus_nonempty"] is True
        assert row["some_finite_index_integral_order_exists"] is True
        assert row["maximal_order_integrality_certified"] is False
        assert row["explicit_order_conductor_computed"] is False
        assert row["status"] == (
            "certified_rational_full_H2_Hodge_coherent_same_P_return"
        )


def test_exact_twists_and_dimensions() -> None:
    result = analyze()
    rows = {row["seed_id"]: row for row in result["rows"]}
    assert rows["triangle_255"]["twist"] == {"a": "1", "b": "0"}
    assert rows["triangle_335"]["twist"] == {"a": "2", "b": "0"}
    assert rows["triangle_355"]["twist"] == {"a": "2", "b": "1"}
    assert rows["triangle_555"]["twist"] == {"a": "2", "b": "1"}
    assert rows["star_5_3_3"]["twist"] == {"a": "1", "b": "0"}
    assert rows["cycle4_one5"]["twist"] == {"a": "2", "b": "1"}
    assert rows["cycle4_opposite5"]["twist"] == {"a": "3", "b": "0"}
    assert rows["fork5_terminal5"]["twist"] == {"a": "1", "b": "0"}
    assert {row["dimension"] for row in result["rows"]} == {2, 3, 4}
    assert {row["picard_rank"] for row in result["rows"]} == {6, 8, 10}
    assert {row["complement_K_rank"] for row in result["rows"]} == {6, 7, 8}


def test_k3_hasse_minkowski_target_is_frozen() -> None:
    target = analyze()["k3_rational_invariants"]
    assert target["rank"] == 22
    assert target["signature"] == [3, 19]
    assert target["determinant"] == "-1"
    assert target["hasse"] == {
        "2": -1,
        "3": 1,
        "5": 1,
        "7": 1,
        "infinity": -1,
    }
    assert target["all_nontrivial_places_enumerated"] is True
    assert target["hilbert_reciprocity_product"] == 1


def test_every_total_form_has_complete_hasse_places_and_reciprocity() -> None:
    for row in analyze()["rows"]:
        for key in (
            "picard_rational_invariants",
            "complement_rational_invariants",
            "total_rational_invariants",
        ):
            invariants = row[key]
            assert invariants["all_nontrivial_places_enumerated"] is True
            assert invariants["hilbert_reciprocity_product"] == 1


def test_bounded_discovery_recomputes_every_frozen_twist() -> None:
    result = analyze(search_bound=5)
    assert result["positive_count"] == 12
    for row in result["rows"]:
        assert row["bounded_discovery_candidates_tested"] > 0
        assert row["bounded_discovery_twist"] == row["twist"]
        assert row["bounded_discovery_matches_frozen"] is True
