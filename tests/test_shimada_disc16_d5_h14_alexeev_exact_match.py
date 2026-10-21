from __future__ import annotations

import json

from scripts.match_shimada_disc16_d5_h14_alexeev_15_7_1 import OUTPUT


def _load() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_alexeev_15_7_1_weighted_diagram_match() -> None:
    certificate = _load()

    assert certificate["intrinsic_lattice"] == {
        "determinant": 128,
        "same_lattice_by_indefinite_two_elementary_uniqueness": True,
        "signature": [1, 14],
        "two_elementary_invariants": {"a": 7, "delta": 1, "r": 15},
    }

    comparison = certificate["weighted_gram_matrix_comparison"]
    assert comparison["size"] == 60
    assert comparison["shimada_diagonal_histogram"] == {"-4": 10, "-2": 50}
    assert comparison["alexeev_diagonal_histogram"] == {"-4": 10, "-2": 50}
    assert comparison["shimada_off_diagonal_histogram"] == {
        "0": 865,
        "1": 20,
        "2": 395,
        "4": 120,
        "6": 310,
        "10": 60,
    }
    assert (
        comparison["shimada_off_diagonal_histogram"]
        == comparison["alexeev_off_diagonal_histogram"]
    )
    assert comparison["complete_edge_coloured_graph_isomorphic"] is True
    assert comparison["entrywise_equal_after_permutation"] is True
    assert sorted(comparison["shimada_facet_to_alexeev_root_permutation_one_based"]) == list(
        range(1, 61)
    )

    main = certificate["alexeev_main_root_basis_comparison"]
    assert main["count"] == 15
    assert abs(main["coordinate_basis_determinant"]) == 1
    assert main["gram_determinant"] == 128
    assert main["subdiagram"] == "K5^(2)"
    assert main["subdiagram_is_edge_subdivision_of_K5"] is True


def test_novelty_decision_and_parabolic_cross_check() -> None:
    certificate = _load()
    cross_check = certificate["independent_finite_volume_cross_check"]
    assert cross_check["shimada_exact_ideal_extreme_rays"] == 1027
    assert cross_check["alexeev_maximal_parabolic_subdiagrams"] == 1027
    assert cross_check["counts_equal"] is True

    decision = certificate["decision"]
    assert decision["same_intrinsic_lattice"] is True
    assert decision["same_weighted_coxeter_diagram"] is True
    assert decision["same_hyperbolic_chamber_up_to_isometry"] is True
    assert decision["new_hyperbolic_example"] is False
