from coxeter_polytopes.exact.golden import Golden, PHI
from scripts.analyze_shimada_disc11 import (
    OPPOSITE_WALL_MINUS,
    OPPOSITE_WALL_PLUS,
    analyze_data,
    eigen_nef_opposite_wall_obstruction,
    expected_shimada_gram,
    field_vector_from_shimada,
)


def test_shimada_basis_matches_our_real_multiplication_trace_lattice() -> None:
    # The last coordinate is a simple root in Shimada's second E8(-1).
    root = [0] * 19 + [1]
    vector = field_vector_from_shimada(root)
    expected_vector = [Golden(0)] * 10
    expected_vector[8] = PHI - 1
    assert vector == expected_vector

    result = analyze_data(expected_shimada_gram(), [root])
    assert result["shimada_gram_matches_published_normal_form"] is True
    assert result["basis_change_recovers_shimada_gram"] is True
    assert result["published_root_representatives"] == 1
    assert result["trace_minus_two_verified"] == 1
    assert result["positive_at_both_places"] == 1
    assert result["positive_shell_distribution"] == {"4 - 2*phi": 1}
    assert result["field_integral_reflection_roots"] == 1
    assert result["field_integral_shell_distribution"] == {"4 - 2*phi": 1}
    assert result["field_span_rank"] == 1
    assert result["pair_sign_distribution"] == {}
    assert result["seed_matches_two_e8_bases_plus_one_coupled_root"] is False
    assert result["field_integral_root_seed"][0]["shimada_coordinates"] == root
    assert result["field_integral_root_seed"][0]["field_norm"] == {"a": "4", "b": "-2"}


def test_two_effective_walls_prove_fixed_eigen_nef_section_empty() -> None:
    obstruction = eigen_nef_opposite_wall_obstruction(
        [OPPOSITE_WALL_PLUS, OPPOSITE_WALL_MINUS]
    )
    assert obstruction["verified"] is True
    assert obstruction["q_e_e"] == {"a": "2/5", "b": "6/5"}
    plus = obstruction["walls"]["d_plus"]
    minus = obstruction["walls"]["d_minus"]
    assert plus["field_coordinates"][:2] == [
        {"a": "1", "b": "1"},
        {"a": "-1", "b": "2"},
    ]
    assert plus["h_norm"] == {"a": "6", "b": "-4"}
    assert plus["q_pairing_with_e"] == {"a": "8/5", "b": "9/5"}
    assert minus["field_coordinates"][:2] == [
        {"a": "-3", "b": "-1"},
        {"a": "1", "b": "-3"},
    ]
    assert minus["h_norm"] == {"a": "8", "b": "-6"}
    assert minus["q_pairing_with_e"] == {"a": "-11/5", "b": "-18/5"}
    assert plus["trace_square"] == minus["trace_square"] == "-2"
