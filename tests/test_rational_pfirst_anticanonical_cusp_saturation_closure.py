from scripts.certify_rational_pfirst_anticanonical_cusp_saturation_closure import (
    payload,
)


def test_exact_anticanonical_cusp_saturation_closure() -> None:
    data = payload()
    assert data["simple_reflection_descent_steps"] == 10
    assert data["saturated_root_squares"] == [-1] + [-2] * 9
    assert data["saturated_canonical_pairings"] == [-1] + [0] * 9
    assert data["saturated_root_basis_determinant"] == -1
    assert data["saturated_coxeter_order_histogram"] == {
        "2": 36,
        "3": 8,
        "4": 1,
    }
    assert data["old_walls_with_nonnegative_simple_expansion"] == 32
    assert data["old_walls_remaining_simple_facets"] == 1
    assert data["vertex_squares"] == [0, 0, 2, 6, 5, 4, 3, 2, 1, 1]
    assert data["characteristic_ideal_vertex_count"] == 1
    assert data["hamming_point_equals_sum_of_two_cusps"] is True
    assert data["anticanonical_affine_marks"] == [2, 4, 6, 5, 4, 3, 2, 1, 3]
    assert data["anticanonical_affine_support_euler"] == 10
    assert data["host_reflection_maximality"][
        "relative_reflection_maximality_in_O_plus_L"
    ] is True
    assert data["explicit_realization"]["pencil"]["total_blowups"] == 9
    assert data["explicit_realization"]["trivial_lattice_basis_determinant"] == -1
    assert data["explicit_realization"]["full_nef_chamber_returned"] is True
    assert data["anticanonical_cusp_rigidity"][
        "twice_minus_four_wall_at_F_cusp_possible"
    ] is False
    assert data["scope"]["hamming_saturation_outcome"] == "MAXIMAL_AND_REALIZABLE"
    assert data["scope"]["all_ranks_or_non_anticanonical_cusps"] == "OPEN"
