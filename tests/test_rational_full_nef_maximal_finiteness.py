import json
from pathlib import Path

from scripts.certify_rational_full_nef_maximal_finiteness import (
    ARCHIVE,
    build_certificate,
    gaussian_binomial,
    number_of_binary_subspaces,
)


def certificate():
    return build_certificate()


def test_direct_rational_host_list_is_finite_and_unimodular():
    cert = certificate()
    assert cert["checks"]["all_direct_hosts_unimodular"]
    assert cert["checks"]["direct_host_count_is_10_in_positive_dimension"]
    assert cert["direct_positive_dimensional_host_count_upper_bound"] == 10
    assert [row["rank"] for row in cert["direct_hosts"]] == list(range(1, 11))
    assert len(cert["direct_hosts"][1]["forms"]) == 2
    assert all(len(row["forms"]) == 1 for row in cert["direct_hosts"][2:])


def test_u_plus_minus_one_has_the_claimed_integral_odd_split():
    row = certificate()["controls"]["U_plus_minus_one_split"]
    assert abs(row["change_of_basis_determinant"]) == 1
    assert row["conjugated_gram"] == row["target_gram"]
    assert row["target_gram"] == [[1, 0, 0], [0, -1, 0], [0, 0, -1]]


def test_binary_subspace_formula_and_finite_coble_host_bounds():
    cert = certificate()
    assert [number_of_binary_subspaces(rank) for rank in range(1, 6)] == [2, 5, 16, 67, 374]
    assert gaussian_binomial(4, 2) == 35
    assert cert["checks"]["all_binary_subspace_bounds_hold"]
    bounds = cert["coble_extension_bounds"]
    assert bounds["complex_ranks"] == list(range(11, 21))
    assert bounds["odd_positive_ranks"] == list(range(11, 23))
    assert bounds["complex_marked_extension_upper_bound"] == 9334239813471146416714351279847
    assert (
        bounds["odd_positive_marked_extension_upper_bound"]
        == 19586680934637969160979021208209691846
    )
    assert bounds["bound_is_before_integrality_and_isometry_quotients"]


def test_standard_half_overlattice_determinant_scaling():
    cert = certificate()
    assert cert["checks"]["all_determinant_scaling_identities_hold"]
    rows = cert["determinant_scaling"]
    rank_20_middle = next(
        row for row in rows if row["rank"] == 20 and row["H_dimension"] == 10
    )
    assert rank_20_middle["base_abs_determinant"] == 2**20
    assert rank_20_middle["overlattice_index"] == 2**10
    assert rank_20_middle["extension_abs_determinant"] == 1


def test_f4_maximal_triangle_is_outside_both_ag_rank_gates():
    cert = certificate()
    row = cert["controls"]["F4_maximal_non_coble"]
    assert row["NS_unimodular"]
    assert row["physical_gram"] == [[-4, 0, 1], [0, -1, 1], [1, 1, -1]]
    assert row["half_extension_abs_determinant"] == 2
    assert row["half_root_gram"] == [[-2, 0, 1], [0, -2, 2], [1, 2, -2]]
    assert row["log_class_pairings"] == {"B/2": 0, "E": -2, "G": -1}
    assert row["maximal_triangle"] == "(2,3,infinity)"
    assert row["area_over_pi"] == "1/6"
    assert not row["passes_direct_stage0_gate"]
    assert not row["passes_coble_log_gate"]
    assert cert["checks"]["F4_is_maximal_but_outside_both_rank_branches"]


def test_nonmaximal_fixed_lattice_ideal_polygon_tower():
    cert = certificate()
    rows = cert["controls"]["nonmaximal_fixed_lattice_ideal_polygon_tower"]
    assert [row["facet_count"] for row in rows] == list(range(3, 27))
    assert [row["area_over_pi"] for row in rows] == list(range(1, 25))
    assert all(set(row["root_squares"]) == {-1} for row in rows)
    assert all(set(row["adjacent_pairings"]) == {1} for row in rows)
    assert all(all(value > 1 for value in row["nonadjacent_pairings"]) for row in rows)
    assert all(row["integral_reflections"] and row["reflection_isometries"] for row in rows)
    assert [
        row["index"] for row in rows if row["proper_reflection_subgroup_of_W1"]
    ] == list(range(2, 25))
    assert cert["checks"]["ideal_tower_rows_above_base_are_certifiably_nonmaximal"]
    assert cert["checks"]["ideal_tower_never_passes_direct_canonical_gate"]


def test_fixed_ambient_lattice_has_unbounded_saturated_eigenslices():
    cert = certificate()
    rows = cert["controls"]["fixed_host_unbounded_eigenslices"]
    assert [row["gram"] for row in rows] == [
        [[2 * m - 1, 0], [0, -1]] for m in range(1, 25)
    ]
    assert [row["abs_determinant"] for row in rows] == list(range(1, 49, 2))
    assert cert["checks"]["eigenslice_discriminants_are_unbounded_sample"]
    assert cert["checks"]["eigenslice_rows_are_saturated"]


def test_theorem_scope_and_field_are_not_overstated():
    theorem = certificate()["theorem"]
    assert theorem["scope"] == "literal full-dimensional rational full-nef carriers"
    assert theorem["saturation_conclusion"] == "W(P)=R(M)"
    assert theorem["fixed_host_chamber_orbits"] == 1
    assert theorem["minimal_vinberg_field"] == "Q"
    assert not theorem["continuous_surface_moduli_classified"]
    assert not theorem["eigenslices_covered"]
    assert not theorem["arbitrary_commensurable_modules_covered"]
    assert not theorem["global_carrier_existence_proved"]


def test_frozen_archive_is_exact_rebuild():
    expected = json.dumps(certificate(), indent=2, sort_keys=True) + "\n"
    assert Path(ARCHIVE).read_text() == expected
