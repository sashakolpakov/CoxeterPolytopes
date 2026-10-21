from __future__ import annotations

import json

from scripts.certify_mobius_q2_cubic_p6 import OUTPUT, build_certificate


CERTIFICATE = build_certificate()


def test_f4_frobenius_centers_and_steiner_graphs() -> None:
    finite = CERTIFICATE["finite_geometry"]
    assert len(finite["centers"]) == 5
    graphs = finite["frobenius_graphs"]
    assert graphs["intersection_number"] == 5
    assert graphs["transverse_at_all_five_centers"] is True
    assert graphs["bezout_proves_no_further_geometric_intersections"] is True
    general = finite["general_position"]
    assert general["triple_O_1_1_evaluation_rank_histogram"] == {"3": 10}
    assert general["four_point_O_1_1_evaluation_rank_histogram"] == {"4": 5}
    assert general["no_four_centers_on_an_O_1_1_curve"] is True
    design = finite["S_3_3_5"]
    assert design["block_count"] == 10
    assert design["point_replication"] == [6]
    assert design["pair_replication"] == [3]
    assert design["triple_replication"] == [1]
    assert design["all_ten_block_curves_are_smooth_irreducible_projectivity_graphs"] is True


def test_anticanonical_model_is_the_smooth_hermitian_cubic() -> None:
    model = CERTIFICATE["finite_geometry"]["anticanonical_model"]
    assert model["O_2_2_evaluation_rank"] == 5
    assert model["section_dimension_after_five_basepoints"] == 4
    assert model["relation_kernel_dimensions_by_degree"] == {"1": 0, "2": 0, "3": 1}
    assert model["equation"] == "z0^2*z3+z0*z3^2+z1^2*z2+z1*z2^2=0"
    assert model["Hermitian_matrix_rank"] == 4
    assert model["gradient"] == ["z3^2", "z2^2", "z1^2", "z0^2"]
    assert model["geometrically_smooth"] is True
    assert model["F4_rational_point_count_from_blowup"] == 45


def test_all_27_exceptional_classes_are_realized_and_exhaustive() -> None:
    surface = CERTIFICATE["surface"]
    assert surface["canonical_square"] == 3
    assert surface["NS_determinant"] == 1
    assert surface["NS_signature"] == [1, 6, 0]
    assert surface["NS_identification"] == "odd unimodular I_(1,6)"
    assert surface["F4_rational_minus_one_curve_count"] == 27
    assert surface["all_27_minus_one_curves_are_geometrically_realized"] is True
    marking = CERTIFICATE["P2_blowdown_and_ERT_identification"]
    assert marking["basis_determinant_absolute"] == 1
    assert marking["exceptional_class_solution_histogram_by_degree"] == {
        "0": 6,
        "1": 15,
        "2": 6,
    }
    assert marking["all_exceptional_lattice_solutions_realized"] is True


def test_complete_nef_chamber_has_exact_p6_census() -> None:
    chamber = CERTIFICATE["line_chamber"]
    assert chamber["wall_count"] == 27
    assert chamber["distinct_wall_pairing_histogram"] == {"0": 216, "1": 135}
    assert chamber["all_distinct_pairings_are_zero_or_one"] is True
    reflections = chamber["reflections"]
    assert reflections["all_27_integral"] is True
    assert reflections["all_preserve_NS_gram"] is True
    assert reflections["all_congruent_to_identity_mod_2"] is True
    cone = chamber["exact_double_description"]
    assert cone["extreme_ray_count"] == 99
    assert cone["finite_vertex_count"] == 72
    assert cone["ideal_vertex_count"] == 27
    assert cone["hyperideal_vertex_count"] == 0
    assert cone["active_facet_count_histogram"] == {"finite_6": 72, "ideal_10": 27}
    assert cone["every_wall_is_irredundant"] is True


def test_literal_vinberg_criterion_distinguishes_finite_volume_from_compactness() -> None:
    vinberg = CERTIFICATE["line_chamber"]["literal_Vinberg"]
    assert vinberg["decision"] == "finite_volume_noncompact"
    assert vinberg["finite_volume"] is True
    assert vinberg["compact"] is False
    assert vinberg["elliptic_rank_5_edge_count"] == 648
    assert vinberg["edge_endpoint_histogram_finite_then_ideal"] == {
        "0,2": 216,
        "1,1": 432,
    }
    assert vinberg["finite_volume_bad_edge_count"] == 0
    assert vinberg["compact_bad_edge_count"] == 648
    assert vinberg["finite_vertex_count"] == 72
    assert vinberg["ideal_vertex_count"] == 27
    assert vinberg["invalid_full_zero_neighborhood_count"] == 0


def test_marking_is_literally_ert_p6_not_only_same_counts() -> None:
    marking = CERTIFICATE["P2_blowdown_and_ERT_identification"]
    assert marking["ERT_Table_2_normal_multiset_exact"] is True
    assert len(marking["line_records"]) == 27
    lineage = CERTIFICATE["lineage"]
    assert lineage["abstract_polytope"] == "Potyagailo--Vinberg P^6"
    assert lineage["Euclidean_dual"] == "Gosset/Coxeter 2_21"
    assert lineage["new_abstract_polytope"] is False
    assert lineage["new_surface"] is False
    assert "no novelty claim" in lineage["worldwide_priority"]


def test_full_ns_circle_reflections_select_exactly_q2_and_q3() -> None:
    dichotomy = CERTIFICATE["Mobius_factory_integral_reflection_dichotomy"]
    assert dichotomy["divisibility"] == 1
    assert dichotomy["integrality_criterion"] == "q-1 divides 2"
    assert dichotomy["prime_power_solutions"] == [2, 3]
    samples = dichotomy["sample_exact_arithmetic"]
    assert samples["2"]["circle_reflection_integral"] is True
    assert samples["3"]["circle_reflection_integral"] is True
    assert all(
        not record["circle_reflection_integral"]
        for q, record in samples.items()
        if q not in {"2", "3"}
    )
    assert "not nonexistence" in dichotomy["scope_warning"]


def test_scope_and_frozen_archive() -> None:
    assert CERTIFICATE["all_checks_pass"] is True
    assert CERTIFICATE["inverse_scope"]["rating"] == "L3 exact structural calibration, not class-wide L4"
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == CERTIFICATE
