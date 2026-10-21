from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path

from scripts.audit_mobius_q2_cubic_p6_hostile import OUTPUT, build_audit


@lru_cache(maxsize=1)
def audit() -> dict[str, object]:
    return build_audit()


def test_no_production_module_is_imported() -> None:
    source = Path("scripts/audit_mobius_q2_cubic_p6_hostile.py").read_text(
        encoding="utf-8"
    )
    assert "from scripts.certify_mobius_q2_cubic_p6 import" not in source
    assert "import scripts.certify_mobius_q2_cubic_p6" not in source
    assert audit()["target_comparison"]["production_module_imported"] is False


def test_f4_graphs_general_position_and_hermitian_cubic() -> None:
    geometry = audit()["finite_geometry"]
    graphs = geometry["Frobenius_graph_intersection"]
    assert graphs["intersection_number"] == 5
    assert graphs["listed_geometric_intersection_count"] == 5
    assert graphs["local_Jacobian_determinant"] == 1
    assert graphs["transverse_at_all_five_points"] is True
    assert graphs["Bezout_exhaustive"] is True

    general = geometry["general_position"]
    assert general["triple_O_1_1_rank_histogram"] == {"3": 10}
    assert general["four_point_O_1_1_rank_histogram"] == {"4": 5}
    assert general["no_four_on_O_1_1"] is True
    assert general["all_ten_triple_curves_smooth_irreducible_projectivity_graphs"] is True

    model = geometry["anticanonical_model"]
    assert model["O_2_2_evaluation_rank"] == 5
    assert model["section_dimension"] == 4
    assert model["relation_dimensions"] == {"1": 0, "2": 0, "3": 1}
    assert model["equation"] == "z0^2*z3+z0*z3^2+z1^2*z2+z1*z2^2=0"
    assert model["Hermitian_matrix_rank"] == 4
    assert model["gradient_has_no_geometric_projective_zero"] is True
    assert model["P3_F4_point_count"] == 85
    assert model["cubic_F4_point_count_by_direct_enumeration"] == 45
    assert model["geometrically_smooth"] is True


def test_all_27_curves_and_literal_ert_marking_are_exhaustive() -> None:
    lattice = audit()["surface_lattice_and_curves"]
    assert lattice["NS_determinant"] == 1
    assert lattice["NS_signature"] == [1, 6, 0]
    assert lattice["canonical_square"] == 3
    assert lattice["curve_species"] == {
        "exceptional": 5,
        "ruling_fibre": 10,
        "O_1_1_projectivity_graph": 10,
        "Frobenius_graph": 2,
    }
    assert lattice["curve_count"] == 27
    assert lattice["pairing_histogram"] == {"0": 216, "1": 135}
    assert lattice["exceptional_class_exhaustion"]["solution_histogram"] == {
        "0": 6,
        "1": 15,
        "2": 6,
    }
    assert lattice["exceptional_class_exhaustion"][
        "every_solution_is_one_of_the_displayed_curves"
    ] is True
    assert lattice["blowdown_marking"]["basis_determinant_absolute"] == 1
    assert lattice["ERT_Table_2"]["normal_multiset_exact"] is True
    assert lattice["full_NS_reflections"]["all_27_integral"] is True
    assert lattice["full_NS_reflections"]["all_identity_modulo_two"] is True
    assert lattice["Mori_nef_gate"]["the_27_halfspaces_are_the_complete_nef_cone"] is True


def test_exact_cdd_has_72_finite_27_ideal_and_no_hyperideal_ray() -> None:
    cone = audit()["exact_projective_cone"]
    assert cone["ray_count"] == 99
    assert (cone["finite"], cone["ideal"], cone["hyperideal"]) == (72, 27, 0)
    assert cone["Lorentz_square_histogram"] == {"-1": 72, "0": 27}
    assert cone["active_facet_histogram"] == {"finite_6": 72, "ideal_10": 27}
    assert cone["every_wall_irredundant"] is True
    assert len(cone["active_sets"]["finite"]) == 72
    assert len(cone["active_sets"]["ideal"]) == 27
    assert cone["active_sets"]["hyperideal"] == []


def test_literal_vinberg_endpoint_and_full_cusp_criteria() -> None:
    vinberg = audit()["literal_Vinberg"]
    assert vinberg["elliptic_rank_5_edges"] == 648
    assert vinberg["finite_rank_6_vertices"] == 72
    assert vinberg["connected_affine_components"] == 135
    assert vinberg["full_zero_neighborhood_size_histogram"] == {"8": 135}
    assert vinberg["maximal_rank_5_parabolic_cusp_supports"] == 27
    assert vinberg["affine_component_sizes_per_cusp"] == [2, 2, 2, 2, 2]
    assert vinberg["connected_affine_zero_neighborhood_audits"] == 135
    assert vinberg["invalid_zero_neighborhood_audits"] == 0
    assert vinberg["edge_endpoint_histogram_finite_then_ideal"] == {
        "0,2": 216,
        "1,1": 432,
    }
    assert vinberg["bad_finite_volume_edges"] == 0
    assert vinberg["bad_compact_edges"] == 648
    assert vinberg["cdd_finite_active_sets_equal_all_elliptic_vertices"] is True
    assert vinberg["cdd_ideal_active_sets_equal_all_full_parabolic_cusps"] is True
    assert vinberg["decision"] == "finite_volume_noncompact"
    assert vinberg["finite_volume"] is True
    assert vinberg["compact"] is False


def test_full_ns_q_dichotomy_scope_lineage_and_frozen_archive() -> None:
    result = audit()
    dichotomy = result["full_NS_factory_dichotomy"]
    assert dichotomy["circle_divisibility_in_full_NS"] == 1
    assert dichotomy["necessary_and_sufficient_criterion"] == "q-1 divides 2"
    assert dichotomy["prime_power_solutions"] == [2, 3]
    assert all(
        record["integral_full_NS_reflection"] == (q in {"2", "3"})
        for q, record in dichotomy["sample_arithmetic"].items()
    )
    assert "not abstract" in dichotomy["scope"]

    assert result["verdict"] == "GREEN_NO_MATHEMATICAL_DEFECT"
    assert result["scope"]["direction"] == "SURFACE_TO_POLYTOPE positive factory calibration"
    assert result["scope"]["inverse_level"] == (
        "L3 exact structural calibration, not arbitrary-P class-wide L4"
    )
    assert result["scope"]["new_abstract_polytope"] is False
    assert result["scope"]["new_surface"] is False
    assert result["primary_lineage"]["lineage_verdict"].startswith(
        "surface, 27-line configuration, and abstract P6 are classical"
    )
    assert all(result["target_comparison"]["checks"].values())
    assert result["all_checks_pass"] is True
    assert all(result["checks"].values())
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == result
