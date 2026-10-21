from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
import subprocess
import sys

from scripts.certify_coble_type_ii_b3_plane_inverse import (
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_coble_type_ii_b3_plane_inverse.py"


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    return build_certificate()


def test_raw_P_has_exactly_two_clean_7_plus_3_plane_branches() -> None:
    raw = certificate()["raw_polytope_extraction"]
    assert raw["input_policy"]["reads"] == [
        "integral lattice Gram",
        "thirty primitive wall vectors",
    ]
    assert raw["all_sequential_P2_endpoint_count"] == 573
    assert raw["endpoint_degree_signature_count"] == 28
    assert raw["selected_endpoint_count"] == 2
    assert raw["selected_signature"] == {
        "half_boundary_degree": 3,
        "root_degree_histogram": {"0": 3, "1": 9},
        "contracted_original_norm_histogram": {"-2": 3, "-1": 7},
    }
    assert [branch["contracted_wall_set_zero_based"] for branch in raw["branches"]] == [
        [1, 5, 9, 13, 15, 17, 19, 22, 26, 28],
        [3, 7, 11, 13, 15, 17, 19, 22, 26, 28],
    ]
    assert [branch["basis_determinant"] for branch in raw["branches"]] == [1, -1]
    assert all(
        branch["basis_gram"]
        == [[1 if i == j == 0 else -1 if i == j else 0 for j in range(11)] for i in range(11)]
        for branch in raw["branches"]
    )


def test_both_unimodular_bases_transport_the_complete_marked_wall_set() -> None:
    branches = certificate()["raw_polytope_extraction"]["branches"]
    assert [branch["plane_class_in_raw_basis"] for branch in branches] == [
        [13, 9, 19, 27, 23, 20, 15, 12, 7, 3, 3],
        [12, 9, 18, 26, 22, 19, 14, 11, 7, 3, 3],
    ]
    expected_roots = {
        "EA_minus_Ea",
        "ED_minus_Ed",
        "EG_minus_Eg",
        "L_DG",
        "L_AG",
        "L_AD",
        "L_EFG",
        "L_BCG",
        "L_ACF",
        "L_ABE",
        "L_CDE",
        "L_BDF",
    }
    for branch in branches:
        assert branch["canonical_class_coordinates_H_E"] == [-3, *([1] * 10)]
        assert branch["primitive_half_boundary_coordinates_H_E"] == [3, *([-1] * 10)]
        assert branch["colored_incidence_isomorphism_count"] == 12
        assert len(branch["all_thirty_raw_wall_coordinates_H_E"]) == 30
        assert {record["canonical_curve"] for record in branch["root_transport"]} == expected_roots
        assert branch["all_twelve_raw_roots_equal_canonical_curve_classes"] is True
        assert branch["all_thirty_wall_vectors_reconstruct_exactly"] is True


def test_rigid_B3_cluster_and_all_twelve_root_curves() -> None:
    geometry = certificate()["P_derived_incidence_realization"]
    rigidity = geometry["incidence_rigidity"]
    assert rigidity["last_collinearity_determinant"] == "-(c - 1)*(c + 1)/c"
    assert rigidity["nondegenerate_parameter_solution"] == "c=-1; c=1 gives C=B"
    assert set(rigidity["already_solved_collinearity_determinants"].values()) == {"0"}
    assert geometry["proper_points"] == [
        {"label": "A", "coordinates": [1, 0, 0]},
        {"label": "D", "coordinates": [0, 1, 0]},
        {"label": "G", "coordinates": [0, 0, 1]},
        {"label": "B", "coordinates": [1, 1, 1]},
        {"label": "C", "coordinates": [1, 1, -1]},
        {"label": "E", "coordinates": [-1, 1, 1]},
        {"label": "F", "coordinates": [-1, 1, -1]},
    ]
    arrangement = geometry["line_arrangement"]
    assert arrangement["name"] == "projectivized B3 reflection arrangement"
    assert arrangement["line_count"] == 9
    assert arrangement["proper_point_line_counts"] == {
        "A": 4,
        "D": 4,
        "G": 4,
        "B": 3,
        "C": 3,
        "E": 3,
        "F": 3,
    }
    assert arrangement["line_pairs_meeting_at_blown_proper_points"] == 30
    assert arrangement["remaining_distinct_transverse_double_intersections"] == 6
    inventory = geometry["root_inventory"]
    assert inventory["line_roots"] == 9
    assert inventory["exceptional_chain_roots"] == 3
    assert inventory["total_root_pair_intersection_one"] == 15
    assert inventory["total_root_pair_intersection_zero"] == 51


def test_Hilbert_system_forces_the_Choi_Lam_sextic_and_strong_Coble_functions() -> None:
    payload = certificate()
    systems = payload["linear_systems"]
    assert systems["anticanonical_cubics"]["ambient_vector_dimension"] == 10
    assert systems["anticanonical_cubics"]["condition_rank"] == 10
    assert systems["anticanonical_cubics"]["kernel_dimension"] == 0
    sextics = systems["antibicanonical_sextics"]
    assert sextics["ambient_vector_dimension"] == 28
    assert sextics["condition_row_count"] == 30
    assert sextics["condition_rank"] == 27
    assert sextics["kernel_dimension"] == 1
    assert sextics["polynomial"] == "x^4*y^2-3*x^2*y^2*z^2+x^2*z^4+y^4*z^2"
    assert payload["blowup_surface"]["boundary_class"] == "B=6H-2sum E_i=-2K"
    assert payload["blowup_surface"]["strong_Coble"] is True
    assert "Choi-Lam" in payload["lineage"]["sextic_equation_lineage"]


def test_absolute_integrality_and_complete_singular_scheme() -> None:
    singular = certificate()["absolute_factor_and_singular_scheme_audit"]
    assert singular["factor_count_with_unit_over_Q"] == 2
    assert singular["nonunit_exponent_over_Q"] == 1
    assert singular["absolute_factor_count_with_unit"] == 2
    assert singular["absolute_nonunit_exponent"] == 1
    assert singular["absolute_component_count"] == 1
    assert singular["absolutely_irreducible"] is True
    assert singular["singular_affine_cone_dimension"] == 1
    assert singular["singular_scheme_degree"] == 13
    assert singular["singular_radical_degree"] == 7
    assert singular["singular_minimal_primes"] == 7
    assert singular["radical_mod_seven_points"] == 0
    assert singular["seven_points_mod_radical"] == 0
    assert singular["singular_support_is_exactly_the_seven_proper_points"] is True


def test_local_types_infinitely_near_resolution_and_boundary_disjointness() -> None:
    payload = certificate()
    local = payload["local_boundary_resolution"]
    assert local["A3_tacnode_count"] == 3
    assert local["A1_node_count"] == 4
    assert local["delta_sum"] == 10
    by_point = {record["point"]: record for record in local["records"]}
    for label in ("A", "D", "G"):
        assert by_point[label]["tangent_discriminant"] == 0
        assert by_point[label]["first_blowup_near_multiplicity"] == 2
        assert by_point[label]["first_blowup_tangent_discriminant"] == -4
    for label in ("B", "C", "E", "F"):
        assert by_point[label]["tangent_discriminant"] == -48
    assert local["all_singularities_resolved_by_the_declared_7_plus_3_blowups"] is True
    assert local["strict_boundary_smooth"] is True
    assert local["strict_boundary_geometrically_integral_rational"] is True
    assert local["strict_boundary_self_intersection"] == -4
    assert payload["checks"]["all_twelve_root_curves_smooth_geometrically_integral"] is True
    assert "total intersection zero" in payload["root_boundary_disjointness"]["geometric"]


def test_same_marked_H10_30_return_scope_and_frozen_replay() -> None:
    payload = certificate()
    returned = payload["same_polytope_return"]
    assert returned["independent_CM_root_cone_gate"] == {
        "finite_rays": 38,
        "ideal_rays": 17,
        "hyperideal_rays": 0,
        "all_seventeen_cusps_equal_the_raw_E_plus_b_classes": True,
        "root_gram_is_Kondo_Type_II": True,
        "root_span_index_in_E10": 2,
        "unique_E10_saturation_up_to_diagram_automorphism": True,
    }
    assert returned["returned_wall_inventory"] == {
        "primitive_half_boundary_b_minus_K": 1,
        "explicit_minus_two_roots": 12,
        "automatic_cusp_lift_minus_one_curves": 17,
        "total_facets": 30,
    }
    assert returned["independent_raw_cone_gate"] == {
        "all_thirty_facets_irredundant": True,
        "finite_vertices": 1333,
        "ideal_vertices": 207,
        "hyperideal_vertices": 0,
        "elliptic_edges": 10030,
        "connected_affine_neighborhoods": 436,
        "bad_edges": 0,
        "bad_parabolic_neighborhoods": 0,
        "reflection_index_in_W_P10": 12142080,
    }
    assert returned["exact_same_marked_H10_30"] is True
    assert returned["finite_volume"] is True
    assert returned["noncompact"] is True
    assert payload["inverse_level"]["rating"] == (
        "L3_STRICT_SINGLE_RAW_P_INVERSE_REALIZATION_NOT_L4"
    )
    assert payload["field_scope"]["equations_and_blowups_defined_over"] == "Q"
    assert "characteristic 3" in payload["field_scope"]["excluded_overstatement"]
    assert "no new equation" in payload["lineage"]["priority_boundary"]
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))

    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--check-output", str(OUTPUT)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0
