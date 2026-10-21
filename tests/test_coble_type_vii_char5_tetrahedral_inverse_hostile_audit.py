from __future__ import annotations

from functools import lru_cache
import json

from scripts.audit_coble_type_vii_char5_tetrahedral_inverse_hostile import (
    OUTPUT,
    build_audit,
)


@lru_cache(maxsize=1)
def payload() -> dict[str, object]:
    return build_audit()


def test_independent_finite_geometry_and_curve_systems() -> None:
    audit = payload()
    geometry = audit["finite_geometry"]
    assert geometry["projective_plane_point_count"] == 31
    assert geometry["perfect_matching_count"] == 15
    assert geometry["concurrent_matching_count"] == 10
    assert geometry["one_factorization_count"] == 6
    assert geometry["one_factorization_concurrency_histogram"] == {"0": 1, "4": 5}
    assert geometry["usable_factorization_count"] == 5
    assert (geometry["line_count"], geometry["conic_count"], geometry["quartic_count"]) == (
        12,
        4,
        4,
    )
    assert all(geometry["conic_determinants"])
    assert geometry["quartic_condition_ranks"] == [14, 14, 14, 14]

    systems = audit["linear_systems"]
    assert systems["anticanonical_cubic"] == {
        "ambient_dimension": 10,
        "evaluation_rank": 10,
        "kernel_dimension": 0,
    }
    sextic = systems["antibicanonical_sextic"]
    assert (sextic["ambient_dimension"], sextic["gradient_rank"], sextic["kernel_dimension"]) == (
        28,
        27,
        1,
    )
    assert sextic["value_augmented_rank"] == 27


def test_exact_quartic_and_boundary_singular_schemes() -> None:
    audit = payload()
    schemes = audit["scheme_theoretic_curve_audit"]
    assert schemes["all_quartics_irreducible_over_F5"] is True
    assert schemes["all_quartics_geometrically_integral"] is True
    for quartic in schemes["quartics"]:
        assert quartic["jacobian_degree"] == quartic["radical_degree"] == 3
        assert quartic["minimal_primes"] == 3
        assert quartic["jacobian_mod_radical"] == 0
        assert quartic["radical_mod_jacobian"] == 0
        assert quartic["jacobian_mod_points"] == 0
        assert quartic["points_mod_jacobian"] == 0

    boundary = schemes["boundary"]
    assert schemes["boundary_irreducible_over_F5"] is True
    assert schemes["boundary_geometrically_integral"] is True
    assert boundary["jacobian_degree"] == 16
    assert boundary["radical_degree"] == 10
    assert boundary["minimal_primes"] == 10
    assert boundary["jacobian_mod_radical"] == 0
    assert boundary["radical_mod_jacobian"] != 0
    assert boundary["radical_mod_points"] == boundary["points_mod_radical"] == 0


def test_local_four_A1_six_A2_and_boundary_disjointness() -> None:
    audit = payload()
    local = audit["local_geometry_and_transversality"]
    assert (local["boundary_A1_count"], local["boundary_A2_count"]) == (4, 6)
    assert local["boundary_node_discriminant_histogram"] == {"2": 2, "3": 2}
    assert local["boundary_cusp_cubic_histogram"] == {
        "1": 1,
        "2": 3,
        "3": 1,
        "4": 1,
    }
    assert local["quartic_node_count"] == 12
    assert local["quartic_node_discriminant_histogram"] == {"2": 6, "3": 6}
    assert local["quartic_boundary_resultant_histogram"] == {"1": 6, "4": 6}
    assert local["smooth_root_incidence_count"] == 84
    assert local["smooth_root_tangent_value_histogram_by_kind"] == {
        "line": {"1": 18, "4": 18},
        "conic": {"1": 12, "4": 12},
        "quartic": {"1": 12, "4": 12},
    }
    assert local["all_root_branches_avoid_boundary_branches"] is True

    strong = audit["blowup_and_strong_Coble"]
    assert strong["minus_K_empty"] is True
    assert strong["minus_2K_unique"] is True
    assert strong["boundary_Jacobian_length"] == 16
    assert strong["boundary_radical_degree"] == 10
    assert strong["boundary_strict_transform_smooth_rational"] is True
    assert strong["strong_Coble_verdict"] == "PASS"


def test_all_frames_markings_and_full_wall_transport() -> None:
    audit = payload()
    lattice = audit["raw_lattice_and_cone"]
    assert lattice["primitive_wall_norm_histogram"] == {"-2": 20, "-1": 52}
    assert lattice["intrinsic_boundary_wall_index_zero_based"] == 0
    assert lattice["maximal_disjoint_clique_size_histogram"] == {
        "7": 80,
        "9": 25,
        "10": 5,
    }
    assert lattice["intrinsic_plane_frame_count"] == 5
    assert lattice["all_walls_primitive"] is True
    assert lattice["CM_gram_determinant"] == -1
    assert lattice["CM_root_span_rank"] == 10
    assert lattice["root_pairing_histogram"] == {"0": 135, "1": 30, "2": 25}
    assert lattice["all_twenty_root_reflections_integral_involutive_isometries"] is True
    assert [frame["frame_set_wall_indices_zero_based"] for frame in lattice["all_five_frames"]] == [
        [21, 22, 23, 24, 50, 51, 56, 59, 62, 65],
        [21, 28, 29, 30, 50, 51, 56, 59, 62, 65],
        [22, 28, 34, 35, 50, 51, 56, 59, 62, 65],
        [23, 29, 34, 39, 50, 51, 56, 59, 62, 65],
        [24, 30, 35, 39, 50, 51, 56, 59, 62, 65],
    ]
    assert all(frame["basis_determinant"] == 1 for frame in lattice["all_five_frames"])
    assert all(frame["boundary_coordinates"] == [3, *([-1] * 10)] for frame in lattice["all_five_frames"])
    assert all(frame["colored_marking_count"] == 24 for frame in lattice["all_five_frames"])
    assert all(frame["exact_root_vector_marking_count"] == 24 for frame in lattice["all_five_frames"])
    assert lattice["all_120_markings_exact"] is True
    assert lattice["aggregate_full_wall_coordinate_sha256"] == (
        "74be27cfc3a7e6d21264b3e85d08b2d1c956e9cfcbaf12fdb1e9f212453d0c49"
    )


def test_direct_root_cone_and_non_circular_same_P_handoff() -> None:
    audit = payload()
    cone = audit["raw_lattice_and_cone"]["twenty_root_cone"]
    assert (cone["finite_ray_count"], cone["ideal_ray_count"], cone["hyperideal_ray_count"]) == (
        450,
        51,
        0,
    )
    assert cone["ray_count"] == 501
    assert cone["primitive_ray_list_sha256"] == (
        "80e1e81adc197bb92e818f52f3e695aa1d4d61fc2f4aecf6a550ea1b99fbb5ec"
    )
    assert cone["primitive_ideal_ray_list_sha256"] == (
        "8621adee95e09a0213356e044f923a0fd902e22aa55ead6ffc518959fbfc73d7"
    )
    assert cone["minimum_root_slack"] > 0
    assert cone["minimum_time_sheet_pairing"] > 0
    assert cone["facet_face_ranks"] == [9] * 20
    assert cone["all_twenty_root_inequalities_irredundant"] is True
    assert len(cone["cusp_records"]) == 51
    assert cone["ideal_rays_equal_all_51_E_plus_b_classes"] is True

    upstream = audit["raw_lattice_and_cone"]["upstream_P_side_vinberg_import"]
    assert upstream["facet_count"] == 72
    assert (upstream["finite_vertices"], upstream["ideal_vertices"], upstream["hyperideal_rays"]) == (
        18761,
        1731,
        0,
    )
    assert upstream["bad_edges"] == upstream["bad_parabolic_neighbourhoods"] == 0
    handoff = audit["AG_same_P_handoff"]
    assert handoff["surface_existence_or_Type_VII_classification_used_to_construct_carrier"] is False
    assert handoff["general_AG_theorems_used_after_explicit_carrier_exists"] is True
    assert handoff["exact_same_marked_H10_72"] is True


def test_scope_target_agreement_and_frozen_archive() -> None:
    audit = payload()
    assert audit["verdict"] == "GREEN_INDEPENDENT_L3_SAME_MARKED_H10_72"
    assert audit["sources"]["production_module_imported"] is False
    assert all(audit["target_certificate_agreement"].values())
    scope = audit["lineage_and_scope"]
    assert scope["new_abstract_hyperbolic_polytope"] is False
    assert scope["new_dimension"] is False
    assert scope["inverse_level"] == "L3 structural one-datum calibration, not class-wide L4"
    assert "not ten-nodal" in audit["hostile_warning"]
    assert audit["all_checks_pass"] is True
    assert all(audit["checks"].values())
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == json.loads(
        json.dumps(audit, sort_keys=True)
    )
