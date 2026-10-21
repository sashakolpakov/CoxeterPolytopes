from __future__ import annotations

from functools import lru_cache
import json

from scripts.certify_coble_type_vii_char5_tetrahedral_inverse import (
    OUTPUT,
    build_certificate,
)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    return build_certificate()


def test_raw_polytope_forces_exactly_five_marked_plane_frames() -> None:
    payload = certificate()
    raw = payload["marked_same_polytope_return"]
    assert raw["primitive_wall_norm_histogram"] == {"-2": 20, "-1": 52}
    assert raw["intrinsic_boundary_wall_index_zero_based"] == 0
    assert raw["maximal_disjoint_clique_size_histogram"] == {
        "7": 80,
        "9": 25,
        "10": 5,
    }
    assert raw["intrinsic_plane_frame_count"] == 5
    assert [frame["frame_set_wall_indices_zero_based"] for frame in raw["all_five_frames"]] == [
        [21, 22, 23, 24, 50, 51, 56, 59, 62, 65],
        [21, 28, 29, 30, 50, 51, 56, 59, 62, 65],
        [22, 28, 34, 35, 50, 51, 56, 59, 62, 65],
        [23, 29, 34, 39, 50, 51, 56, 59, 62, 65],
        [24, 30, 35, 39, 50, 51, 56, 59, 62, 65],
    ]
    assert all(frame["basis_determinant"] == 1 for frame in raw["all_five_frames"])
    assert all(
        frame["boundary_coordinates_in_H_E_basis"] == [3, *([-1] * 10)]
        for frame in raw["all_five_frames"]
    )
    assert all(
        frame["root_kind_histogram"] == {"conic": 4, "line": 12, "quartic": 4}
        for frame in raw["all_five_frames"]
    )
    assert all(
        frame["colored_incidence_marking_count"] == 24
        and frame["all_24_markings_give_exact_twenty_root_vectors"] is True
        and frame["all_72_wall_coordinates_integral"] is True
        for frame in raw["all_five_frames"]
    )
    assert raw["all_120_colored_markings_give_exact_root_vectors"] is True
    assert raw["aggregate_all_five_full_wall_coordinate_sha256"] == (
        "74be27cfc3a7e6d21264b3e85d08b2d1c956e9cfcbaf12fdb1e9f212453d0c49"
    )


def test_F5_conic_factorizations_and_forced_curve_systems() -> None:
    payload = certificate()
    geometry = payload["finite_geometry"]
    assert geometry["ordinary_points"] == [
        {"point": "o0", "coordinates": [0, 0, 1]},
        {"point": "o1", "coordinates": [1, 1, 1]},
        {"point": "o2", "coordinates": [1, 3, 4]},
        {"point": "o3", "coordinates": [1, 2, 4]},
        {"point": "o4", "coordinates": [1, 4, 1]},
        {"point": "o5", "coordinates": [1, 0, 0]},
    ]
    assert geometry["special_points"] == [
        {"point": "s0", "coordinates": [1, 1, 4]},
        {"point": "s1", "coordinates": [1, 3, 1]},
        {"point": "s2", "coordinates": [1, 2, 2]},
        {"point": "s3", "coordinates": [1, 0, 2]},
    ]
    assert geometry["K6_perfect_matching_count"] == 15
    assert geometry["concurrent_chord_matching_count"] == 10
    assert geometry["K6_one_factorization_count"] == 6
    assert geometry["one_factorization_concurrency_histogram"] == {"0": 1, "4": 5}
    assert geometry["usable_tetrahedral_factorization_count"] == 5
    assert (
        geometry["triangle_line_count"],
        geometry["conic_count"],
        geometry["quartic_count"],
    ) == (12, 4, 4)
    assert all(conic["symmetric_matrix_determinant_mod_5"] for conic in geometry["conics"])
    assert all(quartic["condition_rank"] == 14 for quartic in geometry["quartics"])

    systems = payload["linear_systems"]
    assert systems["anticanonical_cubics"] == {
        "ambient_vector_dimension": 10,
        "evaluation_rank": 10,
        "kernel_dimension": 0,
    }
    sextic = systems["antibicanonical_sextics"]
    assert sextic["ambient_vector_dimension"] == 28
    assert sextic["gradient_condition_row_count"] == 30
    assert sextic["gradient_condition_rank"] == 27
    assert sextic["kernel_dimension"] == 1
    assert sextic["value_augmented_rank"] == 27
    assert sextic["coefficients"] == [
        0,
        0,
        0,
        1,
        1,
        4,
        4,
        0,
        2,
        1,
        1,
        4,
        0,
        1,
        1,
        4,
        3,
        2,
        0,
        2,
        0,
        4,
        2,
        4,
        3,
        1,
        0,
        0,
    ]


def test_quartics_and_boundary_have_exact_singular_schemes() -> None:
    payload = certificate()
    audit = payload["exact_factor_and_singular_scheme_audit"]
    for quartic in audit["quartics"]:
        assert quartic["irreducible_over_F5"] is True
        assert quartic["geometrically_integral"] is True
        assert quartic["singular_degree"] == 3
        assert quartic["radical_degree"] == 3
        assert quartic["minimal_primes"] == 3
        assert quartic["j_mod_r"] == quartic["r_mod_j"] == 0
        assert quartic["jacobian_mod_points"] == quartic["points_mod_jacobian"] == 0

    boundary = audit["boundary"]
    assert boundary["irreducible_over_F5"] is True
    assert boundary["geometrically_integral"] is True
    assert boundary["singular_degree"] == 16
    assert boundary["radical_degree"] == 10
    assert boundary["minimal_primes"] == 10
    assert boundary["j_mod_r"] == 0
    assert boundary["r_mod_j"] != 0
    assert boundary["radical_mod_points"] == boundary["points_mod_radical"] == 0
    assert boundary["jacobian_scheme_is_nonreduced"] is True
    assert "does not have ten nodes" in payload["hostile_warning"]


def test_local_A1_A2_types_transversality_and_strong_Coble() -> None:
    payload = certificate()
    local = payload["local_geometry_and_transversality"]
    assert local["boundary_node_count"] == 4
    assert local["boundary_cusp_count"] == 6
    assert local["boundary_node_discriminant_histogram"] == {"2": 2, "3": 2}
    assert local["boundary_cusp_tangent_cubic_coefficient_histogram"] == {
        "1": 1,
        "2": 3,
        "3": 1,
        "4": 1,
    }
    assert local["boundary_delta_sum"] == 10
    assert local["quartic_node_count"] == 12
    assert local["quartic_node_discriminant_histogram"] == {"2": 6, "3": 6}
    assert local["smooth_root_center_incidence_count"] == 84
    assert local["smooth_root_boundary_tangent_value_histogram_by_kind"] == {
        "line": {"1": 18, "4": 18},
        "conic": {"1": 12, "4": 12},
        "quartic": {"1": 12, "4": 12},
    }
    assert local["quartic_node_boundary_node_resultant_count"] == 12
    assert local["quartic_node_boundary_node_resultant_histogram"] == {
        "1": 6,
        "4": 6,
    }
    assert local["all_root_branches_avoid_all_boundary_branches"] is True

    strong = payload["blowup_and_strong_Coble"]
    assert strong["minus_K_empty"] is True
    assert strong["minus_2K_unique"] is True
    assert strong["boundary_strict_transform_smooth_rational"] is True
    assert strong["all_twenty_root_strict_transforms_smooth_geometrically_integral_rational"] is True
    assert strong["all_twenty_roots_disjoint_from_boundary"] is True
    assert strong["strong_Coble_verdict"] == "PASS"


def test_direct_root_cone_cusps_and_same_H10_72_return() -> None:
    payload = certificate()
    raw = payload["marked_same_polytope_return"]
    root_cone = raw["direct_twenty_root_cone"]
    assert (
        root_cone["finite_ray_count"],
        root_cone["ideal_ray_count"],
        root_cone["hyperideal_ray_count"],
    ) == (450, 51, 0)
    assert root_cone["ray_count"] == 501
    assert root_cone["primitive_ray_list_sha256"] == (
        "80e1e81adc197bb92e818f52f3e695aa1d4d61fc2f4aecf6a550ea1b99fbb5ec"
    )
    assert root_cone["primitive_ideal_ray_list_sha256"] == (
        "8621adee95e09a0213356e044f923a0fd902e22aa55ead6ffc518959fbfc73d7"
    )
    assert root_cone["minimum_root_inequality_slack"] > 0
    assert root_cone["minimum_time_sheet_pairing"] > 0
    assert len(root_cone["cusp_records"]) == 51
    assert root_cone["all_51_ideal_rays_equal_the_51_primitive_E_plus_b_classes"] is True

    full = raw["upstream_full_P_vinberg_return"]
    assert full == {
        "polytope_id": "COBLE-VII-p5-b1-H10-72",
        "dimension": 10,
        "facet_count": 72,
        "finite_vertices": 18761,
        "ideal_vertices": 1731,
        "hyperideal_extreme_rays": 0,
        "elliptic_edges": 143565,
        "finite_volume": True,
        "compact": False,
        "all_facets_irredundant": True,
        "literal_vinberg_zero_failures": True,
    }
    assert raw["general_one_boundary_return"]["same_full_wall_set"] is True
    assert raw["exact_same_marked_H10_72"] is True


def test_scope_and_frozen_archive_replay() -> None:
    payload = certificate()
    assert payload["status"] == (
        "PROVED_EXPLICIT_P_FIRST_STRONG_COBLE_AND_SAME_MARKED_H10_72"
    )
    assert payload["inverse_level"]["rating"] == (
        "L3 structural P-first calibration, not class-wide L4"
    )
    assert payload["sources"]["raw_H10_72_archive"][
        "surface_existence_or_classification_imported"
    ] is False
    assert payload["lineage"]["new_abstract_hyperbolic_polytope"] is False
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == json.loads(
        json.dumps(payload, sort_keys=True)
    )
