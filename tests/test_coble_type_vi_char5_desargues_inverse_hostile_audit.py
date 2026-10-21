from __future__ import annotations

from functools import lru_cache
import json

from scripts.audit_coble_type_vi_char5_desargues_inverse_hostile import (
    OUTPUT,
    build_audit,
)


@lru_cache(maxsize=1)
def payload() -> dict[str, object]:
    return build_audit()


def test_independent_quotient_design_and_linear_systems() -> None:
    audit = payload()
    quotient = audit["quotient_model"]
    assert quotient["dimension"] == 3
    assert len(quotient["edge_points"]) == 10
    assert len({tuple(row["point"]) for row in quotient["edge_points"]}) == 10

    design = audit["finite_design"]
    assert design["line_count"] == 10
    assert design["conic_count"] == 10
    assert design["line_replication"] == [3]
    assert design["conic_replication"] == [6]
    assert design["all_supports_exact"] is True
    assert design["all_conics_smooth_geometrically_integral"] is True
    assert 0 not in {
        int(value) for value in design["conic_determinant_histogram"]
    }

    systems = audit["linear_systems"]
    assert systems["anticanonical_cubics"] == {
        "ambient_vector_dimension": 10,
        "evaluation_rank": 10,
        "kernel_dimension": 0,
        "empty_projective_system": True,
    }
    assert systems["antibicanonical_sextics"]["ambient_vector_dimension"] == 28
    assert systems["antibicanonical_sextics"]["gradient_condition_rank"] == 27
    assert systems["antibicanonical_sextics"]["gradient_kernel_dimension"] == 1
    assert systems["antibicanonical_sextics"]["value_augmented_rank"] == 27


def test_exact_singular_scheme_nodes_and_geometric_integrality() -> None:
    audit = payload()
    scheme = audit["singular_scheme_and_integrality"]
    assert scheme["irreducible_over_F5"] is True
    assert scheme["geometrically_integral"] is True
    assert scheme["projective_singular_scheme_degree"] == 10
    assert scheme["radical_degree"] == 10
    assert scheme["minimal_prime_count"] == 10
    assert scheme["saturated_jacobian_equals_its_radical"] is True
    assert scheme["saturated_jacobian_equals_explicit_ten_point_ideal"] is True

    local = audit["local_geometry"]
    assert local["ordinary_node_count"] == 10
    assert local["discriminant_histogram"] == {"2": 9, "3": 1}
    assert local["root_center_incidence_count"] == 90
    assert local["all_root_directions_avoid_both_boundary_branches"] is True
    assert local["tangent_value_histogram_by_root_kind"] == {
        "line": {"1": 24, "4": 6},
        "conic": {"2": 18, "3": 42},
    }


def test_strong_coble_package_and_boundary_disjointness() -> None:
    audit = payload()
    package = audit["blowup_and_strong_Coble"]
    assert package["minus_K_empty"] is True
    assert package["minus_2K_has_unique_member"] is True
    assert package["boundary_smooth_after_one_blowup_at_each_node"] is True
    assert package["boundary_geometric_genus"] == 0
    assert package["twenty_roots_smooth_geometrically_integral_rational"] is True
    assert package["twenty_roots_disjoint_from_boundary_after_blowup"] is True
    assert package["strong_Coble_verdict"] == "PASS"


def test_exact_marked_same_H10_78_bridge_repairs_profile_only_gate() -> None:
    audit = payload()
    bridge = audit["marked_raw_polytope_handoff"]
    assert bridge["primitive_wall_norm_histogram"] == {"-2": 20, "-1": 58}
    assert bridge["intrinsic_boundary_wall_index_zero_based"] == 0
    assert bridge["maximal_disjoint_clique_size_histogram"] == {
        "7": 20,
        "8": 60,
        "9": 55,
        "10": 1,
    }
    assert bridge["unique_ten_curve_frame_wall_indices_zero_based"] == [
        25,
        26,
        29,
        35,
        40,
        43,
        49,
        57,
        64,
        77,
    ]
    assert bridge["marked_basis_determinant"] == 1
    assert bridge["marked_basis_gram"] == [
        [1, *([0] * 10)],
        *[
            [0, *([0] * index), -1, *([0] * (9 - index))]
            for index in range(10)
        ],
    ]
    assert bridge["boundary_coordinates_in_H_E_frame"] == [3, *([-1] * 10)]
    assert bridge["root_plane_profile"] == {
        "degree_1_three_centers": 10,
        "degree_2_six_centers": 10,
    }
    assert bridge["colored_incidence_marking_count"] == 120
    assert bridge["marking_torsor_sha256"] == (
        "d608b7795e0f72cee97fc598ccd067c42135059168f016c4a4c6443eb58ceb6f"
    )
    assert bridge["all_twenty_root_coordinates_match_triangle_cut_classes"] is True
    assert bridge["transported_full_wall_coordinate_count"] == 78
    assert bridge["primitive_isotropic_E_plus_b_count"] == 57
    assert bridge["exact_same_marked_H10_78"] is True
    assert bridge["upstream_P_side_import"]["finite_volume"] is True
    assert bridge["upstream_P_side_import"]["compact"] is False
    assert bridge["upstream_P_side_import"]["hyperideal_extreme_rays"] == 0

    target = audit["target_certificate_audit"]
    assert target["independent_coefficients_and_points_agree"] is True
    assert target["target_same_P_gate_had_explicit_integral_marking"] is False


def test_scope_and_frozen_hostile_archive() -> None:
    audit = payload()
    assert audit["verdict"] == "GREEN_WITH_MARKED_HANDOFF_STRENGTHENING"
    assert audit["lineage_and_scope"]["inverse_level"] == (
        "L3 datum calibration, not a class-wide L4 carrier theorem"
    )
    assert len(audit["lineage_and_scope"]["not_proved"]) == 4
    assert audit["all_checks_pass"] is True
    assert all(audit["checks"].values())
    frozen = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert json.loads(json.dumps(audit, sort_keys=True)) == frozen
