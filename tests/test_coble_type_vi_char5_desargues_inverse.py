from __future__ import annotations

import json

from scripts.certify_coble_type_vi_char5_desargues_inverse import (
    OUTPUT,
    build_certificate,
)


def certificate() -> dict[str, object]:
    return build_certificate()


def test_desargues_and_cut_designs_are_exact() -> None:
    payload = certificate()
    assert payload["point_count"] == 10
    assert payload["triangle_line_count"] == 10
    assert payload["triangle_replication"] == [3]
    assert payload["cut_conic_count"] == 10
    assert payload["cut_replication"] == [6]


def test_anticanonical_and_antibicanonical_systems() -> None:
    payload = certificate()
    assert payload["anticanonical_system"] == {
        "ambient_dimension": 10,
        "evaluation_rank": 10,
        "kernel_dimension": 0,
        "empty": True,
    }
    sextic = payload["sextic_double_point_system"]
    assert sextic["ambient_dimension"] == 28
    assert sextic["condition_rank"] == 27
    assert sextic["kernel_dimension"] == 1
    assert sextic["value_conditions_follow"] is True


def test_sextic_is_geometrically_integral_with_exact_nodes() -> None:
    payload = certificate()
    singular = payload["singular_and_factor_audit"]
    assert singular["irreducible_over_F5"] is True
    assert singular["geometrically_integral"] is True
    assert singular["projective_singular_scheme_dimension"] == 0
    assert singular["projective_singular_degree"] == 10
    assert singular["singular_scheme_is_exactly_the_ten_reduced_F5_points"] is True
    local = payload["sextic_local_audit"]
    assert local["ordinary_node_count"] == 10
    assert local["root_center_incidence_count"] == 90
    assert local["all_root_directions_avoid_both_boundary_branches"] is True


def test_strong_coble_and_same_polytope_return() -> None:
    payload = certificate()
    assert payload["strong_Coble"]["verdict"] == "PASS"
    handoff = payload["raw_polytope_handoff"]
    assert handoff["unique_simultaneous_P2_frame"] is True
    assert handoff["constructed_root_profile"] == handoff["raw_root_profile"]
    marking = handoff["exact_raw_marking"]
    assert marking["basis_determinant_absolute"] == 1
    assert marking["exceptional_frame_wall_indices"] == (
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
    )
    assert marking["colored_incidence_isomorphism_count"] == 120
    assert marking["all_root_vectors_equal_marked_plane_classes"] is True
    assert marking["exact_marked_root_seed_transport"] is True
    assert handoff["full_cone_rays"] == {
        "finite": 21661,
        "ideal": 1797,
        "hyperideal": 0,
    }
    assert handoff["exact_same_P"] is True


def test_inverse_scope_is_not_overclaimed() -> None:
    payload = certificate()
    assert payload["inverse_level"]["rating"] == "L3 structural P-first calibration, not L4"
    assert payload["lineage"]["worldwide_priority"] == "open; no new abstract polytope is claimed"


def test_frozen_archive_matches_live_certificate() -> None:
    live_json = json.loads(json.dumps(certificate(), sort_keys=True))
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == live_json
