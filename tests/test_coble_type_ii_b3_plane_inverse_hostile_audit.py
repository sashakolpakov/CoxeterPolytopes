from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.audit_coble_type_ii_b3_plane_inverse_hostile import (
    OUTPUT,
    build_audit,
)


@pytest.fixture(scope="module")
def live() -> dict:
    return build_audit()


def test_frozen_archive_is_exact_live_replay(live: dict) -> None:
    frozen = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert live == frozen
    assert live["all_checks_pass"] is True
    assert live["verdict"] == "GREEN_WITH_CORRECTED_CHARACTERISTIC_3_BOUNDARY"


def test_raw_573_denominator_and_two_branches(live: dict) -> None:
    raw = live["raw_endpoint_reconstruction"]
    assert raw["state_count_by_contraction_depth_0_through_10"] == [
        1, 17, 87, 329, 948, 2094, 3472, 4417, 4455, 3291, 573
    ]
    assert raw["distinct_degree_signature_count"] == 28
    assert [branch["contracted_wall_set_zero_based"] for branch in raw["branches"]] == [
        [1, 5, 9, 13, 15, 17, 19, 22, 26, 28],
        [3, 7, 11, 13, 15, 17, 19, 22, 26, 28],
    ]
    assert [branch["basis_determinant"] for branch in raw["branches"]] == [1, -1]
    assert all(branch["colored_B3_isomorphism_count"] == 12 for branch in raw["branches"])
    assert all(branch["all_thirty_reconstruct_exactly"] for branch in raw["branches"])


def test_surface_geometry_over_Q_and_C(live: dict) -> None:
    systems = live["B3_incidence_and_linear_systems"]["linear_systems_over_Q"]
    absolute = live["characteristic_zero_absolute_and_singular_audit"]
    local = live["characteristic_zero_local_resolution"]
    roots = live["root_curve_geometry"]
    assert (systems["h0_minus_K"], systems["h0_minus_2K"]) == (0, 1)
    assert absolute["absolute_component_count"] == 1
    assert absolute["singular_scheme_degree"] == 13
    assert absolute["radical_degree"] == 7
    assert (local["A3_count"], local["A1_count"], local["delta_sum"]) == (3, 4, 10)
    assert roots["pair_intersection_histogram"] == {"0": 51, "1": 15}
    assert roots["boundary_intersection_with_every_root"] == 0


def test_complete_cone_and_same_marked_P(live: dict) -> None:
    cone = live["direct_cone_and_same_P_return"]["direct_raw_cdd"]
    returned = live["direct_cone_and_same_P_return"]["same_marked_P_return"]
    assert (cone["extreme_rays"], cone["finite"], cone["ideal"], cone["hyperideal"]) == (1540, 1333, 207, 0)
    assert cone["primitive_ray_digest"] == "f88e2b4397b6e6c27cee784cd1ecad4f534dd39decf219774e3f3f12b83e94f9"
    assert cone["all_thirty_facets_irredundant"] is True
    assert returned["wall_inventory"] == {"boundary": 1, "roots": 12, "cusp_lift_minus_one": 17, "total": 30}
    assert returned["exact_same_marked_H10_30"] is True


def test_characteristic_three_rejection_is_narrow_and_correct(live: dict) -> None:
    char3 = live["characteristic_three_scope_correction"]
    assert char3["affine_cone_dimension"] == 1
    assert char3["radical_affine_cone_dimension"] == 1
    assert char3["singular_scheme_degree"] == 21
    assert char3["radical_degree"] == 7
    assert len(char3["repeated_tangent_points"]) == 4
    assert "not positive-dimensional" in char3["explicit_nonclaim"]


def test_hostile_script_does_not_import_target_builder() -> None:
    source = (Path(__file__).parents[1] / "scripts/audit_coble_type_ii_b3_plane_inverse_hostile.py").read_text(encoding="utf-8")
    assert "from scripts.certify_coble_type_ii_b3_plane_inverse import" not in source
    assert "import scripts.certify_coble_type_ii_b3_plane_inverse" not in source

