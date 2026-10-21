from __future__ import annotations

import json
from pathlib import Path

from scripts.audit_rational_full_nef_maximal_finiteness_hostile import (
    OUTPUT,
    _binary_subspace_table,
    _rref_subspaces,
    _totally_isotropic,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
PRODUCER = ROOT / "data/enumerations/rational_full_nef_maximal_finiteness.json"


def certificate() -> dict:
    return build_certificate()


def test_W_equals_R_M_under_exactly_the_three_valid_maximality_modes() -> None:
    saturation = certificate()["saturation_logic"]
    counts = saturation["missing_reflection_countermodel_counts"]
    assert counts == {
        "M_relative_reflection_maximality": 0,
        "arithmetic_reflection_maximality_with_finite_volume": 0,
        "maximal_integral_lattice_only_hostile_control": 4,
        "strong_discrete_reflection_maximality": 0,
    }
    assert saturation["three_valid_maximality_modes_force_W_equals_R_M"]
    assert saturation["maximal_integral_lattice_does_not_force_saturation"]
    hostile = saturation["insufficient_maximal_integral_lattice_countermodel"]
    assert hostile["missing_integral_reflection"]
    assert hostile["maximal_integral_lattice_only"]


def test_complete_reflection_arrangement_has_one_chamber_orbit() -> None:
    row = certificate()["one_chamber_orbit_audit"]
    assert row["chamber_count"] == 6
    assert row["complete_reflection_group_orbit_of_chamber_zero"] == list(range(6))
    assert row["complete_arrangement_has_one_chamber_orbit"]
    assert row["incomplete_reflection_set_need_not_be_transitive"]
    assert len(row["general_proof_dependencies"]) == 3


def test_direct_host_list_and_U_split_are_independent() -> None:
    row = certificate()["direct_host_audit"]
    assert row["positive_dimensional_host_count"] == 10
    assert row["odd_forms_ranks_3_through_10"] == list(range(3, 11))
    assert row["U_plus_minus_one_change_determinant"] == -1
    assert row["U_plus_minus_one_conjugated_gram"] == [
        [1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
    ]


def test_binary_host_bounds_use_recurrence_and_independent_RREF() -> None:
    cert = certificate()
    producer = json.loads(PRODUCER.read_text(encoding="utf-8"))
    table = _binary_subspace_table(22)
    counts = [sum(table[rank]) for rank in range(23)]
    assert counts[1:6] == [2, 5, 16, 67, 374]
    for rank in range(1, 6):
        subspaces = _rref_subspaces(rank)
        assert len(subspaces) == counts[rank]
        assert len(set(subspaces)) == len(subspaces)
    small = cert["standard_Coble_host_audit"]["small_rank_independent_RREF_counts"]
    assert [row["totally_isotropic_subspaces_for_I_mod_2"] for row in small] == [
        1,
        2,
        4,
        11,
        31,
    ]
    assert all(
        _totally_isotropic(rows)
        for rows in _rref_subspaces(4)
        if _totally_isotropic(rows)
    )
    assert cert["standard_Coble_host_audit"]["complex_rank_11_through_20_bound"] == producer[
        "coble_extension_bounds"
    ]["complex_marked_extension_upper_bound"]
    assert cert["standard_Coble_host_audit"]["odd_positive_rank_11_through_22_bound"] == producer[
        "coble_extension_bounds"
    ]["odd_positive_marked_extension_upper_bound"]


def test_F4_maximality_gap_is_repaired_by_sharp_noncompact_area() -> None:
    row = certificate()["F4_reflection_maximality_repair"]
    assert row["root_norms"] == [-2, -2, -2]
    assert sorted(map(str, row["coxeter_orders"])) == ["2", "3", "infinity"]
    assert row["triangle_area_over_pi"] == "1/6"
    assert row["minimum_positive_one_ideal_triangle_area_over_pi"] == "1/6"
    assert row["minimum_signature_in_search"] == [2, 3, "infinity"]
    assert row["all_three_reflections_integral_isometric_involutions"]
    assert row["producer_archive_maximality_certificate_was_missing"]
    assert row["reflection_maximality_proved"]


def test_nonmaximal_fixed_lattice_control_keeps_maximality_essential() -> None:
    row = certificate()["nonmaximal_fixed_lattice_control"]
    assert row["unbounded_facets_at_fixed_I_1_2_without_maximality"]
    assert row["all_exact_ideal_polygon_checks_pass"]
    assert [item["area_over_pi_from_all_ideal_polygon"] for item in row["rows"]] == list(
        range(1, 25)
    )
    assert [item["proper_index_over_base_for_index_above_one"] for item in row["rows"]][
        1:
    ] == list(range(2, 25))


def test_field_Q_and_effectiveness_scopes_are_literal() -> None:
    cert = certificate()
    field = cert["vinberg_field_audit"]
    assert field["rational_normals_force_Q_cyclic_products"]
    assert field["irrational_normal_hostile_cyclic_product_Q_plus_sqrt2_coefficients"] == [
        "8/3",
        "4/3",
    ]
    assert field["irrational_normal_control_is_not_rational"]
    assert field["trace_module_degree_inequality_is_optional_not_used"]
    assert field["carrier_surface_fields_of_definition_are_not_bounded_by_this_argument"]
    scope = cert["effectiveness_scope"]
    assert scope["finite_host_list_is_a_terminating_mathematical_enumeration"]
    assert scope["finite_abstract_chamber_bound_proved"]
    assert not scope["uniform_reflectivity_or_facet_stopping_bound_proved"]
    assert not scope["terminating_census_of_reflective_hosts_proved"]
    assert not scope["finite_surface_isomorphism_or_definition_field_list_proved"]
    assert not scope["raw_P_to_carrier_essential_surjectivity_proved"]


def test_verdict_is_pass_with_certificate_repair_not_scope_narrowing() -> None:
    cert = certificate()
    assert cert["status"] == "PASS_WITH_ONE_INDEPENDENTLY_REPAIRED_CERTIFICATE_GAP"
    assert cert["verdict"]["core_finiteness_theorem"].startswith("sound")
    assert cert["verdict"]["scope_narrowing_required"] is False
    assert all(cert["claims_present"].values())
    assert all(cert["checks"].values())
    assert cert["all_checks_pass"]


def test_frozen_hostile_archive_is_exact_rebuild() -> None:
    assert json.loads(Path(OUTPUT).read_text(encoding="utf-8")) == certificate()
