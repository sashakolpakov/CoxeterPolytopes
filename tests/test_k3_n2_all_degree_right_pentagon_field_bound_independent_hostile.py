"""Independent hostile tests for the all-degree right-pentagon theorem."""

from __future__ import annotations

from scripts.audit_k3_n2_all_degree_right_pentagon_field_bound_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    cyclic_field_packet_certificate,
    gram_recurrence_certificate,
    norm_and_upper_bound_certificate,
    quadratic_certificate,
    sha256,
)


def test_frozen_owner_and_local_scope_bytes() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_all_five_Gram_relations_are_the_A2_recurrence() -> None:
    record = gram_recurrence_certificate()
    assert record["tree_determinant"] == "1-A_left-A_mid-A_right+A_left*A_right"
    assert record["cyclic_recurrence"] == "x[i-1]*x[i+1]=x[i]+1"
    assert record["all_five_relations_checked"] is True


def test_same_degree_dependent_epsilon_is_used_twice() -> None:
    record = norm_and_upper_bound_certificate()
    assert record["degree_bound_from_3d_le_20"] == 6
    assert record["same_epsilon_used_twice"] is True
    terminal = record["rows"][-1]
    assert terminal["floor_used_for_x_i_plus_3"] == "1/4096"
    assert terminal["floor_used_for_z"] == "1/4096"
    assert terminal["strict_shift_upper"] == 16_785_408
    assert terminal["strict_A_upper"] == 16_785_409
    assert record["invalid_mixed_floor_value_not_used"] == 69_648


def test_degree_gt_two_uses_complete_packet_not_pair_generation() -> None:
    record = cyclic_field_packet_certificate()
    assert record["nonzero_Cartan_graph"] == "C5"
    assert record["closed_walks_factor_into_simple_cycles"] is True
    assert record["individual_directed_entries_assumed_integral"] is False
    assert record["complete_packet_generates_k_is_hypothesis"] is True
    assert record["single_pair_generates_for_degree_gt_2"] is False
    assert record["uniform_packet_height"] == 36_938_546_139_294_728_224
    assert record["Vandermonde_difference_height"] == 1_182_033_476_457_431_303_168
    assert record["Vandermonde_exponent"] == 30
    assert record["finite_field_set"] is True


def test_quadratic_specialization_retains_sharp_bound() -> None:
    record = quadratic_certificate()
    assert record["shift_interval"] == "(1/16,288)"
    assert record["disjoint_intervals_force_quadratic_generation"] is True
    assert record["field_discriminant_upper_bound"] == 1_343_277


def test_combined_audit_scope() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_CONDITIONAL_EXACT_SCOPE"
    assert all(audit["scope"].values())
