"""Independent hostile tests for the finite-area triangle field theorem."""

from __future__ import annotations

from scripts.audit_k3_n2_finite_area_triangle_field_bound_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    closed_walk_certificate,
    height_and_field_certificate,
    sha256,
    trace_integrality_certificate,
)


def test_owner_v2_bytes_are_locked() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_trace_integrality_does_not_need_primitive_normals() -> None:
    record = trace_integrality_certificate()
    assert record["pair_cycle_formula"] == "Tr(T_i*T_j)=p_ij"
    assert record["triple_cycle_formula"] == "Tr(T_0*T_1*T_2)=c"
    assert record["individual_directed_Cartan_entries_assumed_integral"] is False
    assert record["primitive_facet_normals_required"] is False
    assert record["cyclic_products_invariant_under_normal_rescaling"] is True
    assert record["finite_projective_module_trace_lies_in_order"] is True


def test_closed_walks_decompose_beyond_naive_alternation() -> None:
    record = closed_walk_certificate()
    assert record["examples"]["0-1-0-2-0"] == ["p01", "p02"]
    assert record["examples"]["0-1-2-0"] == ["c"]
    assert record["all_closed_walk_coordinates_in_Z[p01,p12,p02,c]"] is True


def test_ideal_endpoint_uses_weak_heights_and_same_discriminant() -> None:
    record = height_and_field_certificate()
    assert record["full_trace_degree_bound"] == 6
    assert record["finite_or_ideal_pair_cycle_weak_height"] == 4
    assert record["triple_cycle_weak_height"] == 8
    assert record["primitive_generator_weak_height"] == 128
    assert record["conjugate_difference_weak_height"] == 256
    assert record["Vandermonde_exponent"] == 30
    assert record["field_discriminant_bound"] == "256^30"
    assert record["Hermite_Minkowski_gives_finite_field_set"] is True


def test_combined_scope_is_finite_area_not_only_compact() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_EXACT_SCOPE"
    assert all(audit["scope"].values())
