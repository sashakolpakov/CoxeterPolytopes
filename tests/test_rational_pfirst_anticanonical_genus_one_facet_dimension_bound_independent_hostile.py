"""Independent hostile tests for the anticanonical genus-one facet theorem."""

from __future__ import annotations

from scripts.audit_rational_pfirst_anticanonical_genus_one_facet_dimension_bound_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    endpoint_certificate,
    intrinsic_theorem_certificate,
    sha256,
)


def test_frozen_owner_and_reservoir_bytes() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_intrinsic_unimodular_reflection_argument() -> None:
    record = intrinsic_theorem_certificate()
    assert record["K_primitive"] is True
    assert record["primitive_in_unimodular_implies_divisibility"] == 1
    assert record["integral_reflection_forces"] == "2/F^2 is an integer"
    assert record["allowed_F_squares"] == [-1, -2]
    assert record["allowed_rho"] == [11, 12]
    assert record["allowed_projective_nef_dimensions"] == [10, 11]
    assert record["characteristic_free"] is True


def test_genuine_hyperbolic_facet_supplies_negativity() -> None:
    record = intrinsic_theorem_certificate()
    assert record["genuine_hyperbolic_facet_forces_negative_square"] is True


def test_p10_endpoint_is_exact_and_full_nef() -> None:
    record = endpoint_certificate(10)
    assert (record["F_square"], record["rho"], record["n"]) == (-1, 11, 10)
    assert (record["finite_vertices"], record["ideal_vertices"], record["hyperideal_vertices"]) == (17, 2, 0)
    assert record["F_primitive"] and record["F_reflection_integral"]
    assert record["all_facets_genuine"] and record["complete_nef_source_replayed"]


def test_p11_is_the_sharp_endpoint() -> None:
    record = endpoint_certificate(11)
    assert (record["F_square"], record["rho"], record["n"]) == (-2, 12, 11)
    assert (record["finite_vertices"], record["ideal_vertices"], record["hyperideal_vertices"]) == (26, 2, 0)
    assert record["arithmetic_genus"] == 1
    assert record["all_facets_genuine"] and record["complete_nef_source_replayed"]


def test_combined_scope() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_EXACT_SCOPE"
    assert all(audit["scope"].values())
