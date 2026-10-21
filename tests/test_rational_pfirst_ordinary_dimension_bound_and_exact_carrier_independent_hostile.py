"""Independent hostile tests for the ordinary P-first combined theorem."""

from __future__ import annotations

from scripts.audit_rational_pfirst_ordinary_dimension_bound_and_exact_carrier_independent_hostile import (
    OWNER_PINS,
    ROOT,
    build_audit,
    characteristic_gap_records,
    dual_rays,
    escape_record,
    flex_record,
    rank_ten_record,
    sha256_file,
    simplex_record,
)


def test_owner_bytes_are_the_refrozen_automatic_escape_package() -> None:
    for item in OWNER_PINS.values():
        assert sha256_file(ROOT / item["path"]) == item["sha256"]


def test_standard_dual_rays_and_sharp_branch_split() -> None:
    records = [simplex_record(n) for n in range(3, 10)]
    assert [record["anticanonical_square"] for record in records] == [
        6,
        5,
        4,
        3,
        2,
        1,
        0,
    ]
    assert all(record["facet_basis_determinant"] in (-1, 1) for record in records)
    assert records[-1]["null_dual_indices"] == [0, 1]
    assert tuple(dual_rays(9)[0]) == (3,) + (-1,) * 9
    assert tuple(dual_rays(9)[1]) == (1, -1) + (0,) * 8


def test_spacelike_escape_crossing_remains_strict() -> None:
    for n in range(10, 16):
        record = escape_record(n)
        assert record["all_strict_wall_slacks"] == [1] * (n + 1)
        assert record["F_wall_slacks"] == [1] + [0] * n
        assert record["escape_quadratic"][2] == 9 - n
        assert record["root_product_is_negative"] is True


def test_reduced_characteristic_gap_is_symbolically_strict() -> None:
    records = characteristic_gap_records()
    assert [record["n"] for record in records] == list(range(3, 9))
    assert all(
        record["all_nonzero_nonnegative_inputs_give_strict_gap"]
        for record in records
    )
    assert [record["unique_reduced_characteristic_vector"][0] for record in records] == [
        3
    ] * 6


def test_ninefold_cusp_relation_and_characteristic_split() -> None:
    record = rank_ten_record()
    assert record["primitive_affine_multiplicities"] == [
        2,
        4,
        6,
        5,
        4,
        3,
        2,
        1,
        3,
    ]
    assert record["reduced_support_euler_number"] == 10
    assert record["F_is_characteristic"] is True
    assert record["other_ideal_ray_is_characteristic"] is False
    assert abs(record["trivial_lattice_determinant"]) == 1


def test_flex_local_length_and_partial_residuals_are_exact() -> None:
    record = flex_record()
    assert record["base_scheme_colength"] == 9
    assert record["contact_with_flex_line"] == 3
    assert record["residual_exponents_after_n_blowups"] == {
        3: 6,
        4: 5,
        5: 4,
        6: 3,
        7: 2,
        8: 1,
        9: 0,
    }
    assert record["finite_discriminant_has_two_simple_roots_over_C"] is True


def test_full_independent_verdict_and_nonclaims() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS"
    assert audit["automatic_canonical_escape"][
        "extra_K_not_in_P_hypothesis_required"
    ] is False
    assert audit["scope"]["dimension_bound_uses_no_surface_or_maximality"] is True
    assert audit["scope"]["exact_carrier_uses_relative_maximality"] is True
    assert audit["scope"]["twice_primitive_minus_four_lane_included"] is False
