"""Independent hostile tests for the compact right-hexagon field theorem."""

from __future__ import annotations

from fractions import Fraction

from scripts.audit_k3_n2_compact_right_hexagon_field_bound_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    closure_record,
    extension_record,
    field_record,
    independent_coordinates,
    sha256_file,
    short_pair_record,
)


def test_all_owner_and_local_scope_bytes_are_pinned() -> None:
    for item in PINS.values():
        assert sha256_file(ROOT / item["path"]) == item["sha256"]


def test_cofactor_coordinates_are_independent_exact_replay() -> None:
    record = closure_record()
    assert record == {
        "exact_rational_samples": 125,
        "common_numerator_identity": True,
        "six_rotated_four_consecutive_relations": True,
        "distinguished_square_roots_only": True,
    }
    coordinates = independent_coordinates(Fraction(2), Fraction(3), Fraction(4))
    assert len(coordinates) == 6
    assert min(coordinates) > 0


def test_short_pair_dichotomy_keeps_the_exact_constants() -> None:
    record = short_pair_record()
    assert record["alternating_order_permutations_checked"] == 6
    assert record["small_alternating_branch_threshold"] == 16
    assert record["large_alternating_branch_odd_upper"] == "(51+34sqrt(17))/256"
    assert record["rational_upper"] == "221/256"
    assert record["universal_some_shifted_coordinate_upper"] == 16
    assert record["universal_some_pair_cycle_upper"] == 68
    assert record["does_not_bound_every_coordinate"] is True


def test_unique_large_conjugate_rejects_every_proper_degree_pattern() -> None:
    record = extension_record()
    assert record["unique_value_above_four_forces_full_field"] is True
    assert record["totally_real_extension_used"] is True
    assert all(
        row["real_extension_count"] >= 2
        for row in record["proper_subfield_patterns"]
    )


def test_K3_degree_and_vandermonde_bound() -> None:
    record = field_record()
    assert record["field_degree_cap"] == 6
    assert record["discriminant_exponent"] == 30
    assert record["Vandermonde_difference_bound"] == 68
    assert record["field_discriminant_bound"] == "68^30"


def test_scope_does_not_expand_to_general_hexagons_or_carrier_existence() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_EXACT_COMPACT_RIGHT_ANGLED_SCOPE"
    scope = audit["scope"]
    assert scope["primitive_exact_full_trace_complex_K3_required"] is True
    assert scope["exact_compact_right_angled_hexagon_required"] is True
    assert scope["carrier_construction_proved"] is False
    assert scope["essential_surjectivity_proved"] is False
    assert scope["nonright_or_noncompact_hexagon_closed"] is False
    assert scope["arbitrary_six_or_more_facet_chamber_closed"] is False
    assert scope["every_dotted_pair_bounded"] is False
    assert scope["mixed_height_bounded"] is False
    assert scope["maximality_used"] is False
