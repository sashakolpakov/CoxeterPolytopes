"""Independent hostile tests for the enriched AG/ABSW transfer."""

from __future__ import annotations

from scripts.audit_enriched_canonical_chamber_ag_absw_transfer_independent_hostile import (
    OWNER_PINS,
    ROOT,
    build_audit,
    dominant_height_certificate,
    field_constant_record,
    primitive_grid_certificate,
    sha256_file,
    uniform_constant_certificate,
    verify_owner_archive,
)


def test_locked_owner_bytes_match() -> None:
    for item in OWNER_PINS.values():
        assert sha256_file(ROOT / item["path"]) == item["sha256"]


def test_all_dimensions_two_through_twenty_five_are_partitioned() -> None:
    records = [field_constant_record(n) for n in range(2, 26)]
    assert [record["n"] for record in records] == list(range(2, 26))
    assert [
        record["n"] for record in records if record["nonrational_possible"]
    ] == list(range(2, 10))
    assert field_constant_record(19)["full_trace_K3_degree_ceiling"] == 1
    assert field_constant_record(20)["full_trace_K3_degree_ceiling"] == 0
    assert field_constant_record(20)["actual_field_degree_bound"] == 1
    assert field_constant_record(25)["field"] == "Q"


def test_8192_factor32_and_exponent30_are_literal_maxima() -> None:
    certificate = uniform_constant_certificate()
    assert certificate["maximum_primitive_factor"] == 32
    assert certificate["maximum_small_coordinate_bound"] == 8192
    assert certificate["maximum_discriminant_exponent"] == 30
    assert certificate["uniform_formula"] == "(32*max(8192,C/s0))^30"
    assert certificate["n2_formula"] == "(32*max(64,C/s0))^30"


def test_n2_sharp_coarse_record_is_exact() -> None:
    assert field_constant_record(2) == {
        "n": 2,
        "full_trace_K3_degree_ceiling": 6,
        "actual_field_degree_bound": 6,
        "chordless_length_bound": 5,
        "selected_generator_count": 2,
        "primitive_grid_side": 16,
        "small_coordinate_bound": 64,
        "primitive_bound_factor": 32,
        "discriminant_exponent": 30,
        "nonrational_possible": True,
    }


def test_dominant_one_positive_place_identity_is_replayed() -> None:
    certificate = dominant_height_certificate()
    assert certificate["exact_samples"] > 20
    assert certificate[
        "dominance_and_negative_definite_other_summands_imply_height"
    ] is True


def test_positive_grid_and_vandermonde_exponents_are_exact() -> None:
    records = primitive_grid_certificate()
    assert [record["degree_ceiling"] for record in records] == list(range(2, 7))
    for record in records:
        assert record["union_bound_for_bad_points"] < record["grid_points"]
        degree = record["degree_ceiling"]
        assert record["vandermonde_exponent"] == degree * (degree - 1)


def test_owner_archive_is_independently_replayed() -> None:
    record = verify_owner_archive()
    assert record["payload_digest_verified"] is True
    assert record["dimensions_replayed"] == [2, 25]
    assert record["n20_trace_ceiling"] == 0
    assert record["n20_actual_field_degree_bound"] == 1


def test_final_verdict_retains_all_nonclaims() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_CONDITIONAL"
    separation = audit["logical_separation"]
    assert separation["finite_field_set_proved_under_uniform_height"] is True
    assert separation["finite_fields_imply_finite_hosts"] is False
    assert separation[
        "finite_hosts_plus_relative_saturation_give_finite_chambers"
    ] is True
    assert separation["gate5_automatic_from_exact_maximality_in_n_ge_2"] is False
    assert separation["raw_essential_surjectivity"] is False
