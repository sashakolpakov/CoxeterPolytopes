"""Independent hostile tests for compact physical-star K3 host finiteness."""

from __future__ import annotations

import json

from scripts.audit_k3_compact_physical_star_trace_host_finiteness_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    energy_and_codifferent_certificate,
    full_picard_descent_certificate,
    order_and_overlattice_certificate,
    sha256,
)


def test_all_owner_and_dependency_bytes_are_frozen() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_energy_budget_closes_the_codifferent_box() -> None:
    record = energy_and_codifferent_certificate()
    assert record["place_energy_identity"] == "sum_tau(-tau(g_ii))=2"
    assert record["codifferent_box_is_finite"] is True
    assert record["bounded_codifferent_box_points"] > 0
    assert record["order_codifferent_basis"] == ["1/2", "sqrt(5)/10"]


def test_owner_trace_bound_is_conservative() -> None:
    record = energy_and_codifferent_certificate()
    assert record["owner_trace_bound"] == "abs(B(delta_i,delta_j))<16"
    assert record["independent_stronger_trace_bound"] == (
        "abs(B(delta_i,delta_j))<14"
    )
    assert record["sampled_weighted_bound_max"] < 14


def test_cycle_order_blocks_conductor_escape() -> None:
    record = order_and_overlattice_certificate()
    assert record["cycle_fraction_field_is_k"] is True
    assert record["coefficient_order_interval"] == "O0 subset O subset O_k"
    assert record["toy_overorder_index"] == 2
    assert record["toy_overorder_count"] == 2


def test_root_frame_duality_captures_every_integral_module() -> None:
    record = order_and_overlattice_certificate()
    assert record["root_frame_containment"] == "F subset M subset F#"
    assert record["root_frame_discriminant_group_order"] == 8
    assert record["toy_overlattice_index"] == 2
    assert record["toy_overlattice_even_and_integral"] is True
    assert record["index_square_identity_checked"] is True


def test_full_picard_real_line_scalar_descends_to_k() -> None:
    record = full_picard_descent_certificate()
    assert record["physical_wall_vector_lies_in_M"] is True
    assert record["distinguished_real_line_proportionality_descends_to_k"] is True
    assert record["primitive_trace_block_in_larger_NS_is_not_enough"] is True


def test_native_and_ordinary_reflection_species_are_not_conflated() -> None:
    record = full_picard_descent_certificate()
    assert record["native_quadratic_reflection_Q_minus_dimension"] == 2
    assert record["ordinary_K3_root_reflection_Q_minus_dimension"] == 1
    assert record["reflection_species_are_distinct_for_degree_gt_one"] is True


def test_same_host_composition_and_nonclaims_are_pinned() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_EXACT_STATED_SCOPE"
    composition = audit["same_host_composition"]
    assert composition["same_marked_native_host_maximality_required"] is True
    assert composition[
        "maximality_implies_W_equals_complete_native_reflection_group"
    ] is True
    assert composition[
        "complete_locally_finite_mirror_arrangement_one_chamber_orbit"
    ] is True
    assert all(audit["scope"].values())


def test_owner_archive_is_self_hashed_and_pins_owner_report() -> None:
    archive_path = ROOT / PINS["owner_archive"][0]
    archive = json.loads(archive_path.read_text(encoding="utf-8"))
    payload = archive["payload"]
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    assert archive["payload_sha256"] == __import__("hashlib").sha256(
        canonical
    ).hexdigest()
    assert payload["report_sha256"] == PINS["owner_report"][1]
