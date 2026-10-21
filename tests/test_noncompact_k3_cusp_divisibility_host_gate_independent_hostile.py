"""Independent hostile tests for the noncompact K3 cusp-divisibility gate."""

from __future__ import annotations

import hashlib
import json

from scripts.audit_noncompact_k3_cusp_divisibility_host_gate_independent_hostile import (
    PINS,
    ROOT,
    bounded_spanning_repair_certificate,
    build_audit,
    lattice_and_embedding_certificate,
    multisection_certificate,
    perpendicular_roots_certificate,
    sha256,
)


def test_owner_bytes_are_terminally_frozen() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_explicit_primitive_K3_embedding_has_a_unit_minor() -> None:
    for index in range(2, 50):
        record = lattice_and_embedding_certificate(index)
        assert record["signature"] == [1, 2]
        assert record["absolute_discriminant"] == 2 * index**2
        assert record["primitive_embedding_unit_minor"] == 1
        assert record["primitive_K3_embedding"] is True
        assert record["generic_period_exact_Picard"] is True


def test_affine_cusp_packet_is_fixed_while_divisibility_escapes() -> None:
    rows = [lattice_and_embedding_certificate(index) for index in range(2, 30)]
    assert {tuple(map(tuple, row["fixed_affine_A1_gram"])) for row in rows} == {
        ((-2, 2), (2, -2))
    }
    assert [row["cusp_divisibility"] for row in rows] == list(range(2, 30))
    assert len({row["absolute_discriminant"] for row in rows}) == len(rows)


def test_all_roots_perpendicular_to_the_cusp_have_the_claimed_form() -> None:
    record = perpendicular_roots_certificate(11)
    assert record["root_formula"] == "m*u plus_or_minus a"
    assert record["enumerated_root_count"] == 2 * (2 * 30 + 1)
    assert record["alpha_beta_are_consecutive_affine_simple_roots"] is True
    assert record["Weyl_reframing_makes_them_physical_curves"] is True


def test_divisibility_is_an_attained_multisection_degree() -> None:
    for index in (2, 3, 5, 11, 29):
        record = multisection_certificate(index, twist=4)
        assert record["test_divisor_square"] == 8 * index
        assert record["test_divisor_fibre_degree"] == index
        assert record["RR_forces_test_divisor_effective"] is True
        assert record["one_horizontal_component_has_degree_exactly_d"] is True
        assert record["minimum_positive_multisection_degree"] == index


def test_bounded_full_spanning_root_packet_repairs_host_finiteness() -> None:
    record = bounded_spanning_repair_certificate()
    assert record["finite_gram_box_size"] == 64
    assert record["nonsingular_gram_count"] > 0
    assert record["full_span_gives_R_subset_NS_subset_R_dual"] is True
    assert record["finite_discriminant_group_has_finitely_many_subgroups"] is True
    assert record["finite_integral_even_K3_overlattices"] is True
    assert record["same_host_maximality_is_downstream"] is True


def test_owner_archive_is_self_hashed() -> None:
    archive_path = ROOT / PINS["owner_archive"][0]
    envelope = json.loads(archive_path.read_text(encoding="utf-8"))
    canonical = json.dumps(
        envelope["payload"], sort_keys=True, separators=(",", ":")
    ).encode()
    assert envelope["payload_sha256"] == hashlib.sha256(canonical).hexdigest()


def test_scope_does_not_promote_the_local_family_to_the_target_theorem() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_EXACT_STATED_SCOPE"
    assert all(audit["scope"].values())
