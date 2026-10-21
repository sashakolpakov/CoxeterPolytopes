"""Regression tests for the fixed-Jacobian K3 cusp-index no-go."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_noncompact_k3_fixed_jacobian_index_blindness import (
    build_certificate,
    digest,
    encoded_certificate,
    family_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/noncompact_k3_fixed_jacobian_index_blindness.md"
SCOPE_AUDIT = ROOT / "results/noncompact_k3_fixed_jacobian_index_blindness_hostile_scope_audit.md"
CERTIFIER = ROOT / "scripts/certify_noncompact_k3_fixed_jacobian_index_blindness.py"
ARCHIVE = ROOT / "data/enumerations/noncompact_k3_fixed_jacobian_index_blindness.json"


def test_fixed_jacobian_family_exact_lattice_invariants():
    for d in range(2, 45):
        record = family_record(d)
        assert record["signature"] == [1, 17]
        assert record["picard_rank"] == 18
        assert record["absolute_discriminant"] == d * d
        assert record["discriminant_group_invariant_factors"] == [d, d]
        assert record["discriminant_length"] == 2
        assert record["primitive_embedding_unit_minor"] == 1
        assert record["transcendental_embedding_unit_minor"] == -1
        assert record["orthogonal_complement"] == "U(d)+U"
        assert record["actual_fibre_packet_marking"] == "F=u and vertical roots=E8(-1)^2"
        assert record["primitive_fibre_divisibility"] == d
        assert record["multisection_index"] == d


def test_affine_packet_is_full_only_in_the_local_parabolic_sense():
    for d in (2, 3, 5, 11, 29):
        record = family_record(d)
        assert record["fixed_fibres"] == ["II*", "II*", "I1", "I1", "I1", "I1"]
        assert record["euler_sum"] == 24
        assert record["jacobian_mordell_weil_rank"] == 0
        assert record["incident_curve_count"] == 18
        assert record["incident_vector_span_rank"] == 17
        assert record["incident_gram_rank"] == 16
        assert record["finite_root_rank"] == record["expected_finite_root_rank"] == 16
        assert record["common_orthogonal_dimension"] == 1


def test_transverse_height_escapes_exactly():
    for d in range(2, 60):
        record = family_record(d)
        assert record["transverse_lattice_root_square"] == -2
        assert record["transverse_lattice_root_fibre_degree"] == d
        assert min(record["transverse_root_incident_pairings"]) >= 0
        assert record["all_positive_divisor_fibre_degrees_divisible_by_d"] is True


def test_hosts_are_pairwise_nonisometric_by_discriminant():
    discriminants = [family_record(d)["absolute_discriminant"] for d in range(2, 60)]
    assert len(discriminants) == len(set(discriminants))


def test_archive_replays_exactly():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_scope_forbids_a_false_finite_volume_counterfamily():
    scope = build_certificate()["scope"]
    assert scope["same_fixed_relative_jacobian_for_all_d"] is True
    assert scope["projective_K3_torsor_for_all_d"] is True
    assert scope["full_local_parabolic_rank_for_selected_cusp"] is True
    assert scope["full_finite_volume_chamber_claimed_for_X_d"] is False
    assert scope["full_automorphism_group_finite_claimed_for_X_d"] is False
    assert scope["displayed_transverse_lattice_root_claimed_irreducible_physical"] is False
    assert scope["finite_volume_counterfamily_claimed"] is False
    assert scope["global_bounded_spanning_root_theorem_proved"] is False
    assert scope["published_118_host_classification_reproved"] is False


def test_report_states_the_exact_gate_and_sources():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "one fixed Jacobian elliptic K3 surface",
        r"\operatorname{NS}(X_d)\simeq U(d)\oplus E_8(-1)^2",
        r"\operatorname{ind}(X_d)=\operatorname{ord}(\alpha_d)=d",
        r"T(X_d)\simeq U(d)\oplus U",
        "full-rank **local parabolic**",
        "not asserted to be the complete walls of a finite-volume chamber",
        "not a counterfamily to the global noncompact host theorem",
        "field-collapse theorem has already reduced",
        "Meinsma--Shinder, Lemma 4.4",
        r"\Sha_{\rm alg}(S)\simeq\operatorname{Br}(S)",
        "Theorem 14.1.12 (Nikulin, Theorem 1.14.4)",
        "actual marked splitting",
        r"J^k(X)_t\simeq X_t",
        "bounded physical wall outside",
    ):
        assert phrase in text


def test_hostile_scope_audit_separates_local_packets_from_global_volume():
    text = " ".join(SCOPE_AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "PASS, with a deliberately local conclusion",
        r"projective K3 torsor for every \(d\)",
        "every closed fibre is isomorphic",
        "full **local parabolic rank**, not full rank in the host",
        "determinant equality alone is not being used",
        "actual vertical roots",
        "Bounded discriminant **length** is blind",
        "complete discriminant form is not blind",
        "not the complete walls of a finite-volume chamber",
        "Neither finite volume",
        r"\(d\le390\) is only a **cusp-adapted conditional**",
        r"\le(4+17\cdot13^2)^9=2877^9",
    ):
        assert phrase in text


def test_transport_controls_are_clean():
    combined = REPORT.read_bytes() + SCOPE_AUDIT.read_bytes() + CERTIFIER.read_bytes()
    assert b"\x08" not in combined
    text = combined.decode("utf-8")
    assert ",quad" not in text
    assert ",qquad" not in text


def test_invalid_small_parameter_is_rejected():
    with pytest.raises(ValueError):
        family_record(1)
