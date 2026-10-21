"""Regression tests for the noncompact K3 cusp-divisibility host gate."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_noncompact_k3_cusp_divisibility_host_gate import (
    build_certificate,
    digest,
    encoded_certificate,
    family_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/noncompact_k3_cusp_divisibility_host_gate.md"
CERTIFIER = ROOT / "scripts/certify_noncompact_k3_cusp_divisibility_host_gate.py"
ARCHIVE = ROOT / "data/enumerations/noncompact_k3_cusp_divisibility_host_gate.json"


def test_unbounded_family_formulas_and_primitive_embeddings():
    for d in range(2, 80):
        record = family_record(d)
        assert record["signature"] == [1, 2]
        assert record["absolute_discriminant"] == 2 * d * d
        assert record["cusp_divisibility"] == d
        assert record["fixed_affine_A1_gram"] == [[-2, 2], [2, -2]]
        assert record["transverse_root_square"] == -2
        assert record["transverse_root_cusp_degree"] == d
        assert record["root_basis_change_determinant"] == 1
        assert record["all_three_root_reflections_integral"] is True
        assert record["primitive_K3_embedding_minor_gcd"] == 1


def test_family_is_pairwise_nonisometric_by_discriminant():
    discriminants = [family_record(d)["absolute_discriminant"] for d in range(2, 80)]
    assert len(discriminants) == len(set(discriminants))


def test_archive_replays_exactly():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_nonclaims_are_literal():
    scope = build_certificate()["scope"]
    assert scope["finite_volume_claimed_for_M_d"] is False
    assert scope["reflection_maximality_claimed_for_M_d"] is False
    assert scope["infinite_family_under_all_target_hypotheses_claimed"] is False
    assert scope["global_bounded_spanning_root_lemma_proved_here"] is False
    assert scope["finite_host_theorem_proved_without_that_lemma"] is False
    assert scope["Property_4_proved"] is False


def test_report_states_the_gate_and_avoids_the_false_counterexample_claim():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "cusp-local host finiteness fails",
        r"M_d=U(d)\oplus\langle-2\rangle",
        r"|\operatorname{disc}M_d|=2d^2",
        r"\operatorname{div}_{M_d}(u)",
        "fixed affine",
        "primitive embedding",
        "bounded spanning physical-root packet",
        r"R\subseteq\operatorname{NS}(X)\subseteq R^*",
        "same-host reflection maximality is downstream",
        "not asserted to have finite volume",
        "not a counterfamily to the full target theorem",
    ):
        assert phrase in text


def test_transport_controls_are_clean():
    combined = REPORT.read_bytes() + CERTIFIER.read_bytes()
    assert b"\x08" not in combined
    text = combined.decode("utf-8")
    assert ",quad" not in text
    assert ",qquad" not in text


def test_invalid_small_parameter_is_rejected():
    with pytest.raises(ValueError):
        family_record(1)
