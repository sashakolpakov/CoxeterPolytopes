"""Regression tests for the finite-volume K3 bounded spanning-root lemma."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import pytest

from scripts.certify_noncompact_k3_bounded_spanning_root_lemma import (
    bounded_gram_count,
    build_certificate,
    canonical_vertex,
    digest,
    encoded_certificate,
    lorentz_pair,
    matrix_rank,
    orthogonal_hyperplane_basis,
    rank_record,
    vertex_pair_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/noncompact_k3_bounded_spanning_root_lemma.md"
SCOPE_AUDIT = ROOT / "results/noncompact_k3_bounded_spanning_root_lemma_hostile_scope_audit.md"
CERTIFIER = ROOT / "scripts/certify_noncompact_k3_bounded_spanning_root_lemma.py"
ARCHIVE = ROOT / "data/enumerations/noncompact_k3_bounded_spanning_root_lemma.json"


@pytest.mark.parametrize("rho", range(3, 21))
@pytest.mark.parametrize(
    ("first_kind", "second_kind"),
    [
        ("finite", "finite"),
        ("finite", "ideal"),
        ("ideal", "finite"),
        ("ideal", "ideal"),
    ],
)
def test_two_distinct_complete_vertex_packets_span(rho, first_kind, second_kind):
    record = vertex_pair_record(rho, first_kind, second_kind)
    assert record["first_packet_rank"] == rho - 1
    assert record["second_packet_rank"] == rho - 1
    assert record["two_distinct_packet_union_rank"] == rho


@pytest.mark.parametrize("rho", range(3, 21))
@pytest.mark.parametrize("kind", ["finite", "ideal"])
def test_one_vertex_packet_is_exactly_the_cusp_local_deficit(rho, kind):
    vertex = canonical_vertex(kind, rho, 0)
    packet = orthogonal_hyperplane_basis(vertex)
    assert matrix_rank(packet) == rho - 1
    assert all(lorentz_pair(vertex, normal) == 0 for normal in packet)
    assert matrix_rank(packet + packet) == rho - 1


@pytest.mark.parametrize("rho", range(3, 21))
def test_integral_gram_box_is_finite_and_uses_bound_13(rho):
    record = rank_record(rho)
    assert record["strict_normalized_bound"] == 14
    assert record["integral_off_diagonal_bound"] == 13
    assert record["off_diagonal_slots"] == rho * (rho - 1) // 2
    assert record["bounded_gram_box_count"] == bounded_gram_count(rho)
    assert record["bounded_gram_box_count"] == 14 ** (rho * (rho - 1) // 2)
    assert record["sample_lorentz_gram_determinant"] != 0


def test_archive_replays_exactly():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_command_line_archive_check():
    result = subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "PASS archive replay" in result.stdout


def test_scope_flags_keep_every_essential_gate():
    scope = build_certificate()["hostile_scope"]
    assert scope["finite_sided_exact_return_required"] is True
    assert scope["intersecting_vs_ultraparallel_visual_split_required"] is True
    assert scope["finite_volume_required_at_ideal_vertices"] is True
    assert scope["projective_closed_star_required"] is True
    assert scope["complete_incident_cusp_packet_required"] is True
    assert scope["two_distinct_vertices_required"] is True
    assert scope["physical_literal_square_minus_two_facets_required"] is True
    assert scope["simplicity_required"] is False
    assert scope["compactness_required"] is False
    assert scope["cusp_classification_required"] is False
    assert scope["reflective_lattice_classification_required"] is False
    assert scope["reflection_maximality_required_for_spanning_lemma"] is False
    assert scope["Property_4_proved"] is False


def test_report_states_every_noncompact_bridge_and_exact_conclusion():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "uniform constant \\(B=13\\)",
        "closed star",
        "strictly below \\(14\\)",
        "It remains to treat ultraparallel supports",
        "not being asserted for overlapping visual cones",
        "two distinct extreme points",
        "boundedness of the horospherical link",
        "nonsimple cusps",
        r"x^\perp+y^\perp=L_{\mathbf R}",
        r"0\le E_i.E_j\le13",
        r"R\subseteq L\subseteq R^*",
        "Finite-sided exact return",
        "Projective closed incidence",
        "Full incident cusp packet",
        "Literal physical roots",
        "No Property 4 conclusion",
    ):
        assert phrase in text


def test_scope_audit_records_the_sharp_one_cusp_failure():
    text = " ".join(SCOPE_AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "PASS",
        "one complete cusp packet has rank \\(\\rho-1\\)",
        "two distinct projective points",
        "closed star",
        "Intersecting or asymptotic supports",
        "ultraparallel supports",
        "same physical square-\\((-2)\\) vector",
        "Property 4",
    ):
        assert phrase in text


def test_invalid_rank_and_vertex_kind_are_rejected():
    with pytest.raises(ValueError):
        rank_record(2)
    with pytest.raises(ValueError):
        rank_record(21)
    with pytest.raises(ValueError):
        canonical_vertex("hyperideal", 3, 0)


def test_transport_controls_are_clean():
    combined = REPORT.read_bytes() + SCOPE_AUDIT.read_bytes() + CERTIFIER.read_bytes()
    assert b"\x08" not in combined
    text = combined.decode("utf-8")
    assert ",quad" not in text
    assert ",qquad" not in text
