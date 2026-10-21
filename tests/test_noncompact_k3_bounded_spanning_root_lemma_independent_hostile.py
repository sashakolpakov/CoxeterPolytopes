"""Independent hostile tests for the cusped K3 spanning-root theorem."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import pytest

from scripts.audit_noncompact_k3_bounded_spanning_root_lemma_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    complete_box_cusp_packet,
    host_and_chamber_certificate,
    sha256,
    two_vertex_certificate,
    visual_formula_certificate,
)


AUDIT = ROOT / "results/noncompact_k3_bounded_spanning_root_lemma_independent_hostile_audit.md"
SCRIPT = ROOT / "scripts/audit_noncompact_k3_bounded_spanning_root_lemma_independent_hostile.py"


def test_all_owner_bytes_are_frozen() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_owner_archive_and_owner_certifier_replay() -> None:
    archive = ROOT / PINS["owner_archive"][0]
    envelope = json.loads(archive.read_text(encoding="utf-8"))
    canonical = json.dumps(
        envelope["payload"], sort_keys=True, separators=(",", ":")
    ).encode()
    assert envelope["payload_sha256"] == hashlib.sha256(canonical).hexdigest()
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / PINS["owner_certifier"][0]),
            "--check",
            str(archive),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "PASS archive replay" in result.stdout


def test_visual_identity_and_strict_bound_are_independent() -> None:
    record = visual_formula_certificate()
    assert record["observer_frame_symbolic_residual"] == "0"
    assert record["positive_gap_required_for_formula"] is True
    assert record["gap_at_most_minimum_aperture_required"] is True
    assert record["intersecting_or_asymptotic_branch_uses_b_at_most_2"] is True
    assert record["ultraparallel_facet_shadows_disjoint_after_irredundancy"] is True
    assert record["nested_ultraparallel_shadows_would_make_one_wall_redundant"] is True
    assert record["sample_count"] > 30_000
    assert 13.9 < record["sample_maximum_normalized_pairing"] < 14
    assert record["sample_maximum_sine_quotient"] < 2


@pytest.mark.parametrize("rho", range(3, 21))
def test_complete_bounded_cusp_link_spans_u_perp_even_when_nonsimple(rho: int) -> None:
    record = complete_box_cusp_packet(rho)
    assert record["packet_rank"] == rho - 1
    assert record["expected_u_perp_rank"] == rho - 1
    assert record["all_lifts_orthogonal_to_u"] is True
    assert record["all_sample_lifts_square_minus_two"] is True
    assert record["positive_horizontal_dependence_lifts_to_nonzero_u"] is True
    assert record["nonsimple_when_rank_at_least_four"] is (rho >= 4)


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
def test_two_distinct_extreme_point_packets_span_full_space(
    rho: int, first_kind: str, second_kind: str
) -> None:
    record = two_vertex_certificate(rho, first_kind, second_kind)
    assert record["one_packet_rank"] == rho - 1
    assert record["duplicated_packet_rank"] == rho - 1
    assert record["two_packet_union_rank"] == rho
    assert record["distinct_projective_lines"] is True


@pytest.mark.parametrize("rho", range(3, 21))
def test_physical_bound_host_finiteness_and_weyl_interface(rho: int) -> None:
    record = host_and_chamber_certificate(rho)
    assert record["strict_normalized_less_than"] == 14
    assert record["integral_physical_intersection_at_most"] == 13
    assert record["bounded_gram_box_count"] == 14 ** (rho * (rho - 1) // 2)
    assert record["sample_lorentz_gram_determinant"] != 0
    assert record["sample_discriminant_group_order"] > 0
    assert record["R_subset_L_subset_R_dual_from_integrality_and_full_rank"] is True
    assert record["intermediate_overlattices_in_finite_R_dual_mod_R"] is True
    assert record["root_reflection_is_integral_involutive_isometry"] is True
    assert record["complete_root_mirrors_are_locally_finite"] is True
    assert record["generic_compact_path_crosses_finitely_many_mirrors"] is True
    assert record["successive_wall_reflections_give_one_Weyl_chamber_orbit"] is True


def test_scope_keeps_property_four_and_weaker_returns_out() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_EXACT_FINITE_SIDED_FULL_PICARD_SCOPE"
    scope = audit["scope"]
    for required in (
        "finite_sidedness_needed_for_minimum_and_uniform_cusp_height",
        "facet_irredundancy_used_for_disjoint_ultraparallel_shadows",
        "projective_closed_star_needed",
        "finite_volume_needed_for_bounded_cusp_link",
        "complete_incident_packet_needed",
        "two_distinct_extreme_points_needed",
        "literal_same_vector_physical_minus_two_roots_needed",
        "full_Picard_complete_root_arrangement_needed_for_one_orbit",
    ):
        assert scope[required] is True
    assert scope["cusp_simplicity_needed"] is False
    assert scope["compactness_needed"] is False
    assert scope["reflection_maximality_needed"] is False
    assert scope["Property_4_proved"] is False


def test_independent_command_line_rechecks_pins() -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--check-pins"],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)
    assert payload["verdict"] == "PASS_EXACT_FINITE_SIDED_FULL_PICARD_SCOPE"
    assert set(payload["pins"]) == set(PINS)


def test_audit_report_contains_every_hostile_gate() -> None:
    text = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "PASS in the exact finite-sided, finite-volume",
        "implicit geometric sublemma, not a defect",
        "nesting would make one of their supporting half-space inequalities redundant",
        "positive gap",
        "intersect in hyperbolic space",
        "if they are asymptotic",
        r"0<c\leq\theta_E\leq\min",
        "complete nonsimple cusp packet",
        "before any vertex assumption",
        r"x^\perp+y^\perp=L_{\mathbf R}",
        "including cross-pairs",
        "literal class of an irreducible smooth rational K3 curve",
        r"R\subseteq L\subseteq R^*",
        "one Weyl group orbit",
        "Property 4",
    ):
        assert phrase in text


def test_transport_and_control_bytes_are_clean() -> None:
    payload = AUDIT.read_bytes() + SCRIPT.read_bytes() + Path(__file__).read_bytes()
    forbidden = set(range(0, 9)) | {11, 12} | set(range(14, 32)) | {127}
    assert not any(byte in forbidden for byte in payload)
    text = payload.decode("utf-8")
    assert "," + "quad" not in text
    assert "," + "qquad" not in text
    assert "|" + "le" not in text
    assert "TO" + "DO" not in text
    assert "FIX" + "ME" not in text


def test_no_inline_tex_delimiters_were_stripped() -> None:
    source = AUDIT.read_text(encoding="utf-8")
    outside_display: list[str] = []
    in_display = False
    for line in source.splitlines():
        if line.strip() == r"\[":
            in_display = True
            continue
        if line.strip() == r"\]":
            in_display = False
            continue
        if not in_display:
            outside_display.append(line)
    prose = "\n".join(outside_display)
    prose = re.sub(r"\\\(.*?\\\)", "", prose)
    prose = re.sub(r"\x60[^\x60]*\x60", "", prose)
    for command in (
        r"\rho",
        r"\theta",
        r"\pi",
        r"\mathbf",
        r"\perp",
        r"\leq",
        r"\geq",
        r"\overline",
        r"\lambda",
        r"\binom",
        r"\det",
    ):
        assert command not in prose
    for stripped in (
        "below (14)",
        "observer (O)",
        "pairing is (2p",
        "Let (E)",
        "spans (u^",
        "square-((-2))",
        "Choose (rho)",
    ):
        assert stripped not in prose
