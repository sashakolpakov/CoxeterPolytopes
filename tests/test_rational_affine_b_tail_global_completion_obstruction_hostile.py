"""Hostile regressions for the affine-B-tail completion obstruction."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.audit_rational_affine_b_tail_global_completion_obstruction_hostile import (
    OUTPUT,
    build_payload,
    digest,
    pairing,
    r_root,
    truncator_record,
    visible_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = (
    ROOT
    / "results"
    / "rational_affine_b_tail_global_completion_obstruction_hostile_audit.md"
)
PRODUCER_ARCHIVE = (
    ROOT
    / "data"
    / "enumerations"
    / "rational_affine_b_tail_global_completion_obstruction.json"
)
AUDITOR = (
    ROOT
    / "scripts"
    / "audit_rational_affine_b_tail_global_completion_obstruction_hostile.py"
)


def test_hostile_archive_is_exact_and_independent():
    archive = json.loads(OUTPUT.read_text())
    payload = build_payload()
    assert archive["payload"] == payload
    assert archive["payload_sha256"] == digest(payload)
    text = AUDITOR.read_text()
    assert "from scripts.certify_rational_affine_b_tail" not in text
    assert payload["verdict"] == "PASS_AFTER_SCOPE_AND_CERTIFICATE_REPAIRS"


def test_visible_simplex_threshold_and_first_escape():
    records = [visible_record(length) for length in range(2, 49)]
    assert [record["m"] for record in records if record["finite_volume"]] == [
        2,
        3,
        4,
    ]
    for record in records:
        length = record["m"]
        assert abs(record["coordinate_determinant"]) == 1
        assert record["dual_squares"][-1] == 4 - length
        assert record["dual_squares"][:3] == [0, 0, 1]

    gate = build_payload()["hidden_wall_gate"]
    assert gate["necessary_inequality"] == "sum(mu_1,...,mu_5) > 2*a+b"
    assert gate["finite_volume_requires_some_new_wall_negative_on_v5"] is True


def test_ruling_swap_is_uniform_full_NS_obstruction():
    payload = build_payload()
    for record in payload["swap_records"]:
        assert record["root_square"] == -2
        assert record["integral"] is True
        assert record["swaps_rulings"] is True
        assert record["fixes_exceptional_coordinates"] is True
        assert record["symmetric_ample_pairing"] == 0

    result = payload["full_NS_uniform_obstruction"]
    assert result["obstructing_root"] == "c-d"
    assert result["absolute_full_NS_reflection_maximality"] is False
    assert result["adjoining_it_is_a_proper_discrete_reflection_overgroup"] is True


def test_T5_survives_and_producer_archive_retains_it():
    for length in (5, 6, 17, 40):
        record = truncator_record(length)
        assert record["square"] == -1
        assert record["canonical_pairing"] == -1
        assert record["arithmetic_genus"] == 0
        assert record["pairing_with_v5"] == -1
        assert record["effective_divisor_class"] is True
        assert record["irreducible_smooth_facet_automatic"] is False
        assert len(record["vector"]) == length + 2

    producer = json.loads(PRODUCER_ARCHIVE.read_text())["payload"]
    samples = producer["first_truncating_root_samples"]
    for length in range(5, 41):
        assert samples[str(length)][0]["label"] == "T_5"
        assert len(samples[str(length)][0]["vector"]) == length + 2
        assert "facet status are not automatic" in samples[str(length)][0]["status"]


def test_first_height_R_family_adjunction_pairings_and_jet_scope():
    payload = build_payload()
    gate = payload["first_height_classification"]
    assert gate["beta_one_negative_root_exists"] is False
    assert gate["minimum_positive_normalized_height"] == 2
    for record in gate["height_two_roots"]:
        assert record["nonzero_coefficients"] == 4 * record["alpha"] + 2
        assert record["square"] == -2

    records = payload["R_family_records"]
    assert [record["arithmetic_genus"] for record in records] == list(range(12))
    for record in records:
        alpha = record["a"]
        assert record["support"] == 4 * alpha + 2
        assert record["canonical_pairing"] == 2 * alpha - 2
        assert record["minimum_row_dependencies_for_effectivity"] == alpha
        assert (
            record[
                "irreducible_exact_multiplicity_locus_codimension_lower_bound"
            ]
            == alpha
        )
        assert record["smooth_rational_adjunction_compatible"] == (alpha == 1)

    length = 45
    for alpha in range(1, 8):
        for beta in range(alpha + 1, 9):
            assert pairing(r_root(alpha, length), r_root(beta, length)) == (
                2 * (beta - alpha - 1)
            )


def test_adapted_host_and_facet_scope_are_not_overclaimed():
    payload = build_payload()
    adapted = payload["adapted_host_boundary"]
    assert adapted["parity_required_for_every_half_boundary"] is True
    assert adapted["full_NS_obstruction_automatically_transfers_to_every_adapted_host"] is False

    boundary = payload["claim_boundary"]
    assert boundary["T5_effective_divisor_class_is_automatic"] is True
    assert boundary["T5_irreducible_smooth_nef_facet_is_automatic"] is False
    assert boundary["Ra_effectivity_is_automatic"] is False
    assert boundary["reflection_maximality_makes_every_root_a_physical_facet"] is False


def test_hostile_report_records_repairs_and_exact_scope():
    text = REPORT.read_text()
    required = [
        "PASS after scope and certificate repairs",
        r"For every \(m\ge2\)",
        r"\(\mathbf P^1\times\mathbf P^1\)",
        "adjoined as a proper discrete reflection",
        r"\boxed{m\le4}",
        "first **square-minus-two,",
        "does not prove irreducibility, smoothness, or",
        "maximality does not turn every lattice root",
        "Neron--Severi host",
        "requires its own parity calculation",
        "no root census or external classification",
        "does not cover arbitrary reducible",
    ]
    for phrase in required:
        assert phrase in text
    assert "K.T_5=-1,quad" not in text
