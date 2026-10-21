"""Regression tests for the compact right-angled K3 hexagon theorem."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import pytest

from scripts.certify_k3_n2_compact_right_hexagon_field_bound import (
    build_certificate,
    closure_record,
    digest,
    encoded_certificate,
    field_record,
    right_hexagon_from_alternating_cosh,
    short_pair_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_compact_right_hexagon_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_compact_right_hexagon_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_compact_right_hexagon_field_bound.json"


def test_rank_three_cofactor_and_four_minor_closure_on_exact_grid():
    tested = 0
    values = (
        Fraction(5, 4),
        Fraction(4, 3),
        Fraction(3, 2),
        Fraction(5, 3),
        Fraction(2),
        Fraction(7, 3),
    )
    for u in values:
        for v in values:
            for w in values:
                record = closure_record(u, v, w)
                assert record["common_numerator_residuals"] == ["0"] * 3
                assert record["four_consecutive_Gram_minor_residuals"] == ["0"] * 6
                coordinates = right_hexagon_from_alternating_cosh(u, v, w)
                assert min(coordinates) > 0
                tested += 1
    assert tested == 216


def test_uniform_short_pair_constant_is_exact():
    record = short_pair_record()
    assert record["alternating_threshold"] == 16
    assert record["large_alternating_case_odd_upper_bound"] == "(51+34sqrt(17))/256"
    assert record["resulting_rational_numerator_upper_bound"] == 221
    assert record["large_alternating_case_odd_strict_upper_bound"] == 1
    assert record["universal_some_shifted_coordinate_strict_upper_bound"] == 16
    assert record["universal_some_pair_cycle_strict_upper_bound"] == 68
    assert 51 + 34 * 5 == 221 < 256


def test_large_alternating_coordinates_force_a_short_odd_coordinate():
    tested = 0
    for u in (Fraction(5), Fraction(7), Fraction(10), Fraction(20)):
        for v in (Fraction(5), Fraction(8), Fraction(11), Fraction(23)):
            for w in (Fraction(6), Fraction(9), Fraction(13), Fraction(29)):
                coordinates = right_hexagon_from_alternating_cosh(u, v, w)
                alternating = coordinates[0], coordinates[2], coordinates[4]
                if min(alternating) < 16:
                    continue
                assert min(coordinates[1], coordinates[3], coordinates[5]) < 1
                assert min(coordinates) < 16
                tested += 1
    assert tested > 40


def test_field_generator_and_discriminant_constants_recompute():
    record = field_record()
    assert record["field_degree_upper_bound"] == 6
    assert record["primitive_generator_conjugate_interval"] == "(0,68)"
    assert record["unique_distinguished_generator_interval"] == "(4,68)"
    assert record["definite_generator_interval"] == "(0,4)"
    assert record["Vandermonde_difference_bound"] == 68
    assert record["field_discriminant_bound"] == "68^30"
    assert record["unique_large_conjugate_forces_full_field_generation"] is True


def test_archive_replays_and_scope_is_exact():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    conclusions = payload["conclusions"]
    assert conclusions["K3_full_trace_field_degree_at_most_6"] is True
    assert conclusions["some_dotted_pair_cycle_below_68"] is True
    assert conclusions["that_dotted_pair_cycle_generates_the_full_field"] is True
    assert conclusions["field_discriminant_at_most_68_power_30"] is True
    assert conclusions["imposed_mixed_height_required"] is False
    assert conclusions["reflection_maximality_used"] is False
    scope = payload["scope"]
    assert scope["minimal_cyclic_Cartan_ground_field_required"] is True
    assert scope["exact_same_compact_right_angled_hexagon_eigen_nef_chamber_required"] is True
    assert scope["non_right_hexagon_closed"] is False
    assert scope["noncompact_hexagon_closed"] is False
    assert scope["arbitrary_six_or_more_facet_chamber_closed"] is False
    assert scope["carrier_existence_proved"] is False


def test_report_states_hypotheses_proof_and_nonclaims():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "proved for every full-trace complex-K3 field degree",
        "minimal complex projective K3 surface",
        "primitive full trace",
        "exactly a compact right-angled Coxeter hexagon",
        "minimal cyclic Cartan ground field",
        r"x_1&={S\over x_0x_2}",
        r"x_i x_{i+1}=x_{i+3}x_{i+4}",
        r"{51+34\sqrt{17}\over256}",
        "Lemma 4.1 (unique-large-conjugate primitivity)",
        r"has \([k:F]\) distinct extensions",
        r"&\le|\operatorname {Disc}(\mathbf Z[p_i])|",
        "Theorem 5.1 (compact right-hexagon field theorem)",
        "does not cover non-right hexagons",
        "short-dotted-pair statement",
    ):
        assert phrase in text


def test_no_hidden_maximality_or_classification_is_used():
    combined = " ".join(
        (REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(encoding="utf-8")).split()
    )
    assert "No imposed height, power coupling, Coxeter classification" in combined
    assert "Saturation and maximality are not needed" in combined
    assert '"uses_reflective_lattice_classification": true' not in combined


def test_invalid_alternating_cosh_is_rejected():
    with pytest.raises(ValueError):
        right_hexagon_from_alternating_cosh(1, 2, 3)
