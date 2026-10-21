"""Tests for the finite/ideal-vertex quadrilateral endpoint addendum."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_k3_n2_finite_area_quadrilateral_field_bound_addendum import (
    ARCHIVE,
    build_certificate,
    diagonal_record,
    digest,
    encoded_certificate,
    endpoint_coefficient_record,
    field_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_finite_area_quadrilateral_field_bound_addendum.md"


def test_ideal_endpoint_coefficient_bounds_are_weak() -> None:
    record = endpoint_coefficient_record()
    assert record["hypercube_endpoints_checked"] == 16
    assert record["endpoint_maxima_P_Q_R_absK"] == [2, 2, 2, 4]
    assert record["weak_coefficient_bounds"] == "0<=P,Q,R<=2 and |K|<=8"
    assert record["weak_determinant_growth"].endswith("<=4U+4V+4UV+8")


def test_strict_norm_floor_recovers_same_strict_diagonal_bound() -> None:
    record = diagonal_record()
    assert record["weak_growth_plus_strict_forbidden_endpoint_recovers_strict_quadratic"] is True
    assert record["remainder_square_difference"] == 17_606_953_271_521
    assert record["both_opposite_magnitudes_strict_upper_bound"] == 32_772
    assert record["both_opposite_squared_ratios_strict_upper_bound"] == 1_074_003_984


def test_field_and_quadratic_constants_are_unchanged() -> None:
    record = field_record()
    assert record["simple_cycle_weak_conjugate_height"] == 18_455_752_922_365_956_096
    assert record["Vandermonde_difference_bound"] == 590_584_093_515_710_595_072
    assert record["Vandermonde_exponent"] == 30
    assert record["quadratic_field_discriminant_upper_bound"] == 4_857_950_597


def test_archive_replays_and_scope_remains_narrow() -> None:
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    conclusions = payload["conclusions"]
    assert conclusions["compact_quadrilateral_theorem_extended_to_finite_or_ideal_vertices"] is True
    assert conclusions["noncompact_infinite_area_or_ultraideal_cases_closed"] is False
    assert conclusions["arbitrary_n2_polygon_closed"] is False


def test_report_has_the_acute_polygon_lemma_and_strict_endpoint_recovery() -> None:
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "two angles supplementary to polygon angles",
        r"0\le a,b,c,d_0\le1",
        r"0\le P,Q,R\le2",
        "strictly decreasing for",
        "restores the strict inequality",
        r"U,V<32772",
        "does not cover ultraideal/truncated or infinite-area cases",
    ):
        assert phrase in text
