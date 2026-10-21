"""Regression tests for compact exact full-trace K3 field finiteness."""

from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
import sys

import sympy as sp

from scripts.certify_k3_compact_narrow_facet_field_finiteness import (
    ARCHIVE,
    REPORT,
    build_certificate,
    digest,
    encoded_certificate,
    rank_rows,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = (
    ROOT / "scripts/certify_k3_compact_narrow_facet_field_finiteness.py"
)


def test_k3_rank_degree_rows_and_dimension_one_exclusion() -> None:
    rows = rank_rows()
    by_rank = {row["module_rank"]: row for row in rows}
    assert [by_rank[m]["field_degree_bound"] for m in range(3, 11)] == [
        6,
        5,
        4,
        3,
        2,
        2,
        2,
        2,
    ]
    assert [by_rank[m]["discriminant_exponent"] for m in range(3, 11)] == [
        30,
        20,
        12,
        6,
        2,
        2,
        2,
        2,
    ]
    assert all(by_rank[m]["rational_only"] for m in range(11, 21))
    assert min(by_rank) == 3


def test_narrow_facet_constant_and_sine_ratio() -> None:
    assert 4 * 2 * 2 - 2 == 14
    for j in range(1, 100):
        a = math.pi * j / 100
        for i in range(j + 1):
            b = a * i / max(1, j)
            ratio = math.sin((a + b) / 2) / math.sin(a / 2)
            assert ratio < 2


def test_cartan_cycle_coboundary_similarity() -> None:
    matrix = sp.Matrix(
        [
            [2, -2, -3],
            [-1, 2, -5],
            [-2, -7, 2],
        ]
    )
    lambdas = [sp.Rational(2), sp.Rational(3), sp.Rational(5)]
    diagonal = sp.diag(*lambdas)
    conjugate = diagonal.inv() * matrix * diagonal
    for i, j in ((0, 1), (1, 2), (2, 0)):
        assert conjugate[i, j] * conjugate[j, i] == matrix[i, j] * matrix[j, i]
    assert (
        conjugate[0, 1] * conjugate[1, 2] * conjugate[2, 0]
        == matrix[0, 1] * matrix[1, 2] * matrix[2, 0]
    )
    assert conjugate.charpoly().as_expr() == matrix.charpoly().as_expr()


def test_primitive_grid_and_uniform_discriminant_constant() -> None:
    payload = build_certificate()
    assert payload["primitive_element"]["grid_size"] == 16
    assert payload["primitive_element"]["embedding_pair_count_ceiling"] == 15
    assert payload["uniform_cycle_conjugate_height"] == 14**10
    assert payload["uniform_discriminant_base"] == 32 * 14**10
    assert payload["uniform_discriminant_exponent"] == 30
    assert payload["uniform_field_discriminant_bound_expression"] == (
        "(32*14^10)^30"
    )


def test_certificate_scope_and_nonclaims() -> None:
    payload = build_certificate()
    descent = payload["field_descent"]
    assert descent["cycle_field_equals_coefficient_field"] is True
    assert descent["separate_minimal_field_assumption_required"] is False
    scope = payload["scope"]
    assert payload["minimum_hyperbolic_dimension"] == 2
    assert scope["hyperbolic_dimension_one_excluded"] is True
    assert scope["reflection_maximality_required"] is False
    assert scope["imposed_correspondence_height_required"] is False
    conclusions = payload["conclusions"]
    assert conclusions["finite_field_set_for_compact_exact_K3_carriers"] is True
    assert conclusions["with_noncompact_cusp_theorem_all_finite_volume_closed"] is True
    assert conclusions["carrier_existence_proved"] is False
    assert conclusions["property4_essential_surjectivity_proved"] is False
    assert conclusions["finite_chamber_set_proved"] is False


def test_report_contains_exact_scope_and_global_corollary() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for needle in (
        "m=n+1\\ge3",
        "No separate minimal-ground-field hypothesis is imposed",
        "32\\cdot14^{10}",
        "every finite-volume exact full-trace",
        "does **not** apply to hyperbolic dimension one",
        "Nor does this theorem construct",
    ):
        assert needle in text


def test_archive_is_exact_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
