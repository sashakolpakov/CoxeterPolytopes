"""Tests for unconditional maximal arithmetic reflection chamber finiteness."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from scripts.certify_unconditional_maximal_arithmetic_reflection_chamber_finiteness import (
    ARCHIVE,
    BLOWUP_MAX,
    CONJUGATE_HEIGHT,
    FIELD_DEGREE_MAX,
    FRAME_CONSTANT,
    RANK_MAX,
    REPORT,
    build_certificate,
    coefficient_bounds,
    digest,
    encoded_certificate,
    simple_directed_cycle_slot_bound,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = (
    ROOT
    / "scripts/certify_unconditional_maximal_arithmetic_reflection_chamber_finiteness.py"
)


def test_terminal_scope_is_unconditional_and_global() -> None:
    theorem = build_certificate()["terminal_theorem"]
    assert theorem["minimum_dimension"] == 2
    assert theorem["maximum_dimension"] == 29
    assert theorem["maximum_ground_field_degree"] == 25
    assert "noncompact_ground_field" not in theorem
    assert "noncompact_dimensions" not in theorem
    assert theorem["finite_native_pair_set"] is True
    assert theorem["finite_unmarked_normalized_chambers"] is True
    assert theorem["AG_carrier_hypothesis"] is False


def test_cycle_height_and_coefficient_box_are_finite() -> None:
    assert RANK_MAX == 30
    assert FIELD_DEGREE_MAX == 25
    assert FRAME_CONSTANT == 62
    assert CONJUGATE_HEIGHT == 62**30
    assert simple_directed_cycle_slot_bound(30) > 0
    bounds = coefficient_bounds(25, 62**30)
    assert len(bounds) == 25
    assert bounds[0] == 25 * 62**30
    assert all(value > 0 for value in bounds)


def test_finite_pair_route_does_not_use_absw_conclusion() -> None:
    routes = build_certificate()["two_terminal_routes"]
    assert routes["uses_ABSW_finiteness_conclusion_in_factorization"] is False
    assert routes["finite_pair_factorization"][-2:] == [
        "reflection maximality W=R(L)",
        "one complete-mirror chamber orbit",
    ]


def test_universal_weak_rational_surface_bound() -> None:
    corollary = build_certificate()["weak_rational_surface_corollary"]
    assert BLOWUP_MAX == 4 * 25 * 30 - 3 == 2997
    assert corollary["maximum_blowups"] == 2997
    assert corollary["requires_reflection_maximality"] is False
    assert corollary["effective_or_nef_return_claimed"] is False
    assert corollary["reflection_action_return_claimed"] is False


def test_report_keeps_imports_and_nonclaims_explicit() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for needle in (
        "unconditional global finiteness theorem outright",
        "does not invoke the ABSW finiteness conclusion",
        "published imported input",
        "not an AG-relative result",
        "Any weak-action upgrade is a",
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
    )
