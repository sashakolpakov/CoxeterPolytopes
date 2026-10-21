"""Regression tests for the genus-three vertical-blowup invariant no-go."""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.certify_rational_genus3_vertical_blowup_invariant_no_go import (
    build_certificate,
    fibration_row,
    tail_row,
    tail_vectors,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "rational_genus3_vertical_blowup_invariant_no_go.md"


def test_vertical_blowups_translate_the_raw_noether_pair_exactly():
    for core_rho in (1, 6, 11, 18, 26, 53):
        zero = fibration_row(core_rho, 0)
        for blowups in (0, 1, 2, 9, 41, 200):
            row = fibration_row(core_rho, blowups)
            assert row["rho"] == zero["rho"] + blowups
            assert row["surface_K_square"] == zero["surface_K_square"] - blowups
            assert row["surface_euler"] == zero["surface_euler"] + blowups
            assert row["relative_K_square"] == zero["relative_K_square"] - blowups
            assert row["relative_euler"] == zero["relative_euler"] + blowups
            assert row["relative_noether_sum"] == 36
            assert row["hodge_degree"] == 3


def test_free_tail_has_exact_Bm_roots_colors_and_orbit_span():
    for length in range(1, 31):
        row = tail_row(length)
        assert row["surface_root_squares"] == [-2] * (length - 1) + [-1]
        assert row["host_root_squares"] == [-4] * (length - 1) + [-2]
        assert row["canonical_host_colors"] == [0] * (length - 1) + [-2]
        assert row["mixed_projection"] == [-1] + [0] * (length - 1)
        assert row["mixed_orbit_size"] == 2 * length
        assert row["mixed_orbit_rank"] == length


def test_tail_gram_is_the_signed_permutation_simple_system():
    for length in range(2, 15):
        row = tail_row(length)
        gram = row["surface_gram"]
        assert len(gram) == length
        for left in range(length):
            for right in range(length):
                if left == right:
                    expected = -1 if left == length - 1 else -2
                elif abs(left - right) == 1:
                    expected = 1
                else:
                    expected = 0
                assert gram[left][right] == expected


def test_certificate_keeps_the_claim_boundary_explicit():
    certificate = build_certificate()
    assert certificate["genus"] == 3
    assert certificate["hodge_degree"] == 3
    assert certificate["relative_noether_constant"] == 36
    assert certificate["stable_data_bounds_tail"] is False
    assert certificate["relative_minimal_data_bounds_tail"] is False
    assert certificate["raw_relative_K_nef_after_vertical_blowup"] is False
    assert certificate["finite_volume_orbit_span_bounds_tail"] is False
    assert certificate["complete_nef_counterfamily_claimed"] is False
    assert certificate["reflection_maximal_counterfamily_claimed"] is False


def test_report_states_exact_invariance_and_exclusions():
    text = REPORT.read_text()
    required = [
        "vertical blowup erasure",
        "relative canonical algebra",
        r"\((K_f^2,e_f)\)",
        r"\((-1,+1)\)",
        "applying a relatively-minimal slope or nefness inequality",
        "full signed-permutation group",
        "finite-volume attachment condition",
        "not prove that the full\nMori cone",
        "does not refute a depth theorem",
        "total marked degree",
        "equivariant contraction orbit/descent",
    ]
    for phrase in required:
        assert phrase in text


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        fibration_row(0, 1)
    with pytest.raises(ValueError):
        fibration_row(2, -1)
    with pytest.raises(ValueError):
        tail_vectors(0)
