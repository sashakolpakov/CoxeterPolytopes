"""Regression tests for the direct full-NS minimal-model reduction."""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.certify_rational_direct_full_ns_minimal_model_reduction import (
    adapted_section_exponent_two_necessary,
    adapted_swap_parity,
    build_certificate,
    f0_swap_record,
    section_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "rational_direct_full_ns_minimal_model_reduction.md"


def test_every_blowup_of_F0_retains_the_ruling_swap_mirror():
    for exceptional_rank in range(0, 51):
        row = f0_swap_record(exceptional_rank)
        assert row["root_square"] == -2
        assert row["c_maps_to_d"] is True
        assert row["d_maps_to_c"] is True
        assert row["all_total_exceptional_coordinates_fixed"] is True
        assert row["c_plus_d_pairing"] == 0
        assert row["mirror_contains_an_ample_class"] is True


def test_negative_section_integrality_is_exactly_e_plus_hits_divides_two():
    for e in range(0, 20):
        for hits in range(0, 12):
            row = section_record(e, hits)
            if e + hits == 0:
                assert row["reflection_integral_on_full_NS"] is None
                continue
            assert row["section_square"] == -(e + hits)
            assert row["fibre_pairing"] == 1
            assert row["reflection_coefficient_on_fibre"] == str(
                __import__("fractions").Fraction(2, e + hits)
            )
            assert row["reflection_integral_on_full_NS"] is (
                (e + hits) in (1, 2)
            )


def test_minimal_model_and_F2_center_conclusions():
    certificate = build_certificate()
    assert certificate["allowed_minimal_models"] == ["P2", "F2"]
    assert certificate["F2_section_hit_count"] == 0
    assert certificate["F0_excluded_by_ruling_swap_through_ample_cone"] is True
    assert certificate["F_e_e_gt_2_excluded_by_section_reflection_integrality"] is True
    assert certificate["F1_is_not_minimal"] is True
    assert certificate["adapted_host_theorem_claimed"] is False
    assert certificate["vertical_depth_bound_claimed"] is False
    assert certificate["reflective_classification_used"] is False


def test_adapted_ruling_swap_parity_is_not_silently_dropped():
    assert adapted_swap_parity(())["reflection_preserves_adapted_host"] is True
    assert adapted_swap_parity((0, 2, -4))[
        "reflection_preserves_adapted_host"
    ] is True
    assert adapted_swap_parity((1,))[
        "reflection_preserves_adapted_host"
    ] is False


def test_adapted_section_exponent_two_condition_is_only_necessary():
    for e in range(1, 12):
        for hits in range(0, 8):
            row = adapted_section_exponent_two_necessary(e, hits)
            assert row["can_preserve_exponent_two_host_only_if_glue_allows"] is (
                e + hits in (1, 2, 4)
            )
            assert row["condition_is_sufficient"] is False
    assert adapted_swap_parity((2, 3))[
        "reflection_preserves_adapted_host"
    ] is False


def test_report_states_scope_and_surviving_gate():
    text = REPORT.read_text()
    required = [
        r"\boxed{\mathbf P^2\quad\hbox{or}\quad\mathbb F_2.}",
        "arbitrary\nsequence of blowups",
        "mirror crosses its interior",
        r"e+t\mid2",
        "no centre of the chosen blowdown",
        "direct/full-lattice theorem",
        "No adapted-host\nversion",
        "Odd boundary parity",
        r"\boxed{n\mid4.}",
        "necessary, not sufficient",
        "not yet a depth bound",
        "plane model or from \\(\\mathbb F_2\\)",
    ]
    for phrase in required:
        assert phrase in text


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        section_record(-1, 0)
    with pytest.raises(ValueError):
        section_record(2, -1)
    with pytest.raises(ValueError):
        f0_swap_record(-1)
    with pytest.raises(ValueError):
        adapted_section_exponent_two_necessary(-1, 0)
