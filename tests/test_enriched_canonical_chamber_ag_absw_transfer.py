"""Regression tests for the enriched canonical-chamber AG/ABSW transfer."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_enriched_canonical_chamber_ag_absw_transfer import (
    ARCHIVE,
    DEPENDENCIES,
    REPORT,
    ROOT,
    build_certificate,
    digest,
    encoded_certificate,
    field_record,
    sha256_file,
)


def test_every_frozen_dependency_pin_matches_live_bytes():
    assert {
        relative: sha256_file(ROOT / relative)
        for relative in sorted(DEPENDENCIES)
    } == DEPENDENCIES


def test_three_dimension_branches_and_ordinary_closure_are_exact():
    branches = build_certificate()["dimension_branches"]
    assert branches["direct_rational_full_NS"] == {
        "hyperbolic_dimension_bound": 11,
        "minimal_Vinberg_field": "Q",
    }
    assert branches["weak_log_rational"] == {
        "hyperbolic_dimension_bound": 25,
        "minimal_Vinberg_field": "Q",
    }
    assert branches["complex_full_trace"]["rank_inequality"] == (
        "[k:Q]*(n+1)<=20"
    )
    assert branches["complex_full_trace"]["nonrational_dimension_bound"] == 9
    assert branches["ordinary_property4"]["exact_carriers_for_dimensions"] == (
        list(range(3, 10))
    )
    assert branches["global_dimension_bound"] == 25


def test_n2_chordless_height_constant_is_literal():
    record = field_record(2, correspondence_height=256, square_lower_bound=2)
    assert record["field_degree_bound"] == 6
    assert record["chordless_cycle_length_bound"] == 5
    assert record["selected_generator_count_bound"] == 2
    assert record["primitive_combination_coefficient_bound"] == 16
    assert record["small_affine_cycle_conjugate_bound"] == 64
    assert record["height_ratio_ceiling"] == 128
    assert record["primitive_generator_conjugate_bound"] == 4096
    assert record["field_discriminant_exponent"] == 30
    assert int(record["field_discriminant_bound"]) == 4096**30


def test_nonrational_dimension_uniform_constant_is_safe():
    certificate = build_certificate()
    coarse = 32 * 8192
    assert certificate["dimension_uniform_nonrational_formula"] == (
        "(32*max(8192,C/s0))^30"
    )
    assert int(
        certificate["dimension_uniform_C256_s2_discriminant_bound"]
    ) == coarse**30
    for dimension in range(2, 10):
        record = field_record(dimension, 256, 2)
        assert record["primitive_generator_conjugate_bound"] <= coarse
        assert record["field_discriminant_exponent"] <= 30
    assert field_record(10)["ground_field"] == "Q"
    assert field_record(19)["full_trace_K3_degree_ceiling"] == 1
    assert field_record(19)["field_degree_bound"] == 1
    assert field_record(20)["full_trace_K3_degree_ceiling"] == 0
    assert field_record(20)["field_degree_bound"] == 1
    assert field_record(25)["ground_field"] == "Q"


def test_archive_replays_byte_for_byte_and_has_no_scope_upgrade():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    assert payload["dependency_pins_match"] is True
    conclusions = payload["conclusions"]
    assert conclusions["conditional_global_dimension_at_most_25"] is True
    assert conclusions["bounded_large_cycle_bidegrees_give_finite_fields"]
    assert conclusions["bounded_complete_Hilbert_Chow_grammar_supplies_gate5"]
    assert conclusions["finite_fields_follow_from_rank_alone"] is False
    assert conclusions["finite_fields_imply_finite_integral_hosts"] is False
    assert conclusions["all_raw_maximal_arithmetic_chambers_have_carriers"] is False
    assert conclusions["maximality_forces_spectral_height"] is False


def test_report_states_all_five_gates_and_the_two_remaining_frontiers():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "**Trace fidelity.**",
        "**Chamber fidelity.**",
        "**Canonical/MMP fidelity.**",
        "**Essential surjectivity.**",
        "**Uniform arithmetic height.**",
        "finite projective and every facet reflection preserves",
        r"q=\operatorname {Tr}_{k/\mathbf Q}(\alpha h)",
        "complete negative hyperbolic facets",
        r"D=-(K+B/2)",
        "canonical-preserving trace action",
        r"d(n+1)\le20",
        r"n\le25",
        r"\ell\le n+3",
        r"\theta_C=(-1)^\ell c(C)+2^\ell",
        r"\left(32\max\{8192,C/s_0\}\right)^{30}",
        "Corollary 3.2 (bounded Hilbert/Chow grammar supplies Gate 5)",
        "complete signed marked Neron--Severi basis",
        "degree of **both** parts",
        "bound every rational denominator",
        "bounded field discriminant does not by itself bound",
        r"twice-primitive \((-4)\) lane",
        r"H\cdot\tau_t(H)=32a-10",
        "does not prove that Gate 5 is logically necessary",
        "not an arbitrary raw-polytope carrier compiler",
    ]
    for phrase in required:
        assert phrase in text


def test_report_and_certificate_pin_one_another():
    certificate = build_certificate()
    assert certificate["report_sha256"] == sha256_file(REPORT)
    assert REPORT == Path(ROOT / "results/enriched_canonical_chamber_ag_absw_transfer.md")


def test_invalid_field_record_inputs_are_rejected():
    with pytest.raises(ValueError):
        field_record(1)
    with pytest.raises(ValueError):
        field_record(2, 0, 2)
    with pytest.raises(ValueError):
        field_record(2, 10, 0)
