"""Regression tests for the affine-B-tail global completion obstruction."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_rational_affine_b_tail_global_completion_obstruction import (
    DEFAULT_OUTPUT,
    build_certificate,
    digest,
    pair,
    truncator_certificate,
    visible_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "rational_affine_b_tail_global_completion_obstruction.md"


def test_visible_simplex_has_exact_dimension_threshold():
    records = [visible_certificate(length) for length in range(2, 31)]
    assert [
        record["m"]
        for record in records
        if record["finite_volume_visible_simplex"]
    ] == [2, 3, 4]
    for record in records:
        length = record["m"]
        assert record["wall_count"] == length + 2
        assert abs(record["gram_determinant"]) == 1
        assert record["dual_squares"][-1] == 4 - length


def test_first_hyperideal_escape_and_truncators_are_exact():
    v5 = (1, 2, -1, -1, -1, -1, -1)
    assert pair(v5, v5) == -1

    t5 = truncator_certificate(5)
    assert t5 == [
        {
            "label": "T_5",
            "vector": list(v5),
            "support_length": 5,
            "square": -1,
            "canonical_pairing": -1,
            "arithmetic_genus": 0,
            "pairing_with_first_escape_v5": -1,
            "status": (
                "effective divisor class; irreducibility, smoothness, "
                "and facet status are not automatic"
            ),
        }
    ]

    records = truncator_certificate(30)
    assert records[0]["label"] == "T_5"
    assert records[0]["vector"] == list(v5) + [0] * 25
    root_records = records[1:]
    assert [record["label"] for record in root_records] == [
        "R_1",
        "R_2",
        "R_3",
        "R_4",
        "R_5",
        "R_6",
        "R_7",
    ]
    for index, record in enumerate(root_records, start=1):
        assert record["square"] == -2
        assert record["canonical_pairing"] == 2 * index - 2
        assert record["arithmetic_genus"] == index - 1
        assert record["support_length"] == 4 * index + 2
        assert record["minimum_jet_row_dependencies_for_effectivity"] == index
        assert (
            record[
                "irreducible_exact_multiplicity_locus_codimension_lower_bound"
            ]
            == index
        )


def test_ruling_swap_archive_and_claim_boundary():
    archive = json.loads(DEFAULT_OUTPUT.read_text())
    payload = build_certificate()
    assert archive["payload"] == payload
    assert archive["payload_sha256"] == digest(payload)
    assert payload["absolute_reflection_maximality_on_full_NS"] is False
    assert payload["obstructing_root"] == "c-d"
    for sample in payload["swap_reflection_samples"]:
        assert sample["surface_square"] == -2
        assert sample["doubled_host_square"] == -4
        assert sample["reflection_integral"] is True
        assert sample["reflection_swaps_c_and_d"] is True
        assert sample["pullback_c_plus_d_pairing"] == 0


def test_affine_h3_semiample_index_and_lattice_no_go():
    gate = build_certificate()["affine_h3_index_gate"]
    assert gate["physical_null_coefficients"] == [4, 1, 1, 1, 1]
    assert gate["canonical_pairing"] == 4
    assert gate["arithmetic_genus"] == 3
    assert gate["semiample_index_over_C"] == 1
    assert gate["relatively_minimal_core_picard_bound"] == 26
    assert gate["genus_two_branch"] is False
    samples = gate["lattice_saturation_countermodel"]
    assert [sample["N"] for sample in samples] == [1, 2, 3, 5, 11, 37]
    for sample in samples:
        assert sample["ambient_discriminant_exponent"] == 2
        assert sample["affine_gram_independent_of_N"] is True
        assert sample["affine_null_relation"] == [sample["N"], 0, 0, 0, 0, 0]


def test_report_states_uniform_result_and_exclusions():
    text = REPORT.read_text()
    required = [
        "Uniform completion obstruction",
        r"For every \(m\ge2\)",
        "c-d",
        "proper discrete reflection overgroup",
        r"finite volume exactly for \(m\le4\)",
        "does not say that these surfaces cannot have rational-polyhedral",
        r"proper sublattice chosen to omit \(c-d\)",
        "positive-genus wall systems",
        "semiampleness has index one",
        "no genus-two/\\(18\\) branch",
        "exponent-two ambient discriminant does not bound",
        "before any appeal to a classification theorem",
        "does **not** assert that the member is irreducible",
        "does not exclude the higher-height square-minus-one candidate",
        "not asserted for arbitrary reducible",
    ]
    for phrase in required:
        assert phrase in text
    assert "not the seed of an unbounded reflection-maximal" in text
