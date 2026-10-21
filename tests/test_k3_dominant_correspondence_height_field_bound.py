"""Regression tests for the dominant-correspondence field-height gate."""

from __future__ import annotations

from fractions import Fraction
import json
from math import comb
from pathlib import Path

import pytest

from scripts.certify_k3_dominant_correspondence_height_field_bound import (
    build_certificate,
    coefficient_bounds,
    digest,
    encoded_certificate,
    k3_rank_record,
    polynomial_search_count,
    primitive_combination_bound,
    quadratic_escape_record,
    spectral_identity_record,
    squarefree_part,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_dominant_correspondence_height_field_bound.md"
CERTIFIER = (
    ROOT / "scripts/certify_k3_dominant_correspondence_height_field_bound.py"
)
ARCHIVE = (
    ROOT / "data/enumerations/k3_dominant_correspondence_height_field_bound.json"
)
AUDIT = (
    ROOT
    / "results/k3_dominant_correspondence_height_field_bound_hostile_scope_audit.md"
)


def test_signed_place_identity_is_exact_on_a_large_fraction_grid():
    for degree in range(2, 9):
        for square in (2, 7, 14):
            for scale in (1, 3, 11):
                weights = tuple(
                    Fraction(-scale * index, index + 1)
                    for index in range(1, degree)
                )
                distinguished = Fraction(5 * degree + scale, 2)
                eigenvalues = (distinguished,) + tuple(
                    distinguished - Fraction(index, degree + 1)
                    for index in range(1, degree)
                )
                record = spectral_identity_record(
                    square, weights, eigenvalues
                )
                assert record["dominant_spectral_inequality"] is True
                assert Fraction(record["spectral_gap"]) >= 0
                assert record["spectral_gap"] == record["expanded_gap"]


def test_zero_negative_projection_gives_equality():
    record = spectral_identity_record(14, (0, 0, 0), (31, 9, 2, 1))
    assert record["mixed_height"] == str(14 * 31)
    assert record["spectral_gap"] == "0"


def test_primitive_combination_grid_has_strict_room():
    for degree in range(2, 13):
        for count in range(1, 6):
            record = primitive_combination_bound(degree, count)
            assert record["good_grid_point_exists"] is True
            assert record["bad_grid_union_bound"] < record["grid_size"]
            assert record["positive_coefficient_bound"] == (
                degree * (degree - 1) // 2 + 1
            )


def test_coefficient_box_is_the_elementary_symmetric_bound():
    for degree in range(1, 9):
        for root_bound in (Fraction(1), Fraction(7, 2), Fraction(19)):
            bounds = coefficient_bounds(degree, root_bound)
            assert bounds == [
                int(Fraction(comb(degree, index)) * root_bound**index)
                for index in range(1, degree + 1)
            ]
    assert polynomial_search_count(3, 2) > 0
    assert polynomial_search_count(4, 2) > polynomial_search_count(3, 2)


def test_rank_four_family_passes_only_with_unbounded_mixed_height():
    heights = []
    fields = set()
    for a in range(2, 513):
        record = quadratic_escape_record(a)
        assert record["polarization_square"] == 14
        assert record["mixed_correspondence_height"] == 32 * a - 10
        assert record["dominance_square_margin_numerator_over_49"] > 0
        assert (
            record["height_ratio_squared_discriminant_margin_numerator"] > 0
        )
        assert record["fixed_square_but_unbounded_mixed_height"] is True
        heights.append(record["mixed_correspondence_height"])
        fields.add(record["field_squarefree_parameter"])
    assert heights == sorted(heights)
    assert heights[-1] > 100 * heights[0]
    assert len(fields) > 200


def test_squarefree_decomposition_is_exact():
    for value in range(1, 5000):
        squarefree, square = squarefree_part(value)
        assert squarefree * square * square == value
        prime = 2
        residual = squarefree
        while prime * prime <= residual:
            assert residual % (prime * prime) != 0
            prime += 1


def test_k3_n2_record_has_degree_six_and_explicit_discriminant_bound():
    record = k3_rank_record(
        hyperbolic_dimension=2,
        generator_count=2,
        correspondence_height_bound=32,
    )
    assert record["field_degree_bound"] == 6
    assert record["positive_square_lower_bound"] == 2
    assert record["individual_conjugate_bound"] == "16"
    assert record["positive_combination_coefficient_bound"] == 16
    assert record["primitive_generator_conjugate_bound"] == "512"
    assert int(record["field_discriminant_bound"]) == 512**30
    assert record["finite_monic_polynomial_search_count"] > 0


def test_archive_replays_exactly_and_scope_is_literal():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate(max_a=128)
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate(128)
    conclusions = payload["conclusions"]
    assert conclusions["one_positive_place_spectral_inequality"] is True
    assert (
        conclusions[
            "bounded_rank_and_generator_packet_give_finite_field_set"
        ]
        is True
    )
    assert conclusions["bounded_polarization_square_alone_suffices"] is False
    assert conclusions["finite_field_set_not_only_bounded_degree"] is True
    scope = payload["scope"]
    assert scope["requires_exact_full_Picard_or_positive_H_in_trace_block"]
    assert scope["arbitrary_field_generators_assumed_dominant"] is False
    assert scope["pair_coordinate_integrality_from_diagram_alone"] is False
    assert scope["proves_maximality_forces_the_height_gate"] is False
    assert scope["constructs_nonorthogonal_rank_six_counterfamily"] is False
    assert scope["uses_reflective_lattice_classification"] is False


def test_report_states_the_theorem_and_the_open_implication():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "Theorem 1.1 (dominant spectral height)",
        r"H\cdot\tau _\theta(H)",
        r"\sigma _0(\theta)H^2",
        "Lemma 2.1 (positive primitive combination)",
        "Theorem 3.1 (pair-generated K3 field gate)",
        "exact full-Picard trace",
        r"\operatorname {NS}(X)=N",
        "an embedding of a field is injective",
        r"d\le6",
        r"\bigl(32\max\{4,C/s_0\}\bigr)^{30}",
        r"Z_\theta\in CH^2(X\times X)_{\mathbf Q}",
        r"H_a\cdot\tau_{t_a}(H_a)=32a-10",
        "does not prove that reflection maximality supplies the bound",
        "finite field set, not merely bounded field degree",
        "an ample or big-and-nef class need not lie in",
        "still-open global implication",
    ]
    for phrase in required:
        assert phrase in text


def test_report_and_certifier_import_no_reflective_classification():
    combined = REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(
        encoding="utf-8"
    )
    assert "uses_reflective_lattice_classification\": True" not in combined
    assert "No dimension or field" in combined
    assert "finiteness theorem for reflection groups" in combined


def test_hostile_audit_freezes_every_scope_boundary():
    text = " ".join(AUDIT.read_text(encoding="utf-8").split())
    required = [
        "PASS after scope repairs",
        r"\operatorname {NS}(X)=N",
        "Dominant generators are not arbitrary generators",
        "The theorem therefore requires the field to be pair-generated",
        "Cartan integrality and zero conjugates",
        r"2^r\le d",
        r"R^{d(d-1)}",
        "not merely a degree bound",
        r"(32R_0)^{30}",
        r"H\cdot\tau_t(H)=32a-10",
        "derive the scalar height bound from reflection maximality",
    ]
    for phrase in required:
        assert phrase in text


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        spectral_identity_record(0, (-1,), (2, 1))
    with pytest.raises(ValueError):
        spectral_identity_record(2, (1,), (2, 1))
    with pytest.raises(ValueError):
        spectral_identity_record(2, (-1,), (1, 2))
    with pytest.raises(ValueError):
        primitive_combination_bound(1, 1)
    with pytest.raises(ValueError):
        coefficient_bounds(2, Fraction(1, 2))
    with pytest.raises(ValueError):
        quadratic_escape_record(1)
    with pytest.raises(ValueError):
        k3_rank_record(1, 1, 10)
    with pytest.raises(ValueError):
        build_certificate(1)
