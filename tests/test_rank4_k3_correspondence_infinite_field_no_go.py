"""Regression tests for the rank-four K3 infinite-field scoped no-go."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_rank4_k3_correspondence_infinite_field_no_go import (
    add,
    b_negative,
    b_zero,
    build_certificate,
    companion,
    determinant,
    digest,
    embedding_record,
    field_record,
    gram,
    identity,
    multiply,
    operator_record,
    scale,
    squarefree_decomposition,
    tau,
    transpose,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rank4_k3_correspondence_infinite_field_no_go.md"
ARCHIVE = (
    ROOT
    / "data/enumerations/rank4_k3_correspondence_infinite_field_no_go.json"
)
CERTIFIER = (
    ROOT
    / "scripts/certify_rank4_k3_correspondence_infinite_field_no_go.py"
)


def zero(rank: int):
    return tuple(tuple(0 for _ in range(rank)) for _ in range(rank))


def test_companion_action_is_self_adjoint_and_quadratic_for_a_long_range():
    for a in range(2, 401):
        c = companion(a)
        assert multiply(transpose(c), b_zero()) == multiply(b_zero(), c)
        assert multiply(transpose(c), b_negative(a)) == multiply(
            b_negative(a), c
        )
        operator = tau(a)
        g = gram(a)
        assert multiply(transpose(operator), g) == multiply(g, operator)
        assert (
            add(
                add(multiply(operator, operator), scale(-2 * a, operator)),
                identity(4),
            )
            == zero(4)
        )


def test_lattice_determinant_signature_criteria_and_eigenline_signs():
    for a in range(2, 301):
        ba = b_negative(a)
        assert ba[0][0] == 1 - 2 * a < 0
        assert determinant(ba) == 2 * (a - 1) > 0
        assert determinant(gram(a)) == -32 * (a - 1)
        # For the two roots lambda of x^2-2a*x+1, their product is one.
        # The large root is >1 and the small root lies in (0,1), so the
        # B0 eigenline square lambda^2-1 has the required opposite signs.
        assert a * a - 1 > 0
        row = operator_record(a)
        assert row["lattice_signature"] == [1, 3]
        assert row["t_eigenspace_signature"] == [1, 1]
        assert row["conjugate_eigenspace_signature"] == [0, 2]


def test_explicit_K3_embedding_is_isometric_and_primitive_uniformly():
    for a in range(2, 301):
        row = embedding_record(a)
        assert row["pullback_gram"] == [list(entries) for entries in gram(a)]
        assert abs(row["primitive_minor_determinant"]) == 1
        assert row["primitive_embedding"] is True
        assert row["K3_orthogonal_complement_signature"] == [2, 16]


def test_squarefree_field_parameter_and_Pell_identity():
    for a in range(2, 1001):
        row = field_record(a)
        d = row["squarefree_field_parameter"]
        b = row["Pell_y"]
        assert a * a - d * b * b == 1
        assert d > 1
        for prime in range(2, isqrt_for_test(d) + 1):
            assert d % (prime * prime) != 0


def isqrt_for_test(number: int) -> int:
    # Deliberately independent of the producer's factorization loop.
    candidate = 0
    while (candidate + 1) * (candidate + 1) <= number:
        candidate += 1
    return candidate


def test_finite_sample_exhibits_field_repetition_but_rapidly_growing_image():
    rows = [field_record(a) for a in range(2, 257)]
    fields = {row["squarefree_field_parameter"] for row in rows}
    assert len(fields) > 100
    # Repetition really occurs, so the report does not confuse parameters
    # with fields: a=2 and a=7 both generate Q(sqrt(3)).
    assert field_record(2)["squarefree_field_parameter"] == 3
    assert field_record(7)["squarefree_field_parameter"] == 3


def test_archive_is_an_exact_replay_with_a_digest():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    expected = build_certificate(max_a=128)
    assert envelope["payload"] == expected
    assert envelope["payload_sha256"] == digest(expected)
    assert envelope["payload"]["distinct_sample_field_count"] > 50


def test_certificate_keeps_the_strong_field_gate_open():
    payload = build_certificate(max_a=32)
    family = payload["family"]
    assert family["NS_rank"] == 4
    assert family["primitive_full_trace_block"] is True
    assert family["integral_self_adjoint_operator"] is True
    assert family["algebraic_cycle_coefficients"].startswith("rational")
    scope = payload["scope"]
    assert scope["exact_complete_restricted_root_component_returned"] is True
    assert scope["prescribed_input_chamber_returned"] is False
    assert scope["restricted_component_Coxeter_claimed"] is False
    assert scope["restricted_component_finite_volume_claimed"] is False
    assert scope["reflection_saturation_claimed"] is False
    assert scope["reflection_maximality_claimed"] is False
    assert scope["exact_maximal_same_P_infinite_field_family_claimed"] is False
    assert scope["strong_exact_maximal_field_set_gate"] == "open"


def test_report_states_the_exact_construction_and_nonclaims():
    text = REPORT.read_text(encoding="utf-8")
    normalized = " ".join(text.split())
    required = [
        r"G_a=2B_0\oplus2B_a",
        r"\det G_a=-32(a-1)",
        r"\tau_a^2-2a\tau_a+I=0",
        "literal full trace package",
        "primitive isometric embedding",
        r"\operatorname {NS}(X_a)\simeq L_a",
        "rational algebraic correspondence",
        "complete restricted root arrangement",
        "finite-union contradiction",
        "Hermite--Minkowski",
        "strongest exact-maximal field-set gate remains open",
    ]
    for phrase in required:
        assert phrase in normalized
    exclusions = [
        "does **not** prove that the restricted eigen-nef component is finite-volume",
        "that the restricted component is a Coxeter polytope",
        "that its facet group is saturated",
        "that it is reflection-maximal",
        "that it returns a prescribed input chamber",
    ]
    for phrase in exclusions:
        assert phrase in normalized


def test_no_reflective_classification_is_imported():
    combined = REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(
        encoding="utf-8"
    )
    for forbidden in ("Nikulin", "Vinberg", "ABSW"):
        assert forbidden not in combined


def test_invalid_parameters_are_rejected():
    with pytest.raises(ValueError):
        companion(1)
    with pytest.raises(ValueError):
        b_negative(0)
    with pytest.raises(ValueError):
        embedding_record(1)
    with pytest.raises(ValueError):
        field_record(1)
    with pytest.raises(ValueError):
        squarefree_decomposition(0)
    with pytest.raises(ValueError):
        build_certificate(max_a=1)
