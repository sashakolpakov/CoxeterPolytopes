import json
from fractions import Fraction
from pathlib import Path

import pytest

from scripts.certify_rational_finite_vertex_discrepancy_tail_gate import (
    ARCHIVE,
    build_certificate,
    discrepancy_vector,
    exceptional_tail,
    weak_barrier_residual,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_finite_vertex_discrepancy_tail_gate.md"


def test_general_barrier_identity_on_asymmetric_safe_packet():
    species = ["B", "E", "R", "R"]
    intersections = {
        (0, 1): 1,
        (1, 2): 1,
        (2, 3): 1,
    }
    _, discrepancies = discrepancy_vector(species, intersections)
    residual = weak_barrier_residual(species, intersections)
    assert all(value >= 0 for value in residual)
    assert all(
        value >= (Fraction(-1, 2) if item == "B" else 0)
        for item, value in zip(species, discrepancies, strict=True)
    )


def test_three_hostile_prototypes_are_sharp():
    payload = build_certificate()
    rows = {row["id"]: row for row in payload["sharp_hostile_controls"]}
    assert rows["H1_BB_intersection_2"]["discrepancies"] == ["-1", "-1"]
    assert rows["H2_BR_intersection_2"]["discrepancies"] == ["-1", "-1"]
    assert rows["H3_E_three_B_leaves"]["discrepancies"] == [
        "-2",
        "-1",
        "-1",
        "-1",
    ]


def test_arbitrary_exceptional_tail_formula_and_input_rejection():
    for length in (1, 2, 3, 8, 32, 64):
        row = exceptional_tail(length)
        assert row["leading_principal_determinants"] == [
            str(((-1) ** size) * (size + 1))
            for size in range(1, length)
        ] + [str((-1) ** length)]
        assert row["determinant"] == str((-1) ** length)
        assert row["discrepancies"] == [str(i) for i in range(1, length + 1)]
        assert row["canonical_square_loss_before_contraction"] == length
        assert row["discrepancy_correction"] == length
    for bad in (0, -1, True, 1.5, "3"):
        with pytest.raises(ValueError):
            exceptional_tail(bad)  # type: ignore[arg-type]


def test_archive_and_report_scope():
    payload = build_certificate()
    assert json.loads(ARCHIVE.read_text(encoding="utf-8")) == payload
    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)
    text = raw.decode("utf-8")
    assert "a_B\\geq-\\frac12" in text
    assert "epsilon<\\frac12" in text
    assert "does not bound" in text
    assert "reflection maximality" in text
    assert "H1--H3" in text
