"""Regression checks for the exact silver-Lambda22 H5000 frontier."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h5000_replacement.json"


def test_h5000_artifact_has_complete_scoped_frontier() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["new_control_count"] == 69
    assert payload["retained_gate_compatible_type_count"] == 635
    assert payload["nonnegative_relaxed_row_count"] == 7
    assert len(payload["rows"]) == 7
    assert all(row["resolution"] for row in payload["rows"])
    assert payload["discovery"] is None
    assert "bounded prefix only" in payload["scope"]


def test_replacement_is_hyperideal_and_every_pinned_threat_is_short() -> None:
    payload = json.loads(ARTIFACT.read_text())
    norm = Silver.from_json(payload["replacement_norm"]["exact"])
    assert norm == Silver(Fraction(1361, 961), -1)
    assert norm.sign(0) > 0

    # For q=2 rows, the stated combinations u-v-kd have square two.
    q2_k = (6, 6, 7, 7)
    q2_rows = [row for row in payload["rows"] if row["norm"]["display"] == "2"]
    assert len(q2_rows) == 4
    for row, k in zip(q2_rows, q2_k, strict=True):
        U, V, C = row["gram"]
        p, q = row["maximizing_d_pairings"]
        assert U + V - 2 * C + 4 * k * k - 2 * k * (p - q) == 2

    # For q=4 rows, u-v-kd lies in S2 and has square four.
    q4_k = (8, 8, 10)
    q4_rows = [row for row in payload["rows"] if row["norm"]["display"] == "4"]
    assert len(q4_rows) == 3
    for row, k in zip(q4_rows, q4_k, strict=True):
        U, V, C = row["gram"]
        p, q = row["maximizing_d_pairings"]
        assert U + V - 2 * C + 4 * k * k - 2 * k * (p - q) == 4
        assert k % 2 == 0


def test_h5000_audit_source_retains_exhaustive_scan() -> None:
    source = (ROOT / "scripts/audit_silver_lambda22_h5000_replacement.py").read_text()
    assert "enumerate_lambda22_controls(BOUND)" in source
    assert "candidate_leech_gram_triples(level.target)" in source
    assert "gated_relaxed_witness_bound" in source
    assert "Vinberg termination" in source
