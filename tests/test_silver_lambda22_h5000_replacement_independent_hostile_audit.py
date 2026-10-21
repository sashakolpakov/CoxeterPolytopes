"""Regression tests for the clean-room silver Lambda22 H5000 audit."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/silver_lambda22_h5000_replacement_independent_hostile_audit.json"
)
SOURCE = ROOT / "scripts/independent_hostile_audit_silver_lambda22_h5000_replacement.py"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def sign(value: dict[str, str]) -> int:
    a = Fraction(value["a"])
    b = Fraction(value["b"])
    if not a:
        return (b > 0) - (b < 0)
    if not b or (a > 0) == (b > 0):
        return (a > 0) - (a < 0)
    return ((a if a * a - 2 * b * b > 0 else b) > 0) * 2 - 1


def test_clean_room_frontier_counts_and_cutoff_semantics() -> None:
    frozen = payload()
    interval = frozen["new_interval"]
    assert frozen["verdict"] == "PASS"
    assert interval["control_count"] == 69
    assert interval["retained_type_count"] == 635
    assert len(interval["threats"]) == 7
    assert interval["controls_by_norm"] == {
        "1": 8,
        "2": 20,
        "3": 4,
        "4": 8,
        "4 + 2*sqrt(2)": 14,
        "6": 6,
        "12": 4,
        "12 + 6*sqrt(2)": 5,
    }
    assert "strictly greater" in frozen["cutoff_semantics"]["lower"]
    assert "less than or equal" in frozen["cutoff_semantics"]["upper"]


def test_every_nonnegative_extension_is_unique_and_short_vector_impossible() -> None:
    frozen = payload()
    assert frozen["logical_check"] == {
        "all_nonnegative_extensions_enumerated": True,
        "all_nonnegative_extensions_short_vector_impossible": True,
        "nonnegative_extensions_per_threat": [1] * 7,
    }
    for row in frozen["new_interval"]["threats"]:
        assert row["nonnegative_extension_count"] == 1
        extension = row["nonnegative_extensions"][0]
        assert sign(extension["margin"]) >= 0
        U, V, M = row["gram"]
        p, q = extension["d_dot_u"], extension["d_dot_v"]
        resolution = extension["resolution"]
        k = resolution["k"]
        square = U + V - 2 * M + 4 * k * k - 2 * k * (p - q)
        assert square == resolution["square"]
        assert 0 < square < resolution["certified_minimum"]
        if resolution["lattice"] == "S2":
            assert k % 2 == 0


def test_audit_is_structurally_independent_and_checks_continuity() -> None:
    source = SOURCE.read_text()
    assert "from coxeter_polytopes" not in source
    assert "def enumerate_controls" in source
    assert "def candidate_grams" in source
    assert "def extensions" in source
    assert "def continuity_audit" in source
    continuity = payload()["continuity"]
    assert continuity["prior_ray_strict_through_H1598"] is True
    assert continuity["ambient_minimum"] == 4
    assert continuity["S2_minimum"] == 12
