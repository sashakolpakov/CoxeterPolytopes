"""Regression tests for the independent golden triangle packet replay."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/golden_agcat_triangle_packets_independent_hostile_audit.json"
)
SOURCE = ROOT / "scripts/independent_hostile_audit_golden_agcat_triangle_packets.py"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_independent_counts_and_stream_hashes() -> None:
    frozen = payload()
    replay = frozen["replay"]
    assert frozen["verdict"] == "PASS"
    assert replay["alphabet_count"] == 6_875
    assert replay["nonzero_count"] == 6_874
    assert replay["bucket_count"] == 3_661
    assert replay["maximum_bucket_size"] == 109
    assert replay["compatible_bucket_triples"] == 27_088
    assert replay["raw_unordered_nonzero_multisets"] == 54_158_527_500
    assert replay["formula_square_product_packets"] == 4_768_141
    assert replay["square_product_packets"] == 4_768_141
    assert replay["conjugate_positive_packets"] == 2_633_885
    assert replay["conjugate_nonpositive_packets"] == 2_134_256
    assert replay["canonical_sha256"] == (
        "41386c1bfb64c497acc2ddb967c9a0d646f82e67dfc8907b8bc7d3b50e85bbe3"
    )
    assert replay["retained_sha256"] == (
        "8be12b4371f49e1a0dd35f2c290c4e5735e64da56dc209655d1e0af5b49476b6"
    )


def test_selected_sign_partition_and_raw_retained_zero_split() -> None:
    replay = payload()["replay"]
    assert replay["retained_selected_positive"] == 0
    assert replay["retained_selected_lorentzian"] == 2_633_885
    assert replay["raw_selected_degenerate"] == 1
    assert replay["retained_selected_degenerate"] == 0
    assert replay["zero_records"] == [
        {
            "edges": [[1, 0], [1, 0], [1, 0]],
            "cycle": [-1, 0],
            "determinant": [0, 0],
            "conjugate_sign": 0,
        }
    ]
    assert sum(row["count"] for row in replay["strata"]) == 4_768_141


def test_replay_is_core_independent_and_scope_is_local() -> None:
    source = SOURCE.read_text()
    assert "from coxeter_polytopes.classification.golden_triangle_packets" not in source
    assert "def squareclass" in source
    assert "def identity_positive_square_root" in source
    assert "def expected_multiset_count" in source
    assert "def replay" in source
    frozen = payload()
    assert "no high-rank compatibility or emptiness claim" in frozen["scope"]
    assert "remain open" in frozen["logic"]["compact_scope"]
