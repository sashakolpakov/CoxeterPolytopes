import json
from pathlib import Path

import pytest

from scripts.certify_silver_alexeev_h9_folding_obstruction import verify


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/candidates/silver_alexeev_h9_folding_obstruction.json"


def test_frozen_silver_alexeev_folding_obstruction_replays_exactly() -> None:
    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))

    result = verify(manifest)

    assert result["framing_steps"] == 108
    assert result["ridge_zero_walls"] == [2, 4, 5]
    assert result["cosine_squared"] == {"a": "1/4", "b": "1/8"}
    assert all(manifest["checks"].values())


def test_silver_alexeev_folding_obstruction_rejects_mutated_witness() -> None:
    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    manifest["ridge_witness_square"] = {"a": "0", "b": "0"}

    with pytest.raises(ValueError, match="ridge_witness_square"):
        verify(manifest)
