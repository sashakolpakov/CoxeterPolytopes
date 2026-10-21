import copy
import json
from pathlib import Path

import pytest

from scripts.certify_silver_alexeev_h9_local_star import verify


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = ROOT / "data/candidates/silver_alexeev_h9_local_star.json"


def test_frozen_silver_alexeev_local_star_replays_exactly() -> None:
    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))

    result = verify(manifest)

    assert result["local_group_order"] == 12
    assert result["full_sector_count"] == 6
    assert result["nonfull_sector_count"] == 6
    assert result["angle_histogram"] == [
        {"cosine_squared": {"a": "1/4", "b": "1/8"}, "count": 4},
        {"cosine_squared": {"a": "3/8", "b": "-1/4"}, "count": 2},
    ]
    assert all(manifest["checks"].values())


def test_silver_alexeev_local_star_rejects_mutated_sector() -> None:
    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    mutated = copy.deepcopy(manifest)
    mutated["sectors"][0]["full_dimensional_silver_sector"] = False

    with pytest.raises(ValueError, match="sectors"):
        verify(mutated)
