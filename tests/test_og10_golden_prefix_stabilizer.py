from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_og10_golden_prefix_stabilizer.py"
ARTIFACT = ROOT / "data/candidates/og10_golden_prefix_stabilizer.json"


def _module():
    spec = importlib.util.spec_from_file_location("og10_prefix_stabilizer", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_frozen_prefix_stabilizer_certificate() -> None:
    module = _module()
    fresh = module.build()
    frozen = json.loads(ARTIFACT.read_text())
    assert fresh == frozen
    assert fresh["verified"] is True
    assert [record["wall_count"] for record in fresh["prefixes"]] == [171, 321, 457]
    assert all(
        record["oriented_wall_set_stabilizer"] == "identity"
        for record in fresh["prefixes"]
    )
    assert fresh["infinite_order_candidate"]["stabilizes_457_wall_prefix"] is False
