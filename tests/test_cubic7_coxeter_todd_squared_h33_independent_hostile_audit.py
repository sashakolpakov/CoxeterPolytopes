import importlib.util
import json
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/independent_hostile_audit_cubic7_coxeter_todd_squared_h33.py"
ARTIFACT = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h33_independent_hostile_audit.json"


@lru_cache(maxsize=1)
def load_module():
    spec = importlib.util.spec_from_file_location("cubic7_k12_h33_hostile", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_clean_room_hostile_replay_is_terminal_pass():
    module = load_module()
    payload = module.build()
    assert payload["imports_production_search_code"] is False
    assert payload["control_count"] > 0
    assert [row["types"] for row in payload["frontier"]] == [1, 1, 1, 2]
    assert payload["local_arithmetic"]["canonical_norms"] == [2, 6]
    assert payload["chronological_cross"]["other_H33_type"] == {
        "maximum": {"a": "-2", "b": "-8", "c": "-4"},
        "relaxed_schur_count": 1225,
    }
    assert payload["chronological_cross"]["same_cutter_type"] == {
        "schur_count": 8173,
        "minimum_four_count": 6960,
        "maximum": {"a": "0", "b": "-4", "c": "-2"},
    }
    assert payload["hostile_finding"]["kind"] == "REPAIRED_PROOF_GAP"
    assert payload["replacement"]["scale"] == "7/5"
    assert all(
        module.C7(Fraction(value["a"]), Fraction(value["b"]), Fraction(value["c"])).sign() < 0
        for value in payload["replacement"]["strict_pairing_bounds"]
    )
    assert payload["terminal_verdict"] == "PASS"


def test_frozen_hostile_artifact_matches_replay():
    module = load_module()
    assert json.loads(ARTIFACT.read_text()) == module.build()
