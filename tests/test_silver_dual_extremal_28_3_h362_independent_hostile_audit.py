import importlib.util
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/independent_hostile_audit_silver_dual_extremal_28_3_h362_replacement.py"
ARTIFACT = ROOT / "data/enumerations/silver_dual_extremal_28_3_h362_replacement_independent_hostile_audit.json"


def payload():
    return json.loads(ARTIFACT.read_text())


def test_hostile_replay_is_owner_import_free() -> None:
    source = SCRIPT.read_text()
    assert "from coxeter_polytopes.search.silver_dual_extremal_28_3 import" not in source
    assert "import coxeter_polytopes.search.silver_dual_extremal_28_3" not in source
    completed = subprocess.run([sys.executable, str(SCRIPT), "--verify"], cwd=ROOT, text=True, capture_output=True, check=True)
    assert "PASS" in completed.stdout


def test_hostile_replay_exact_conclusions() -> None:
    p = payload()
    assert p["verdict"] == "PASS"
    assert p["control_count"] == 103
    assert p["height_count"] == 70
    assert p["first_threat"]["gram"] == [362, 182, 256]
    assert p["first_threat"]["maximum"] == "53 + 37*sqrt(2)"
    assert p["first_threat"]["old_witness_pairing"] == "1/2 + 1/2*sqrt(2)"
    assert p["replacement"]["scale"] == "12/25"
    assert p["replacement"]["strict"] is True
