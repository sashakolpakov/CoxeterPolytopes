from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_nikulin_336_quasilanner_counterexample.py"
ARCHIVE = ROOT / "data/enumerations/nikulin_336_quasilanner_counterexample.json"
REPORT = ROOT / "results/nikulin_conceptual_dependency_hostile_audit.md"


def load_module():
    spec = importlib.util.spec_from_file_location("nikulin_336_quasilanner", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_frozen_archive_is_exact_replay() -> None:
    module = load_module()
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived == module.build_certificate()
    assert archived["all_checks_pass"] is True


def test_exact_inertia_deletions_and_affine_null_root() -> None:
    cert = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    diagram = cert["diagram"]
    assert diagram["determinant"] == -1
    assert diagram["ldl_pivots"] == ["-2", "-3/2", "-4/3", "1/4"]
    assert diagram["inertia"] == {"positive": 1, "negative": 3, "zero": 0}
    deletions = diagram["vertex_deletions"]
    assert deletions["1"]["inertia"] == {"positive": 0, "negative": 2, "zero": 1}
    assert deletions["1"]["primitive_null_vector"] == ["1", "2", "sqrt(3)"]
    assert deletions["1"]["negative_definite"] is False
    assert all(deletions[str(node)]["negative_definite"] for node in (2, 3, 4))


def test_marked_finite_edge_and_ground_field() -> None:
    cert = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert cert["marked_finite_edge"] == {
        "endpoint_exclusive_nodes": [2, 3],
        "u": 1,
        "u_less_than_2": True,
        "both_endpoint_deletions_negative_definite": True,
    }
    packet = cert["arithmetic_packet"]
    assert packet["nonzero_two_cycle_products"] == [1, 1, 3]
    assert packet["ground_field"] == "Q"
    assert cert["geometric_verdict"]["lanner"] is False
    assert cert["geometric_verdict"]["quasi_lanner"] is True
    assert cert["repair"]["downstream_field_sets_change"] is False


def test_report_pins_certificate_and_scope() -> None:
    raw = REPORT.read_bytes()
    assert all(byte not in raw for byte in (0, 8, 11, 13))
    text = raw.decode("utf-8")
    assert "certify_nikulin_336_quasilanner_counterexample.py" in text
    assert "nikulin_336_quasilanner_counterexample.json" in text
    assert "field-set conclusion is nevertheless repairable" in text
    assert "finite field list or bounded [k:Q]  -/->  bounded hyperbolic dimension" in text
