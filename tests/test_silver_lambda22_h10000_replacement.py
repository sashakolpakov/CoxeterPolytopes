import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/audit_silver_lambda22_h10000_replacement.py"
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h10000_replacement.json"


def test_frozen_h10000_replacement_certificate() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["height_bound_inclusive"]["display"] == "10000"
    assert payload["replacement_norm"]["signs"][0] == 1
    assert payload["new_control_count"] == 129
    assert payload["retained_gate_compatible_type_count"] == 1555
    assert payload["nonnegative_relaxed_row_count"] == 1
    threat = payload["unique_threat"]
    assert threat["gram"] == [3194, 1600, 2260]
    assert threat["maximizing_d_pairings"] == [113, 80]
    assert threat["locally_admissible_extension_count"] == 665
    assert threat["nonnegative_extension_count"] == 1
    assert threat["all_nonnegative_extensions"][0]["extension_determinant"] == 0
    assert threat["short_vector_norm"] == 2


def test_unique_threat_reduction_replays_exactly() -> None:
    module = runpy.run_path(str(SCRIPT))
    threat = module["threat_resolution"]()
    assert threat["nonnegative_extension_count"] == 1
    assert threat["all_nonnegative_extensions"][0]["d_pairings"] == [113, 80]
    assert threat["v_Schur_numerator"] == 0
    assert threat["cross_Schur_numerator"] == 0
    assert threat["short_vector"] == "e=u-28d"
    assert threat["short_vector_norm"] == 2
