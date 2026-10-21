import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/cubic7_coxeter_todd_squared_h169_independent_hostile_audit.json"


def test_h169_hostile_artifact_records_clean_room_bounded_pass():
    payload = json.loads(ARTIFACT.read_text())
    assert payload["imports_production_search_code"] is False
    assert [row["type_count"] for row in payload["levels"]] == [4, 1, 8, 13, 16, 19, 37]
    assert payload["dangerous_shell_case_count"] == 8
    assert payload["first_realized_threat"]["type_index"] == 35
    assert payload["first_realized_threat"]["status"] == "not_chronologically_accepted"
    assert payload["replacement"]["scale"] == "11/8"
    assert payload["replacement"]["all_necessary_bounds_strict"] is True
    assert payload["terminal_verdict"] == "PASS_BOUNDED_THREAT"
