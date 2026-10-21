from __future__ import annotations

from pathlib import Path

from scripts.independent_h9_14_record12_coxeter_probe import (
    EXPECTED_A0_SUMMARY_SHA256,
    EXPECTED_ALL_SUMMARY_SHA256,
    replay,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/independent_h9_14_record12_coxeter_probe.py"


def test_secondary_probe_is_no_import_and_freezes_both_scopes() -> None:
    source = SOURCE.read_text(encoding="utf-8")
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    assert EXPECTED_ALL_SUMMARY_SHA256 == (
        "06a0d502d2f739190762f5b7e2e6442a2cadb35a4c4a96577ee3885d9d822258"
    )

    quick = replay(False)
    assert quick["orbit_count"] == 48
    assert quick["schedule_histogram"] == {"0": 48}
    assert quick["all_selected_orbits_excluded"] is True
    assert quick["summary_sha256"] == EXPECTED_A0_SUMMARY_SHA256
    assert {
        tuple(stage["elliptic_actual_vertex_survivor_count"] for stage in orbit["stages"])
        for orbit in quick["orbits"]
    } == {(420, 24, 0)}
