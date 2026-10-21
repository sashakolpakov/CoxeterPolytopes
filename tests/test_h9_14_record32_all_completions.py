from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

import scripts.certify_h9_14_record32_all_completions as aggregate_certifier
from scripts.certify_h9_14_record32_all_completions import (
    EXPECTED_ARCHIVE_FILE_SHA256,
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_HOSTILE_FILE_SHA256,
    OUTPUT,
    _canonical_digest,
    _normalized_certifier_source,
)
from scripts.certify_h9_14_record32_all_completions_candidate import (
    OUTPUT as CANDIDATE_OUTPUT,
    verify as verify_candidate,
)
from scripts.certify_h9_14_record32_whole_shadow_hostile import (
    DEFAULT_ARCHIVE as HOSTILE_OUTPUT,
    sha256 as file_sha256,
    validate_archive as validate_hostile,
)


ROOT = Path(__file__).resolve().parents[1]


def test_certificate_digest_is_json_roundtrip_invariant() -> None:
    document: dict[str, object] = {
        "histogram": {2: 7, 10: 3},
        "nested": [{5: 1, 12: 4}],
    }
    document["certificate_sha256"] = _canonical_digest(document)
    parsed = json.loads(json.dumps(document))
    assert parsed["certificate_sha256"] == _canonical_digest(parsed)


def test_certifier_self_hash_normalizes_only_outer_freezes() -> None:
    source = (
        'EXPECTED_HOSTILE_FILE_SHA256 = "first"\n'
        'preserved = "EXPECTED_HOSTILE_FILE_SHA256 = trailing"\n'
        'EXPECTED_CERTIFICATE_SHA256 = "second"\n'
        'EXPECTED_ARCHIVE_FILE_SHA256 = "third"\n'
        'def preserved_function():\n'
        '    return "all other bytes"\n'
    )
    assert _normalized_certifier_source(source) == (
        'EXPECTED_HOSTILE_FILE_SHA256 = "<OUTER_FREEZE>"\n'
        'preserved = "EXPECTED_HOSTILE_FILE_SHA256 = trailing"\n'
        'EXPECTED_CERTIFICATE_SHA256 = "<OUTER_FREEZE>"\n'
        'EXPECTED_ARCHIVE_FILE_SHA256 = "<OUTER_FREEZE>"\n'
        'def preserved_function():\n'
        '    return "all other bytes"\n'
    )


def test_historical_candidate_remains_unpromoted() -> None:
    candidate = json.loads(CANDIDATE_OUTPUT.read_text(encoding="utf-8"))
    summary = verify_candidate(candidate)
    assert summary["production_candidate_complete"] is True
    assert summary["independent_replay_complete"] is False
    assert summary["whole_graph_shadow_compact_coxeter_theorem"] is False
    assert summary["finite_volume_noncompact_chamber_excluded"] is False


def test_build_cli_writes_canonical_archive(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    output = tmp_path / "record32.json"
    document = {"z": [2, 1], "a": {"proof": True}}
    monkeypatch.setattr(aggregate_certifier, "OUTPUT", output)
    monkeypatch.setattr(aggregate_certifier, "build", lambda: document)
    monkeypatch.setattr(
        aggregate_certifier,
        "verify",
        lambda value: {"rebuilt": value == document},
    )
    monkeypatch.setattr(
        aggregate_certifier,
        "EXPECTED_ARCHIVE_FILE_SHA256",
        "TO_BE_FROZEN",
    )
    aggregate_certifier.main(["--build"])
    assert output.read_text(encoding="utf-8") == (
        json.dumps(document, indent=2, sort_keys=True) + "\n"
    )
    assert json.loads(capsys.readouterr().out) == {"rebuilt": True}


def test_no_project_import_hostile_sources_are_frozen_and_complete() -> None:
    hostile = json.loads(HOSTILE_OUTPUT.read_text(encoding="utf-8"))
    validate_hostile(hostile)
    assert file_sha256(HOSTILE_OUTPUT) == EXPECTED_HOSTILE_FILE_SHA256
    assert hostile["status"] == {
        "independent_replay_complete": True,
        "whole_graph_shadow_compact_coxeter_no_go": True,
    }
    sources = hostile["independent_sources"]
    assert set(sources) == {
        "positive_orbit_replay",
        "a0_replay",
        "a0_orbit_comparator",
        "additions_replay",
        "downstream_replay",
    }
    for source in sources.values():
        path = ROOT / source["path"]
        text = path.read_text(encoding="utf-8")
        assert sha256(path.read_bytes()).hexdigest() == source["sha256"]
        assert "from coxeter_polytopes" not in text
        assert "import coxeter_polytopes" not in text

    comparison = hostile["a0"]["full_action_orbit_comparison"]
    assert comparison["production_orbit_count"] == 292
    assert comparison["independent_orbit_count"] == 292
    assert comparison["target_orbit_sets_equal"] is True
    assert comparison["production_only_orbit_count"] == 0
    assert comparison["independent_only_orbit_count"] == 0


def test_frozen_whole_shadow_certificate_rebuilds_compact_only() -> None:
    assert EXPECTED_HOSTILE_FILE_SHA256 != "TO_BE_FROZEN"
    assert EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
    assert EXPECTED_ARCHIVE_FILE_SHA256 != "TO_BE_FROZEN"
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert sha256(OUTPUT.read_bytes()).hexdigest() == EXPECTED_ARCHIVE_FILE_SHA256
    result = subprocess.run(
        [sys.executable, str(aggregate_certifier.CERTIFIER_SOURCE)],
        cwd=ROOT,
        env=os.environ | {"PYTHONPATH": "src:."},
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    assert json.loads(result.stdout) == {
        "record_index": 32,
        "a0_full_action_target_orbits": 292,
        "a1_ridge_valid_target_orbits": 2,
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert document["branch_closure"][
        "dehn_sommerville_schedules_a_b_facets_link_g4_sum"
    ] == [
        [0, 6, 114, 36],
        [1, 12, 112, 30],
        [2, 18, 110, 24],
        [3, 24, 108, 18],
        [4, 30, 106, 12],
        [5, 36, 104, 6],
    ]
    assert document["scope"]["compact_only"] is True
    assert document["scope"]["ideal_or_parabolic_faces_enumerated"] is False
    assert document["scope"]["finite_volume_noncompact_chamber_excluded"] is False
    assert document["theorem"] == {
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
    }
