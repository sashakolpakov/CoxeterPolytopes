from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess

import pytest

from scripts.certify_h9_14_record18_positive_A4_scout import (
    ARCHIVE,
    fresh_replay,
    verify,
)
from scripts.certify_h9_14_record18_a3_additions import (
    ARCHIVE as A3_ARCHIVE,
    verify as verify_a3,
)


ROOT = Path(__file__).resolve().parents[1]
EXCEPTION_SOURCE = ROOT / "scripts/explore_h9_14_record18_a4_exception.cpp"


def test_record18_positive_archive_is_exact_but_not_a_theorem() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 18,
        "production_positive_census_complete": True,
        "positive_orbit_count": 70,
        "a4_quotient_weighting_verified": True,
        "independent_replay_complete": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def test_record18_positive_source_freshly_reproduces_all_four_branches() -> None:
    replay = fresh_replay()
    assert replay["fresh_source_replay_matches_archive"] is True


@pytest.fixture(scope="session")
def exceptional_a4_audit(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        pytest.skip("no C++20 compiler is available")
    executable = tmp_path_factory.mktemp("record18-a4") / "exception"
    subprocess.run(
        [
            compiler,
            "-std=c++20",
            "-O3",
            "-Wall",
            "-Wextra",
            "-pedantic",
            str(EXCEPTION_SOURCE),
            "-o",
            str(executable),
        ],
        check=True,
        cwd=ROOT,
    )
    completed = subprocess.run(
        [str(executable)],
        check=True,
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    return json.loads(completed.stdout)


def test_record18_exceptional_a4_orbit_has_no_scheduled_completion(
    exceptional_a4_audit: dict[str, object],
) -> None:
    assert exceptional_a4_audit["forced_A5_shadow_count"] == 28
    assert exceptional_a4_audit["additional_A5_count"] == 4
    assert exceptional_a4_audit["candidate_count"] == 764
    assert exceptional_a4_audit["raw_unordered_addition_count"] == 14_084_657_251
    assert exceptional_a4_audit["base_small_transversal_count"] == 47
    assert exceptional_a4_audit["no_small_transversal_addition_count"] == 3_960
    assert exceptional_a4_audit["facet_count_histogram"] == {
        "82": 160,
        "83": 1_920,
        "85": 240,
        "87": 120,
        "88": 900,
        "92": 260,
        "94": 360,
        "102": 0,
    }
    assert exceptional_a4_audit["scheduled_102_facet_target_count"] == 0


def test_record18_all_a3_orbits_are_production_closed_but_not_promoted() -> None:
    document = json.loads(A3_ARCHIVE.read_text(encoding="utf-8"))
    assert verify_a3(document) == {
        "record_index": 18,
        "a3_positive_orbit_count": 40,
        "production_a3_branch_empty": True,
        "independent_replay_complete": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }
