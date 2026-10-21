from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess

import pytest

from coxeter_polytopes.combinatorics.h9_14_record34_a0 import verify_manifest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/enumerations/h9_14_record34_a0_census.json"
SOURCE = ROOT / "scripts/certify_h9_14_record34_a0.cpp"


def test_record34_a0_manifest_replays_independently() -> None:
    document = verify_manifest(MANIFEST)
    assert document["totals"] == {
        "genbg_representatives": 1_288_153,
        "active_colourings": 26_598_204,
        "distinct_row_states": 25_501_701,
        "no_small_transversal_states": 3_101_794,
        "target_114_states": 18,
        "pure_target_states": 18,
        "ridge_valid_states": 0,
    }
    assert document["canonical_orbit_count"] == 9


def test_record34_a0_cpp_replays_all_archived_target_sources(tmp_path: Path) -> None:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        pytest.skip("a C++20 compiler is required for the enumerator smoke test")
    executable = tmp_path / "record34_a0"
    subprocess.run(
        [compiler, "-O2", "-std=c++20", str(SOURCE), "-o", str(executable)],
        check=True,
        cwd=ROOT,
    )
    document = json.loads(MANIFEST.read_text(encoding="utf-8"))
    sources = tuple(
        dict.fromkeys(record["source_graph6"] for record in document["target_states"])
    )
    completed = subprocess.run(
        [str(executable)],
        input="".join(f"{source}\n" for source in sources),
        text=True,
        capture_output=True,
        check=True,
        cwd=ROOT,
    )
    assert (
        "graphs 6 degree_graphs 6 colourings 132 distinct 90 "
        "no_small 69 target 18 pure 18 ridge 0"
    ) in completed.stdout
    assert completed.stderr.count("TARGET profile=") == 18
    assert "RIDGE profile=" not in completed.stderr
    assert "degree_two_dual_components=" in completed.stderr

