from __future__ import annotations

from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts" / "audit_h9_14_record34_ridge_csp.cpp"

EXPECTED_LINES = (
    "0,0,0 binary=UNSAT nodes=19189069",
    "0,0,1 binary=UNSAT nodes=9790583",
    "0,0,2 binary=UNSAT nodes=4573694",
    "0,1,1 binary=UNSAT nodes=14114074",
    "0,1,2 binary=UNSAT nodes=7363401",
    "0,2,2 binary=UNSAT nodes=11272993",
    "1,0,0 binary=UNSAT nodes=25658567",
    "1,0,1 binary=UNSAT nodes=13437573",
    "1,0,2 binary=UNSAT nodes=6196706",
    "1,1,1 binary=UNSAT nodes=19465238",
    "1,1,2 binary=UNSAT nodes=10088675",
    "1,2,2 binary=UNSAT nodes=17687387",
    "2,0,0 binary=UNSAT nodes=17327195",
    "2,0,1 binary=UNSAT nodes=7609697",
    "2,0,2 binary=UNSAT nodes=2571830",
    "2,1,1 binary=UNSAT nodes=12607068",
    "2,1,2 binary=UNSAT nodes=5612096",
    "2,2,2 binary=UNSAT nodes=13239748",
    (
        "3,0,0 binary=SAT nodes=6648352 "
        "binary+ternary=UNSAT nodes=2076565"
    ),
    "3,0,1 binary=UNSAT nodes=2737208",
    "3,0,2 binary=UNSAT nodes=487629",
    "3,1,1 binary=UNSAT nodes=6005109",
    "3,1,2 binary=UNSAT nodes=2517144",
    "profiles 23 binary_unsat 22 binary_sat 1 final_unsat 23",
)


def test_record34_ridge_csp_has_no_profile_survivor(tmp_path: Path) -> None:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        pytest.skip("a C++20 compiler is required for the ridge-CSP audit")
    executable = tmp_path / "record34_ridge_csp"
    subprocess.run(
        [compiler, "-O3", "-std=c++20", str(SOURCE), "-o", str(executable)],
        check=True,
        cwd=ROOT,
    )
    completed = subprocess.run(
        [str(executable)],
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
        cwd=ROOT,
    )
    assert tuple(completed.stdout.splitlines()) == EXPECTED_LINES
    assert completed.stderr == ""
