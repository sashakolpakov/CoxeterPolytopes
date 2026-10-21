"""Independent hostile tests for exact-trace cusp field collapse."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from scripts.audit_exact_trace_noncompact_ground_field_collapse_independent_hostile import (
    build_audit,
    cube_cusp_record,
)


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/exact_trace_noncompact_ground_field_collapse.md"
OWNER_CERT = ROOT / "scripts/certify_exact_trace_noncompact_ground_field_collapse.py"
OWNER_ARCHIVE = ROOT / "data/enumerations/exact_trace_noncompact_ground_field_collapse.json"
OWNER_TEST = ROOT / "tests/test_exact_trace_noncompact_ground_field_collapse.py"
AUDIT_REPORT = (
    ROOT / "results/exact_trace_noncompact_ground_field_collapse_independent_hostile_audit.md"
)
AUDIT_SCRIPT = (
    ROOT / "scripts/audit_exact_trace_noncompact_ground_field_collapse_independent_hostile.py"
)
PINS = {
    OWNER_REPORT: "d738018bacc3803f78fa725da375bb8045ecce92a1c879b03b30a9998e980255",
    OWNER_CERT: "61f7d0a21cab3b80975849d015c8f2cf13386b3d5f35890917dd0b672e6ebfcb",
    OWNER_ARCHIVE: "fb87928a12165c226272078ce57f8c98418ce02f6629a80fecb3e5769e7cfcd7",
    OWNER_TEST: "d2fe708d1c368ed2ee02c673c4e456e5a497b8dc7d8ce8e318ebbc9c2d00d40f",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_frozen_owner_bytes_are_exact():
    assert {path: sha256(path) for path in PINS} == PINS


def test_cube_cusps_include_nonsimple_fixed_line_pressure():
    for dimension in range(2, 31):
        record = cube_cusp_record(dimension)
        assert record["normal_rank"] == dimension
        assert record["fixed_kernel_dimension"] == 1
        assert record["incident_walls"] == 2 * (dimension - 1)
        assert record["nonsimple"] == (dimension > 2)


def test_dimension_one_is_outside_finite_volume_noncompact_scope():
    with pytest.raises(ValueError):
        cube_cusp_record(1)


def test_logical_scope_flags_are_exact():
    payload = build_audit(12)
    assert payload["verdict"] == "PASS_EXACT_FULL_CUSP_SCOPE"
    assert all(payload["logical_gates"].values())
    assert payload["scope"]["full_incident_packet_required"]
    assert not payload["scope"]["simple_cusp_required"]
    assert not payload["scope"]["free_order_module_required"]
    assert not payload["scope"]["primitive_normals_required"]
    assert not payload["scope"]["compact_case_closed"]


def test_independent_script_replays_without_owner_import():
    completed = subprocess.run(
        [sys.executable, str(AUDIT_SCRIPT), "--max-dimension", "14"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    assert payload["frames"][-1]["dimension"] == 14
    assert "certify_exact_trace_noncompact" not in AUDIT_SCRIPT.read_text(
        encoding="utf-8"
    )


def test_report_states_each_fixed_line_dependency_and_nonclaim():
    text = " ".join(AUDIT_REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "PASS in the exact finite-sided",
        "zero recession cone",
        r"\sum_i\lambda_i a_i>0",
        r"E+\mathbf R u_0=u_0^\perp",
        "kernel formation commutes",
        r"h(u,u)=0",
        r"k=\mathbf Q",
        "partial affine diagram does not suffice",
        "says nothing about compact chambers",
    ):
        assert phrase in text


def test_transport_controls_are_clean():
    combined = AUDIT_REPORT.read_bytes() + AUDIT_SCRIPT.read_bytes()
    assert not any(byte < 32 and byte not in (9, 10, 13) for byte in combined)
    text = combined.decode("utf-8")
    assert "|le" not in text
    assert ",quad" not in text
    assert ",qquad" not in text
