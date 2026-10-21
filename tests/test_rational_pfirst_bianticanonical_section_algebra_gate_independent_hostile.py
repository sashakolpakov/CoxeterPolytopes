"""Independent hostile replay of the bianticanonical section-algebra gate."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_pfirst_bianticanonical_section_algebra_gate.md"
CERTIFIER = ROOT / "scripts/certify_rational_pfirst_bianticanonical_section_algebra_gate.py"
ARCHIVE = ROOT / "data/enumerations/rational_pfirst_bianticanonical_section_algebra_gate.json"
OWNER_TESTS = ROOT / "tests/test_rational_pfirst_bianticanonical_section_algebra_gate.py"
AUDIT = ROOT / "results/rational_pfirst_bianticanonical_section_algebra_gate_independent_hostile_audit.md"
PINS = {
    REPORT: "f68d21dd80d6963467a062cc988d3464104fd6ad737602c929001193e730e6c0",
    CERTIFIER: "f3e3c96464d6eb6f9a762bddb8fe8fcc003c41bb44910f0bd2d602821081b54d",
    ARCHIVE: "6055b14f762e26d2e86d9fcda08c2fb3e20856094f460e6a21e565ef2e310db7",
    OWNER_TESTS: "27922c4f163d1c42aa5616c22faa387e0fe4ac3134f81870ad0cd2c3eb94af45",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lorentz_dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return left[0] * right[0] - sum(a * b for a, b in zip(left[1:], right[1:]))


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_rank_eleven_classes_and_adjunction_rederive_independently():
    canonical = (-3,) + (1,) * 10
    anticanonical = tuple(-entry for entry in canonical)
    boundary = tuple(2 * entry for entry in anticanonical)
    assert lorentz_dot(canonical, canonical) == -1
    assert lorentz_dot(anticanonical, anticanonical) == -1
    assert lorentz_dot(boundary, boundary) == -4
    assert lorentz_dot(canonical, boundary) == 2
    arithmetic_genus = 1 + (
        lorentz_dot(boundary, boundary) + lorentz_dot(canonical, boundary)
    ) // 2
    assert arithmetic_genus == 0


def test_all_bezout_budgets_are_strict_and_have_the_right_local_orders():
    budgets = {
        "cubic_cubic": (3 * 3, 10 * 1 * 1),
        "cubic_double_sextic": (3 * 6, 10 * 1 * 2),
        "nodal_sextic_cubic": (6 * 3, 10 * 2 * 1),
        "nodal_sextic_double_sextic": (6 * 6, 10 * 2 * 2),
    }
    assert budgets == {
        "cubic_cubic": (9, 10),
        "cubic_double_sextic": (18, 20),
        "nodal_sextic_cubic": (18, 20),
        "nodal_sextic_double_sextic": (36, 40),
    }
    assert all(local > global_degree for global_degree, local in budgets.values())


def test_formal_matrix_shapes_and_labelled_exceptional_fan_are_exact():
    assert math.comb(3 + 2, 2) == 10
    assert math.comb(6 + 2, 2) == 28
    assert 10 * math.comb(1 + 1, 2) == 10
    assert 10 * math.comb(2 + 1, 2) == 30
    assert 2**10 == 1024
    assert math.factorial(10) == 3_628_800


def test_winger_projection_supplies_the_exact_nodal_positive_control():
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    projected = payload["explicit_ten_node_existence_dependency"][
        "projected_certificate"
    ]
    assert projected["center_count"] == 10
    assert projected["all_ten_nodes"] is True
    assert projected["absolutely_irreducible"] is True
    assert projected["normalization_genus"] == 0
    assert (projected["cubic_rank"], projected["sextic_rank"]) == (10, 27)
    assert projected["strong_coble"] is True


def test_owner_archive_replays_and_keeps_global_nonclaims():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-archive"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    scope = json.loads(ARCHIVE.read_text(encoding="utf-8"))["scope"]
    assert scope["same_full_nef_cone_claimed"] is False
    assert scope["same_raw_coxeter_chamber_claimed"] is False
    assert scope["global_property4_claimed"] is False
    assert scope["universal_incidence_nonemptiness_claimed"] is False
    assert scope["finite_generation_of_Cox_ring_required"] is False


def test_repaired_corollary_states_every_maximal_return_hypothesis():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "finite-sided finite-volume Coxeter chamber",
        "distinct oriented physical facet classes",
        r"b\ge1",
        "pairwise orthogonal",
        r"\sum_{i=1}^bB_i=-2K",
        "displayed irredundant halfspaces",
        "Every side reflection preserves the adapted host",
        "no proper discrete reflection overgroup",
        "on this same",
        "Intersecting any other member",
        "corresponding prime-section locus",
    ):
        assert phrase in report


def test_report_and_hostile_audit_keep_labelled_fan_and_prime_scope_exact():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "complete labelled-exceptional contraction fan",
        "does not exclude a product of two different divisor classes",
        "divisor primality, equivalently geometric integrality",
        "universal incidence nonemptiness",
        "same full nef cone or the same raw Coxeter chamber",
    ):
        assert phrase in report
    for phrase in (
        "PASS after making the exact-return hypotheses literal",
        "not claimed to enumerate every Mori contraction",
        "positive branch does not rest on generic-projection language",
        "residual divisor is effective and linearly trivial",
        "one and the same marked surface",
        "componentwise uniqueness",
        "does not replace coordinate-valued evaluation matrices",
    ):
        assert phrase in audit


def test_tex_transport_controls_are_clean():
    combined = REPORT.read_bytes() + AUDIT.read_bytes()
    assert b"\x08" not in combined
    text = combined.decode("utf-8")
    assert ",quad" not in text
    assert ",qquad" not in text
