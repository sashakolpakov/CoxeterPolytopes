"""Independent hostile tests for the compact K3 Coxeter-polygon theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from scripts.audit_k3_n2_compact_polygon_field_bound_independent_hostile import (
    build_certificate,
    constants_record,
    fermi_strip_linear_forms,
    fermi_wedge_terms,
    incident_facets,
    nonadjacent,
    norm_angle_record,
    proper_subfield_extension_record,
    separation_record,
)


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = ROOT / "results/k3_n2_compact_polygon_field_bound.md"
OWNER_CERT = ROOT / "scripts/certify_k3_n2_compact_polygon_field_bound.py"
OWNER_ARCHIVE = ROOT / "data/enumerations/k3_n2_compact_polygon_field_bound.json"
OWNER_TEST = ROOT / "tests/test_k3_n2_compact_polygon_field_bound.py"
LOCAL_HOSTILE_REPORT = (
    ROOT / "results/k3_n2_compact_polygon_field_bound_hostile_scope_audit.md"
)
LOCAL_HOSTILE_TEST = (
    ROOT / "tests/test_k3_n2_compact_polygon_field_bound_hostile_scope.py"
)
INDEPENDENT_REPORT = (
    ROOT / "results/k3_n2_compact_polygon_field_bound_independent_hostile_audit.md"
)
INDEPENDENT_SCRIPT = (
    ROOT / "scripts/audit_k3_n2_compact_polygon_field_bound_independent_hostile.py"
)

PINS = {
    OWNER_REPORT: "45d6bd99d54ce8e8ff3ca60f90f9b9f6f854dba8f8a38e4a953a64f722083c85",
    OWNER_CERT: "9bc53b8f95cf2eccc0f2ddb55cff594313d89b07a17e77a30cf651d0f8143aa5",
    OWNER_ARCHIVE: "e0ff5ebdd592f99f4b96a5f7febd043d8b9b30bef676b80746f49bc75dae0819",
    OWNER_TEST: "52ebc589cee65ee3efe61a9aaf841ade4b3f573e0e9f0a848c141435f5a8cfff",
    LOCAL_HOSTILE_REPORT: "4474ba6c5565d5a0ec65cb04704efcbde9f9a906fa1f5573d8c6f827fd9bce19",
    LOCAL_HOSTILE_TEST: "6e9d95def280f11c662e4c2fdc9b930d660b0a73bb02a84ba8c5efd635bc4d4b",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_frozen_owner_and_local_hostile_bytes_are_exact():
    assert {path: sha256(path) for path in PINS} == PINS


def test_nonobtuse_fermi_wedge_forbids_a_foot_beyond_the_endpoint():
    for alpha in (math.pi / 64, math.pi / 7, math.pi / 3, math.pi / 2):
        for t in (-0.01, -0.5, -4.0):
            for r in (0.0, 0.25, 3.0):
                first, second = fermi_wedge_terms(alpha, t, r)
                assert first > 0
                assert second >= -1e-14
                assert first + second > 0


def test_perpendicular_strip_is_cut_out_by_the_stated_linear_halfspaces():
    left, right = -0.7, 1.2
    for r in (-3.0, -0.1, 0.0, 0.8, 4.0):
        for t in (-2.0, left, 0.3, right, 2.5):
            left_form, right_form = fermi_strip_linear_forms(t, r, left, right)
            in_strip = left - 1e-12 <= t <= right + 1e-12
            in_halfspaces = left_form >= -1e-10 and right_form <= 1e-10
            assert in_strip == in_halfspaces


def test_vertex_wall_and_all_vertex_pair_selectors_cover_n_four_through_256():
    record = separation_record(256)
    assert record["facet_range"] == [4, 256]
    assert record["all_selectors_pass"] is True
    assert record["vertex_wall_pairs"] > 5_000_000
    assert record["vertex_pairs"] > 2_000_000


def test_named_quadrilateral_and_even_wraparound_witnesses():
    for size in range(4, 258, 2):
        for vertex in range(size):
            opposite = (vertex + size // 2) % size
            assert any(
                nonadjacent(a, b, size)
                for a in incident_facets(vertex, size)
                for b in incident_facets(opposite, size)
            )
    assert nonadjacent(-1, 1, 4)
    assert nonadjacent(0, 2, 4)


def test_projective_trace_norm_floor_has_exact_division_by_four():
    rows = norm_angle_record()
    assert len(rows) == 6
    for row in rows:
        degree = int(row["degree"])
        assert Fraction(row["z0_floor"]) == Fraction(1, 4 ** (degree - 1))
        assert Fraction(row["sine_square_floor"]) == Fraction(1, 4**degree)
        assert Fraction(row["sine_floor"]) == Fraction(1, 2**degree)
    assert Fraction(rows[-1]["sine_floor"]) == Fraction(1, 64)


def test_sector_area_and_pair_cycle_constants_rederive_exactly():
    record = constants_record()
    assert record == {
        "degree_bound": 6,
        "angle_floor": "1/64",
        "cosh_half_distance_bound": 257,
        "cosh_distance_bound": 132_097,
        "pair_cycle_bound": 69_798_469_636,
        "discriminant_exponent": 30,
    }
    assert 2 * 257**2 - 1 == 132_097
    assert 4 * 132_097**2 == 69_798_469_636


def test_unique_large_value_cannot_come_from_a_proper_subfield():
    rows = proper_subfield_extension_record()
    assert rows
    assert all(row["real_extensions"] >= 2 for row in rows)
    assert {row["total_degree"] for row in rows} == {2, 3, 4, 5, 6}


def test_independent_script_replays_without_owner_imports():
    completed = subprocess.run(
        [sys.executable, str(INDEPENDENT_SCRIPT), "--max-facets", "128"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    assert payload["verdict"] == "PASS_EXACT_COMPACT_FULL_TRACE_SCOPE"
    assert payload["projection"]["uses_intersection_of_both_inward_halfplanes"]
    assert payload["constants"]["pair_cycle_bound"] == 69_798_469_636
    assert not payload["scope"]["finite_area_noncompact_closed"]
    source = INDEPENDENT_SCRIPT.read_text(encoding="utf-8")
    assert "certify_k3_n2_compact_polygon_field_bound" not in source


def test_report_states_every_pressure_point_and_scope_boundary():
    text = " ".join(INDEPENDENT_REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "needs both half-planes at an endpoint",
        "intersection of two geodesic half-planes",
        "every nonadjacent support pair is ultraparallel",
        "quadrilateral wraparound",
        "untruncated sector",
        r"\sin^2\alpha_i\ge 4^{-d}",
        "The exponent is",
        r"69\,798\,469\,636",
        "exactly one conjugate larger",
        r"\mathbf Q(p)=k",
        "literal finite field set",
        "minimal cyclic Cartan ground-field hypothesis",
        "Finite-area noncompact",
    ):
        assert phrase in text


def test_tex_transport_and_forbidden_dependency_controls_are_clean():
    combined = INDEPENDENT_REPORT.read_bytes() + INDEPENDENT_SCRIPT.read_bytes()
    assert not any(byte < 32 and byte not in (9, 10, 13) for byte in combined)
    text = combined.decode("utf-8")
    assert "|le" not in text
    assert ",quad" not in text
    assert ",qquad" not in text
    for forbidden in ("Nikulin", "Vinberg", "ABSW"):
        assert forbidden not in text


def test_certificate_scope_and_geometry_flags_are_exact():
    payload = build_certificate(64)
    assert payload["projection"] == {
        "uses_intersection_of_both_inward_halfplanes": True,
        "fermi_strip_is_intersection_of_geodesic_halfplanes": True,
        "support_ultraparallelity_requires_compact_facet_minimizers": True,
    }
    assert payload["scope"]["compact_polygon_required"]
    assert payload["scope"]["primitive_exact_full_trace_K3_required"]
    assert payload["scope"]["minimal_cyclic_ground_field_retained_by_owner"]
    assert not payload["scope"]["minimality_used_by_unique_large_generator_argument"]
    assert not payload["scope"]["carrier_construction_proved"]
    assert not payload["scope"]["reflection_maximality_used"]
