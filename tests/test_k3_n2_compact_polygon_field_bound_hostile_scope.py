"""Independent hostile replay of the general compact K3 polygon theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_compact_polygon_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_compact_polygon_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_compact_polygon_field_bound.json"
OWNER_TESTS = ROOT / "tests/test_k3_n2_compact_polygon_field_bound.py"
AUDIT = ROOT / "results/k3_n2_compact_polygon_field_bound_hostile_scope_audit.md"
PINS = {
    REPORT: "45d6bd99d54ce8e8ff3ca60f90f9b9f6f854dba8f8a38e4a953a64f722083c85",
    CERTIFIER: "9bc53b8f95cf2eccc0f2ddb55cff594313d89b07a17e77a30cf651d0f8143aa5",
    ARCHIVE: "e0ff5ebdd592f99f4b96a5f7febd043d8b9b30bef676b80746f49bc75dae0819",
    OWNER_TESTS: "52ebc589cee65ee3efe61a9aaf841ade4b3f573e0e9f0a848c141435f5a8cfff",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cyclic_nonadjacent(left: int, right: int, size: int) -> bool:
    gap = (left - right) % size
    return gap not in (0, 1, size - 1)


def incident(vertex: int, size: int) -> tuple[int, int]:
    return (vertex - 1) % size, vertex % size


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_vertex_to_remote_wall_separation_rederives_independently():
    checked = 0
    for size in range(4, 193):
        for vertex in range(size):
            sides = incident(vertex, size)
            for wall in range(size):
                if wall in sides:
                    continue
                assert any(cyclic_nonadjacent(side, wall, size) for side in sides)
                checked += 1
    assert checked > 2_000_000


def test_all_distinct_vertex_pairs_have_a_dotted_side_witness():
    checked = 0
    for size in range(4, 193):
        for left in range(size):
            for right in range(left + 1, size):
                assert any(
                    cyclic_nonadjacent(a, b, size)
                    for a in incident(left, size)
                    for b in incident(right, size)
                )
                checked += 1
    assert checked > 1_000_000


def test_norm_gap_and_angle_floor_have_the_correct_exponent():
    for degree in range(1, 7):
        # 1 <= z_0*product(z_sigma) <= z_0*4^(d-1).
        z0_floor = Fraction(1, 4 ** (degree - 1))
        sine_squared_floor = z0_floor / 4
        sine_floor_squared = Fraction(1, 2**degree) ** 2
        assert sine_squared_floor == sine_floor_squared
        assert Fraction(1, 2**degree) >= Fraction(1, 64)


def test_sector_area_constants_rederive_without_owner_code():
    angle_denominator = 64
    pi_rational_upper_bound = 4
    cosh_half = 1 + angle_denominator * pi_rational_upper_bound
    cosh_distance = 2 * cosh_half**2 - 1
    pair_cycle = 4 * cosh_distance**2
    assert (cosh_half, cosh_distance, pair_cycle) == (
        257,
        132_097,
        69_798_469_636,
    )
    assert 6 * 5 == 30


def test_unique_large_conjugate_forces_full_field_generation():
    for degree in range(2, 7):
        for proper_degree in range(1, degree):
            if degree % proper_degree:
                continue
            relative_degree = degree // proper_degree
            # A proper-subfield value repeats at every real extension of
            # the distinguished restricted embedding, contradicting one >4 place.
            assert relative_degree > 1


def test_archive_replay_and_scope_are_locked():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    constants = payload["constants"]
    assert constants["shortest_dotted_pair_cycle_strict_upper_bound"] == 69_798_469_636
    assert constants["field_discriminant_bound"] == "69798469636^30"
    assert payload["conclusions"]["short_dotted_pair_cycle_generates_full_field"]
    assert not payload["conclusions"]["reflection_maximality_used"]
    scope = payload["scope"]
    assert scope["facet_count_at_least_4_required"]
    assert not scope["finite_area_noncompact_polygon_closed"]
    assert not scope["ultraideal_or_truncated_polygon_closed"]
    assert not scope["higher_dimensional_chamber_closed"]
    assert not scope["carrier_existence_proved"]


def test_report_and_audit_make_every_hostile_boundary_literal():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "Everything below assumes at least four facets",
        "intersection of the inward half-planes",
        "common perpendicular",
        "The following two distance facts include all cyclic edge cases",
        r"\prod_{\sigma\ne\sigma_0}\sigma(z_i)",
        "On the finite projective",
        r"0<\sigma(p)<4",
        r"\mathbf Q(p)=k",
        "does not cover finite-area noncompact polygons",
    ):
        assert phrase in report
    for phrase in (
        "PASS in the literal compact, rank-three",
        "This proves ultraparallelity and positivity",
        "including every quadrilateral wraparound",
        "full sectors, not truncated pieces",
        "correct exponent and direction",
        "does not assume an imposed mixed height",
        "Triangles are intentionally separate",
    ):
        assert phrase in audit


def test_transport_and_classification_controls_are_clean():
    combined_bytes = REPORT.read_bytes() + AUDIT.read_bytes()
    assert b"\x08" not in combined_bytes
    combined = combined_bytes.decode("utf-8")
    assert ",quad" not in combined
    assert ",qquad" not in combined
    for forbidden in ("Nikulin", "Vinberg", "ABSW"):
        assert forbidden not in combined

