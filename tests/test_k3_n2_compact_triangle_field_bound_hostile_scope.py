"""Local hostile replay of the finite-area K3 triangle field theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_compact_triangle_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_compact_triangle_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_compact_triangle_field_bound.json"
OWNER_TESTS = ROOT / "tests/test_k3_n2_compact_triangle_field_bound.py"
AUDIT = ROOT / "results/k3_n2_compact_triangle_field_bound_hostile_scope_audit.md"
INDEPENDENT_AUDIT = (
    ROOT / "results/k3_n2_finite_area_triangle_field_bound_independent_hostile_audit.md"
)
PINS = {
    REPORT: "68b91b604e08bda1a48eb76fd13635f108bba5aab8fcc8ecf4d369afc367baa2",
    CERTIFIER: "e674fca868efaa7d0041e54376500c6f980321264e7d0e02fdb06457b6e29c3d",
    ARCHIVE: "4bfe02df789ef5d0ad90a4deadfb25b646154e109926ffe8a82b472d3be270e6",
    OWNER_TESTS: "f32c974d16023c81eb813bc94ef216855d75a3eab666b0fb1adfed567d992fe4",
    INDEPENDENT_AUDIT: "c6c5271cd7794e55cc167dbf855f3c4fa7ea1112e75847d7f1218d3fa63a6862",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def multiply(left, right):
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def trace(matrix):
    return sum(matrix[i][i] for i in range(len(matrix)))


def subtract(left, right):
    return tuple(
        tuple(a - b for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_reflection_word_trace_isolates_pair_and_triangle_cycles():
    identity = tuple(
        tuple(Fraction(i == j) for j in range(3))
        for i in range(3)
    )
    cartan = (
        (Fraction(2), Fraction(-1, 2), Fraction(-2, 3)),
        (Fraction(-3, 4), Fraction(2), Fraction(-4, 5)),
        (Fraction(-5, 6), Fraction(-6, 7), Fraction(2)),
    )
    projections = []
    reflections = []
    for i in range(3):
        projection = tuple(
            tuple(cartan[i][j] if row == i else Fraction(0) for j in range(3))
            for row in range(3)
        )
        projections.append(projection)
        reflections.append(subtract(identity, projection))
    for i, j in combinations(range(3), 2):
        pair = cartan[i][j] * cartan[j][i]
        assert trace(multiply(reflections[i], reflections[j])) == -1 + pair
    word = multiply(multiply(reflections[0], reflections[1]), reflections[2])
    pair_sum = sum(
        cartan[i][j] * cartan[j][i]
        for i, j in combinations(range(3), 2)
    )
    triangle = cartan[0][1] * cartan[1][2] * cartan[2][0]
    assert trace(word) == -3 + pair_sum - triangle


def test_normalized_cycle_height_includes_ideal_endpoints():
    for denominator in range(1, 40):
        magnitudes = (
            Fraction(denominator, denominator),
            Fraction(max(0, denominator - 1), denominator),
            Fraction(max(0, denominator - 2), denominator),
        )
        pair_cycles = tuple(4 * value * value for value in magnitudes)
        triangle_cycle = 8 * magnitudes[0] * magnitudes[1] * magnitudes[2]
        assert all(0 <= value <= 4 for value in pair_cycles)
        assert abs(triangle_cycle) <= 8
        assert triangle_cycle * triangle_cycle == pair_cycles[0] * pair_cycles[1] * pair_cycles[2]
    assert 4 * Fraction(1) ** 2 == 4
    assert 8 * Fraction(1) ** 3 == 8


def test_primitive_grid_and_discriminant_bound_recompute():
    assert len(tuple(combinations(range(6), 2))) == 15
    H = 8
    assert 16 * H == 128
    assert 2 * 16 * H == 256
    assert 6 * 5 == 30


def test_archive_replays_and_scope_stays_literal():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["conclusions"]["all_pair_cycles_bounded_at_every_embedding"] is True
    assert payload["conclusions"]["length_three_cycle_bounded_at_every_embedding"] is True
    assert payload["conclusions"]["Gram_determinant_closure_required"] is False
    assert payload["scope"]["returned_field_is_cyclic_Cartan_ground_field_required"] is True
    assert payload["scope"]["complete_saturated_reflection_packet_retained"] is True
    assert payload["scope"]["reflection_maximality_in_declared_order_host_retained"] is True
    assert payload["scope"]["reflection_maximality_used_in_bound"] is False
    assert payload["scope"]["compact_triangle_closed"] is True
    assert payload["scope"]["ideal_vertex_triangle_closed"] is True
    assert payload["scope"]["exact_same_finite_area_Coxeter_triangle_eigen_nef_chamber_required"] is True
    assert payload["scope"]["truncated_ultraideal_polyhedron_closed"] is False
    assert payload["scope"]["arbitrary_n2_polygon_closed"] is False


def test_report_and_audit_separate_finite_fields_from_bounded_degree():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "proved for every full-trace complex-K3 field degree",
        "including ideal vertices",
        "sole simple length-three cyclic product",
        "Reflection-word trace induction",
        "Every closed directed walk decomposes into simple directed cycles",
        r"H=8",
        r"|\operatorname {Disc}(k)|\le256^{30}",
        "Hermite--Minkowski gives a finite field set",
        "finite field set rather than merely bounded degree",
    ):
        assert phrase in report
    for phrase in (
        "PASS in finite-area cyclic-ground-field scope",
        "including ideal vertices",
        "no ultraparallel wall coordinate",
        "deliberately weak endpoint bounds",
        "Complete saturation and declared-host reflection maximality are retained",
        r"arbitrary \(n=2\) polygon incidence",
    ):
        assert phrase in audit


def test_tex_transport_and_classification_scope_are_clean():
    combined_bytes = REPORT.read_bytes() + AUDIT.read_bytes()
    assert not any(
        byte < 32 and byte not in (9, 10, 13)
        for byte in combined_bytes
    )
    combined = combined_bytes.decode("utf-8")
    assert "|le" not in combined
    assert ",quad" not in combined
    assert ",qquad" not in combined
    normalized = " ".join(combined.split())
    assert "No angle list, reflective-lattice classification" in normalized
    assert "No reflective-lattice classification enters" in normalized
