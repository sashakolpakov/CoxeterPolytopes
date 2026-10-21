"""Independent hostile replay of the all-degree K3 right-pentagon theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_all_degree_right_pentagon_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_all_degree_right_pentagon_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_all_degree_right_pentagon_field_bound.json"
OWNER_TESTS = ROOT / "tests/test_k3_n2_all_degree_right_pentagon_field_bound.py"
AUDIT = ROOT / "results/k3_n2_all_degree_right_pentagon_field_bound_hostile_scope_audit.md"
PINS = {
    REPORT: "4177f65e13e4fdc19333d5b5ab6097799cd293d961a76aca685630d62cbf2bfc",
    CERTIFIER: "499a1ff95ea7f40ad789e83b8533c2fffcc064c0927a08d7fd6464698190907c",
    ARCHIVE: "5c8c05c10a5f1c5c79dac0abc4df2be18447f743986a4bd7783f347c4b1d669f",
    OWNER_TESTS: "fd16417cbfbce39bd467e1ea488cbdded18f435bba23569650b2e347ee12c88e",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_coordinates(r: Fraction, s: Fraction) -> tuple[Fraction, ...]:
    delta = r * s - 1
    return r, delta, s, (s + 1) / delta, (r + 1) / delta


def residuals(coordinates: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(
        coordinates[(i - 1) % 5] * coordinates[(i + 1) % 5]
        - coordinates[i]
        - 1
        for i in range(5)
    )


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_five_equation_closure_rederives_independently():
    tested = 0
    for denominator_r in range(1, 12):
        for denominator_s in range(1, 12):
            for numerator_r in range(1, 7 * denominator_r):
                for numerator_s in range(1, 7 * denominator_s):
                    r = Fraction(numerator_r, denominator_r)
                    s = Fraction(numerator_s, denominator_s)
                    if r * s <= 1:
                        continue
                    coordinates = independent_coordinates(r, s)
                    assert residuals(coordinates) == (0, 0, 0, 0, 0)
                    A, B, C, D, E = (coordinate + 1 for coordinate in coordinates)
                    assert (
                        C * E - C - D - E + 1,
                        -A + B * E - B - E + 1,
                        B * D - B - C - D + 1,
                        A * D - A - D - E + 1,
                        A * C - A - B - C + 1,
                    ) == (0, 0, 0, 0, 0)
                    tested += 1
    assert tested > 100_000


def test_all_place_norm_floor_is_degree_dependent():
    for degree in range(2, 7):
        # |N(4x)|>=1 and d-1 other factors below four imply
        # 4x>4^{-(d-1)}, hence x>4^{-d}.
        y_floor = Fraction(1, 4 ** (degree - 1))
        x_floor = y_floor / 4
        assert x_floor == Fraction(1, 4**degree)
    assert Fraction(1, 4**6) < Fraction(1, 16)


def test_corrected_recurrence_bound_uses_the_floor_twice():
    for degree in range(2, 7):
        L = 4**degree
        epsilon = Fraction(1, L)
        upper = Fraction(1, epsilon**2) + Fraction(2, epsilon)
        assert upper == L**2 + 2 * L == (L + 1) ** 2 - 1
    assert 4096**2 + 2 * 4096 == 16_785_408
    assert 16 + 17 * 4096 == 69_648
    assert 16_785_408 > 69_648


def test_degree_six_bound_is_approached_by_positive_A2_packets():
    # The order epsilon^{-2} is real, not an artifact of loose algebra.
    epsilon = Fraction(1, 4096)
    z = epsilon + Fraction(1, 4096**3)
    target_x3 = epsilon + Fraction(1, 4096**3)
    x0 = ((z + 1) / target_x3 + 1) / z
    coordinates = independent_coordinates(x0, z)
    assert residuals(coordinates) == (0, 0, 0, 0, 0)
    assert min(coordinates) > epsilon
    assert x0 < 16_785_408
    assert x0 > 16_000_000


def test_nonzero_Cartan_graph_has_only_pairs_and_one_simple_five_cycle():
    vertices = range(5)
    edges = {
        frozenset((i, j))
        for i in vertices
        for j in vertices
        if i < j and (i - j) % 5 not in (1, 4)
    }
    assert len(edges) == 5
    assert all(sum(i in edge for edge in edges) == 2 for i in vertices)
    # No triangles or quadrilaterals occur in a five-cycle.
    for size in (3, 4):
        for subset in combinations(vertices, size):
            induced_edges = [edge for edge in edges if edge <= frozenset(subset)]
            assert len(induced_edges) < size


def test_length_five_and_primitive_grid_constants_are_exact():
    L = 4096
    pair_ratio_bound = (L + 1) ** 2
    H = 32 * (L + 1) ** 5
    assert pair_ratio_bound == 16_785_409
    assert H == 36_938_546_139_294_728_224
    assert len(tuple(combinations(range(6), 2))) == 15
    assert 32 * H == 1_182_033_476_457_431_303_168


def test_quadratic_specialization_remains_sharp_and_pair_generated():
    epsilon = Fraction(1, 16)
    assert Fraction(1, epsilon**2) + Fraction(2, epsilon) == 288
    assert 4 * (288 + 1) == 1156
    assert 1159**2 - 4 == 1_343_277
    # Distinguished pair-cycle values exceed four; definite values are below four.
    assert Fraction(17, 4) > 4 > Fraction(15, 4)


def test_archive_replays_and_scope_does_not_expand():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["conclusions"]["automatic_shift_floor_is_4_to_minus_d"] is True
    assert payload["conclusions"]["suggested_69648_bound_is_not_claimed"] is True
    assert payload["conclusions"]["bounded_cyclic_packet_gives_finite_field_set"] is True
    assert payload["conclusions"]["bounded_mixed_height_premise_required"] is False
    scope = payload["scope"]
    assert scope["returned_field_is_cyclic_Cartan_ground_field_required"] is True
    assert scope["complete_saturated_reflection_packet_retained"] is True
    assert scope["reflection_maximality_in_declared_order_host_retained"] is True
    assert scope["reflection_maximality_used_in_bound"] is False
    assert scope["non_right_or_noncompact_pentagon_closed"] is False
    assert scope["arbitrary_n2_polygon_closed"] is False


def test_report_and_audit_pin_the_correction_and_nonclaims():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "proved for all full-trace complex-K3 degrees",
        r"&<\varepsilon^{-2}+2\varepsilon^{-1}",
        r"16+17\cdot4^6=69,648",
        "is not automatic",
        "cyclic Cartan packet",
        "Hermite--Minkowski proves the finite field set",
        "gratuitous declared field extension",
        r"does not: - prove that every exact \(n=2\) chamber is a right pentagon",
    ):
        assert phrase in report
    for phrase in (
        "PASS, subject to the literal cyclic-ground-field scope",
        "There is no all-degree justification for the quadratic floor",
        r"provisional value \(69,648\) is invalid without an extra theorem",
        r"no single pair coordinate is claimed to generate \(k\) when \(d>2\)",
        "Maximality is retained but unused",
    ):
        assert phrase in audit


def test_tex_transport_and_classification_scope_are_clean():
    combined = REPORT.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
    assert ",quad" not in combined
    assert ",qquad" not in combined
    normalized = " ".join(combined.split())
    assert "angle classification, or reflective-lattice classification is used" in normalized
    assert "No angle or reflective-lattice classification enters" in normalized
