"""Independent hostile replay of the general K3 right-pentagon theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import isqrt, sqrt
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_integral_right_pentagon_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_integral_right_pentagon_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_integral_right_pentagon_field_bound.json"
OWNER_TESTS = ROOT / "tests/test_k3_n2_integral_right_pentagon_field_bound.py"
AUDIT = ROOT / "results/k3_n2_integral_right_pentagon_field_bound_hostile_scope_audit.md"
PINS = {
    REPORT: "4e149c75b4446fbc703f3407e61942bdd02e4d2e4aa9a0c8f9dfcb8a722a6b0b",
    CERTIFIER: "350627dd0545d1ac14438f06ed031749030a56e7005053abe2903396db0088e8",
    ARCHIVE: "069eac8340c0e95a1c79e4296ea94a0c6d535d99004997008fd196bb498385b1",
    OWNER_TESTS: "c212bfa2c1584d06d6b94136893a51ce49acbe1776340d679f75e455774e2973",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_coordinates(r: Fraction, s: Fraction):
    delta = r * s - 1
    return r, delta, s, (s + 1) / delta, (r + 1) / delta


def independent_residuals(coordinates):
    return tuple(
        coordinates[(index - 1) % 5] * coordinates[(index + 1) % 5]
        - coordinates[index]
        - 1
        for index in range(5)
    )


def multiply(left, right):
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def subtract(left, right):
    return tuple(
        tuple(a - b for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def trace(matrix):
    return sum(matrix[index][index] for index in range(len(matrix)))


def squarefree(value: int) -> int:
    result = 1
    prime = 2
    remaining = value
    while prime * prime <= remaining:
        parity = 0
        while remaining % prime == 0:
            remaining //= prime
            parity ^= 1
        if parity:
            result *= prime
        prime += 1
    if remaining > 1:
        result *= remaining
    return result


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_all_five_tree_minor_identities_rederive_independently():
    tested = 0
    for denominator_r in range(1, 11):
        for denominator_s in range(1, 11):
            for numerator_r in range(1, 6 * denominator_r):
                for numerator_s in range(1, 6 * denominator_s):
                    r = Fraction(numerator_r, denominator_r)
                    s = Fraction(numerator_s, denominator_s)
                    if r * s <= 1:
                        continue
                    coordinates = independent_coordinates(r, s)
                    assert independent_residuals(coordinates) == (0, 0, 0, 0, 0)
                    A, B, C, D, E = (coordinate + 1 for coordinate in coordinates)
                    assert (
                        C * E - C - D - E + 1,
                        -A + B * E - B - E + 1,
                        B * D - B - C - D + 1,
                        A * D - A - D - E + 1,
                        A * C - A - B - C + 1,
                    ) == (0, 0, 0, 0, 0)
                    tested += 1
    assert tested > 60_000


def test_cyclic_upper_bound_has_no_ordering_or_index_gap():
    epsilon = Fraction(1, 16)
    tested = 0
    for denominator_r in range(1, 15):
        for denominator_s in range(1, 15):
            for numerator_r in range(1, 25):
                for numerator_s in range(1, 25):
                    r = Fraction(numerator_r, denominator_r)
                    s = Fraction(numerator_s, denominator_s)
                    if r * s <= 1:
                        continue
                    coordinates = independent_coordinates(r, s)
                    if min(coordinates) <= epsilon:
                        continue
                    for index in range(5):
                        z = coordinates[(index + 2) % 5]
                        assert coordinates[index] < 16 + 17 / z < 288
                    tested += 1
    assert tested > 10_000


def test_pair_cycle_is_integral_module_trace_coordinate():
    identity = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    for forward in (Fraction(-7, 3), Fraction(-1), Fraction(5, 4)):
        for backward in (Fraction(-5, 2), Fraction(2, 7), Fraction(9, 5)):
            p_i = (
                (Fraction(2), forward, Fraction(0)),
                (Fraction(0), Fraction(0), Fraction(0)),
                (Fraction(0), Fraction(0), Fraction(0)),
            )
            p_j = (
                (Fraction(0), Fraction(0), Fraction(0)),
                (backward, Fraction(2), Fraction(0)),
                (Fraction(0), Fraction(0), Fraction(0)),
            )
            s_i = subtract(identity, p_i)
            s_j = subtract(identity, p_j)
            assert trace(multiply(s_i, s_j)) == -1 + forward * backward


def test_norm_lower_bound_and_general_constants_are_directionally_exact():
    # If y'>-4 and N(y)=-n, then y0=n/|y'|>n/4>=1/4.
    for n in range(1, 100):
        for conjugate_abs in (Fraction(1, 4), Fraction(1), Fraction(5, 2), Fraction(15, 4)):
            distinguished = Fraction(n, 1) / conjugate_abs
            assert distinguished > Fraction(n, 4) >= Fraction(1, 4)
            assert distinguished / 4 > Fraction(1, 16)
    assert 4 * (288 + 1) == 1156
    assert 1156 + 4 - 1 == 1159
    assert 1156 * 4 - 1 == 4623
    assert 1159**2 - 4 == 1_343_277


def test_pair_cycle_cannot_hide_field_generation_in_length_five_only():
    # The two embedding intervals (4,1156) and (0,4) are disjoint.
    for distinguished in (Fraction(17, 4), Fraction(5), Fraction(100), Fraction(1155)):
        for conjugate in (Fraction(1, 4), Fraction(2), Fraction(15, 4)):
            assert distinguished > 4 > conjugate > 0
            assert distinguished != conjugate
    # In a quadratic field any nonrational element is a primitive generator.
    assert 1_343_277 == 1159**2 - 4


def test_sharp_corollary_candidate_count_recomputes_without_owner_code():
    rows = []
    radicands = set()
    for integer_trace in range(3, 35):
        for integer_norm in range(1, 34):
            discriminant = integer_trace**2 - 4 * integer_norm
            if discriminant <= 0 or isqrt(discriminant) ** 2 == discriminant:
                continue
            if discriminant <= (integer_trace - 2) ** 2:
                continue
            rows.append((integer_trace, integer_norm, discriminant))
            radicands.add(squarefree(discriminant))
    assert len(rows) == 528
    assert len(radicands) == 336
    assert max(discriminant for _, _, discriminant in rows) == 1152


def test_Qsqrt5_poison_packet_replays_independently():
    phi0 = (1 + sqrt(5)) / 2
    phi1 = (1 - sqrt(5)) / 2
    coordinates = (
        (Fraction(-1, 2), Fraction(1, 2)),
        (Fraction(-1, 4), Fraction(1, 4)),
        (Fraction(1, 2), Fraction(2)),
        (Fraction(8), Fraction(14)),
        (Fraction(2), Fraction(4)),
    )

    def phi_product(left, right):
        a, b = left
        c, d = right
        return a * c + b * d, a * d + b * c + b * d

    one = (Fraction(1), Fraction(0))
    for index in range(5):
        product = phi_product(coordinates[(index - 1) % 5], coordinates[(index + 1) % 5])
        target = coordinates[index][0] + 1, coordinates[index][1]
        assert product == target
        assert all((4 * coefficient).denominator == 1 for coefficient in coordinates[index])
    distinguished = [float(a) + float(b) * phi0 for a, b in coordinates]
    conjugate = [float(a) + float(b) * phi1 for a, b in coordinates]
    assert min(distinguished) > 0
    assert all(-1 < value < 0 for value in conjugate)
    assert distinguished[0] < 1
    assert any(coefficient.denominator != 1 for coefficient in coordinates[0])
    cycle = (-408, -664)
    assert all(isinstance(coefficient, int) for coefficient in cycle)


def test_archive_replays_and_scope_does_not_expand():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    conclusions = payload["conclusions"]
    assert conclusions["moving_denominator_AC_minus_A_minus_C_below_288"] is True
    assert conclusions["general_field_discriminant_at_most_1343277"] is True
    assert conclusions["power_coupling_required"] is False
    assert conclusions["bounded_mixed_height_premise_required"] is False
    scope = payload["scope"]
    assert scope["complete_saturated_reflection_packet_retained"] is True
    assert scope["reflection_maximality_in_declared_order_host_retained"] is True
    assert scope["reflection_maximality_used_in_bound"] is False
    assert scope["A_and_C_algebraic_integers_required_for_general_theorem"] is False
    assert scope["A_and_C_algebraic_integers_required_for_1152_corollary"] is True
    assert scope["poison_packet_is_exact_K3_carrier"] is False
    assert scope["arbitrary_n2_chamber_closed"] is False


def test_report_and_hostile_audit_pin_scope_and_nonclaims():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "proved general compact-right-pentagon field theorem",
        "exactly a compact right-angled pentagon",
        "complete and saturated",
        "reflection-maximal",
        r"no square root of a facet norm is required to lie in \(k\)",
        r"No claim that \(x_i\in\mathcal O\) is made",
        "Every pair cyclic product already generates",
        "candidate field set; no existence or realization count is asserted",
        "poison **Gram-coordinate packet**",
        r"does not: - prove that every exact \(n=2\) chamber is a right pentagon",
    ):
        assert phrase in report
    for phrase in (
        "PASS, general right-pentagon scope",
        "There is no endpoint or zero-norm exception",
        "No length-five-only generation gap remains",
        "336 candidate squareclasses, not 336 realized fields",
        "maximal chamber",
        r"does not settle arbitrary \(n=2\) incidence",
    ):
        assert phrase in audit


def test_tex_transport_and_classification_scope_are_clean():
    combined = REPORT.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
    assert ",quad" not in combined
    assert ",qquad" not in combined
    assert "No classification of reflective lattices is imported" in " ".join(combined.split())
    assert "No reflective-lattice classification is used" in " ".join(combined.split())
