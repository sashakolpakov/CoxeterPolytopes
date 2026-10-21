"""Independent hostile replay of the scoped power-coupled K3 obstruction."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_power_coupled_right_pentagon_field_obstruction.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_power_coupled_right_pentagon_field_obstruction.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_power_coupled_right_pentagon_field_obstruction.json"
OWNER_TESTS = ROOT / "tests/test_k3_n2_power_coupled_right_pentagon_field_obstruction.py"
AUDIT = ROOT / "results/k3_n2_power_coupled_right_pentagon_field_obstruction_hostile_scope_audit.md"
PINS = {
    REPORT: "29109ff20a0edfd327f1eb1c2cb20eeff5ebd0c39e512d2880ff357813ce8f01",
    CERTIFIER: "cb506082db0fe031e051b7769101448250f3473626fa3542ca443bac91ef3b6c",
    ARCHIVE: "d6748868987fa800d0bfcce98dc7995e9b95b8c269a625e089ec2718dcd2d98f",
    OWNER_TESTS: "42273f5a266f07a6c1af547fc47560267d08d8c23c547135b2360679c312572e",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_multiply(left, right):
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def matrix_subtract(left, right):
    return tuple(
        tuple(a - b for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def matrix_trace(matrix):
    return sum(matrix[index][index] for index in range(len(matrix)))


def independent_closure(A: Fraction, C: Fraction):
    delta = A * C - A - C
    B = (A - 1) * (C - 1)
    D = A * (C - 1) / delta
    E = C * (A - 1) / delta
    minors = (
        C * E - C - D - E + 1,
        -A + B * E - B - E + 1,
        B * D - B - C - D + 1,
        A * D - A - D - E + 1,
        A * C - A - B - C + 1,
    )
    return B, D, E, minors


def exact_envelope(trace: int, norm: int) -> bool:
    discriminant = trace * trace - 4 * norm
    r = trace - 1
    remaining = 69 - r * r - discriminant
    return remaining > 0 and remaining * remaining > 4 * r * r * discriminant


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_all_closure_equations_rederive_without_owner_code():
    tested = 0
    for a_denominator in range(1, 9):
        for c_denominator in range(1, 9):
            for a_numerator in range(a_denominator + 1, 7 * a_denominator):
                for c_numerator in range(c_denominator + 1, 7 * c_denominator):
                    A = Fraction(a_numerator, a_denominator)
                    C = Fraction(c_numerator, c_denominator)
                    if A * C - A - C <= 0:
                        continue
                    B, D, E, minors = independent_closure(A, C)
                    assert minors == (0, 0, 0, 0, 0)
                    assert min(B, D, E) > 1
                    assert E - 1 == A / (A * C - A - C)
                    assert D - 1 == C / (A * C - A - C)
                    tested += 1
    assert tested > 30_000


def test_power_locus_fixed_divisor_and_conjugate_signs_independently():
    for A in (Fraction(3, 2), Fraction(5, 3), Fraction(2), Fraction(7, 3), Fraction(3)):
        for q in range(2, 25):
            C = A**q
            u = A**q - A ** (q - 1) - 1
            if u <= 0:
                continue
            _, D, E, _ = independent_closure(A, C)
            assert A * C - A - C == A * u
            assert E == 1 + 1 / u
            assert D == 1 + A ** (q - 1) / u
            assert u * (4 * E - 4) == 4
    for denominator in range(2, 70):
        for numerator in range(1, denominator):
            conjugate = Fraction(numerator, denominator)
            for q in (2, 3, 8, 31):
                u1 = conjugate**q - conjugate ** (q - 1) - 1
                v1 = 4 / u1
                assert -2 < u1 < -1
                assert -4 < v1 < -2


def test_pair_cycle_isolated_by_reflection_product_trace():
    identity = (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    )
    for a_ij in (Fraction(-1, 3), Fraction(-2), Fraction(5, 7)):
        for a_ji in (Fraction(-3, 5), Fraction(4), Fraction(7, 2)):
            # P_i=e_i alpha_i and P_j=e_j alpha_j with alpha_i(e_i)=2.
            p_i = (
                (Fraction(2), a_ij, Fraction(0)),
                (Fraction(0), Fraction(0), Fraction(0)),
                (Fraction(0), Fraction(0), Fraction(0)),
            )
            p_j = (
                (Fraction(0), Fraction(0), Fraction(0)),
                (a_ji, Fraction(2), Fraction(0)),
                (Fraction(0), Fraction(0), Fraction(0)),
            )
            s_i = matrix_subtract(identity, p_i)
            s_j = matrix_subtract(identity, p_j)
            assert matrix_trace(matrix_multiply(s_i, s_j)) == -1 + a_ij * a_ji


def test_norm_divisor_and_height_direction_have_no_sign_reversal():
    for n in (1, 2, 4, 8, 16):
        for v1_abs in (Fraction(9, 4), Fraction(5, 2), Fraction(3), Fraction(15, 4)):
            v0 = Fraction(n, 1) / v1_abs
            u0 = 4 / v0
            assert 2 < v1_abs < 4
            assert u0 < Fraction(16, n)
    assert [n for n in range(1, 17) if 16 % n == 0] == [1, 2, 4, 8, 16]


def test_five_field_box_is_recomputed_independently():
    rows = []
    for trace in range(2, 6):
        for norm in range(1, 5):
            discriminant = trace * trace - 4 * norm
            if discriminant <= 0 or isqrt(discriminant) ** 2 == discriminant:
                continue
            if discriminant <= (trace - 2) ** 2:
                continue
            if not exact_envelope(trace, norm):
                continue
            radicand = discriminant
            square = 2
            while square * square <= radicand:
                while radicand % (square * square) == 0:
                    radicand //= square * square
                square += 1
            rows.append((trace, norm, discriminant, radicand))
    assert rows == [
        (3, 1, 5, 5),
        (4, 1, 12, 3),
        (4, 2, 8, 2),
        (5, 2, 17, 17),
        (5, 3, 13, 13),
    ]


def test_archive_replays_and_does_not_expand_scope():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["conclusions"]["possible_field_radicands"] == [2, 3, 5, 13, 17]
    scope = payload["scope"]
    assert scope["reflection_maximality_in_declared_order_host_retained"] is True
    assert scope["reflection_maximality_used_to_derive_height_bound"] is False
    assert scope["additional_power_coupling_C_equals_A_to_q_required"] is True
    assert scope["additional_A_algebraic_integer_field_generator_required"] is True
    assert scope["arbitrary_independent_A_C_coupling_closed"] is False
    assert scope["existence_of_a_carrier_in_the_ansatz_asserted"] is False
    assert scope["uniform_mixed_correspondence_degree_derived"] is False


def test_report_and_audit_keep_exact_hypotheses_and_counter_scope():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "The exact carrier package retained",
        "The additional power-coupled ansatz",
        "minimal complex projective K3 surface",
        "exactly a compact right-angled pentagon",
        "complete and saturated",
        "reflection-maximal",
        "Only order integrality is needed",
        r"algebraic integrality of \(A\) is an actual extra hypothesis",
        "no existence of such a K3 pentagon is asserted",
        "Reflection maximality is not used",
        r"independently, without the relation \(C=A^q\)",
        "The exact surviving problem is to control this moving denominator",
    ):
        assert phrase in report
    for phrase in (
        "PASS, scoped",
        r"no square root of a facet norm is assumed to lie in \(k\)",
        "identity between algebraic integers",
        "No inequality direction reverses",
        "necessary five-field list",
        r"does not solve the general \(n=2\) automatic-height problem",
    ):
        assert phrase in audit


def test_tex_transport_and_classification_scope_are_clean():
    combined = REPORT.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
    assert ",quad" not in combined
    assert ",qquad" not in combined
    assert "No reflective-lattice classification is used" in combined
    assert "No reflective-lattice classification enters" in combined
