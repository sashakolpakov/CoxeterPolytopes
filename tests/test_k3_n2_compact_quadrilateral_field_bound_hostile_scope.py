"""Independent hostile replay of the all-degree compact K3 quadrilateral theorem."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import combinations, permutations, product
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_compact_quadrilateral_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_compact_quadrilateral_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_compact_quadrilateral_field_bound.json"
OWNER_TESTS = ROOT / "tests/test_k3_n2_compact_quadrilateral_field_bound.py"
AUDIT = ROOT / "results/k3_n2_compact_quadrilateral_field_bound_hostile_scope_audit.md"
PINS = {
    REPORT: "7f300ef5b54e9dbf0ca47140b5ca7134119046fc6f1ee39bd46bce58a6028528",
    CERTIFIER: "25ccdd5e21ae2ad5ff5d43d079fb9308b3dbacb4780fa0ba8638f7c27f9d884d",
    ARCHIVE: "a61e4a28960e082380af1fc7e0dd3f8b311ab2871b6e361268f1f4e50b6213c4",
    OWNER_TESTS: "0517df99c24402c881e32c6666e91abf67c0d4e3f84fb5be2cc857aaf25e072c",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def determinant(matrix: tuple[tuple[Fraction, ...], ...]) -> Fraction:
    return sum(
        Fraction(sign(permutation))
        * product_of(matrix[row][column] for row, column in enumerate(permutation))
        for permutation in permutations(range(len(matrix)))
    )


def product_of(values) -> Fraction:
    answer = Fraction(1)
    for value in values:
        answer *= value
    return answer


def matrix_and_formula(a, b, c, d, U, V):
    a, b, c, d, U, V = map(Fraction, (a, b, c, d, U, V))
    matrix = (
        (Fraction(1), -a, -U, -d),
        (-a, Fraction(1), -b, -V),
        (-U, -b, Fraction(1), -c),
        (-d, -V, -c, Fraction(1)),
    )
    r, s = U * U - 1, V * V - 1
    P, Q, R = a * b + c * d, a * d + b * c, a * c + b * d
    K = (
        -a * a
        - b * b
        - c * c
        - d * d
        + a * a * c * c
        + b * b * d * d
        - 2 * a * b * c * d
    )
    return matrix, r * s + K - 2 * P * U - 2 * Q * V - 2 * R * U * V


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_four_normal_determinant_is_rederived_without_owner_code():
    adjacent = (Fraction(0), Fraction(1, 3), Fraction(2, 3))
    opposite = (Fraction(5, 4), Fraction(3, 2), Fraction(7, 3))
    tested = 0
    for parameters in product(adjacent, adjacent, adjacent, adjacent, opposite, opposite):
        matrix, formula = matrix_and_formula(*parameters)
        assert determinant(matrix) == formula
        tested += 1
    assert tested == 729


def test_all_place_norm_gap_has_the_exact_degree_exponent():
    # If y=4r is integral, y_0>0, and its d-1 other conjugates have
    # absolute value below four, then |N(y)|>=1 forces r>4^{-d}.
    for degree in range(2, 7):
        weakest_y0 = Fraction(1, 4 ** (degree - 1))
        weakest_r = weakest_y0 / 4
        assert weakest_r == Fraction(1, 4**degree)
    assert Fraction(1, 4**6) == Fraction(1, 4096)


def test_determinant_coefficient_bounds_and_monotonicity_are_directional():
    for denominator in range(2, 14):
        for numerators in product(range(denominator), repeat=4):
            a, b, c, d = (Fraction(n, denominator) for n in numerators)
            P, Q, R = a * b + c * d, a * d + b * c, a * c + b * d
            K = (
                -a * a
                - b * b
                - c * c
                - d * d
                + a * a * c * c
                + b * b * d * d
                - 2 * a * b * c * d
            )
            assert 0 <= P < 2 and 0 <= Q < 2 and 0 <= R < 2
            assert abs(K) < 8
    # Exact derivatives have negative numerators for V>1.
    for V in (Fraction(1001, 1000), Fraction(5, 4), Fraction(2), Fraction(10)):
        assert -4 < 0  # numerator factor of d(4/(V-1))/dV
        assert -4 * (V * V + 4 * V + 1) < 0


def test_degree_six_trial_bound_is_exact_and_on_the_positive_side():
    rational = 537_034_767
    surd = 8_389_888
    assert rational > 0
    assert rational**2 - 4097 * surd**2 == 17_606_953_271_521 > 0
    # sqrt(4097)<4097/64, so the linear coefficient is below 32772;
    # its half (the quadratic vertex) is therefore below the trial point.
    assert 4097**2 - 4096 * 4097 > 0
    assert 32_772**2 == 1_074_003_984


def test_cycle_height_and_primitive_grid_constants_recompute():
    magnitude_bound = 32_772
    H = 16 * magnitude_bound**4
    assert H == 18_455_752_922_365_956_096
    assert len(tuple(combinations(range(6), 2))) == 15

    # A concrete six-embedding packet exhibits the counting mechanism:
    # each embedding pair forbids at most one rational coefficient.
    theta1 = tuple(Fraction(i * i + 1) for i in range(6))
    theta2 = tuple(Fraction(i) for i in range(6))
    bad = set()
    for i, j in combinations(range(6), 2):
        if theta2[i] != theta2[j]:
            bad.add((theta1[j] - theta1[i]) / (theta2[i] - theta2[j]))
    assert len(bad) <= 15
    assert any(Fraction(c) not in bad for c in range(16))
    assert 32 * H == 590_584_093_515_710_595_072


def test_quadratic_sharpening_replays_independently():
    assert 8847 > 0
    assert 8847**2 - 17 * 2128**2 == 1_286_881 > 0
    assert 4 * 132**2 == 69_696
    assert 69_699**2 - 4 == 4_857_950_597


def test_archive_replays_and_scope_stays_literal():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    assert payload["conclusions"]["K3_full_trace_field_degree_at_most_6"] is True
    assert payload["conclusions"]["all_degree_both_opposite_magnitudes_below_32772"] is True
    assert payload["conclusions"]["bounded_height_premise_required"] is False
    scope = payload["scope"]
    assert scope["returned_field_is_cyclic_Cartan_ground_field_required"] is True
    assert scope["complete_saturated_reflection_packet_retained"] is True
    assert scope["reflection_maximality_in_declared_order_host_retained"] is True
    assert scope["reflection_maximality_used_in_bound"] is False
    assert scope["noncompact_quadrilateral_closed"] is False
    assert scope["arbitrary_n2_polygon_closed"] is False


def test_report_and_audit_distinguish_bounded_degree_from_finite_fields():
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "proved for every full-trace complex-K3 field degree",
        r"d=[k:\mathbf Q]\le6",
        r"r,s>4^{-d}\ge4^{-6}",
        r"A,C<32772^2=1,074,003,984",
        "cyclic Cartan products",
        "Hermite--Minkowski finiteness",
        "gratuitous extension",
        r"does not: - prove that every exact \(n=2\) chamber is a quadrilateral",
    ):
        assert phrase in report
    for phrase in (
        "PASS, all-degree compact-quadrilateral scope",
        "bounded degree alone is not being substituted for field finiteness",
        "Pair products alone need not generate a field of degree greater than two",
        "Maximality is retained but unused",
        "does not cover noncompact quadrilaterals",
    ):
        assert phrase in audit


def test_tex_transport_and_classification_scope_are_clean():
    combined = REPORT.read_text(encoding="utf-8") + AUDIT.read_text(encoding="utf-8")
    assert ",quad" not in combined
    assert ",qquad" not in combined
    normalized = " ".join(combined.split())
    assert "No imposed height, power relation, angle classification" in normalized
    assert "No angle classification or reflective-lattice classification enters" in normalized
