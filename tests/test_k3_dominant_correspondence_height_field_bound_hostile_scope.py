"""Standalone hostile replay of the K3 mixed-correspondence height gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
from itertools import product
from math import comb
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_dominant_correspondence_height_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_dominant_correspondence_height_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_dominant_correspondence_height_field_bound.json"
OWNER_TESTS = ROOT / "tests/test_k3_dominant_correspondence_height_field_bound.py"
AUDIT = ROOT / "results/k3_dominant_correspondence_height_field_bound_hostile_scope_audit.md"
PINS = {
    REPORT: "3c1f779bbb6a765c368ff928c41287ddb221b1c3c56b60a4f8c920dec39c01b8",
    CERTIFIER: "552b598393e4f5f17616d5a04e75bb841d2cf838c5175e672d96b59fdd6b1ee6",
    ARCHIVE: "938496ca23ad8da741d3bfa7bff9735f3a3d37e97492d3da1c045003d2b9e423",
    OWNER_TESTS: "237389fbc61c40ebaa378125e38de34bd35f65522c69b5ed23667fa18628ae15",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def multiply(left, right):
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(x * y for x, y in zip(row, column)) for column in columns)
        for row in left
    )


def test_owner_bytes_are_exactly_pinned():
    assert {path: sha256(path) for path in PINS} == PINS


def test_signed_place_identity_and_sign_are_exact():
    for degree in range(2, 11):
        for square in (Fraction(2), Fraction(9, 2), Fraction(14)):
            weights = tuple(
                -Fraction(3 * index + 1, index + 2)
                for index in range(1, degree)
            )
            lambda0 = Fraction(7 * degree + 3, 2)
            lambdas = tuple(
                lambda0 - Fraction(index, degree + 1)
                for index in range(1, degree)
            )
            distinguished_weight = square - sum(weights)
            height = lambda0 * distinguished_weight + sum(
                value * weight for value, weight in zip(lambdas, weights)
            )
            gap = height - lambda0 * square
            assert gap == sum(
                (value - lambda0) * weight
                for value, weight in zip(lambdas, weights)
            )
            assert gap >= 0

    # Dominance is indispensable: reversing it reverses the gap.
    assert (Fraction(5) - Fraction(3)) * Fraction(-2) < 0


def test_positive_primitive_grid_has_strict_counting_room():
    for degree in range(2, 13):
        hyperplanes = degree * (degree - 1) // 2
        bound = hyperplanes + 1
        for generator_count in range(1, 6):
            assert hyperplanes * bound ** (generator_count - 1) < bound**generator_count

    # Independently count proper linear hyperplanes on a small positive grid.
    bound = 5
    grid = tuple(product(range(1, bound + 1), repeat=3))
    for coefficients in product(range(-2, 3), repeat=3):
        if coefficients == (0, 0, 0):
            continue
        count = sum(
            sum(c * x for c, x in zip(coefficients, point)) == 0
            for point in grid
        )
        assert count <= bound**2


def test_pair_cauchy_and_lorentzian_dominance_are_exact():
    positive_vectors = ((1, 0, 2), (2, 3, 1), (-1, 4, 2), (5, -2, 7))
    for left in positive_vectors:
        for right in positive_vectors:
            dot = sum(x * y for x, y in zip(left, right))
            left_square = sum(x * x for x in left)
            right_square = sum(x * x for x in right)
            xi = Fraction(4 * dot * dot, left_square * right_square)
            assert 0 <= xi <= 4

    # Minkowski form diag(-1,1,1), with e=(0,1,0).
    e = (0, 1, 0)

    def lorentz(left, right):
        return -left[0] * right[0] + left[1] * right[1] + left[2] * right[2]

    values = []
    for other in ((0, 3, 4), (1, 1, 1), (2, 3, 0)):
        values.append(
            Fraction(
                4 * lorentz(e, other) ** 2,
                lorentz(e, e) * lorentz(other, other),
            )
        )
    assert values[0] < 4
    assert values[1] == 4
    assert values[2] > 4


def test_cyclic_integrality_uses_trace_not_nonprimitive_coefficients():
    # A2 Gram.  The same integral reflection may be displayed with e=2*r.
    gram = ((2, -1), (-1, 2))
    reflection_r = ((-1, 1), (0, 1))
    reflection_f = ((1, 0), (1, -1))
    product_matrix = multiply(reflection_r, reflection_f)
    assert sum(product_matrix[index][index] for index in range(2)) == -1

    e = (2, 0)
    f = (0, 1)

    def pair(left, right):
        return sum(
            left[row] * gram[row][column] * right[column]
            for row in range(2)
            for column in range(2)
        )

    a_ef = Fraction(2 * pair(e, f), pair(e, e))
    a_fe = Fraction(2 * pair(f, e), pair(f, f))
    assert a_ef == Fraction(-1, 2)  # not integral
    assert a_fe == -2
    xi = a_ef * a_fe
    assert xi == 1
    # Here n+1=2, hence n-3=-2 and trace=(n-3)+xi=-1.
    assert -2 + xi == -1


def test_k3_degree_generator_and_discriminant_constants():
    for dimension in range(2, 20):
        degree = 20 // (dimension + 1)
        if degree:
            assert degree <= 6
    for degree in range(2, 7):
        selected = 0
        product_degree = 1
        while product_degree * 2 <= degree:
            product_degree *= 2
            selected += 1
        assert selected <= 2
        assert 2**selected <= degree

    coefficient_bound = 6 * 5 // 2 + 1
    assert coefficient_bound == 16
    root_bound = 2 * coefficient_bound  # times R0
    assert root_bound == 32
    assert 6 * 5 == 30


def test_root_boxes_are_finite_and_bound_the_discriminant():
    for degree in range(1, 7):
        root_bound = 128
        coefficient_box = 1
        for index in range(1, degree + 1):
            coefficient_bound = comb(degree, index) * root_bound**index
            coefficient_box *= 2 * coefficient_bound + 1
        assert coefficient_box > 0
        assert degree * (degree - 1) <= 30
        assert root_bound ** (degree * (degree - 1)) <= root_bound**30


def test_rank_four_escape_forces_unbounded_mixed_height():
    previous = 0
    for a in range(2, 2001):
        height = 32 * a - 10
        assert height > previous
        previous = height
        # ((height/14)-a)^2 >= a^2-1, with positive left base.
        assert height > 14 * a
        assert 32 * a * a - 90 * a + 74 > 0
        assert height * height > 14 * 14 * 4 * (a * a - 1)


def test_archive_replay_and_scope_are_literal():
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(ARCHIVE.read_text(encoding="utf-8"))["payload"]
    scope = payload["scope"]
    assert scope["requires_exact_full_Picard_or_positive_H_in_trace_block"]
    assert scope["pair_coordinate_integrality_from_diagram_alone"] is False
    assert scope["proves_maximality_forces_the_height_gate"] is False
    assert scope["rules_out_counterfamilies_with_unbounded_correspondence_degree"] is False
    assert payload["conclusions"]["finite_field_set_not_only_bounded_degree"] is True


def test_report_and_audit_state_every_boundary():
    theorem = " ".join(REPORT.read_text(encoding="utf-8").split())
    audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "field-generating packet",
        r"\operatorname {tr}_k(s_is_j)=n-3+a_{ij}a_{ji}",
        "nonprimitive displayed normals",
        r"\operatorname {NS}(X)=N",
        r"If \(d=1\), the field is \(\mathbf Q\)",
        r"d\le6",
        r"r\le2",
        r"^{30}",
        "does not prove that reflection maximality supplies the bound",
        "an ample or big-and-nef class need not lie in",
    ):
        assert phrase in theorem
    for phrase in (
        "PASS after scope repairs",
        "trace identity",
        "finite field set, not merely a degree bound",
        "does not derive the scalar height bound from reflection maximality",
    ):
        assert phrase in audit
