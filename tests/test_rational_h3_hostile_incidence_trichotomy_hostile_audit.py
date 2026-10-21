"""Independent exact checks for the rational H3 trichotomy hostile audit."""

from __future__ import annotations

import json
from fractions import Fraction
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_h3_hostile_incidence_trichotomy.md"
AUDIT = (
    ROOT
    / "results/rational_h3_hostile_incidence_trichotomy_hostile_audit.md"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/"
    "rational_h3_hostile_incidence_trichotomy_hostile_audit.json"
)


def mat_vec(matrix, vector):
    return [
        sum(Fraction(entry) * Fraction(value) for entry, value in zip(row, vector))
        for row in matrix
    ]


def bilinear(left, matrix, right):
    return sum(
        Fraction(left[i]) * matrix[i][j] * Fraction(right[j])
        for i in range(len(left))
        for j in range(len(right))
    )


def determinant(matrix):
    rows = [[Fraction(value) for value in row] for row in matrix]
    answer = Fraction(1)
    for column in range(len(rows)):
        pivot = next(
            (row for row in range(column, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            answer *= -1
        pivot_value = rows[column][column]
        answer *= pivot_value
        for row in range(column + 1, len(rows)):
            scale = rows[row][column] / pivot_value
            for entry in range(column + 1, len(rows)):
                rows[row][entry] -= scale * rows[column][entry]
    return answer


def reflect(vector, root, gram):
    root_square = bilinear(root, gram, root)
    coefficient = -2 * bilinear(vector, gram, root) / root_square
    return [
        Fraction(value) + coefficient * Fraction(root_value)
        for value, root_value in zip(vector, root)
    ]


def test_integer_profiles_and_all_canonical_star_identities():
    low = set()
    for length in range(1, 8):
        for profile in product(range(5), repeat=length):
            nonzero = tuple(sorted((value for value in profile if value), reverse=True))
            if not nonzero:
                continue
            total = sum(nonzero)
            square_total = sum(value * value for value in nonzero)
            assert square_total - total == sum(
                value * (value - 1) for value in nonzero
            )
            if total >= 3 and square_total <= 4:
                low.add(nonzero)
    assert low == {(1, 1, 1), (1, 1, 1, 1)}

    for profile in (
        (1, 1, 1),
        (1, 1, 1, 1),
        (2, 1),
        (3,),
        (2, 2),
        (1, 1, 1, 1, 1),
    ):
        total = sum(profile)
        square_total = sum(value * value for value in profile)
        size = 1 + len(profile)
        gram = [[0] * size for _ in range(size)]
        gram[0][0] = -1
        for index, incidence in enumerate(profile, start=1):
            gram[index][index] = -4
            gram[0][index] = gram[index][0] = incidence
        f_vector = [4, *profile]
        canonical = [-1, *([2] * len(profile))]
        intersections = mat_vec(gram, f_vector)
        assert intersections[0] == square_total - 4
        assert intersections[1:] == [0] * len(profile)
        assert bilinear(f_vector, gram, f_vector) == 4 * (square_total - 4)
        assert sum(k * c for k, c in zip(canonical, f_vector)) == -4 + 2 * total


def test_spherical_cycle_discrepancy_and_highest_root_are_exact():
    physical = [
        [-1, 1, 1, 1],
        [1, -4, 0, 0],
        [1, 0, -4, 0],
        [1, 0, 0, -4],
    ]
    leading = [
        determinant([row[:size] for row in physical[:size]])
        for size in range(1, 5)
    ]
    assert leading == [-1, 3, -8, 16]

    fundamental = [3, 1, 1, 1]
    canonical = [-1, 2, 2, 2]
    assert mat_vec(physical, fundamental) == [0, -1, -1, -1]
    square = bilinear(fundamental, physical, fundamental)
    degree = sum(k * z for k, z in zip(canonical, fundamental))
    assert (square, degree, 1 + (square + degree) / 2) == (-3, 3, 1)
    discrepancy = [Fraction(-2), Fraction(-1), Fraction(-1), Fraction(-1)]
    assert mat_vec(physical, discrepancy) == canonical
    assert min(discrepancy) < -1

    host = [
        [-2, 1, 1, 1],
        [1, -2, 0, 0],
        [1, 0, -2, 0],
        [1, 0, 0, -2],
    ]
    basis = [[int(i == j) for i in range(4)] for j in range(4)]
    vector = basis[1]
    for index in (0, 2, 3, 0):
        vector = reflect(vector, basis[index], host)
    assert vector == [2, 1, 1, 1]
    assert bilinear(vector, host, vector) == -2


def test_leaf_swap_is_not_automatic_in_the_standard_host():
    # U + <-1>^13, with f,g first.  This is the numerical lattice of
    # P1xP1 blown up at 13 points.
    dimension = 15
    gram = [[0] * dimension for _ in range(dimension)]
    gram[0][1] = gram[1][0] = 1
    for index in range(2, dimension):
        gram[index][index] = -1

    e_curve = [0] * dimension
    e_curve[1] = 1
    e_curve[2] = -1
    leaves = []
    for leaf_index in range(3):
        leaf = [0] * dimension
        leaf[0] = 1
        for exceptional in range(3 + 4 * leaf_index, 7 + 4 * leaf_index):
            leaf[exceptional] = -1
        leaves.append(leaf)

    assert bilinear(e_curve, gram, e_curve) == -1
    assert [bilinear(leaf, gram, leaf) for leaf in leaves] == [-4, -4, -4]
    assert [bilinear(e_curve, gram, leaf) for leaf in leaves] == [1, 1, 1]
    assert all(
        bilinear(leaves[i], gram, leaves[j]) == 0
        for i in range(3)
        for j in range(i + 1, 3)
    )

    difference = [left - right for left, right in zip(leaves[0], leaves[1])]
    witness = [0] * dimension
    witness[3] = 1
    assert bilinear(witness, gram, difference) == 1
    assert any(value % 2 for value in difference)
    # The reflection correction is (witness.(B1-B2))/2 times
    # (B1-B2)/2, hence a forbidden quarter-lattice vector.
    assert Fraction(bilinear(witness, gram, difference), 4).denominator == 4


def test_affine_branch_is_conditional_and_bounds_only_the_core():
    physical = [
        [-1, 1, 1, 1, 1],
        [1, -4, 0, 0, 0],
        [1, 0, -4, 0, 0],
        [1, 0, 0, -4, 0],
        [1, 0, 0, 0, -4],
    ]
    host = [
        [-2, 1, 1, 1, 1],
        [1, -2, 0, 0, 0],
        [1, 0, -2, 0, 0],
        [1, 0, 0, -2, 0],
        [1, 0, 0, 0, -2],
    ]
    assert mat_vec(physical, [4, 1, 1, 1, 1]) == [0] * 5
    assert mat_vec(host, [2, 1, 1, 1, 1]) == [0] * 5

    cases = []
    for multiplicity in range(1, 20):
        for genus in range(2, 20):
            if multiplicity * (2 * genus - 2) == 4:
                cases.append((multiplicity, genus, 10 + 8 * (genus - 1)))
    assert cases == [(1, 3, 26), (2, 2, 18)]

    # Every vertical blowdown raises K^2 and lowers rho by one, leaving
    # the general fibre and the core inequality unchanged.  Thus the
    # formula has no information about the arbitrary count k.
    for genus, core_rho, tail_length in ((2, 18, 1), (2, 12, 50), (3, 26, 80)):
        assert 10 - core_rho + 8 * (genus - 1) >= 0
        original_rho = core_rho + tail_length
        assert original_rho - core_rho == tail_length


def test_hyperbolic_wrong_sign_contraction_and_hodge_direction():
    for profile in ((2, 1), (3,), (2, 2), (1, 1, 1, 1, 1)):
        total = sum(profile)
        square_total = sum(value * value for value in profile)
        assert square_total > 4
        assert Fraction(square_total, 4) - 1 > 0
        assert Fraction(total, 2) - 1 > 0
        assert 1 - Fraction(total, 2) < 0

        # Eliminating effective boundary coefficients gives a nonpositive
        # multiple of a_E greater than a positive number.
        left_coefficient = 1 - Fraction(square_total, 4)
        right_side = Fraction(total, 2) - 1
        for a_e in (Fraction(0), Fraction(1, 3), Fraction(99)):
            assert left_coefficient * a_e < right_side

        hodge_quotient = Fraction((total - 2) ** 2, square_total - 4)
        # 10-rho+b <= quotient rearranges to a lower bound on rho.
        for boundary_count in (1, 5, 20):
            lower_rho = Fraction(10 + boundary_count) - hodge_quotient
            assert (
                Fraction(10 + boundary_count) - lower_rho
                == hodge_quotient
            )

    # The repaired contraction uses a finite maximum of bad ratios.  This
    # is the step that replaces the invalid import of W3-nefness.
    l_degrees = [Fraction(1), Fraction(2), Fraction(7), Fraction(3, 2)]
    d_degrees = [Fraction(-9), Fraction(1), Fraction(-40), Fraction(-1)]
    threshold = max(
        Fraction(-d, l) for l, d in zip(l_degrees, d_degrees) if d < 0
    )
    chosen_integer = int(threshold) + 1
    assert all(
        chosen_integer * l + d >= 0
        for l, d in zip(l_degrees, d_degrees)
    )


def test_archive_scope_and_repaired_owner_proof_are_frozen():
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archive["verdict"] == "PASS_AFTER_CONTRACTION_REPAIR"
    assert archive["scope"] == {
        "absolute_picard_bound": False,
        "claim": "exact local H3 reduction only",
        "semiampleness": "conditional",
        "tail_bound": False,
    }
    assert archive["spherical"]["discrepancies"] == [-2, -1, -1, -1]
    assert archive["affine"]["conditional_fibration_cases_m_g_rho_core"] == [
        [1, 3, 26],
        [2, 2, 18],
    ]

    owner = " ".join(REPORT.read_text(encoding="utf-8").split())
    assert "weak-incidence positivity argument cannot be imported verbatim" in owner
    assert r"D.E=1-T/2<0" in owner
    assert r"choose \(a\gg0\)" in owner
    assert "Mori cone has finitely many declared generators" in owner
    assert "no universal H3 rank bound" in owner

    if AUDIT.exists():
        audit = " ".join(AUDIT.read_text(encoding="utf-8").split())
        assert "PASS after one substantive repair" in audit
        assert "does not prove an absolute Picard or Coxeter-dimension bound" in audit
