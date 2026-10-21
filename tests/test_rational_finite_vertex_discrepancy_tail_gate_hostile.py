from fractions import Fraction
from pathlib import Path

from scripts.certify_rational_finite_vertex_discrepancy_tail_gate import (
    _determinant,
    _solve,
    build_certificate,
    discrepancy_vector,
    exceptional_tail,
    hostile_controls,
    weak_barrier_residual,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_finite_vertex_discrepancy_tail_gate.md"
SCRIPT = ROOT / "scripts/certify_rational_finite_vertex_discrepancy_tail_gate.py"


def _negative_definite_by_sylvester(matrix: list[list[int]]) -> bool:
    return all(
        ((-1) ** size)
        * _determinant([row[:size] for row in matrix[:size]])
        > 0
        for size in range(1, len(matrix) + 1)
    )


def test_m_matrix_inverse_sign_and_barrier_equation_are_not_reversed():
    species = ["B", "E", "R", "R"]
    intersections = {(0, 1): 1, (1, 2): 1, (2, 3): 1}
    q_matrix, discrepancy = discrepancy_vector(species, intersections)
    a_matrix = [[-value for value in row] for row in q_matrix]
    assert _negative_definite_by_sylvester(q_matrix)

    inverse_columns = [
        _solve(a_matrix, [int(i == j) for i in range(len(species))])
        for j in range(len(species))
    ]
    assert all(value >= 0 for column in inverse_columns for value in column)

    residual = weak_barrier_residual(species, intersections)
    barrier = [Fraction(1, 2) if item == "B" else Fraction(0) for item in species]
    assert [a + u for a, u in zip(discrepancy, barrier, strict=True)] == _solve(
        a_matrix, residual
    )
    assert all(value >= 0 for value in residual)


def test_every_hostile_control_is_negative_definite_and_has_exact_discrepancy():
    expected = {
        "H1_BB_intersection_2": [Fraction(-1), Fraction(-1)],
        "H2_BR_intersection_2": [Fraction(-1), Fraction(-1)],
        "H3_E_three_B_leaves": [
            Fraction(-2),
            Fraction(-1),
            Fraction(-1),
            Fraction(-1),
        ],
    }
    for row in hostile_controls():
        assert _negative_definite_by_sylvester(row["intersection_matrix"])
        assert list(map(Fraction, row["discrepancies"])) == expected[row["id"]]
        assert any(value < 0 for value in map(Fraction, row["weak_barrier_residual"]))


def test_snc_higher_blowup_bound_and_non_snc_tangency_poison():
    # For an SNC point the new log discrepancy is 2 (no component), A+1
    # (one component), or A+B (two components).  Starting at A,B >= 1/2
    # cannot lower the minimum.
    half = Fraction(1, 2)
    assert Fraction(2) >= half
    assert half + 1 >= half
    assert half + half >= half

    # Two curves of discrepancy -1 tangent to order two are not an SNC
    # exceptional divisor.  Resolving first the tangency and then the resulting
    # triple point creates discrepancies -1 and -2, respectively.
    first = 1 + Fraction(-1) + Fraction(-1)
    second = 1 + Fraction(-1) + Fraction(-1) + first
    assert first == -1
    assert second == -2

    text = REPORT.read_text(encoding="utf-8")
    assert "ordinary tangency" in text
    assert "That non-SNC realization is not log canonical" in text
    assert "0<\\epsilon<\\frac12" in text
    assert "not uniformly \\(1/2\\)-klt" in text


def test_arbitrary_tail_quadratic_correction_and_integral_reflections():
    for length in (1, 2, 3, 17, 37):
        row = exceptional_tail(length)
        matrix = row["intersection_matrix"]
        discrepancy = list(map(Fraction, row["discrepancies"]))
        assert _negative_definite_by_sylvester(matrix)
        correction = -sum(
            discrepancy[i] * matrix[i][j] * discrepancy[j]
            for i in range(length)
            for j in range(length)
        )
        assert correction == length
        assert row["discrepancy_correction"] == correction
        assert row["canonical_square_is_preserved_by_contraction"] is True

        doubled = [[2 * value for value in values] for values in matrix]
        for root in range(length):
            for basis_vector in range(length):
                coefficient = Fraction(
                    2 * doubled[basis_vector][root], doubled[root][root]
                )
                assert coefficient.denominator == 1

        if length >= 3:
            rr_product = Fraction(
                2 * doubled[0][1], doubled[0][0]
            ) * Fraction(2 * doubled[1][0], doubled[1][1])
            assert rr_product == 1  # Coxeter order three.
        if length >= 2:
            last = length - 1
            re_product = Fraction(
                2 * doubled[last - 1][last], doubled[last - 1][last - 1]
            ) * Fraction(
                2 * doubled[last][last - 1], doubled[last][last]
            )
            assert re_product == 2  # Coxeter order four.


def test_tail_certificate_computes_rather_than_hardcodes_k_square_correction():
    source = SCRIPT.read_text(encoding="utf-8")
    assert "correction = -sum(" in source
    assert 'row["discrepancy_correction"] = correction.numerator' in source
    payload = build_certificate()
    assert all(
        row["discrepancy_correction"]
        == row["canonical_square_loss_before_contraction"]
        for row in payload["exceptional_tail_controls"]
    )


def test_scope_is_local_q_gorenstein_and_does_not_claim_hostile_core_bounds():
    text = REPORT.read_text(encoding="utf-8")
    assert "projective birational morphism" in text
    assert "normal\n\\(\\mathbb Q\\)-Gorenstein surface" in text
    assert "exceptional locus is exactly" in text
    assert "This note proves no log-MMP core bound in the hostile lanes" in text
    assert "entire* nef cone is a Coxeter chamber" in text
    assert "finite-volume or reflection-maximal" in text
    assert "bounds the resulting\nlog-minimal core" not in text
