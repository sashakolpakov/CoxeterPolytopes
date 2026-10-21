import json
from fractions import Fraction
from pathlib import Path

from scripts.certify_fekete_interval_minus2_1 import (
    ARCHIVE,
    P6_EXPONENTS,
    P6_SUP,
    P20_EXPONENTS,
    P20_SUP,
    build_archive,
    expanded_coefficients,
)


ROOT = Path(__file__).resolve().parents[1]


def test_exact_minimal_polynomial_data():
    assert P6_SUP == Fraction(45, 64)
    assert expanded_coefficients(P6_EXPONENTS) == [1, 3, 0, -5, -1, 2, 0]
    assert Fraction(4, 27) < P6_SUP**2
    assert P6_SUP < 1 < 2 * P6_SUP


def test_exact_degree_twenty_supremum_data():
    assert P20_EXPONENTS == (4, 4, 3, 3, 3)
    assert P20_SUP == Fraction(884_736, 9_765_625)
    assert max(
        Fraction(91_125, 1_048_576),
        Fraction(27, 1_024),
        P20_SUP,
    ) == P20_SUP


def test_one_exceptional_embedding_integer_bound():
    maximum = 210**4 * 209**3 * 208**3
    assert maximum * P20_SUP**22 > 1
    assert maximum * P20_SUP**23 < 1


def test_full_bounded_family_replay_and_archive():
    expected = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert build_archive() == expected
    assert expected["bounded_factor_family"]["row_count"] == 1_280
    assert expected["bounded_factor_family"][
        "unique_winner_exponents_a_b_c_d_e"
    ] == [4, 4, 3, 3, 3]
    assert expected["one_exceptional_embedding"]["degree_bound_for_Q_alpha"] == 23


def test_report_keeps_global_optimality_and_generator_caveats():
    report = (ROOT / "results/fekete_interval_minus2_1_probe_theorem.md").read_text(
        encoding="utf-8"
    )
    assert "optimality **only within (3)**" in report
    assert "The generator hypothesis is essential" in report
    assert "does not assert a global Coxeter" in report
