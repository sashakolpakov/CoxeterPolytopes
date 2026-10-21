from fractions import Fraction

from scripts.certify_og10_whittaker_nonvanishing import (
    central_bessel_bounds,
    certify,
    ideal_count,
    nonunit_ideal_sum_bound,
    nonunit_unit_orbit_bound,
)


def test_exact_og10_whittaker_identity_term_domination():
    result = certify()

    assert result["verified"]
    assert all(result["checks"].values())
    assert result["normalization"]["diagonal_identity_coefficient"] == "1"
    assert (
        result["normalization"][
            "unit_c_weil_matrix_coefficient_absolute_value"
        ]
        == "1/sqrt(176)"
    )
    assert Fraction(result["total_off_identity_error_upper"]) < 1
    assert result["consequence"] == {
        "diagonal_coefficient_nonzero": True,
        "poincare_series_nonzero": True,
        "exists_cusp_form_with_nonzero_target_coefficient": True,
        "single_whittaker_basis_element_weakly_holomorphic": False,
    }


def test_rational_bessel_and_ideal_tail_subcertificates():
    bessel_lower, bessel_upper, first_term_upper = central_bessel_bounds()
    assert 0 < bessel_lower < bessel_upper < Fraction(1, 3)
    assert first_term_upper < Fraction(73, 50)

    assert [ideal_count(norm) for norm in range(1, 12)] == [
        1,
        0,
        0,
        1,
        1,
        0,
        0,
        0,
        1,
        0,
        2,
    ]
    assert nonunit_unit_orbit_bound(4) < Fraction(2, 125)

    finite, tail, total = nonunit_ideal_sum_bound()
    assert finite > 0
    assert 0 < tail < Fraction(1, 50_000_000)
    assert total == finite + tail < Fraction(11, 500)
