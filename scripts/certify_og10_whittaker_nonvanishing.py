#!/usr/bin/env python3
"""Exact identity-term domination for the OG10 Hilbert Poincare series.

Let ``F=Q(sqrt(5))``, ``m0=phi/sqrt(5)``, and let ``rho`` be the
176-dimensional Hilbert--Weil representation in Bruinier's totally-real
Whittaker construction.  The normalized diagonal coefficient of the
Poincare series ``P_(m0,0)`` is ``1 + E``.  This certificate proves
``|E| < 1`` using rational arithmetic only.

The analytic inputs are the standard inequalities

    |J_5(x)| <= 1,
    |J_5(x)| <= (x/2)^5 / 5!,

and the scalar Hilbert Petersson coefficient formula, with its scalar
Kloosterman sum replaced by the corresponding matrix-coefficient sum for
the unitary representation ``rho``.  The report
``results/og10_whittaker_nonvanishing.md`` gives the derivation and source
locators.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial, isqrt


PI_LOWER = Fraction(333, 106)
PI_UPPER = Fraction(355, 113)
SQRT5_LOWER = Fraction(2_236_067, 1_000_000)
SQRT5_UPPER = Fraction(2_236_068, 1_000_000)
SQRT176_LOWER = Fraction(6_633, 500)
PHI_LOWER = Fraction(8, 5)
R = 1 / PHI_LOWER

# A finite cutoff chosen so that log(B)=10 log(2), simplifying the exact
# integral tail estimate.
IDEAL_NORM_CUTOFF = 1024


def _decimal(value: Fraction, digits: int = 15) -> str:
    """Return a deterministic decimal truncation of a nonnegative fraction."""

    scale = 10**digits
    integer = value.numerator * scale // value.denominator
    return f"{integer // scale}.{integer % scale:0{digits}d}"


def _chi5(value: int) -> int:
    residue = value % 5
    if residue == 0:
        return 0
    return 1 if residue in (1, 4) else -1


def _arctan_reciprocal_bounds(
    denominator: int, last_index: int
) -> tuple[Fraction, Fraction]:
    """Alternating-series bounds for ``arctan(1/denominator)``."""

    if last_index % 2 == 0:
        raise ValueError("use an odd last index so the partial sum is a lower bound")
    lower = sum(
        (
            Fraction((-1) ** index, (2 * index + 1) * denominator ** (2 * index + 1))
            for index in range(last_index + 1)
        ),
        Fraction(0),
    )
    # The next (positive) partial sum is an upper bound.
    next_index = last_index + 1
    upper = lower + Fraction(
        1, (2 * next_index + 1) * denominator ** (2 * next_index + 1)
    )
    return lower, upper


def machin_pi_bounds() -> tuple[Fraction, Fraction]:
    """Rationally enclose pi using Machin's identity.

    ``pi = 16 atan(1/5) - 4 atan(1/239)``.  Eleven alternating-series
    terms are already far more accurate than the classical bounds used by
    the Bessel certificate.
    """

    atan5_lower, atan5_upper = _arctan_reciprocal_bounds(5, 11)
    atan239_lower, atan239_upper = _arctan_reciprocal_bounds(239, 3)
    return (
        16 * atan5_lower - 4 * atan239_upper,
        16 * atan5_upper - 4 * atan239_lower,
    )


def ideal_count(norm: int) -> int:
    """Number of ideals of norm ``norm`` in ``Q(sqrt(5))``.

    This is the coefficient ``sum_{d|norm} chi_5(d)`` of
    ``zeta_F(s)=zeta(s)L(s,chi_5)``.
    """

    total = 0
    for divisor in range(1, isqrt(norm) + 1):
        if norm % divisor:
            continue
        total += _chi5(divisor)
        quotient = norm // divisor
        if quotient != divisor:
            total += _chi5(quotient)
    return total


def _normalized_bessel_term(z: Fraction, index: int) -> Fraction:
    """The positive term after factoring ``(x/2)^5/5!`` from J_5(x)."""

    return z**index * Fraction(
        factorial(5), factorial(index) * factorial(index + 5)
    )


def central_bessel_bounds() -> tuple[Fraction, Fraction, Fraction]:
    """Certify ``0 < J_5(4*pi/sqrt(5)) < 1/3``.

    Put ``z=x^2/4=4*pi^2/5``.  After its first term, the defining
    alternating series for J_5 has decreasing terms because
    ``z/((j+1)(j+6)) < 1`` for ``j>=1``.  An odd partial sum is a lower
    bound and an even partial sum is an upper bound.  Termwise rational
    pi-enclosures then give the bounds below.
    """

    z_lower = 4 * PI_LOWER**2 / 5
    z_upper = 4 * PI_UPPER**2 / 5

    upper_normalized = Fraction(0)
    for index in range(9):  # even partial sum through j=8
        z = z_upper if index % 2 == 0 else z_lower
        sign = 1 if index % 2 == 0 else -1
        upper_normalized += sign * _normalized_bessel_term(z, index)

    lower_normalized = Fraction(0)
    for index in range(10):  # odd partial sum through j=9
        z = z_lower if index % 2 == 0 else z_upper
        sign = 1 if index % 2 == 0 else -1
        lower_normalized += sign * _normalized_bessel_term(z, index)

    prefactor_lower = (2 * PI_LOWER / SQRT5_UPPER) ** 5 / factorial(5)
    prefactor_upper = (2 * PI_UPPER / SQRT5_LOWER) ** 5 / factorial(5)
    lower = prefactor_lower * lower_normalized
    upper = prefactor_upper * upper_normalized
    return lower, upper, prefactor_upper


def nonunit_unit_orbit_bound(norm: int) -> Fraction:
    """Uniformly bound the Bessel unit orbit for an ideal of given norm.

    Write the two Bessel arguments as

        s*n^(-1/2)*phi^(t+delta), s*n^(-1/2)*phi^(-t-delta),

    with ``|delta|<=1/2``.  If ``A=(s/(2 sqrt(n)))^5/5! < 1``, the product
    of the two elementary Bessel bounds is at most ``A^2`` and also at
    most ``A phi^(-5|t+delta|)``.  We use the first bound for
    ``|t|<=K=bit_length(n)`` and the second on both geometric tails.

    The returned expression replaces ``A<3/(2 n^(5/2))``, ``phi>8/5``,
    and ``n^(-5/2)<=n^(-2)`` by rational upper bounds.
    """

    if norm < 4:
        raise ValueError("the nonunit bound is only used for norm >= 4")
    cutoff = norm.bit_length()
    central = Fraction((2 * cutoff + 1) * 9, 4 * norm**5)
    tail = (
        Fraction(3, norm**2)
        * R ** (5 * cutoff + 2)
        / (1 - R**5)
    )
    return central + tail


def nonunit_ideal_sum_bound() -> tuple[Fraction, Fraction, Fraction]:
    """Bound all nonunit lower-left rows, including the infinite tail."""

    finite = sum(
        (
            ideal_count(norm) * nonunit_unit_orbit_bound(norm)
            for norm in range(4, IDEAL_NORM_CUTOFF + 1)
        ),
        Fraction(0),
    )

    # For n>B, K=bit_length(n)<=log_2(n)+1 and (R^5)^K<n^-3 because
    # R^5=(5/8)^5<1/8 and 2^K>n.  Also a_F(n)<=d(n)<=n.  Therefore the
    # remaining sum is bounded by the integral from B to infinity of
    #
    # x^-4 [(9/4)(2 log_2 x+3) + 3 R^2/(1-R^5)].
    #
    # Since B=2^10 and log(2)>2/3, its exact upper bound is the expression
    # below: B^-3 [63/4 + C/3], C=27/4+3R^2/(1-R^5).
    constant = Fraction(27, 4) + 3 * R**2 / (1 - R**5)
    tail = Fraction(1, IDEAL_NORM_CUTOFF**3) * (
        Fraction(63, 4) + constant / 3
    )
    return finite, tail, finite + tail


def certify() -> dict[str, object]:
    central_lower, central_upper, a_upper = central_bessel_bounds()
    finite, nonunit_tail, nonunit_sum = nonunit_ideal_sum_bound()
    machin_lower, machin_upper = machin_pi_bounds()

    # For c a unit, rho(S)_(0,0) has absolute value 1/sqrt(176).  The
    # balanced t=0 Bessel product is <1/9.  For |t|>=1, one of the two
    # Bessel arguments gives a geometric bound A phi^(-5|t|), with
    # A<73/50.
    unit_orbit = Fraction(1, 9) + 2 * Fraction(73, 50) / (PHI_LOWER**5 - 1)

    scalar_fourier_prefactor = 4 * PI_UPPER**2 / SQRT5_LOWER
    unit_error = (
        scalar_fourier_prefactor / SQRT176_LOWER * unit_orbit
    )
    nonunit_error = scalar_fourier_prefactor * nonunit_sum
    total_error = unit_error + nonunit_error

    elementary_checks = {
        "machin_pi_enclosure_inside_declared_bounds": (
            PI_LOWER < machin_lower < machin_upper < PI_UPPER
        ),
        "sqrt5_interval": SQRT5_LOWER**2 < 5 < SQRT5_UPPER**2,
        "sqrt176_lower": SQRT176_LOWER**2 < 176,
        "phi_lower": PHI_LOWER**2 < PHI_LOWER + 1,
        "bessel_prefactor_below_73_over_50": a_upper < Fraction(73, 50),
        "central_bessel_positive": central_lower > 0,
        "central_bessel_below_one_third": central_upper < Fraction(1, 3),
        "no_ideals_norm_two_or_three": ideal_count(2) == ideal_count(3) == 0,
        "nonunit_norm_four_present": ideal_count(4) == 1,
        "total_off_identity_below_one": total_error < 1,
    }

    return {
        "verified": all(elementary_checks.values()),
        "theorem": (
            "P_(phi/sqrt(5),0) is nonzero; equivalently the pure norm-two "
            "Whittaker principal part is obstructed"
        ),
        "checks": elementary_checks,
        "normalization": {
            "diagonal_identity_coefficient": "1",
            "scalar_hilbert_fourier_prefactor": "(2*pi)^2/sqrt(5)",
            "unit_c_weil_matrix_coefficient_absolute_value": "1/sqrt(176)",
            "nonunit_generalized_kloosterman_over_norm_bound": "1",
        },
        "central_bessel": {
            "argument": "4*pi/sqrt(5)",
            "lower": str(central_lower),
            "upper": str(central_upper),
            "lower_decimal": _decimal(central_lower),
            "upper_decimal": _decimal(central_upper),
        },
        "unit_lower_left_rows": {
            "bessel_unit_orbit_upper": str(unit_orbit),
            "bessel_unit_orbit_upper_decimal": _decimal(unit_orbit),
            "fourier_error_upper": str(unit_error),
            "fourier_error_upper_decimal": _decimal(unit_error),
        },
        "nonunit_lower_left_rows": {
            "ideal_norm_cutoff": IDEAL_NORM_CUTOFF,
            "finite_bessel_sum_upper": str(finite),
            "finite_bessel_sum_upper_decimal": _decimal(finite),
            "integral_tail_upper": str(nonunit_tail),
            "integral_tail_upper_decimal": _decimal(nonunit_tail),
            "fourier_error_upper": str(nonunit_error),
            "fourier_error_upper_decimal": _decimal(nonunit_error),
        },
        "total_off_identity_error_upper": str(total_error),
        "total_off_identity_error_upper_decimal": _decimal(total_error),
        "consequence": {
            "diagonal_coefficient_nonzero": True,
            "poincare_series_nonzero": True,
            "exists_cusp_form_with_nonzero_target_coefficient": True,
            "single_whittaker_basis_element_weakly_holomorphic": False,
        },
        "scope": (
            "rules out the single pure norm-two Whittaker principal part; "
            "does not rule out compensating principal parts and does not "
            "decide the hyperbolic reflection index"
        ),
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2, sort_keys=True))
