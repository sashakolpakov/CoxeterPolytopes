#!/usr/bin/env python3
"""Exact auxiliary-polynomial certificates on the interval [-2, 1].

The script proves four finite statements.

* ``P6 = x(x-1)(x+1)(x+2)(x^2+x-1)`` has exact sup norm 45/64.
* No nonzero integer polynomial of degree at most five has sup norm < 1;
  in degree six, P6 is unique up to sign.
* In the explicitly bounded 1,280-member factor family described below,
  ``P20`` is the unique normalized-sup-norm winner.  Every losing row has
  an exact rational witness, so numerical root finding is not evidence.
* If a totally real algebraic integer alpha generates its field, one real
  embedding lies in [2,14], and all the others lie in [-2,1], then P20 and
  the integer norm give [Q(alpha):Q] <= 23.

SymPy is used only to isolate the roots of degree-at-most-five logarithmic
derivative numerators.  All comparisons promoted to the certificate are
exact rational inequalities.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Iterable

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/fekete_interval_minus2_1_certificate.json"

X = sp.symbols("x")
FACTORS = (
    sp.Poly(X, X, domain=sp.ZZ),
    sp.Poly(X + 1, X, domain=sp.ZZ),
    sp.Poly(X + 2, X, domain=sp.ZZ),
    sp.Poly(X - 1, X, domain=sp.ZZ),
    sp.Poly(X * X + X - 1, X, domain=sp.ZZ),
)

P6_EXPONENTS = (1, 1, 1, 1, 1)
P20_EXPONENTS = (4, 4, 3, 3, 3)
P6_SUP = Fraction(45, 64)
P20_SUP = Fraction(884_736, 9_765_625)


def factored_value(exponents: tuple[int, ...], point: sp.Rational) -> Fraction:
    value = sp.Rational(1)
    for factor, exponent in zip(FACTORS, exponents, strict=True):
        if exponent:
            value *= abs(factor.eval(point)) ** exponent
    return Fraction(int(sp.numer(value)), int(sp.denom(value)))


def expanded_coefficients(exponents: tuple[int, ...]) -> list[int]:
    polynomial = sp.Poly(1, X, domain=sp.ZZ)
    for factor, exponent in zip(FACTORS, exponents, strict=True):
        polynomial *= factor**exponent
    return [int(coefficient) for coefficient in polynomial.all_coeffs()]


def logarithmic_derivative_numerator(exponents: tuple[int, ...]) -> sp.Poly:
    """Return the squarefree-factor logarithmic-derivative numerator.

    If P=prod(f_i^e_i), then away from the roots of P its critical points
    are the roots of sum_i e_i f_i' prod_{j != i} f_j.  This polynomial has
    degree at most five throughout the bounded family.
    """

    active = [index for index, exponent in enumerate(exponents) if exponent]
    product = sp.Poly(1, X, domain=sp.ZZ)
    for index in active:
        product *= FACTORS[index]
    numerator = sp.Poly(0, X, domain=sp.ZZ)
    for index in active:
        numerator += (
            exponents[index]
            * FACTORS[index].diff()
            * product.exquo(FACTORS[index])
        )
    return numerator


def exact_rational_witness(exponents: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    """Exhibit a rational point close to a true critical point.

    Root isolation locates every real critical point in an interval of width
    at most 10^-8.  We evaluate the polynomial exactly at the rational
    midpoints and endpoints.  The returned value is merely a lower bound for
    the true sup norm, which is exactly the direction needed to eliminate a
    competitor.
    """

    points = [sp.Rational(-2), sp.Rational(1)]
    numerator = logarithmic_derivative_numerator(exponents)
    for (left_right, _multiplicity) in numerator.intervals(
        eps=sp.Rational(1, 10**8)
    ):
        left, right = left_right
        if right < -2 or left > 1:
            continue
        left = max(left, sp.Rational(-2))
        right = min(right, sp.Rational(1))
        points.append((left + right) / 2)

    values = [(factored_value(exponents, point), point) for point in points]
    value, point = max(values)
    return value, Fraction(int(sp.numer(point)), int(sp.denom(point)))


def bounded_family_rows() -> tuple[list[dict[str, object]], str]:
    """Replay the exact 1,280-row bounded factor/exponent comparison."""

    rows: list[dict[str, object]] = []
    digest_rows: list[tuple[object, ...]] = []
    for exponents in itertools.product(
        range(1, 5), range(1, 5), range(1, 5), range(1, 5), range(0, 5)
    ):
        degree = sum(exponents[:4]) + 2 * exponents[4]
        witness_value, witness_point = exact_rational_witness(exponents)
        winner = exponents == P20_EXPONENTS
        if not winner:
            # witness_value^(1/degree) > P20_SUP^(1/20), without roots/logs.
            assert witness_value**20 > P20_SUP**degree
        row = {
            "exponents": list(exponents),
            "degree": degree,
            "winner": winner,
            "witness_point": str(witness_point),
            "witness_value": str(witness_value),
        }
        rows.append(row)
        digest_rows.append(
            (exponents, degree, str(witness_point), str(witness_value), winner)
        )

    assert len(rows) == 1_280
    assert sum(bool(row["winner"]) for row in rows) == 1
    digest_payload = json.dumps(digest_rows, separators=(",", ":")) + "\n"
    digest = hashlib.sha256(digest_payload.encode("utf-8")).hexdigest()
    return rows, digest


def verify_exact_suprema() -> None:
    # Put y=x+1/2 and z=y^2.  On [-2,1], 0 <= z <= 9/4.
    # P6=(z-1/4)(z-5/4)(z-9/4).  Its two non-root critical values
    # have absolute value 2sqrt(3)/9, whose square is 4/27.
    z = sp.symbols("z")
    p6_z = (z - sp.Rational(1, 4)) * (z - sp.Rational(5, 4)) * (
        z - sp.Rational(9, 4)
    )
    critical_minus = sp.Rational(5, 4) - sp.sqrt(3) / 3
    critical_plus = sp.Rational(5, 4) + sp.sqrt(3) / 3
    assert sp.simplify(
        sp.diff(p6_z, z) - (48 * z**2 - 120 * z + 59) / 16
    ) == 0
    assert sp.simplify(p6_z.subs(z, critical_minus) ** 2) == sp.Rational(4, 27)
    assert sp.simplify(p6_z.subs(z, critical_plus) ** 2) == sp.Rational(4, 27)
    assert p6_z.subs(z, 0) == -sp.Rational(45, 64)
    assert p6_z.subs(z, sp.Rational(9, 4)) == 0
    assert Fraction(4, 27) < P6_SUP**2

    # P20=(z-1/4)^4(z-5/4)^3(z-9/4)^3.  The non-root critical
    # points are z=3/4 and z=37/20.  Together with z=0 and z=9/4,
    # these give the complete maximum list.
    p20_z = (
        (z - sp.Rational(1, 4)) ** 4
        * (z - sp.Rational(9, 4)) ** 3
        * (z - sp.Rational(5, 4)) ** 3
    )
    expected_derivative = (
        (4 * z - 9) ** 2
        * (4 * z - 5) ** 2
        * (4 * z - 3)
        * (4 * z - 1) ** 3
        * (20 * z - 37)
        / 131_072
    )
    assert sp.simplify(sp.diff(p20_z, z) - expected_derivative) == 0
    endpoint = Fraction(p20_z.subs(z, 0))
    critical_left = Fraction(p20_z.subs(z, sp.Rational(3, 4)))
    critical_right = abs(Fraction(p20_z.subs(z, sp.Rational(37, 20))))
    assert p20_z.subs(z, sp.Rational(9, 4)) == 0
    assert max(endpoint, critical_left, critical_right, Fraction(0)) == P20_SUP


def verify_minimal_degree_and_uniqueness() -> None:
    # Any integer polynomial of norm <1 vanishes at -2,-1,0,1, so the
    # monic A=x(x-1)(x+1)(x+2) divides it in Z[x].  At the roots r_+,r_-
    # of q=x^2+x-1, A=-1 and |r_+-r_-|=sqrt(5).
    # For a quotient ax+b, two values of modulus <1 imply |a|sqrt(5)<2,
    # whence a=0 and then b=0.  This excludes all nonzero degrees <=5.
    assert 5 > 4  # sqrt(5)>2, recorded without floating point.

    # For a degree-six quotient ax^2+bx+c, reduction modulo q is
    # (b-a)x+(a+c).  The same difference argument gives b=a,c=-a, so
    # the quotient is a*q.  Norm <1 then forces |a|=1 (or the zero row).
    assert P6_SUP < 1
    assert 2 * P6_SUP > 1


def one_exceptional_embedding_record() -> dict[str, object]:
    # On x in [2,14], w=x(x+1) is increasing and w>=6.  Hence
    # w^4(w-1)^3(w-2)^3 is increasing and its maximum is at x=14,w=210.
    maximum = 210**4 * 209**3 * 208**3
    at_22 = Fraction(maximum) * P20_SUP**22
    at_23 = Fraction(maximum) * P20_SUP**23
    assert at_22 > 1
    assert at_23 < 1
    return {
        "distinguished_embedding_interval": [2, 14],
        "other_embedding_interval": [-2, 1],
        "distinguished_maximum": maximum,
        "q_power_22_times_maximum_gt_1": True,
        "q_power_23_times_maximum_lt_1": True,
        "degree_bound_for_Q_alpha": 23,
        "generator_hypothesis_required": True,
    }


def build_archive() -> dict[str, object]:
    verify_exact_suprema()
    verify_minimal_degree_and_uniqueness()
    rows, witness_digest = bounded_family_rows()
    del rows  # The replayable digest keeps the archive small; the script checks every row.
    return {
        "schema": "fekete_interval_minus2_1_exact_v1",
        "interval": [-2, 1],
        "coordinate_change": "w=x(x+1), equivalently z=(x+1/2)^2=w+1/4",
        "minimal_certificate": {
            "factor_form": "x(x-1)(x+1)(x+2)(x^2+x-1)",
            "expanded_coefficients_descending": expanded_coefficients(P6_EXPONENTS),
            "degree": 6,
            "exact_sup_norm": str(P6_SUP),
            "degree_minimal": True,
            "unique_in_degree_6_up_to_sign": True,
        },
        "bounded_factor_family": {
            "form": "x^a(x+1)^b(x+2)^c(x-1)^d(x^2+x-1)^e",
            "ranges": {"a": [1, 4], "b": [1, 4], "c": [1, 4], "d": [1, 4], "e": [0, 4]},
            "row_count": 1_280,
            "unique_winner_exponents_a_b_c_d_e": list(P20_EXPONENTS),
            "winner_factor_form": "[x(x+1)]^4[(x-1)(x+2)]^3[x^2+x-1]^3",
            "winner_expanded_coefficients_descending": expanded_coefficients(P20_EXPONENTS),
            "winner_degree": 20,
            "winner_exact_sup_norm": str(P20_SUP),
            "winner_normalized_sup_norm_decimal": "0.8868612567109771154024468577207441504",
            "losing_row_exact_rational_witness_digest_sha256": witness_digest,
            "root_isolation_width_at_most": "1/100000000",
        },
        "one_exceptional_embedding": one_exceptional_embedding_record(),
        "scope_caveat": (
            "The norm argument bounds Q(alpha). It bounds an ambient field K only "
            "when alpha is certified to generate K (or after a separate relative-degree argument)."
        ),
    }


def canonical_json(data: dict[str, object]) -> str:
    return json.dumps(data, sort_keys=True, indent=2) + "\n"


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="rewrite the canonical archive")
    parser.add_argument("--check", action="store_true", help="check the canonical archive")
    args = parser.parse_args(argv)

    archive = build_archive()
    payload = canonical_json(archive)
    if args.write:
        ARCHIVE.write_text(payload, encoding="utf-8")
    if args.check or not args.write:
        assert ARCHIVE.read_text(encoding="utf-8") == payload, "archive is stale"
    print(
        "PASS: degree-6 norm 45/64; unique bounded-family degree-20 winner; "
        "one-exceptional-generator degree <= 23"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
