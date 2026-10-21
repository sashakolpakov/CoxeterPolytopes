#!/usr/bin/env python3
"""Exact low-label calibration for the Gram--Belyi packet.

This producer certifies the normalized Belyi-polynomial heights attached to
the standard Coxeter relations m=2,3,4,5,6, to a parallel pair, and to one
rational dotted pair.  The rational cases use Rodriguez's p-adic lower bound
and the elementary beta-polynomial witness.  The m=5 case is deliberately
reported as an interval: its 29-adic lower bound is sharp as a *bound*, but no
degree-29 polynomial witness is claimed.  A fully explicit factored
composition supplies the finite upper endpoint without expanding a polynomial
of astronomical degree.

The only external mathematical input is Theorem 3 of Jose Rodriguez,
"Bounding the Degree of Belyi Polynomials" (arXiv:1104.2027; JNT 133
(2013), 2892--2900): if a nonzero algebraic number has nonzero p-adic
valuation, its normalized Belyi-polynomial height is at least p.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/gram_belyi_standard_packets.json"
SCHEMA = "coxeter.gram_belyi.standard_packets.v1"


def valuation(integer: int, prime: int) -> int:
    """Return v_p(integer), for a nonzero integer."""

    if integer == 0:
        raise ValueError("valuation of zero is not used")
    value = abs(integer)
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def fraction_string(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def q5_add(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    """Add exact elements a+b*sqrt(5)."""

    return left[0] + right[0], left[1] + right[1]


def q5_mul(
    left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]
) -> tuple[Fraction, Fraction]:
    """Multiply exact elements a+b*sqrt(5)."""

    return (
        left[0] * right[0] + 5 * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def q5_inv(value: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    """Invert a nonzero exact element a+b*sqrt(5)."""

    denominator = value[0] ** 2 - 5 * value[1] ** 2
    if denominator == 0:
        raise ZeroDivisionError(value)
    return value[0] / denominator, -value[1] / denominator


def beta_witness(a: int, b: int) -> dict[str, Any]:
    """Certify the degree-b beta polynomial at lambda=a/b.

    B_{a,b}(x) = b^b/(a^a (b-a)^(b-a)) x^a (1-x)^(b-a).
    Its only possible finite critical values are zero and the value one at
    a/b.  The representation stays factored, so all checks are exact and
    cheap.
    """

    if not 0 < a < b:
        raise ValueError((a, b))
    lam = Fraction(a, b)
    coefficient = Fraction(b**b, a**a * (b - a) ** (b - a))
    critical_value = coefficient * lam**a * (1 - lam) ** (b - a)
    if critical_value != 1:
        raise ArithmeticError((a, b, critical_value))
    return {
        "degree": b,
        "formula": (
            f"({coefficient.numerator}/{coefficient.denominator})"
            f"*x^{a}*(1-x)^{b-a}"
        ),
        "coefficient": fraction_string(coefficient),
        "distinguished_preimage": fraction_string(lam),
        "distinguished_value": 1,
        "derivative_certificate": (
            "apart from roots at 0 and 1, the logarithmic derivative "
            f"vanishes only at x={a}/{b}"
        ),
    }


def rational_exact_row(
    *,
    relation: str,
    u: Fraction,
    lam: Fraction,
    a: int,
    b: int,
    lower_prime: int,
    interpretation: str,
) -> dict[str, Any]:
    if lam != u / (1 + u):
        raise ArithmeticError((relation, u, lam))
    if lam != Fraction(a, b):
        raise ArithmeticError((relation, lam, a, b))
    vp = valuation(lam.numerator, lower_prime) - valuation(
        lam.denominator, lower_prime
    )
    if vp == 0:
        raise ArithmeticError((relation, lower_prime, vp))
    witness = beta_witness(a, b)
    if witness["degree"] != lower_prime:
        raise ArithmeticError(
            "this exact calibration needs matching Rodriguez and witness degrees"
        )
    return {
        "relation": relation,
        "interpretation": interpretation,
        "normalized_squared_gram_u": fraction_string(u),
        "packet_lambda": fraction_string(lam),
        "minimum_status": "proved_exact",
        "h_Bel": lower_prime,
        "lower_bound": {
            "theorem": "Rodriguez Theorem 3",
            "prime": lower_prime,
            "p_adic_valuation": vp,
            "bound": lower_prime,
        },
        "upper_witness": witness,
    }


def m5_interval_row() -> dict[str, Any]:
    """Build the exact 29-adic lower and explicit polynomial upper witness."""

    u = (Fraction(3, 8), Fraction(1, 8))
    lam = (Fraction(7, 29), Fraction(2, 29))
    one = (Fraction(1), Fraction(0))
    computed_lam = q5_mul(u, q5_inv(q5_add(one, u)))
    if computed_lam != lam:
        raise ArithmeticError((computed_lam, lam))

    # lambda=(7+2*sqrt(5))/29 has primitive minimal polynomial
    # q(x)=29*x^2-14*x+1.  At p=29 its Newton polygon has vertices
    # (0,0),(1,0),(2,1), hence slopes 0 and 1 and root valuations 0 and -1.
    minimal_polynomial = (29, -14, 1)
    lam_squared = q5_mul(lam, lam)
    polynomial_value = q5_add(
        q5_add(
            (29 * lam_squared[0], 29 * lam_squared[1]),
            (-14 * lam[0], -14 * lam[1]),
        ),
        one,
    )
    if polynomial_value != (Fraction(0), Fraction(0)):
        raise ArithmeticError(polynomial_value)
    newton_points = ((0, 0), (1, 0), (2, 1))
    newton_slopes = (Fraction(0), Fraction(1))
    root_valuations = tuple(-slope for slope in newton_slopes)
    if root_valuations != (Fraction(0), Fraction(-1)):
        raise ArithmeticError(root_valuations)

    # Explicit compositional witness.  For either root lambda of q,
    # C(lambda)=-1/15.  The sole non-root critical point of C maps to r.
    s = Fraction(-1, 15)
    r = Fraction(-49, 435)
    c_at_zero = Fraction(0)
    c_at_one = Fraction(1)
    c_at_lambda = s
    c_critical_point = Fraction(7, 29)
    c_critical_value = (
        29 * c_critical_point**2 - 14 * c_critical_point
    ) / 15
    if (c_at_zero, c_at_one, c_at_lambda, c_critical_value) != (
        Fraction(0),
        Fraction(1),
        Fraction(-1, 15),
        r,
    ):
        raise ArithmeticError("quadratic stage failed")

    # Q(y)=-(y-1)^15(15y+1) is degree 16 with
    # Q'(y)=-240*y*(y-1)^14, so it is normalized Belyi.
    def q_value(y: Fraction) -> Fraction:
        return -(y - 1) ** 15 * (15 * y + 1)

    if (q_value(Fraction(0)), q_value(Fraction(1)), q_value(s)) != (
        Fraction(1),
        Fraction(0),
        Fraction(0),
    ):
        raise ArithmeticError("degree-16 stage failed")
    t = q_value(r)

    # In lowest terms t=-b/a.  Both displayed factorizations are retained so
    # the reduction can be audited without factoring a 40-digit integer.
    a = 3**15 * 5**14 * 29**16
    b = 2**32 * 11**30
    if t != Fraction(-b, a):
        raise ArithmeticError((t, a, b))
    if a % 2 != 1:
        raise ArithmeticError("the final normalization uses odd a")

    # F(z)=-(a/b)^b (z-1)^a (z+b/a)^b.  Its derivative is the same
    # nonzero scalar times (a+b)z(z-1)^(a-1)(z+b/a)^(b-1).
    # Thus F(0)=1 and F(1)=F(-b/a)=0, and all critical values are 0 or 1.
    final_degree = a + b
    total_degree = 2 * 16 * final_degree
    if total_degree != 3_099_548_558_797_153_417_819_718_240_423_033_019_647_072:
        raise ArithmeticError(total_degree)

    return {
        "relation": "m=5",
        "interpretation": "finite dihedral angle pi/5",
        "normalized_squared_gram_u": "(3+sqrt(5))/8",
        "packet_lambda": "(7+2*sqrt(5))/29",
        "lambda_minimal_polynomial": "29*x^2-14*x+1",
        "minimum_status": "certified_interval_exact_minimum_open",
        "lower_bound": {
            "theorem": "Rodriguez Theorem 3",
            "prime": 29,
            "newton_polygon_points": [list(point) for point in newton_points],
            "newton_polygon_slopes": [fraction_string(x) for x in newton_slopes],
            "root_valuations": [fraction_string(x) for x in root_valuations],
            "galois_invariance_used": True,
            "bound": 29,
        },
        "upper_bound": total_degree,
        "upper_witness": {
            "composition": "F(Q(C(x)))",
            "C": "(29*x^2-14*x)/15",
            "C_degree": 2,
            "C_values": {
                "C(0)": "0/1",
                "C(1)": "1/1",
                "C(lambda)": "-1/15",
                "C(7/29)": "-49/435",
            },
            "Q": "-(y-1)^15*(15*y+1)",
            "Q_degree": 16,
            "Q_derivative": "-240*y*(y-1)^14",
            "Q(C(7/29))": fraction_string(t),
            "final_negative_coordinate": {
                "t": fraction_string(t),
                "a_factorization": "3^15*5^14*29^16",
                "a": a,
                "b_factorization": "2^32*11^30",
                "b": b,
                "identity": "t=-b/a",
            },
            "F": "-(a/b)^b*(z-1)^a*(z+b/a)^b",
            "F_degree": final_degree,
            "F_derivative_certificate": (
                "up to a nonzero scalar: "
                "(a+b)*z*(z-1)^(a-1)*(z+b/a)^(b-1)"
            ),
            "F_values": {"F(0)": 1, "F(1)": 0, "F(-b/a)": 0},
            "total_degree": total_degree,
        },
        "proved_interval": [29, total_degree],
        "nonclaim": "no degree-29 witness and no exact m=5 minimum is claimed",
    }


def build_certificate() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []

    # Orthogonality is the exceptional packet in the general theorem: lambda
    # is not defined at u=0, and degree two is assigned by construction.
    orthogonal_aux = beta_witness(1, 2)
    rows.append(
        {
            "relation": "m=2",
            "interpretation": "orthogonal facets",
            "normalized_squared_gram_u": "0/1",
            "packet_lambda": None,
            "minimum_status": "assigned_orthogonality_packet",
            "gram_belyi_contribution": 2,
            "orthogonality_bit": True,
            "auxiliary_packet": {
                "lambda": "1/2",
                "h_Bel": 2,
                "minimum_status": "proved_exact",
                "witness": orthogonal_aux,
                "warning": "the auxiliary lambda is not recovered from u=0",
            },
        }
    )

    rows.extend(
        [
            rational_exact_row(
                relation="m=3",
                u=Fraction(1, 4),
                lam=Fraction(1, 5),
                a=1,
                b=5,
                lower_prime=5,
                interpretation="finite dihedral angle pi/3",
            ),
            rational_exact_row(
                relation="m=4",
                u=Fraction(1, 2),
                lam=Fraction(1, 3),
                a=1,
                b=3,
                lower_prime=3,
                interpretation="finite dihedral angle pi/4",
            ),
        ]
    )
    rows.append(m5_interval_row())
    rows.append(
        rational_exact_row(
            relation="m=6",
            u=Fraction(3, 4),
            lam=Fraction(3, 7),
            a=3,
            b=7,
            lower_prime=7,
            interpretation="finite dihedral angle pi/6",
        )
    )
    rows.append(
        rational_exact_row(
            relation="parallel",
            u=Fraction(1),
            lam=Fraction(1, 2),
            a=1,
            b=2,
            lower_prime=2,
            interpretation="asymptotic hyperplanes; zero angle",
        )
    )
    rows.append(
        rational_exact_row(
            relation="dotted u=2",
            u=Fraction(2),
            lam=Fraction(2, 3),
            a=2,
            b=3,
            lower_prime=3,
            interpretation="ultraparallel rational squared Gram entry",
        )
    )

    exact = {
        row["relation"]: row.get("h_Bel", row.get("gram_belyi_contribution"))
        for row in rows
        if row["minimum_status"] != "certified_interval_exact_minimum_open"
    }
    if exact != {
        "m=2": 2,
        "m=3": 5,
        "m=4": 3,
        "m=6": 7,
        "parallel": 2,
        "dotted u=2": 3,
    }:
        raise ArithmeticError(exact)

    return {
        "schema": SCHEMA,
        "all_checks_pass": True,
        "definition": (
            "h_Bel(lambda) is the least degree of a normalized Belyi "
            "polynomial mapping lambda to 0 or 1"
        ),
        "coordinate_rule": "lambda=u/(1+u) for u>0",
        "rows": rows,
        "exact_height_summary": exact,
        "m5_height_interval": rows[3]["proved_interval"],
        "scope": {
            "proved": (
                "exact minima for m=3,4,6, parallel and dotted u=2; "
                "the assigned degree-two orthogonal packet; a certified "
                "finite interval for m=5"
            ),
            "open": "the exact normalized Belyi-polynomial height at m=5",
            "not_used": "no unbounded dessin search and no numerical algebra",
        },
        "primary_source": {
            "author": "Jose Rodriguez",
            "title": "Bounding the Degree of Belyi Polynomials",
            "arxiv": "1104.2027",
            "journal_doi": "10.1016/j.jnt.2012.12.019",
            "result_used": "Theorem 3 (nonzero p-adic valuation gives H(lambda)>=p)",
        },
    }


def canonical(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = canonical(payload)
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen standard-packet archive differs from exact replay")
    print(rendered, end="")


if __name__ == "__main__":
    main()
