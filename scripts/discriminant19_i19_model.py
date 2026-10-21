"""Exact checks for Shioda's elliptic K3 surface with an I19 fiber."""

from __future__ import annotations

import json

import sympy as sp


T19 = [[2, 1], [1, 10]]


def _order_at_zero(polynomial: sp.Expr, variable: sp.Symbol) -> int:
    poly = sp.Poly(sp.cancel(polynomial), variable)
    return min(monomial[0] for monomial, coefficient in poly.terms() if coefficient)


def analyze() -> dict[str, object]:
    t, s = sp.symbols("t s")
    f = (
        t**8
        + 6 * t**7
        + 21 * t**6
        + 50 * t**5
        + 86 * t**4
        + 114 * t**3
        + 109 * t**2
        + 74 * t
        + 28
    )
    g = sp.Rational(1, 2) * (
        2 * t**12
        + 18 * t**11
        + 90 * t**10
        + 312 * t**9
        + 816 * t**8
        + 1692 * t**7
        + 2832 * t**6
        + 3864 * t**5
        + 4272 * t**4
        + 3746 * t**3
        + 2517 * t**2
        + 1167 * t
        + 299
    )
    quintic = 4 * t**5 + 15 * t**4 + 38 * t**3 + 61 * t**2 + 62 * t + 59
    h = -sp.Rational(27, 4) * quintic

    a4 = -3 * f
    a6 = -2 * g
    c4 = sp.expand(-48 * a4)
    discriminant = sp.factor(-16 * (4 * a4**3 + 27 * a6**2))
    expected_discriminant = sp.factor(1728 * h)

    local_c4_at_infinity = sp.cancel(s**8 * c4.subs(t, 1 / s))
    local_discriminant_at_infinity = sp.cancel(s**24 * discriminant.subs(t, 1 / s))
    infinity_c4_order = _order_at_zero(local_c4_at_infinity, s)
    infinity_discriminant_order = _order_at_zero(local_discriminant_at_infinity, s)

    h_squarefree = sp.degree(sp.gcd(sp.Poly(quintic, t), sp.Poly(sp.diff(quintic, t), t))) == 0
    finite_c4_coprime = sp.degree(sp.gcd(sp.Poly(f, t), sp.Poly(quintic, t))) == 0

    return {
        "surface_id": "shioda_i19",
        "base_field": "Q",
        "source": {
            "author": "Tetsuji Shioda",
            "title": "The elliptic K3 surfaces with a maximal singular fibre",
            "theorem": "Theorem 1.1 and Section 3",
            "doi": "10.1016/j.crma.2003.07.007",
            "url": "https://www.numdam.org/articles/10.1016/j.crma.2003.07.007/",
        },
        "published_uniqueness_in_characteristic_zero": True,
        "published_transcendental_lattice": T19,
        "published_neron_severi_lattice": "U + A18(-1)",
        "weierstrass_shape": "y^2 = x^3 - 3*f(t)*x - 2*g(t)",
        "f_coefficients_descending": [1, 6, 21, 50, 86, 114, 109, 74, 28],
        "two_g_coefficients_descending": [
            2,
            18,
            90,
            312,
            816,
            1692,
            2832,
            3864,
            4272,
            3746,
            2517,
            1167,
            299,
        ],
        "h_identity": "f(t)^3-g(t)^2=-(27/4)*(4*t^5+15*t^4+38*t^3+61*t^2+62*t+59)",
        "davenport_stothers_identity_verified": sp.expand(f**3 - g**2 - h) == 0,
        "standard_weierstrass_discriminant": str(discriminant),
        "discriminant_identity_verified": sp.expand(discriminant - expected_discriminant) == 0,
        "finite_discriminant_quintic_squarefree": bool(h_squarefree),
        "finite_c4_coprime_to_discriminant": bool(finite_c4_coprime),
        "infinity_c4_order": infinity_c4_order,
        "infinity_discriminant_order": infinity_discriminant_order,
        "infinity_fiber": "I19",
        "other_singular_fibers": [{"kodaira_type": "I1", "multiplicity": 5}],
        "fiber_euler_sum": 19 + 5,
        "published_picard_number": 20,
        "shioda_tate_mordell_weil_rank": 0,
        "geometric_status": "explicit_singular_k3_model_verified",
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
