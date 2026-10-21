#!/usr/bin/env python3
"""Exact checks for the Vinberg integral-point/arrangement-height report."""

from __future__ import annotations

import json

import sympy as sp


def lorentz_pair(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand(-left[0] * right[0] + left[1] * right[1] + left[2] * right[2])


def vertical_root(c: sp.Expr) -> sp.Matrix:
    return sp.Matrix([c, 1, c])


def circle_root(j: sp.Expr) -> sp.Matrix:
    return sp.Matrix(
        [-(j**2 - j + 1), 1 - 2 * j, -j**2 + j + 1]
    )


def euler_phi(n: int) -> int:
    result = n
    divisor = 2
    remaining = n
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            while remaining % divisor == 0:
                remaining //= divisor
            result -= result // divisor
        divisor += 1
    if remaining > 1:
        result -= result // remaining
    return result


def main() -> None:
    # The Lambert contracting polynomial.
    x, y, z = sp.symbols("x y z", real=True)
    p = sp.expand((x + 1) * (x + 2) * (x + 3) * (x + 4))
    f = sp.expand(p * (p + 1))
    centered = sp.expand(p.subs(x, y - sp.Rational(5, 2)))
    expected_centered = y**4 - sp.Rational(5, 2) * y**2 + sp.Rational(9, 16)
    assert sp.expand(centered - expected_centered) == 0

    p_as_z = z**2 - sp.Rational(5, 2) * z + sp.Rational(9, 16)
    assert p_as_z.subs(z, 0) == sp.Rational(9, 16)
    assert p_as_z.subs(z, sp.Rational(5, 4)) == -1
    assert p_as_z.subs(z, sp.Rational(9, 4)) == 0

    # On p in [-1,9/16], p(p+1) ranges from -1/4 to 225/256.
    u = sp.symbols("u", real=True)
    f_as_u = u * (u + 1)
    assert f_as_u.subs(u, -sp.Rational(1, 2)) == -sp.Rational(1, 4)
    assert f_as_u.subs(u, sp.Rational(9, 16)) == sp.Rational(225, 256)
    contracting_sup = sp.Rational(225, 256)
    exceptional_sup = f.subs(x, 2)
    assert exceptional_sup == 129960

    # d=92 is the last degree not ruled out by the norm inequality.
    assert 129960 * 225**91 > 256**91
    assert 129960 * 225**92 < 256**92
    degree_bound = 92

    # Lehmer's reciprocal polynomial and its degree-five trace polynomial.
    w = sp.symbols("w")
    lehmer = w**10 + w**9 - w**7 - w**6 - w**5 - w**4 - w**3 + w + 1
    trace_poly = x**5 + x**4 - 5 * x**3 - 5 * x**2 + 4 * x + 3
    assert sp.factor(w**5 * trace_poly.subs(x, w + 1 / w) - lehmer) == 0
    trace_as_poly = sp.Poly(trace_poly, x, domain=sp.QQ)
    assert trace_as_poly.count_roots(-2, 2) == 4
    assert trace_as_poly.count_roots(2, 3) == 1

    # Prime cyclotomic trace degrees grow while every conjugate stays in (0,4).
    primes = [5, 7, 11, 13, 17, 19, 23, 29]
    cyclotomic_trace_degrees = {prime: euler_phi(prime) // 2 for prime in primes}
    assert cyclotomic_trace_degrees == {
        5: 2,
        7: 3,
        11: 5,
        13: 6,
        17: 8,
        19: 9,
        23: 11,
        29: 14,
    }

    # The integral rank-two Salem Cartan packet.
    t = sp.symbols("t", positive=True)
    cartan = sp.Matrix([[2, -t], [-1, 2]])
    symmetrizer = sp.diag(1, t)
    salem_gram = symmetrizer * cartan
    assert salem_gram == salem_gram.T
    assert sp.simplify(salem_gram.det() - t * (4 - t)) == 0

    # Farey ideal polygons.
    i, j, k, n = sp.symbols("i j k n", integer=True)
    ri = circle_root(i)
    rj = circle_root(j)
    rk = circle_root(k)
    assert lorentz_pair(ri, ri) == 1
    assert sp.simplify(lorentz_pair(ri, rj) - (1 - 2 * (i - j) ** 2)) == 0
    assert sp.simplify(
        sp.Matrix.hstack(ri, rj, rk).det()
        - 4 * (i - j) * (i - k) * (j - k)
    ) == 0

    v0 = vertical_root(0)
    vn = vertical_root(n)
    wn = -vn
    assert lorentz_pair(v0, v0) == 1
    assert lorentz_pair(vn, vn) == 1
    assert lorentz_pair(v0, circle_root(1)) == -1
    assert sp.factor(lorentz_pair(circle_root(n), wn)) == -1
    assert sp.factor(lorentz_pair(wn, v0)) == -1
    assert sp.factor(sp.Matrix.hstack(v0, vn, circle_root(i)).det()) == -2 * n
    assert sp.simplify(
        sp.Matrix.hstack(v0, ri, rj).det()
        - 2 * (i - j) * (i + j - 1)
    ) == 0
    assert sp.simplify(
        sp.Matrix.hstack(wn, ri, rj).det()
        + 2 * (i - j) * (i + j - 2 * n - 1)
    ) == 0
    assert sp.simplify(sp.Matrix.hstack(v0, wn, ri).det() - 2 * n) == 0

    farey_cycle = sp.expand(4 * (2 * (n - 1) ** 2 - 1) ** 2)
    assert farey_cycle.subs(n, 3) == 196
    facet_count = n + 2
    log_class_coefficient = sp.simplify(facet_count / 2 - 1)
    assert log_class_coefficient == n / 2
    log_volume = sp.expand(log_class_coefficient**2)
    assert log_volume == n**2 / 4
    # [infinity, 0; 1, n] = (0-n)/(0-1) = n.
    farey_cross_ratio = sp.simplify((0 - n) / (0 - 1))
    assert farey_cross_ratio == n

    # At a prime divisor of N, v_0 and -v_N give the same projective root line.
    for prime in primes:
        reduced_v0 = tuple(int(value) % prime for value in v0)
        reduced_wn = tuple(int(value.subs(n, prime)) % prime for value in wn)
        assert reduced_v0 == tuple((-value) % prime for value in reduced_wn)

    certificate = {
        "lambert": {
            "p": str(p),
            "F": str(f),
            "contracting_sup": str(contracting_sup),
            "exceptional_sup": int(exceptional_sup),
            "degree_bound": degree_bound,
            "degree_92_survives_norm_test": True,
            "degree_93_fails_norm_test": True,
        },
        "lehmer": {
            "polynomial": str(lehmer),
            "trace_polynomial": str(trace_poly),
            "trace_roots_in_minus2_2": 4,
            "trace_roots_in_2_3": 1,
        },
        "cyclotomic_trace_degrees": cyclotomic_trace_degrees,
        "salem_cartan": {
            "cartan": str(cartan),
            "symmetric_gram": str(salem_gram),
            "determinant": str(sp.factor(salem_gram.det())),
        },
        "farey": {
            "circle_pairing": "1 - 2*(i-j)^2",
            "triple_determinant": "4*(i-j)*(i-k)*(j-k)",
            "dotted_cycle": str(farey_cycle),
            "facet_count": str(facet_count),
            "log_class_coefficient": str(log_class_coefficient),
            "log_volume": str(log_volume),
            "ordered_cusp_cross_ratio": str(farey_cross_ratio),
            "new_bad_prime_check": primes,
        },
    }
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
