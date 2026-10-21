#!/usr/bin/env python3
"""Exact interval proof of the nine-direction structured subset estimate."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/mixed_eleven_nine_direction_subset_20261020.json'
PINS = [
    'results/right_angled_mixed_eleven_nine_direction_subset_20261020.md',
    'results/right_angled_mixed_eleven_structured_small_and_full_subsets_20261020.md',
    'scripts/check_mixed_eleven_nine_direction_subset_20261020.py',
]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def sqrt_up(q):
    assert q >= 0
    n = isqrt(q.numerator*10**12//q.denominator)
    if n*n*q.denominator < q.numerator*10**12:
        n += 1
    assert n*n*q.denominator >= q.numerator*10**12
    return F(n, 10**6)


def check():
    assert F(1, 75) < F(1, 8)**2
    assert F(288, 49) < 6
    assert F(50, 1323) < F(1, 25)
    assert F(-1, 6)+F(1, 50)+F(7, 900) == F(-5, 36)
    assert F(301, 300)+sqrt_up(F(1, 75)) < F(18, 11)
    least, worst, checked = None, None, 0
    for j in range(100, 10910):
        lo, hi = F(j, 10000), min(F(j+1, 10000), F(12, 11))
        if lo >= hi:
            continue
        vertex = min(hi, max(lo, F(6, 11)))
        v = F(2, 3)*vertex-F(11, 18)*vertex**2
        assert 0 < v <= F(2, 11)
        mu, mu_hi = 1+lo/3, 1+hi/3
        a = max(F(11, 20), mu-sqrt_up(2*v))
        assert 0 < a <= mu
        K = (4*a+mu)/(mu**5*a*a)
        A, B = 3/(7*mu**5), 3*K/(28*v)
        limit = F(5, 6)*v*sqrt_up(v)
        z = min(limit, A/(2*B))
        assert 0 <= z <= limit and A-2*B*z >= 0
        J = 1/mu**2+F(3, 4)*v/mu**4+K*v*v/6+A*z-B*z*z
        correction = F(27, 400)
        if mu_hi < F(18, 11):
            correction = min(correction, F(113, 500)*v/(8*(F(18, 11)-mu_hi)))
        if mu_hi+sqrt_up(2*v) <= F(18, 11):
            correction = F(0)
        gap = 1-J-correction
        assert gap > 0, (j, str(gap))
        if least is None or gap < least:
            least, worst = gap, j
        checked += 1
    return dict(schema=1, status='PASS_EXACT_NINE_DIRECTION_CONTRACT',
                pins={p: digest(ROOT/p) for p in PINS},
                equality='Nine entries 1/3 and two entries 0',
                analytic_neighborhood='0<delta<=1/100: J<=1-5delta/36',
                intervals_checked=checked, worst_left_numerator=worst,
                denominator=10000, least_gap=str(least),
                scope='Analytic subset lemma; no dimension theorem asserted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    assert not (args.write and args.verify)
    result = check()
    if args.write:
        with CERT.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    else:
        assert json.loads(CERT.read_text()) == result, 'Certificate mismatch'
    print(json.dumps({k: v for k, v in result.items() if k != 'pins'}))


if __name__ == '__main__':
    main()
