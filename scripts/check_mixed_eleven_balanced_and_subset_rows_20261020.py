#!/usr/bin/env python3
"""Exact rational proof checks; run inside the repository bounded wrapper."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT / 'results/mixed_eleven_balanced_and_subset_rows_20261020.json'
SOURCES = [
    'results/right_angled_one_cusp_ten_facet_exclusion_research_20261020.md',
    'results/right_angled_mixed_eleven_three_face_corner_balance_20261020.md',
    'results/right_angled_mixed_eleven_structured_small_and_full_subsets_20261020.md',
    'results/right_angled_mixed_nu10_edge_transition_research_20261018.md',
    'results/right_angled_mixed_nu12_carrier_error_sharpening_research_20261019.md',
    'scripts/check_mixed_eleven_balanced_and_subset_rows_20261020.py',
]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def sqrt_up(q):
    assert q >= 0
    scale = 10**6
    n = isqrt(q.numerator * scale**2 // q.denominator)
    if n*n*q.denominator < q.numerator*scale**2:
        n += 1
    assert n*n*q.denominator >= q.numerator*scale**2
    return F(n, scale)


def small_subsets(k):
    start, stop = {7: (1350, 2850), 8: (1700, 3100)}[k]
    alpha2 = F(3*(k-3), k*(k-1))
    alpha3 = F(3*(k-3)*(k-6), k*(k-1)*(k-2))
    alpha4 = F(3*(k*k-17*k+54), k*(k-1)*(k-2))
    b4 = F(18*(k-4), k*(k-1)*(k-2))
    assert alpha3 > 0 and alpha4 < 0 and b4+alpha4/k > 0
    a, t, threshold = F(11, 20), F(37, 100), F(18, 11)
    goal = t + (F(756, 25) + F(27, 100)*comb(k, 3))/comb(k, 3)
    least, worst, checked, skipped = None, None, 0, 0
    for j in range(start, stop):
        lo, hi = F(j, 1100), F(j+1, 1100)
        vertex = min(hi, max(lo, F(3*k, 11)))
        v = 1-vertex**2/k-(3-vertex)**2/(11-k)
        if v < 0:
            skipped += 1
            continue
        assert v <= F(2, 11)
        mu, mu_hi = 2-3*hi/k, 2-3*lo/k
        assert mu > a
        K = (4*a+mu)/(mu**5*a*a)
        J = (1/mu**2+3*alpha2*v/mu**4
             +4*alpha3*v*sqrt_up(v)/mu**5
             +K*(b4+alpha4/k)*v*v)
        correction = F(27, 400)
        if mu_hi < threshold:
            correction = min(correction, F(113, 500)*alpha2*v/(2*(threshold-mu_hi)))
        if mu_hi+sqrt_up(F(3*(k-3), k)*v) <= threshold:
            correction = F(0)
        gap = goal-J-correction
        assert gap > 0, (k, j, str(gap))
        if least is None or gap < least:
            least, worst = gap, j
        checked += 1
    return dict(intervals_checked=checked, intervals_empty=skipped,
                worst_left_numerator=worst, denominator=1100,
                least_gap=str(least), strict=True)


def full_subset():
    mu, v, variance = F(13, 11), F(2, 11), F(24, 605)
    a, t = F(11, 20), F(37, 100)
    assert F(139, 80)**2 > 3
    assert (4*a+mu)/(mu**5*a*a) < 5
    assert F(18**2, 121**2*5) < F(1, 225)
    C = 16/(33*mu**5)
    assert C > F(2, 15)
    assert 37**3 < 226**2 and F(100, 37) > F(18, 11)**2
    J = 1/mu**2+3*variance/mu**4+5*F(15, 121)*v*v+C/15-F(1, 225)
    correction = F(113, 500)*variance/(2*F(5, 11))
    target = (F(756, 25)+F(27, 100)*165)/165
    gap = target-(J-t+correction)
    assert gap > 0
    return dict(positive_part_mean_upper=str(J-t+correction),
                target=str(target), strict_gap=str(gap))


def balanced_rows():
    rows = []
    for k in range(12):
        m = 11-k
        U, V, W = comb(m, 3), comb(k, 2)*m, comb(k, 3)
        value = 40*U+28*V-16*(165-W)
        assert value >= 0
        rows.append(value)
    # Variables in every coefficient vector are B,f,E,c,a2,a3.
    E = [1, 0, 0, 11520, 9, -10]
    R = [-2, 0, 0, 0, -9, 8]
    euler = [0, F(-165, 4), 0, -5760, F(9, 2), 0]
    def eliminate_a3(v):
        return [v[i]+v[5]*euler[i] for i in range(5)]
    e, r = eliminate_a3(E), eliminate_a3(R)
    # Eliminate a2 using e=E, independently in each remaining expression.
    assert e[4] == -36
    r_E = r[4]/e[4]
    r_rest = [r[i]-r_E*e[i] for i in range(4)]
    assert r_E == F(-3, 4)
    assert r_rest == [F(-5, 4), F(-165, 8), 0, F(5760)]
    # 680E+124R-2640f+32B>=0, then the strict positive cusp substitution.
    row = [32+124*r_rest[0], -2640+124*r_rest[1], 680+124*r_E, 124*r_rest[3]]
    assert row == [-123, F(-10395, 2), 587, 714240]
    factor = row[3]/360
    final = [row[0]+factor, row[1]-factor*F(75, 2), row[2]-factor]
    assert final == [1861, F(-159195, 2), -1397]
    assert 159195*1379-117825*1861 == 257580
    assert 2794*1379-2070*1861 == 656
    assert F(2794, 1861)-F(3, 2) == F(5, 3722)
    return dict(integer_corner_rows=rows, strict_lower='3722B>159195f+2794E',
                proposed_upper_excess_gap='5/3722')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    assert not (args.write and args.verify)
    result = dict(schema=1, status='PASS_EXACT_RATIONAL_ROWS',
                  pins={p: digest(ROOT/p) for p in SOURCES},
                  balanced=balanced_rows(),
                  subsets={str(k): small_subsets(k) for k in (7, 8)},
                  full_eleven=full_subset(),
                  scope='Partial structured-subset estimates; no dimension-ten theorem')
    if args.write:
        with CERT.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    else:
        assert json.loads(CERT.read_text()) == result, 'Certificate mismatch'
    print(json.dumps({k: v for k, v in result.items() if k != 'pins'}))


if __name__ == '__main__':
    main()
