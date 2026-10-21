#!/usr/bin/env python3
"""Finite8 budget with whole D11 and two-position whole A17 half-coupons."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_finite8_a17_d11_coupon_budget_20261017'
BASE = 'fresh_finite8_leaf9_half_coupon_budget_20261017'
LOCAL = 'fresh_finite8_ideal5_local_replacement_20261017'


def run():
    prior = json.loads((ROOT / f'results/{BASE}.json').read_text())
    assert prior['dimension_maximum'] == 80 and prior['finite_vertex_budget'] == '8n'
    for path, h in prior['sources'].items():
        assert sha256((ROOT / path).read_bytes()).hexdigest() == h, path
    small = []
    for r in range(1, 6):
        slack = 8*r-F(r*(r-1), 2)
        assert slack >= 6*r
        small.append(dict(rank=r, slack=str(slack), minimum_slack=str(6*r)))
    large = []
    for kind, r, marks in [('D', r, 3) for r in range(11, 17)] + [('A', 13, 1), ('A', 15, 1)]:
        slack = 8*r-F(r*(r-1), 2)
        assert r >= 11 and slack >= 8 and marks <= 3
        large.append(dict(type=kind+str(r), rank=r, eligible_marks=marks, slack=str(slack)))
    assert F(17*16, 2) == 8*17  # A17 has no base slack; charge it separately.
    endpoints = []
    positive_coefficient_cases = []
    for n in range(1, 81):
        for t in range(n//11+1):
            for u in range((n-11*t)//17+1):
                maximum_a = n-11*t-17*u
                coefficient = F(3*t+2*u-12, 2)
                for a in sorted({0, maximum_a}):
                    excess = a*coefficient-8*t
                    assert excess <= 0, (n, a, t, u, excess)
                    endpoints.append(dict(dimension=n, small_rank=a, old_large_count=t,
                        A17_count=u, excess=str(excess)))
                if n == 80 and coefficient > 0:
                    positive_coefficient_cases.append(dict(old_large_count=t, A17_count=u,
                        maximum_small_rank=maximum_a, coefficient=str(coefficient),
                        maximum_excess=str(maximum_a*coefficient-8*t)))
    assert [(r['old_large_count'], r['A17_count']) for r in positive_coefficient_cases] == [
        (3, 2), (4, 1), (4, 2), (5, 0), (5, 1), (6, 0), (7, 0)]
    paths = [Path(__file__)] + [ROOT / f'{folder}/{stem}.{ext}'
        for stem in (BASE, LOCAL) for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    return dict(status='PASS_FINITE8_A17_D11_HALF_COUPON_BUDGET',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        dimension_maximum=80, finite_vertex_budget='8n',
        small_components=small, old_large_components=large,
        new_A17_component=dict(rank=17, eligible_path_positions=[5, 11], slack='0'),
        surplus_bound='a*(3t/2+u-6)-8t, where a+11t+17u<=n<=80',
        dimension_endpoints=endpoints, positive_coefficient_cases_at_80=positive_coefficient_cases,
        rule='Retain unit on same-component distance<=8, distance9 with an actual whole-component leaf mark, or whole component rank<=17. Add half to each pair of marks in separate whole components: one has rank<=5; the other is D11..D16 marked at its fork center or short leaf, A13/A15 marked at its middle, or A17 marked at path position5 or11. Count a pair once.',
        scope='Universal finite-vertex counting theorem. Includes all old coupons and the sealed cycle18 coupons. D11 and small<=5 are safe counting enlargements; their geometric application requires a separate local proof. Full finite8 face coverage remains a separate obligation. No new ambient bound.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    d = run()
    path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert d == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], len(d['dimension_endpoints']), 'endpoint checks', flush=True)
