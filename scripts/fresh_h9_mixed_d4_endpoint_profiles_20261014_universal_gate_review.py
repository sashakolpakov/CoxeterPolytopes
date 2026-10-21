#!/usr/bin/env python3
"""Proof/coverage audit and separate rational cosine bounds; no producer imports.

The full box replay is run separately with the primary checker. This does not
claim a second implementation of its quartic coefficients or box subdivision.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_h9_mixed_d4_endpoint_profiles_20261014'
PRIMARY = BASE + '_universal_gate'
TAG = PRIMARY + '_review'
S = 1 << 256


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def atan_bounds(x):
    # Independent identity uses1/2 and1/3, rather than Machin's1/5 and1/239.
    total = sum(((-1)**j * x**(2*j+1) / (2*j+1) for j in range(180)), F(0))
    following = x**361 / 361
    assert 0 < following < F(1, 1 << 350)
    return total, total + following


def outward(x, lower):
    t = x * S
    return F(t.numerator // t.denominator if lower else -((-t.numerator) // t.denominator), S)


def sqrt_bounds(n):
    x = isqrt(n * S * S)
    return F(x, S), F(x + (x*x != n*S*S), S)


def cos_series(x, degree_index):
    term = F(1)
    total = term
    for j in range(1, degree_index + 1):
        term *= -x*x / ((2*j-1)*(2*j))
        total += term
    return total


def decode_interval(raw, scale):
    return F(int(raw[0]), scale), F(int(raw[1]), scale)


def contained(small, large):
    assert large[0] <= small[0] <= small[1] <= large[1], (small, large)


def run():
    paths = [f'results/{PRIMARY}.json', f'results/{BASE}.json']
    d, profiles = [json.loads((ROOT / p).read_text()) for p in paths]
    seen = set()
    def check(p):
        if p in seen:
            return
        seen.add(p)
        obj = json.loads((ROOT / p).read_text())
        for q, h in obj.get('sources', {}).items():
            assert sha(q) == h, q
            if q.endswith('.json'):
                check(q)
    for p in paths:
        check(p)
    assert d['status'] == 'PASS_ALL_ENDPOINT_PROFILES_AND_CONTINUOUS_TAILS_MISS_GAP'
    assert len(d['edges']) == 52
    assert {str(r['edge_id']) for r in d['edges']} == set(profiles['profiles'])
    totals = Counter()
    for e in d['edges']:
        p = profiles['profiles'][str(e['edge_id'])]
        expected = [('finite', tuple(map(tuple, r['contacts']))) for r in p['finite_profiles']]
        expected += [('tail', r['wall']) for r in p['dihedral_tail_profiles']]
        actual = [(r['kind'], tuple(map(tuple, r['contacts'])) if r['kind'] == 'finite' else r['wall']) for r in e['profiles']]
        assert actual == expected
        outside = set(range(12)) - set(p['support'])
        for r in e['profiles']:
            counts = r['counts']
            leaves = sum(v for k, v in counts.items() if k != 'boxes')
            assert counts['boxes'] == 2 * leaves - 1
            assert isinstance(r['h_bound'], int) and r['h_bound'] >= 1
            assert r['maximum_depth'] < 180
            assert all(k in ('boxes', 'escape') or (k.startswith('wall:') and int(k[5:]) in outside) for k in counts)
            totals.update(counts)
            totals[r['kind'] + '_profiles'] += 1
    assert dict(totals) == d['counts']
    assert totals['finite_profiles'] == 1825 and totals['tail_profiles'] == 72
    # Alternative pi identity: atan(1/2)+atan(1/3)=pi/4.
    a, b = atan_bounds(F(1, 2)), atan_bounds(F(1, 3))
    pi = outward(4*(a[0]+b[0]), True), outward(4*(a[1]+b[1]), False)
    assert 3 < pi[0] < pi[1] < F(22, 7)
    primary_scale = 1 << d['dyadic_precision_bits']
    enc = d['cosine_enclosures']
    contained(pi, decode_interval(enc['pi'], primary_scale))
    roots = {n: sqrt_bounds(n) for n in (2, 3, 5)}
    exact_angles = [(F(0), F(0)), (F(-1), F(-1)),
        (-roots[2][1], -roots[2][0]),
        (-(1+roots[5][1])/2, -(1+roots[5][0])/2),
        (-roots[3][1], -roots[3][0])]
    assert len(enc['isolated']) == 14
    for interval, raw in zip(exact_angles, enc['isolated'][:5]):
        contained(interval, decode_interval(raw, primary_scale))
    checks = []
    for m in range(7, 17):
        # cos is decreasing on these positive arguments. The even-index
        # alternating partial sum is an upper bound; the next odd one lower.
        lo = cos_series(pi[1] / m, 65)
        hi = cos_series(pi[0] / m, 64)
        interval = (-2*hi, -2*lo)
        raw = enc['isolated'][m-2] if m < 16 else enc['tail_threshold']
        contained(interval, decode_interval(raw, primary_scale))
        checks.append(dict(label=m, primary_enclosure=raw,
                           independent_dyadic_enclosure=[str(outward(interval[0], True)), str(outward(interval[1], False))]))
    paths += [f'scripts/{PRIMARY}.py', f'scripts/{TAG}.py']
    return dict(status='PASS_UNIVERSAL_ENDPOINT_GATE_PROOF_COVERAGE_COSINE_AUDIT',
        sources={p: sha(p) for p in paths}, recursive_manifests_checked=len(seen),
        complete_profile_counts=dict(edges=52, finite=1825, tails=72),
        closed_box_count=totals['boxes'], independent_pi_bounds=[str(v) for v in pi],
        independent_cosine_enclosures=checks,
        reviewed_boundaries=['h=0 retained as a superset', 'tail k=2 zero residual retained',
            'joint k-to2 and unbounded-height limit controlled by uniform escape bound',
            'closed Coxeter points and high-label ray', 'overlapping closed dyadic children'],
        scope='Independent source/proof/coverage review and alternative rational cosine enclosures. Full primary saved replay is recorded separately. Quartic coefficients and all individual box decisions are not recomputed by a second implementation.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    result = run()
    out = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['complete_profile_counts'], result['closed_box_count'], 'boxes')
