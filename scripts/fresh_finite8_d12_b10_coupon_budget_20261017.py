#!/usr/bin/env python3
"""Exact finite8 budget for the D12-short-leaf and A2/B10 coupons."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_finite8_d12_b10_coupon_budget_20261017'
BASE = 'fresh_finite8_a4_d13_neighbor_coupon_budget_20261017'


def run():
    prior = json.loads((ROOT / f'results/{BASE}.json').read_text())
    assert prior['dimension_maximum'] == 80 and prior['finite_vertex_budget'] == '8n'
    for path, h in prior['sources'].items():
        assert sha256((ROOT / path).read_bytes()).hexdigest() == h, path
    # These are refinements of the inherited same-component unit budget.
    slacks = {r: 8*r-F(r*(r-1), 2) for r in (10, 12, 13)}
    assert slacks == {10: F(35), 12: F(30), 13: F(26)}
    assert slacks[12]-8 == 22 and slacks[13]-8 == 18
    rows = []
    counts = {}
    for n in range(1, 76):
        before = len(rows)
        for t in range(n//11+1):
            for b in range(t+1):
                for c in range(t-b+1):
                    residue = n-11*t-b-2*c
                    if residue < 0:
                        continue
                    for u in range(residue//17+1):
                        for v in range((residue-17*u)//10+1):
                            maximum_a = residue-17*u-10*v
                            coefficient = F(3*t, 2)+u-6+b+F(c, 4)+F(v, 2)
                            for a in sorted({0, maximum_a}):
                                excess = a*coefficient-8*t-22*b-18*c-35*v
                                assert excess <= 0, (n, a, t, b, c, u, v, excess)
                                assert b+c <= t and a+11*t+b+2*c+17*u+10*v <= n
                                rows.append([n, a, t, b, c, u, v, str(excess)])
        counts[str(n)] = len(rows)-before
    assert len(rows) == 13034
    paths = [Path(__file__)] + [ROOT / f'{folder}/{BASE}.{ext}'
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    return dict(status='PASS_FINITE8_D12_SHORT_LEAF_AND_A2_B10_BUDGET',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        dimension_maximum=75, finite_vertex_budget='8n',
        retained_rules=[prior['retained_rule'], prior['additional_rule']],
        additional_rules=[
            'Upgrade the existing half to one for a marked pair in distinct whole finite components A1 and D12, with the D12 mark a short fork leaf. The fork center is not upgraded. Count each pair once.',
            'Add one half for a marked pair in distinct whole finite components A2 and B10, with the B10 mark its terminal node incident to the doubled bond. Either A2 node is eligible. Count each pair once.'],
        component_slacks={f'B10' if r == 10 else f'D{r}': str(q) for r, q in slacks.items()},
        variable_definitions=dict(a='total rank of components of rank at most5',
            t='number of whole D11..D16, A13 or A15 components', b='number of whole D12 components',
            c='number of whole D13 components', u='number of whole A17 components',
            v='number of whole B10 components'),
        rank_constraint='a+11t+b+2c+17u+10v<=n<=75, b+c<=t',
        coupon_bounds=dict(old='a*(3t/2+u)', D13_neighbor='a*c/4',
            D12_upgrade='a*b', A2_B10='a*v/2'),
        surplus_bound='a*(3t/2+u-6+b+c/4+v/2)-8t-22b-18c-35v',
        endpoint_columns=['dimension', 'small_rank', 'old_large_count', 'D12_count',
            'D13_count', 'A17_count', 'B10_count', 'excess'],
        dimension_endpoint_counts=counts, dimension_endpoints=rows,
        counts=dict(distinct_endpoints_all_dimensions=len(rows), positive_endpoints=0),
        scope='Universal finite-vertex budget through dimension75. This counting theorem does not prove that a deficient face reaches these whole components, change any ideal cost, or establish an ambient dimension bound.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    d = run()
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert d == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['counts'], flush=True)
