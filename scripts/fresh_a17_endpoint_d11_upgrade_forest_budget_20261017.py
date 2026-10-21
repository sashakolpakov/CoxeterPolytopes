#!/usr/bin/env python3
"""Exact finite-forest budgets with an A17 deletion and conditional D11 upgrade.

Keep two rank17 types: coupon-bearing A17 of cost135 and D17 of cost136.
The D11 geometric transfer is a separate obligation; the second variant
here is a counting theorem for its proposed rule, not a face theorem.
"""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_a17_endpoint_d11_upgrade_forest_budget_20261017'
MAXIMUM = 75
BASE = 'fresh_exact_finite8_forest_budget_20261017'
LOCAL = 'fresh_finite_a17_endpoint_payment_support_20261017'
INPUTS = (BASE, BASE+'_review', LOCAL,
    'fresh_finite8_leaf9_half_coupon_budget_20261017',
    'fresh_finite8_a17_d11_coupon_budget_20261017',
    'fresh_finite8_a4_d13_neighbor_coupon_budget_20261017',
    'fresh_finite8_d12_b10_coupon_budget_20261017')


def base_twice(rank):
    return rank*(rank-1) if rank <= 17 else 16*rank-64


def residual_type(rank):
    return ('A' if rank <= 9 else 'D')+str(rank)


def large_type(rank):
    return 'A17' if rank == 17 else 'B10' if rank == 10 else 'D'+str(rank)


def direct_cost_twice(components, upgrade):
    counts = Counter(components)
    small_rank = sum(i*counts['A'+str(i)] for i in range(1, 6))
    t = sum(counts['D'+str(i)] for i in range(11, 17))
    cost = sum(base_twice(int(c[1:]))-2*(c == 'A17') for c in components)
    cost += small_rank*(3*t+2*counts['A17'])
    cost += 2*(counts['A1']*counts['D12']+counts['A2']*counts['B10']+
        counts['A4']*counts['D13'])
    if upgrade:
        cost += 2*counts['A1']*counts['D11']
    return cost


def variant(upgrade):
    # D17 has no cross coupon and belongs in the residual domain. A17 is
    # separately enumerated among the large coupon-bearing components.
    residual_ranks = list(range(1, 10))+list(range(17, MAXIMUM+1))
    tables = {}
    @lru_cache(None)
    def residual(key):
        t, u, v, b, c, w = key
        def weight(r):
            value = base_twice(r)
            if r <= 5:
                value += r*(3*t+2*u)
                value += 2*(b+w) if r == 1 else 2*v if r == 2 else 2*c if r == 4 else 0
            return value
        values, witnesses, last = [0], [()], [0]
        for total in range(1, MAXIMUM+1):
            cost, witness, r = max((values[total-r]+weight(r),
                tuple(sorted(witnesses[total-r]+(r,))), r)
                for r in residual_ranks if r <= total)
            values.append(cost)
            witnesses.append(witness)
            last.append(r)
        tables[key] = dict(coupon_counts=list(key), doubled_optima=values, last_rank=last)
        return values, witnesses

    best = [None]*(MAXIMUM+1)
    large_rows = []
    def visit(rank, remaining, counts):
        if rank <= 17:
            for count in range(remaining//rank+1):
                visit(rank+1, remaining-count*rank, counts+[count])
            return
        used = MAXIMUM-remaining
        key = (sum(counts[1:7]), counts[7], counts[0], counts[2], counts[3],
            counts[1] if upgrade else 0)
        values, witnesses = residual(key)
        base = sum(count*(base_twice(r)-2*(r == 17)) for r, count in zip(range(10, 18), counts))
        large = tuple(large_type(r) for r, count in zip(range(10, 18), counts) for _ in range(count))
        large_rows.append(dict(component_counts=counts, used_rank=used,
            coupon_counts=list(key), doubled_base_cost=base))
        for rest in range(remaining+1):
            n = used+rest
            witness = tuple(sorted(large+tuple(residual_type(r) for r in witnesses[rest]),
                key=lambda c: (int(c[1:]), c)))
            candidate = (base+values[rest], witness)
            if best[n] is None or candidate > best[n]:
                best[n] = candidate
    visit(10, MAXIMUM, [])
    assert len(large_rows) == 1608
    maxima = []
    for n, (cost, witness) in enumerate(best):
        assert sum(int(c[1:]) for c in witness) == n
        assert direct_cost_twice(witness, upgrade) == cost
        beta2 = (n*(n-6) if n % 2 == 0 else (n-1)*(n-5))//4 if n >= 3 else None
        maxima.append(dict(dimension=n, doubled_cost=cost, witness_components=list(witness),
            doubled_slack_below_8n=16*n-cost,
            doubled_curvature_slack=None if beta2 is None else beta2-cost))
    assert maxima[68]['doubled_cost'] == maxima[69]['doubled_cost'] == 1088
    assert maxima[69]['doubled_curvature_slack'] == 0
    return dict(rule='A17_endpoint_deleted_plus_A1_D11_short_leaf_full_upgrade' if upgrade else 'A17_endpoint_deleted',
        D11_upgrade=upgrade, geometric_D11_coverage_proved_here=False,
        residual_ranks=residual_ranks, residual_representatives={str(r): residual_type(r) for r in residual_ranks},
        large_types=[large_type(r) for r in range(10, 18)],
        coupon_count_columns=['D11_through_D16', 'A17', 'B10', 'D12', 'D13', 'upgraded_D11'],
        large_count_vectors=large_rows, residual_tables=[tables[k] for k in sorted(tables)],
        dimension_maxima=maxima,
        counts=dict(large_count_vectors=len(large_rows), residual_tables=len(tables), dimensions=75))


def run():
    inputs = {}
    for stem in INPUTS:
        data = json.loads((ROOT/f'results/{stem}.json').read_text())
        assert data['status'].startswith('PASS'), stem
        for name, expected in data['sources'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest() == expected, (stem, name)
        inputs[stem] = data
    assert inputs[LOCAL]['candidate_removed_pair'] == [0, 16]
    assert inputs[LOCAL]['candidate_same_component_cost'] == 135
    variants = [variant(False), variant(True)]
    prior = inputs[BASE]['dimension_maxima']
    for a, b, old in zip(variants[0]['dimension_maxima'], variants[1]['dimension_maxima'], prior):
        assert a['dimension'] == b['dimension'] == old['dimension']
        assert a['doubled_cost'] <= old['doubled_cost'] and a['doubled_cost'] <= b['doubled_cost']
        assert a['doubled_cost'] == b['doubled_cost']
    paths = [Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}' for stem in INPUTS
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    return dict(status='PASS_EXACT_A17_ENDPOINT_AND_CONDITIONAL_D11_FOREST_BUDGETS',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        dimension_maximum=75, variants=variants, new_ambient_dimension_bound=None,
        exact_cost_formula='sum base(type) + a*(3t/2+u) + a1*b12 + a2*b10 + a4*b13 + optional a1*b11',
        rank17_distinction='A17 base135 with two broad large coupon positions; D17 base136 with no cross coupons. Neither dominates the other uniformly.',
        scope='Exact abstract finite-forest maxima for two explicit rules. D11 full-coupon face coverage is conditional on a separate local theorem. No ideal budget or ambient dimension conclusion is asserted.')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    out = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], [(r['D11_upgrade'], r['counts'],
        [(n, r['dimension_maxima'][n]['doubled_cost']) for n in range(68, 76)]) for r in out['variants']], flush=True)
