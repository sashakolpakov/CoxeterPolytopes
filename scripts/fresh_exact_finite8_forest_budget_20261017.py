#!/usr/bin/env python3
"""Exact current finite-corner maxima over every spherical forest through75.

Finite types are dominated by one coupon-preserving representative per rank.
Enumerate all large coupon counts, then solve the remaining additive problem
by exact integer Bellman recursion. Costs are doubled to avoid fractions.
This is a budget theorem, independent of any new geometric transfer.
"""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_exact_finite8_forest_budget_20261017'
MAXIMUM = 75
INPUTS = ('fresh_finite8_leaf9_half_coupon_budget_20261017',
          'fresh_finite8_a17_d11_coupon_budget_20261017',
          'fresh_finite8_a4_d13_neighbor_coupon_budget_20261017',
          'fresh_finite8_d12_b10_coupon_budget_20261017')


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def representative(rank):
    return 'A'+str(rank) if rank <= 9 or rank == 17 else 'B10' if rank == 10 else 'D'+str(rank)


def base_twice(rank):
    return rank*(rank-1) if rank <= 17 else 16*rank-64


def residual_cost_twice(rank, key):
    t, u, v, b, c = key
    value = base_twice(rank)
    if rank <= 5:
        value += rank*(3*t+2*u)
        value += 2*b if rank == 1 else 2*v if rank == 2 else 2*c if rank == 4 else 0
    return value


def direct_cost_twice(ranks):
    counts = Counter(ranks)
    a = sum(r*counts[r] for r in range(1, 6))
    t = sum(counts[r] for r in range(11, 17))
    return (sum(base_twice(r) for r in ranks) + a*(3*t+2*counts[17]) +
            2*(counts[1]*counts[12]+counts[2]*counts[10]+counts[4]*counts[13]))


def run():
    for stem in INPUTS:
        data = json.loads((ROOT/f'results/{stem}.json').read_text())
        assert data['status'].startswith('PASS')
        for name, checksum in data['sources'].items():
            assert digest(ROOT/name) == checksum, (stem, name)
    residual_ranks = list(range(1, 10))+list(range(18, MAXIMUM+1))
    cache = {}
    @lru_cache(None)
    def residual(key):
        values = [0]
        witnesses = [()]
        last = [0]
        for total in range(1, MAXIMUM+1):
            choices = [(values[total-r]+residual_cost_twice(r, key),
                        tuple(sorted(witnesses[total-r]+(r,))), r)
                       for r in residual_ranks if r <= total]
            value, witness, r = max(choices)
            values.append(value); witnesses.append(witness); last.append(r)
        cache[key] = dict(coupon_counts=list(key), doubled_optima=values, last_rank=last)
        return values, witnesses

    best = [None]*(MAXIMUM+1)
    large_rows = []
    def visit(rank, remaining, counts):
        if rank < 18:
            for count in range(remaining//rank+1):
                visit(rank+1, remaining-count*rank, counts+[count])
            return
        used = MAXIMUM-remaining
        key = (sum(counts[1:7]), counts[7], counts[0], counts[2], counts[3])
        values, witnesses = residual(key)
        base = sum(count*base_twice(r) for count, r in zip(counts, range(10, 18)))
        large = tuple(r for r, count in zip(range(10, 18), counts) for _ in range(count))
        large_rows.append(dict(component_counts=counts, used_rank=used,
            coupon_counts=list(key), doubled_base_cost=base))
        for rest in range(remaining+1):
            n = used+rest
            candidate = (base+values[rest], tuple(sorted(large+witnesses[rest])))
            if best[n] is None or candidate > best[n]:
                best[n] = candidate
    visit(10, MAXIMUM, [])
    assert len(large_rows) == 1608 and len(cache) == 377
    maxima = []
    for n, (cost, ranks) in enumerate(best):
        assert sum(ranks) == n and direct_cost_twice(ranks) == cost <= 16*n
        beta_twice = (n*(n-6) if n % 2 == 0 else (n-1)*(n-5))//4 if n >= 3 else None
        maxima.append(dict(dimension=n, doubled_cost=cost, witness_ranks=list(ranks),
            witness_components=[representative(r) for r in ranks], doubled_slack_below_8n=16*n-cost,
            doubled_curvature_slack=None if beta_twice is None else beta_twice-cost))
    assert [maxima[n]['doubled_cost'] for n in range(68, 76)] == [1088, 1096, 1106, 1118, 1132, 1148, 1156, 1166]
    assert maxima[69]['witness_components'] == ['A1', 'A17', 'A17', 'A17', 'A17']
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in INPUTS
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    return dict(status='PASS_EXACT_CURRENT_FINITE8_FOREST_MAXIMA', dimension_maximum=MAXIMUM,
        sources={str(path.relative_to(ROOT)): digest(path) for path in paths},
        representative_types={str(r): representative(r) for r in range(1, MAXIMUM+1)},
        large_coupon_ranks=list(range(10, 18)), residual_ranks=residual_ranks,
        coupon_count_columns=['D11_through_D16', 'A17', 'B10', 'D12', 'D13'],
        large_count_vectors=large_rows,
        residual_tables=[cache[key] for key in sorted(cache)], dimension_maxima=maxima,
        counts=dict(large_count_vectors=1608, residual_tables=377, dimensions=75),
        exact_cost_formula='sum base(r) + a*(3t/2+u) + a1*b12 + a2*b10 + a4*b13',
        scope='Exact maximum of the unchanged corner-payment rule over every abstract finite Coxeter forest through75. Type domination is justified in the report. Maximizing forests need not occur in an actual polytope; no new ideal transfer or ambient dimension bound is asserted.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'],
          [(n, result['dimension_maxima'][n]['doubled_cost']) for n in range(68, 76)], flush=True)
