#!/usr/bin/env python3
"""Exact rank65 finite budget, conditional on deep distance8 half payments.

The local geometric payment theorem is deliberately not asserted here.
No optimizer is needed: strip one large component and use the frozen budget.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import fresh_d16_d17_d14_coupon_forest_budget_20261017 as S

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_depth17_rank65_conditional_budget_20261018'
DEPTH, RANK = 17, 65


def graph_count(kind, rank):
    edges = [(i, i+1, 4 if kind == 'B' and i == 0 else 3)
             for i in range(rank-1)]
    if kind == 'D':
        edges = [(i, i+1, 3) for i in range(rank-3)]
        edges += [(rank-3, rank-2, 3), (rank-3, rank-1, 3)]
    assert len(edges) == rank-1
    adj = {i: set() for i in range(rank)}
    features = set()
    for a, b, label in edges:
        adj[a].add(b)
        adj[b].add(a)
        if label > 3:
            features.update((a, b))
    features.update(i for i in adj if len(adj[i]) != 2)

    def distances(starts):
        d = {i: 0 for i in starts}
        queue = deque(starts)
        while queue:
            i = queue.popleft()
            for j in sorted(adj[i]):
                if j not in d:
                    d[j] = d[i]+1
                    queue.append(j)
        assert len(d) == rank
        return d

    depth = distances(sorted(features))
    pairs = [[i, j] for i in range(rank) for j, distance in distances([i]).items()
             if i < j and distance == 8 and min(depth[i], depth[j]) >= DEPTH]
    pairs.sort()
    count = max(0, rank - {'A': 42, 'B': 43, 'D': 44}[kind])
    assert len(pairs) == count
    return dict(kind=kind, rank=rank, features=sorted(features), pairs=pairs, count=count)


def run():
    old = json.loads((ROOT/f'results/{S.TAG}.json').read_text())
    for name, digest in old['sources'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    assert [r['small_components_only'] for r in old['results']] == [False, True]
    full, small = [r['dimension_maxima'] for r in old['results']]
    assert full[65]['doubled_cost'] == 974 and small[65]['doubled_cost'] == 958
    assert small[65]['forest'] == ['D16', 'D16', 'D16', 'D17']

    partners = ([('A', r) for r in range(1, 66)] +
                [('B', r) for r in range(2, 66)] +
                [('D', r) for r in range(4, 66)] +
                [('E', 6), ('E', 7), ('E', 8), ('F', 4), ('H', 3), ('H', 4), ('I2*', 2)])
    counts, domination, bounds = [], [], []
    cross_checks = 0
    for rank in range(18, 66):
        cost = {}
        for kind in ('A', 'B', 'D'):
            rec = graph_count(kind, rank)
            counts.append(rec)
            # A unit becomes a half, so its doubled cost decreases by one.
            cost[kind] = S.base((kind, rank)) - rec['count']
            for other in partners:
                assert S.cross((kind, rank), other) == 0
                cross_checks += 1
        assert cost['D'] >= cost['A'] and cost['D'] >= cost['B']
        domination.append(dict(rank=rank, new_doubled_costs=cost,
                               representative='D', no_cross_coupons=True))
        remainder = full[65-rank]
        bound = cost['D'] + remainder['doubled_cost']
        bounds.append(dict(large_rank=rank, doubled_large_cost=cost['D'],
                           remaining_rank=65-rank,
                           inherited_doubled_remainder_maximum=remainder['doubled_cost'],
                           remainder_attainer=remainder['forest'], doubled_bound=bound))
    maximum = max(r['doubled_bound'] for r in bounds)
    assert maximum == 953
    assert [r['large_rank'] for r in bounds if r['doubled_bound'] == maximum] == [65]
    attaining_forest = [('D', 16)]*3 + [('D', 17)]
    # No component of rank<=17 can contain a pair of feature depth17.
    assert S.actual(attaining_forest) == 958
    assert maximum < 958 < 2*S.B.curvature(65) == 960
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}'
        for stem in (S.TAG, S.Q.TAG, S.B.TAG)
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    return dict(status='PASS_CONDITIONAL_DEPTH17_RANK65_FINITE_MAXIMUM479',
        sources={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in paths},
        feature_depth=17, dimension=65,
        rule='Halve precisely the current finite unit pairs at distance8 in whole A/B/D components when both marks have feature distance>=17; retain all other current weights and coupons.',
        geometric_payment_coverage_proved_here=False,
        graph_counts=counts, type_domination=domination, large_component_bounds=bounds,
        counts=dict(graphs=len(counts), large_ranks=len(bounds), zero_cross_checks=cross_checks),
        maximum_with_a_component_at_least18='953/2',
        maximum_with_all_components_at_most17='479',
        exact_modified_maximum='479', attaining_forest=['D16', 'D16', 'D16', 'D17'],
        curvature='480', finite_slack='1',
        proof='If all components have rank<=17 the unchanged frozen bound479 is attained. Otherwise select one component of rankr>=18: finite classification leaves A/B/D, D dominates its modified same-component cost, and every cross coupon involving it is zero. The remaining modified forest costs no more than the frozen old maximum at rank65-r. All48 exact upper bounds are at most953/2<479.',
        scope='Complete exact conditional typed finite-forest maximum at rank65 only. Frozen lower-rank and small-component certificates are inherited, not reoptimized. No local face-coverage theorem, new global dimension bound or finite-volume example is asserted.')


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
    print(result['status'], result['counts'], 'maximum', result['exact_modified_maximum'])
