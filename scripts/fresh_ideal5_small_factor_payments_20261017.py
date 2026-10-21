#!/usr/bin/env python3
"""Complete inherited preferred-pattern maps and ideal5/R5 factor costs.

This is a local census/counting input. Face coverage and dimension comparison
are composed separately; the distance8 overlap domain is not imported.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_ideal7_finite_cross_payments_20261017 as P

ROOT = P.ROOT
TAG = 'fresh_ideal5_small_factor_payments_20261017'
D13 = 'fresh_affine_d13_long_pair_transfer_20261017'
D14 = 'fresh_affine_d14_long_pair_transfer_20261017'
OLD16 = 'fresh_order16_finite_cross_payments_20261017'


def read(stem):
    return json.loads((ROOT/f'results/{stem}.json').read_text())


def basic_pairs(g):
    adj, ds, features = P.distance_data(g)
    near = {i for i in range(len(g)) if any(ds[i][f] <= 5 for f in features)}
    return {p for p in combinations(range(len(g)), 2)
            if ds[p[0]][p[1]] <= 5 or (ds[p[0]][p[1]] <= 9 and set(p) & near)}


def fork_pairs(g):
    adj, ds, _ = P.distance_data(g)
    return {p for p in combinations(range(len(g)), 2)
            if ds[p[0]][p[1]] >= len(g)-4 or len(adj[p[0]]) == len(adj[p[1]]) == 3}


def census():
    old = read(P.OLD)
    diagrams = {r['index']: r for r in read(P.CENSUS)['diagrams']}
    patterns = []
    assert len(old['isolation']['cases']) == 949
    for index, row in enumerate(old['isolation']['cases']):
        assert row['component_order_bound'] <= 17
        if row['component_order_bound'] < 14:
            continue
        di, _ = row['left']
        dj, _ = row['right']
        gram, _ = P.R.unite(P.R.bond_graph(diagrams[di]), P.R.bond_graph(diagrams[dj]),
                            tuple(map(tuple, row['match'])))
        ids = row['marked_component']
        labels = [[int(2-gram[i][j]) if i != j and gram[i][j] else 0 for j in ids] for i in ids]
        patterns.append(dict(index=index, nodes=ids, labels=labels,
            marks=[ids.index(x) for x in row['preferred_corner']],
            order_bound=row['component_order_bound'], method=row['method']))
    assert len(patterns) == 571
    assert Counter(p['order_bound'] for p in patterns) == {14:134, 15:177, 16:180, 17:80}
    equality_pairs = {
        14: {tuple(sorted(r['marks'])) for r in read(D13)['extreme_equality']['rows']},
        15: {tuple(sorted(r['marks'])) for r in read(D14)['equality']['rows']}}
    old16 = read(OLD16)['exceptional_D15']
    old17 = read(P.TAG)['order17']
    records = []
    expected = {14:[70,87,84,82], 15:[75,100,95,96],
                16:[80,107,103,111], 17:[85,111,108,127]}
    for m in (14, 15, 16, 17):
        for kind in 'ABCD':
            target = P.labelled(m, P.A.P.graph(kind, m)[0])
            basic = basic_pairs(target)
            extra = set()
            rows = []
            for pattern in patterns:
                if pattern['order_bound'] < m or len(pattern['nodes']) > m:
                    continue
                maps = P.embeddings(pattern['labels'], target)
                pairs = {tuple(sorted(mp[i] for i in pattern['marks'])) for mp in maps}
                extra |= pairs-basic
                if maps:
                    rows.append(dict(case=pattern['index'], embeddings=maps,
                        marked_images=[list(p) for p in sorted(pairs)],
                        extra_images=[list(p) for p in sorted(pairs-basic)]))
            distances = P.distance_data(target)[1]
            dropped = (fork_pairs(target) if kind == 'D' else
                       {p for p in combinations(range(m), 2) if distances[p[0]][p[1]] > 5}
                       if kind == 'A' else set())
            affected = sorted(r['case'] for r in rows if set(map(tuple, r['marked_images'])) & dropped)
            if kind == 'D':
                assert len(dropped) == 9
                if m in equality_pairs:
                    assert dropped == equality_pairs[m]
                elif m == 16:
                    assert dropped == set(map(tuple, old16['dropped_pairs']))
                    assert affected == old16['affected_cases']
                else:
                    assert dropped == set(map(tuple, old17['targets']['D']['dropped_pairs']))
                    assert affected == old17['dropped_case_indices']
            paid = (basic | extra)-dropped
            assert basic <= paid | dropped
            assert len(paid) == expected[m]['ABCD'.index(kind)]
            record = dict(order=m, kind=kind, target=target,
                basic_pairs=[list(p) for p in sorted(basic)],
                extra_pairs=[list(p) for p in sorted(extra)],
                dropped_pairs=[list(p) for p in sorted(dropped)],
                paid_pairs=[list(p) for p in sorted(paid)], cost=len(paid),
                normalized_cost=str(F(len(paid), m)), case_embeddings=rows,
                embedding_count=sum(len(r['embeddings']) for r in rows),
                dropped_affected_cases=affected)
            records.append(record)
            print('IDEAL5_SMALL', m, kind, len(paid), record['embedding_count'], flush=True)
    return dict(patterns=patterns, records=records,
                order_bound_counts=dict(Counter(p['order_bound'] for p in patterns)))


def formula(kind, m):
    if kind == 'A':
        return 5*m
    offset, cap, first = {'C':(-15,14,6), 'B':(-6,13,8), 'D':(3,12,10)}[kind]
    return 5*m+offset+sum(min(m-k, cap) for k in range(first, first+4))


def large_costs():
    rows = []
    for m in range(18, 80):
        for kind in 'ABCD':
            target = P.labelled(m, P.A.P.graph(kind, m)[0])
            cost = len(basic_pairs(target))
            assert cost == formula(kind, m)
            assert F(cost, m) <= F(163, 23)
            rows.append(dict(kind=kind, order=m, cost=cost))
    branches = []
    for kind, first, intercepts in [('C',20,[-45,-25,-4,18,41]),
                                     ('B',21,[-44,-23,-1,22,46]),
                                     ('D',22,[-43,-21,2,26,51])]:
        for j, intercept in enumerate(intercepts):
            lo = 18 if j == 0 else first+j-1
            hi = first+j if j < 4 else None
            slope = 9-j
            assert formula(kind, lo) == slope*lo+intercept
            if hi is not None:
                assert formula(kind, hi) == slope*hi+intercept
            at = hi if F(slope) > F(163,23) else lo
            assert at is not None
            slack = F(163,23)*at-(slope*at+intercept)
            assert slack >= 0
            branches.append(dict(kind=kind, lower_order=lo, upper_order=hi,
                slope=slope, intercept=intercept, common_bound='163/23',
                minimum_at=at, minimum_slack=str(slack)))
    assert len(branches) == 15
    assert formula('D',23) == 163
    return dict(classical_checks=rows, all_rank_branches=branches,
        cycle_cost='5m', normalized_common_bound='163/23',
        all_rank_formulas=dict(A='5m',
            C='5m-15+sum(min(m-k,14), k=6..9)',
            B='5m-6+sum(min(m-k,13), k=8..11)',
            D='5m+3+sum(min(m-k,12), k=10..13)'),
        scope='Exact BFS checks corroborate the written all-rank spine and interval proof; the final slope5 branches cover arbitrarily large orders.')


def run():
    out = census()
    out['basic_affine_cost'] = large_costs()
    stems = (P.OLD, P.CENSUS, P.TAG, OLD16, D13, D14,
        'fresh_cycle14_ideal5_transfer_20261017', 'fresh_finite8_ideal5_local_replacement_20261017',
        'fresh_cycle16_half_coupon_transfer_20261017')
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in stems
        for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    out.update(status='PASS_IDEAL5_SMALL_FACTOR_CENSUS_AND_ALL_RANK_COSTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        scope='Complete induced labelled maps of all571 inherited preferred patterns with bound>=14, exact pair unions and validated deletion domains, plus all-rank ideal5/R5 costs. Local face transfer and dimension comparison are separate composition obligations. No distance8 overlap coverage is asserted.')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    data = json.loads(json.dumps(run()))
    target = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], flush=True)
