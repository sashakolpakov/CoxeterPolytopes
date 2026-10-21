#!/usr/bin/env python3
"""Complete inherited preferred-pattern maps for the ideal6 small-factor frontier.

No universal replacement of the additional order15 obligations is asserted.
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
TAG = 'fresh_ideal6_small_factor_frontier_20261017'


def basic_pairs(g):
    adj, ds, features = P.distance_data(g)
    near = {i for i in range(len(g)) if any(ds[i][f] <= 6 for f in features)}
    return {p for p in combinations(range(len(g)),2)
            if ds[p[0]][p[1]] <= 6 or (ds[p[0]][p[1]] <= 9 and set(p) & near)}


def fork_pairs(g):
    adj, ds, _ = P.distance_data(g)
    m = len(g)
    return {p for p in combinations(range(m),2)
            if ds[p[0]][p[1]] >= m-4 or len(adj[p[0]]) == len(adj[p[1]]) == 3}


def run():
    old = json.loads((ROOT/f'results/{P.OLD}.json').read_text())
    ds = {r['index']:r for r in json.loads((ROOT/f'results/{P.CENSUS}.json').read_text())['diagrams']}
    patterns = []
    for index, r in enumerate(old['isolation']['cases']):
        if r['component_order_bound'] < 15:
            continue
        di,u = r['left']
        dj,v = r['right']
        g, mp = P.R.unite(P.R.bond_graph(ds[di]), P.R.bond_graph(ds[dj]),
                         tuple(map(tuple,r['match'])))
        ids = r['marked_component']
        h = [[int(2-g[i][j]) if i != j and g[i][j] else 0 for j in ids] for i in ids]
        patterns.append(dict(index=index,nodes=ids,labels=h,
                             marks=[ids.index(x) for x in r['preferred_corner']],
                             order_bound=r['component_order_bound'],method=r['method']))
    assert len(patterns) == 437
    assert Counter(p['order_bound'] for p in patterns) == {15:177,16:180,17:80}
    records = []
    expected = {15:[90,100,95,105],16:[96,107,103,111],17:[102,111,108,127]}
    for m in (15,16,17):
        for kind in 'ABCD':
            target = P.labelled(m,P.A.P.graph(kind,m)[0])
            basic = basic_pairs(target)
            extra = set()
            rows = []
            for pattern in patterns:
                if pattern['order_bound'] < m or len(pattern['nodes']) > m:
                    continue
                maps = P.embeddings(pattern['labels'],target)
                pairs = {tuple(sorted(mp[i] for i in pattern['marks'])) for mp in maps}
                extra |= pairs-basic
                if maps:
                    rows.append(dict(case=pattern['index'],embeddings=maps,
                                     marked_images=[list(p) for p in sorted(pairs)],
                                     extra_images=[list(p) for p in sorted(pairs-basic)]))
            dropped = fork_pairs(target) if kind == 'D' and m >= 16 else set()
            paid = (basic | extra)-dropped
            assert len(paid) == expected[m]['ABCD'.index(kind)]
            possible_drop = fork_pairs(target) if kind == 'D' and m == 15 else set()
            affected = [r['case'] for r in rows if set(map(tuple,r['marked_images'])) & possible_drop]
            if possible_drop:
                assert len(possible_drop) == 9 and len(affected) == 63
            records.append(dict(order=m,kind=kind,target=target,
                                basic_pairs=[list(p) for p in sorted(basic)],
                                extra_pairs=[list(p) for p in sorted(extra)],
                                inherited_dropped_pairs=[list(p) for p in sorted(dropped)],
                                paid_pairs=[list(p) for p in sorted(paid)],cost=len(paid),
                                normalized_cost=str(F(len(paid),m)),case_embeddings=rows,
                                embedding_count=sum(len(r['embeddings']) for r in rows),
                                proposed_unproved_order15_drop=[list(p) for p in sorted(possible_drop)],
                                proposed_drop_affected_cases=affected))
            print('SMALL_IDEAL6',m,kind,'cost',len(paid),'maps',records[-1]['embedding_count'],flush=True)
    stems = [P.OLD,P.CENSUS,P.TAG,'fresh_order16_finite_cross_payments_20261017',
             'fresh_order16_preferred_embeddings_20261017','fresh_ideal6_triangle_payments_20261017']
    sources = [Path(__file__)]
    sources += [ROOT/f'{folder}/{stem}.{ext}' for stem in stems
                for folder,ext in [('scripts','py'),('results','json'),('results','md')]]
    return dict(status='PASS_IDEAL6_SMALL_FACTOR_EMBEDDING_FRONTIER',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources},
                patterns=patterns,records=records,
                order_bound_counts=dict(Counter(p['order_bound'] for p in patterns)),
                scope='Complete induced labelled maps of the437 inherited preferred-component patterns with bound>=15 into classical affine factors of orders15,16,17. Only the already proved order16/17 fork deletions are applied. The analogous order15 deletion is explicitly unproved; no complete new payment theorem or dimension improvement is asserted.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    target = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],flush=True)
