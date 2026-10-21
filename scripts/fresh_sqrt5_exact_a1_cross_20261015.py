#!/usr/bin/env python3
"""Partial rank27 optimization with exact A1 counts and normalized E7 crosses.

At a vertex with a>0 isolated A1 components, assign u/a to an A1--E7
marked pair. Summing all such pairs contributes u per E7 root, so the
vertex separation is additive once the exact A1 count is fixed.
"""
import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_parameter_weight_stress_20261015_cross as X

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_exact_a1_cross_20261015'
BOUND = 27
G, P, C = X.G, X.P, X.C


def cross_corners(row, orbits):
    es = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
    answer = []
    for a, b in product((0, 2), (1, 3)):
        parts = G.pieces(G.graph(row['white'] + [a, b], es))
        left = next(g for g in parts if a in g)
        right = next(g for g in parts if b in g)
        key = None
        if b not in left:
            comp, root = ((right, b) if len(left) == 1 else
                          (left, a) if len(right) == 1 else (None, None))
            if comp is not None and G.family(comp) == 'E7':
                key = 'A1xE7:' + json.dumps(G.code(comp, (root,)), separators=(',', ':'))
                ids = sorted(comp)
                orbits.setdefault(key, dict(key=key, type='E7', nodes=ids,
                    labels=[comp[u].get(v, 2) for u, v in combinations(ids, 2)],
                    marked_root=root))
        answer.append(key)
    return answer


def build():
    old = P.build()
    faces = [f for f in old['faces'] if f['row']['nodes'] - 2 <= BOUND]
    cross = {}
    for f in faces:
        f['cross_corners'] = cross_corners(f['row'], cross)
    connected = {r['key']: r for r in old['orbits'] if r['rank'] <= BOUND}
    variables = [(k, a) for k, r in sorted(connected.items())
                 for a in range(BOUND - r['rank'] + 1)]
    variables += [(k, a) for k in sorted(cross) for a in range(1, BOUND - 7 + 1)]
    index = {v: i for i, v in enumerate(variables)}
    cases, unique = [], {}
    for fi, f in enumerate(faces):
        residual = BOUND - (f['row']['nodes'] - 2)
        for external in range(residual + 1):
            if residual - external == 1:
                continue
            co, baseline = Counter(), F(0)
            for c, ck in zip(f['corners'], f['cross_corners']):
                a = c['core_forest'].count('A1') + external
                baseline += F(c['baseline'])
                if c['orbit'] is not None:
                    co[index[c['orbit'], a]] += 1
                if ck is not None:
                    assert a > 0
                    co[index[ck, a]] += F(1, a)
            key = tuple(sorted(co.items())), baseline - 1
            rec = dict(face=fi, origins=f['origins'], external_A1=external)
            unique.setdefault(key, []).append(rec)
            cases.append(dict(**rec, baseline=str(baseline),
                              coefficients=[[j, str(v)] for j, v in sorted(co.items())]))
    rows, rhs, meta = [], [], []
    for (co, b), witnesses in sorted(unique.items()):
        rows.append({j: -v for j, v in co})
        rhs.append(b)
        meta.append(dict(kind='face', witnesses=witnesses))
    components = []
    for name, g in G.models():
        rank = len(g)
        if rank > BOUND:
            continue
        coef = [Counter() for _ in range(BOUND + 1)]
        baseline = sum((G.base(g, a, b) for a, b in combinations(g, 2)), F(0))
        for a, b in combinations(g, 2):
            key = P.B.code(g, (a, b))
            if key in connected:
                for count in range(BOUND - rank + 1):
                    coef[count][index[key, count]] += 1
        if name == 'E7':
            for root in g:
                key = 'A1xE7:' + json.dumps(G.code(g, (root,)), separators=(',', ':'))
                if key in cross:
                    for count in range(1, BOUND - rank + 1):
                        coef[count][index[key, count]] += 1
        components.append(dict(type=name, rank=rank, baseline=baseline, coeffs=coef))
    return dict(faces=faces, faceconstraints=cases, variables=variables,
                orbits=list(connected.values()), cross_orbits=list(cross.values()),
                components=components, rows=rows, rhs=rhs, meta=meta, zi=len(variables))


def maxforests(data, x, bound=BOUND):
    answer = []
    a1 = next(i for i, c in enumerate(data['components']) if c['type'] == 'A1')
    for count in range(bound + 1):
        remaining = bound - count
        if remaining == 1:
            continue
        totals = [c['baseline'] + sum(v*x[j] for j, v in c['coeffs'][count].items())
                  for c in data['components']]
        dp = {0: (F(0), ())}
        for rank in range(1, remaining + 1):
            for i, comp in enumerate(data['components']):
                if i == a1 or comp['rank'] > rank:
                    continue
                old = dp.get(rank - comp['rank'])
                if old is None:
                    continue
                value = old[0] + totals[i], tuple(sorted(old[1] + (i,)))
                if rank not in dp or value > dp[rank]:
                    dp[rank] = value
        maximum, forest = dp[remaining]
        answer.append(dict(context=count, maximum=maximum,
            forest=tuple(sorted(forest + (a1,)*count)), component_totals=totals,
            rank_values=[None if n < count or n-count not in dp else dp[n-count][0]
                         for n in range(bound + 1)]))
    return answer


def forestrow(data, count, forest):
    assert sum(data['components'][i]['rank'] for i in forest) == BOUND
    assert sum(data['components'][i]['type'] == 'A1' for i in forest) == count
    row, rhs = Counter({data['zi']: -1}), F(0)
    for i in forest:
        comp = data['components'][i]
        row.update(comp['coeffs'][count])
        rhs -= comp['baseline']
    return dict(row), rhs


def run(saved=None):
    data = build()
    print('faces', len(data['faces']), 'cases', len(data['faceconstraints']),
          'variables', data['zi'], flush=True)
    original = C.maxforests, C.forestrow
    C.maxforests, C.forestrow = maxforests, forestrow
    try:
        solution = C.solve(data) if saved is None else saved['solution']
        proof = C.verify(data, solution)
    finally:
        C.maxforests, C.forestrow = original
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_sqrt5_a1_count_weights_20261015.py',
             'scripts/fresh_sqrt5_parameter_weight_stress_20261015.py',
             'scripts/fresh_sqrt5_parameter_weight_stress_20261015_cross.py',
             'results/fresh_sqrt5_parameter_weight_stress_20261015_cross.json']
    return dict(status='PASS_EXACT_PARTIAL_EXACT_A1_NORMALIZED_CROSS_OPTIMUM',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        dimension=BOUND, faces=len(data['faces']), face_cases=len(data['faceconstraints']),
        variables=data['zi'], variable_order=[dict(orbit=k, exact_A1_count=a)
                                            for k, a in data['variables']],
        connected_orbits=data['orbits'], cross_orbits=data['cross_orbits'],
        solution=solution, exact_certificate=proof, target='143/2',
        meets_target=F(solution['optimum']) <= F(143, 2),
        rule='Connected-pair bonuses depend on exact entire A1 count. For a>0 '
             'A1 components each A1/E7 marked pair receives its rooted bonus u_a/a. '
             'Thus each E7 root contributes exactly u_a to the total vertex weight.',
        scope='Conservative partial frozen face catalogue; finite attachment '
              'coverage not asserted. No dimension theorem or independent review.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    path = ROOT / f'results/{TAG}.json'
    saved = json.loads(path.read_text()) if args.verify else None
    result = json.loads(json.dumps(run(saved)))
    if args.verify:
        assert result == saved
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['solution']['optimum'], 'target', result['target'])
