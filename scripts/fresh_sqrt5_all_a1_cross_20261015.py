#!/usr/bin/env python3
"""Partial rank27 optimization with exact A1 counts and normalized A1/component crosses.

At a vertex with a>0 isolated A1 components, assign u/a to an A1--component
marked pair. Summing all such pairs contributes u per component root, so the
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
TAG = 'fresh_sqrt5_all_a1_cross_20261015'
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
            if comp is not None and len(comp) > 1:
                family = G.family(comp)
                key = 'A1x' + family + ':' + json.dumps(G.code(comp, (root,)), separators=(',', ':'))
                ids = sorted(comp)
                orbits.setdefault(key, dict(key=key, type=family, rank=len(comp), nodes=ids,
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
    variables += [(k, a) for k, r in sorted(cross.items()) for a in range(1, BOUND - r['rank'] + 1)]
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
        if rank > 1:
            for root in g:
                key = 'A1x' + name + ':' + json.dumps(G.code(g, (root,)), separators=(',', ':'))
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


def verify_feasible(data, solution):
    x = list(map(F, solution['exact_primal']))
    assert len(x) == data['zi'] + 1 and min(x) >= 0
    slack = [b - sum(v*x[j] for j, v in r.items())
             for r, b in zip(data['rows'], data['rhs'])]
    assert min(slack) >= 0
    fronts = maxforests(data, x)
    maximum = max(r['maximum'] for r in fronts)
    assert maximum == x[-1] <= F(143, 2)
    return dict(maximum=str(maximum), minimum_face_slack=str(min(slack)),
        fronts=[dict(context=r['context'], maximum=str(r['maximum']),
                     maximizing_forest=[data['components'][i]['type'] for i in r['forest']])
                for r in fronts])


def solve_feasible(data):
    rows, rhs = list(data['rows']), list(data['rhs'])
    zi, cuts, iterations = data['zi'], set(), []
    def add(front):
        key = front['context'], tuple(front['forest'])
        if key not in cuts:
            cuts.add(key)
            row, b = forestrow(data, *key)
            rows.append(row); rhs.append(b)
    for front in maxforests(data, [F(0)]*zi):
        add(front)
    for iteration in range(500):
        rr, cc, vv = [], [], []
        for i, row in enumerate(rows):
            for j, v in row.items():
                rr.append(i); cc.append(j); vv.append(float(v))
        matrix = C.coo_matrix((vv, (rr, cc)), shape=(len(rows), zi+1)).tocsr()
        objective = C.np.zeros(zi+1); objective[-1] = 1
        lp = C.linprog(objective, A_ub=matrix,
            b_ub=C.np.array([float(b) for b in rhs]), bounds=(0, None), method='highs')
        assert lp.success, lp.message
        x = [max(F(0), F(float(v)).limit_denominator(10**8)) for v in lp.x]
        # Repair floating-point proposal errors exactly. Increasing every bonus
        # by a common factor preserves all nonnegative face coefficients.
        factor = F(1)
        for row, b in zip(data['rows'], data['rhs']):
            if b >= 0:
                continue
            total = -sum(v*x[j] for j, v in row.items())
            assert total > 0
            factor = max(factor, -b/total)
        x = [v*factor for v in x]
        fronts = maxforests(data, x)
        maximum = max(r['maximum'] for r in fronts)
        iterations.append(dict(iteration=iteration, numerical_lower_bound=float(lp.fun),
                               repair_factor=str(factor), exact_upper_bound=str(maximum)))
        print('LP', iteration, 'lower', float(lp.fun), 'exact feasible upper', float(maximum), flush=True)
        if maximum <= F(143, 2):
            # A readable exact rule: rounding bonuses upward preserves every
            # face inequality. Recompute the complete vertex maximum afterward.
            grid = 1000
            rounded = [F(-((-v.numerator*grid)//v.denominator), grid) for v in x[:-1]]
            clean_maximum = max(r['maximum'] for r in maxforests(data, rounded))
            if clean_maximum <= F(143, 2):
                rounded.append(clean_maximum)
                return dict(exact_primal=list(map(str, rounded)),
                            feasible_maximum=str(clean_maximum), bonus_grid=grid,
                            iterations=iterations, forest_cuts=len(cuts),
                            optimality_claim=False)
        for front in fronts:
            if float(front['maximum']) > float(lp.fun) + 1e-7:
                add(front)
    raise AssertionError('Feasibility search unfinished')


def run(saved=None):
    data = build()
    print('faces', len(data['faces']), 'cases', len(data['faceconstraints']),
          'variables', data['zi'], flush=True)
    solution = solve_feasible(data) if saved is None else saved['solution']
    proof = verify_feasible(data, solution)
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_sqrt5_a1_count_weights_20261015.py',
             'scripts/fresh_sqrt5_parameter_weight_stress_20261015.py',
             'scripts/fresh_sqrt5_parameter_weight_stress_20261015_cross.py',
             'results/fresh_sqrt5_parameter_weight_stress_20261015_cross.json']
    return dict(status='PASS_EXACT_PARTIAL_EXACT_A1_ALL_COMPONENT_NORMALIZED_CROSS_FEASIBILITY',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        dimension=BOUND, faces=len(data['faces']), face_cases=len(data['faceconstraints']),
        variables=data['zi'], variable_order=[dict(orbit=k, exact_A1_count=a)
                                            for k, a in data['variables']],
        connected_orbits=data['orbits'], cross_orbits=data['cross_orbits'],
        solution=solution, exact_certificate=proof, target='143/2',
        meets_target=F(solution['feasible_maximum']) <= F(143, 2),
        rule='Connected-pair bonuses depend on exact entire A1 count. For a>0 '
             'A1 components each A1/component marked pair receives its rooted bonus u_a/a. '
             'Thus each component root contributes exactly u_a to the total vertex weight.',
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
    print(result['status'], result['solution']['feasible_maximum'], 'target', result['target'])
