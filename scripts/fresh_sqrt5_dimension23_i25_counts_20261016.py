#!/usr/bin/env python3
"""Complete-face rank23 scout with normalized entire-I2(5) cross bonuses.

All bonus variables depend on exact A1 count and capped or exact I2(5) count.
Every common spectator A1/I2(5) count and every finite forest is accounted for.
The original baseline is unchanged and every added bonus is nonnegative.
"""
import argparse, json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_dimension24_weight_barrier_20261015 as B
B.W.BOUND = 23
B.W.TARGET = F(99, 2)
W = B.W
G = W.G
ROOT = B.ROOT
TAG = 'fresh_sqrt5_dimension23_i25_counts_20261016'
OUT = ROOT / f'results/{TAG}.json'
CAP = 2

def hub_context(b):
    return min(b, CAP) if CAP else b

def hub_contexts():
    return range(CAP + 1) if CAP else range(12)

def cross_corners(row, orbits):
    es = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
    out = []
    for a, b in product((0, 2), (1, 3)):
        parts = G.pieces(G.graph(row['white'] + [a, b], es))
        left = next((c for c in parts if a in c))
        right = next((c for c in parts if b in c))
        key = None
        if b not in left:
            lf, rf = (G.family(left), G.family(right))
            comp, root = (right, b) if lf == 'I2(5)' and rf not in ('A1', 'I2(5)') else (left, a) if rf == 'I2(5)' and lf not in ('A1', 'I2(5)') else (None, None)
            if comp is not None:
                fam = G.family(comp)
                key = 'I25x' + fam + ':' + json.dumps(G.code(comp, (root,)), separators=(',', ':'))
                ids = sorted(comp)
                orbits.setdefault(key, dict(key=key, type=fam, rank=len(comp), nodes=ids, labels=[comp[u].get(v, 2) for u, v in combinations(ids, 2)], marked_root=root))
        out.append(key)
    return out

def build():
    data = W.build()
    oldzi = data['zi']
    orbits = {}
    for f in data['faces']:
        f['i25_cross_corners'] = cross_corners(f['row'], orbits)
    original_variables = list(data['variables'])

    def named(k, count):
        return k + '@I25context=' + str(count)
    variables = [(named(k, p), a) for k, a in original_variables for p in hub_contexts()] + [(named(k, p), a) for k, r in sorted(orbits.items()) for a in range(24 - r['rank'] - 2) for p in hub_contexts() if p]
    index = {v: i for i, v in enumerate(variables)}
    unique = {}
    cases = 0
    for fi, f in enumerate(data['faces']):
        residual = 25 - f['row']['nodes']
        for outside_a in range(residual + 1):
            for outside_b in range((residual - outside_a) // 2 + 1):
                if residual - outside_a - 2 * outside_b == 1:
                    continue
                co = Counter()
                baseline = F(0)
                for c, ak, bk in zip(f['corners'], f['cross_corners'], f['i25_cross_corners']):
                    a = c['core_forest'].count('A1') + outside_a
                    b = c['core_forest'].count('I2(5)') + outside_b
                    baseline += F(c['baseline'])
                    if c['orbit'] is not None:
                        co[index[named(c['orbit'], hub_context(b)), a]] += 1
                    if ak is not None:
                        assert a > 0
                        co[index[named(ak, hub_context(b)), a]] += F(1, a)
                    if bk is not None:
                        assert b > 0
                        co[index[named(bk, hub_context(b)), a]] += F(1, 2 * b)
                unique.setdefault((tuple(sorted(co.items())), baseline - 1), []).append(dict(face=fi, external_A1=outside_a, external_I25=outside_b))
                cases += 1
    data.update(variables=variables, zi=len(variables), rows=[], rhs=[], meta=[], face_case_count=cases, i25_cross_orbits=list(orbits.values()), old_variables=oldzi)
    for (co, rhs), witnesses in sorted(unique.items()):
        data['rows'].append({j: -v for j, v in co})
        data['rhs'].append(rhs)
        data['meta'].append(dict(kind='face', witnesses=witnesses))
    graphs = dict(G.models())
    for comp in data['components']:
        comp['old_context_coeffs'] = {p: [Counter({index[named(original_variables[j][0], p), original_variables[j][1]]: v for j, v in co.items()}) for co in comp['coeffs']] for p in hub_contexts()}
        comp['i25_coeffs'] = {p: [Counter() for _ in range(24)] for p in hub_contexts()}
        if comp['type'] in ('A1', 'I2(5)'):
            continue
        g = graphs[comp['type']]
        for root in g:
            key = 'I25x' + G.family(g) + ':' + json.dumps(G.code(g, (root,)), separators=(',', ':'))
            if key in orbits:
                for a in range(24 - comp['rank'] - 2):
                    for p in hub_contexts():
                        if p:
                            comp['i25_coeffs'][p][a][index[named(key, p), a]] += 1
    return data

def coeff(comp, a, context):
    out = Counter(comp['old_context_coeffs'][context][a])
    if context:
        out.update(comp['i25_coeffs'][context][a])
    return out

def maxforests(data, x):
    cs = data['components']
    ia = next((i for i, c in enumerate(cs) if c['type'] == 'A1'))
    ib = next((i for i, c in enumerate(cs) if c['type'] == 'I2(5)'))
    out = []
    for a in range(24):
        for context in hub_contexts():
            totals = [c['baseline'] + sum((v * x[j] for j, v in coeff(c, a, context).items())) for c in cs]
            dp = {0: (F(0), ())}
            for r in range(1, 24 - a):
                for i, c in enumerate(cs):
                    if i in (ia, ib) or c['rank'] > r:
                        continue
                    old = dp.get(r - c['rank'])
                    if old is None:
                        continue
                    candidate = (old[0] + totals[i], tuple(sorted(old[1] + (i,))))
                    if r not in dp or candidate > dp[r]:
                        dp[r] = candidate
            for b in range((23 - a) // 2 + 1):
                if hub_context(b) != context or 23 - a - 2 * b not in dp:
                    continue
                v, forest = dp[23 - a - 2 * b]
                out.append(dict(context=[a, b], maximum=v + a * totals[ia] + b * totals[ib], forest=tuple(sorted(forest + (ia,) * a + (ib,) * b))))
    return out

def forestrow(data, context, forest):
    a, b = context
    cs = data['components']
    assert sum((cs[i]['rank'] for i in forest)) == 23
    assert sum((cs[i]['type'] == 'A1' for i in forest)) == a and sum((cs[i]['type'] == 'I2(5)' for i in forest)) == b
    row = Counter({data['zi']: -1})
    rhs = F(0)
    for i in forest:
        row.update(coeff(cs[i], a, hub_context(b)))
        rhs -= cs[i]['baseline']
    return (dict(row), rhs)

def solve(data):
    rows = list(data['rows'])
    rhs = list(data['rhs'])
    meta = list(data['meta'])
    cuts = set()
    steps = []

    def add(front):
        key = (tuple(front['context']), tuple(front['forest']))
        if key in cuts:
            return
        cuts.add(key)
        r, b = forestrow(data, *key)
        rows.append(r)
        rhs.append(b)
        meta.append(dict(kind='forest', context=list(key[0]), forest=list(key[1]), types=[data['components'][i]['type'] for i in key[1]]))
    for f in maxforests(data, [F(0)] * data['zi']):
        add(f)
    for iteration in range(150):
        rr = []
        cc = []
        vv = []
        for i, r in enumerate(rows):
            for j, v in r.items():
                rr.append(i)
                cc.append(j)
                vv.append(float(v))
        A = W.C.coo_matrix((vv, (rr, cc)), shape=(len(rows), data['zi'] + 1)).tocsr()
        objective = W.C.np.zeros(data['zi'] + 1)
        objective[-1] = 1
        lp = W.C.linprog(objective, A_ub=A, b_ub=W.C.np.array(list(map(float, rhs))), bounds=(0, None), method='highs')
        assert lp.success, lp.message
        print('I25 cross LP', iteration, 'lower', lp.fun, 'cuts', len(cuts), flush=True)
        steps.append(dict(iteration=iteration, numerical_lower=float(lp.fun), forest_cuts=len(cuts)))
        if lp.fun > F(99, 2) + 1e-07:
            for denominator in (1000, 10 ** 6, 10 ** 9):
                y = [F(float(v)).limit_denominator(denominator) if abs(v) > 1e-09 else F(0) for v in lp.ineqlin.marginals]
                forest_mass = -sum((y[i] for i, m in enumerate(meta) if m['kind'] == 'forest'))
                if forest_mass > 1:
                    y = [v / forest_mass if m['kind'] == 'forest' else v for v, m in zip(y, meta)]
                capacity = Counter()
                demand = Counter()
                for r, v, m in zip(rows, y, meta):
                    for j, c in r.items():
                        if j == data['zi']:
                            continue
                        if m['kind'] == 'forest':
                            capacity[j] -= v * c
                        else:
                            demand[j] += v * c
                scale = min([F(1)] + [capacity[j] / v for j, v in demand.items() if v > 0])
                y = [v * scale if m['kind'] == 'face' else v for v, m in zip(y, meta)]
                try:
                    lower, slack = B.dual_check(data, rows, rhs, y)
                except AssertionError:
                    continue
                if lower > F(99, 2):
                    return dict(kind='exact_model_obstruction', lower_bound=str(lower), terms=[dict(multiplier=str(y[i]), rhs=str(rhs[i]), coefficients=[[j, str(v)] for j, v in sorted(rows[i].items())], metadata=meta[i]) for i in range(len(rows)) if y[i]], numerical_steps=steps)
            raise AssertionError('Exact dual recovery unfinished')
        x = [max(F(0), F(float(v)).limit_denominator(10 ** 8)) for v in lp.x]
        factor = F(1)
        for r, b in zip(data['rows'], data['rhs']):
            if b < 0:
                value = -sum((v * x[j] for j, v in r.items()))
                assert value > 0
                factor = max(factor, -b / value)
        x = [v * factor for v in x]
        fronts = maxforests(data, x)
        maximum = max((f['maximum'] for f in fronts))
        print('exact repaired upper', float(maximum), flush=True)
        if maximum <= F(99, 2):
            rounded = [F(-(-v.numerator * 10000 // v.denominator), 10000) for v in x[:-1]]
            upper = max((f['maximum'] for f in maxforests(data, rounded)))
            if upper <= F(99, 2):
                rounded.append(upper)
                return dict(kind='exact_feasible_rule', exact_primal=list(map(str, rounded)), maximum=str(upper), numerical_steps=steps)
        old = len(cuts)
        for f in fronts:
            if float(f['maximum']) > lp.fun + 1e-07:
                add(f)
        assert len(cuts) > old, 'No new separator; exact repair unresolved'
    raise AssertionError('Bounded solve unfinished')

def verify(data, cert):
    if cert['kind'] == 'exact_feasible_rule':
        x = list(map(F, cert['exact_primal']))
        assert len(x) == data['zi'] + 1 and min(x) >= 0
        assert min((b - sum((v * x[j] for j, v in r.items())) for r, b in zip(data['rows'], data['rhs']))) >= 0
        fronts = maxforests(data, x)
        maximum = max((f['maximum'] for f in fronts))
        assert maximum == F(cert['maximum']) == x[-1] <= F(99, 2)
        return dict(maximum=str(maximum), contexts=len(fronts), complete_exact_face_and_forest_check=True, fronts=[dict(context=f['context'], maximum=str(f['maximum']), types=[data['components'][i]['type'] for i in f['forest']]) for f in fronts])
    rows = []
    rhs = []
    ys = []
    lookup = {(tuple(sorted(r.items())), b): i for i, (r, b) in enumerate(zip(data['rows'], data['rhs']))}
    for t in cert['terms']:
        r = {j: F(v) for j, v in t['coefficients']}
        b = F(t['rhs'])
        m = t['metadata']
        if m['kind'] == 'forest':
            er, eb = forestrow(data, m['context'], m['forest'])
            assert r == er and b == eb
        else:
            assert m == data['meta'][lookup[tuple(sorted(r.items())), b]]
        rows.append(r)
        rhs.append(b)
        ys.append(F(t['multiplier']))
    lower, slack = B.dual_check(data, rows, rhs, ys)
    assert lower == F(cert['lower_bound']) > F(99, 2)
    return dict(lower_bound=str(lower), readable_lower_bound=str(F(lower.numerator * 100 // lower.denominator, 100)), terms=len(rows), exact_weak_duality=True)

def run(saved=None):
    data = build()
    print('I25 model', len(data['faces']), 'faces', data['face_case_count'], 'cases', len(data['rows']), 'unique', data['zi'], 'variables', flush=True)
    cert = solve(data) if saved is None else saved['certificate']
    proof = verify(data, cert)
    manifests = B.sources()
    manifests[str(Path(__file__).relative_to(ROOT))] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(status='PASS_EXACT_DIMENSION23_I25_COUNT_CONTEXT_MODEL', sources=manifests, dimension=23, target='99/2', I25_count_context_cap=CAP, I25_context_mode=('exact' if CAP == 0 else 'capped'), faces=len(data['faces']), face_cases=data['face_case_count'], unique_face_constraints=len(data['rows']), variables=data['zi'], old_variables=data['old_variables'], connected_orbits=data['orbits'], cross_orbits=data['cross_orbits'], new_cross_orbits=data['i25_cross_orbits'], variable_order=[dict(orbit=k.rsplit('@I25context=',1)[0], A1_count=a, I25_context=int(k.rsplit('@I25context=',1)[1])) for k, a in data['variables']], certificate=cert, verification=proof, rule='Inherited nonnegative q2/H4 baseline; all connected, normalized A1-to-component and normalized I2(5)-to-component bonuses depend on exact A1 count a and hub context min(b,cap), or exact b if cap=0. I2(5) crosses have weight u/(2b), destinations exclude A1/I2(5); A1 crosses have weight v/a. Every common spectator count(a,b) is included; remaining rank0 or>=2 is filled by non-A1/non-I2(5) components. Complete finite-forest separation fixes exact(a,b).', scope='Prototype exact weight-model calculation only. No independent review or global bound claimed; any feasible rule still needs independent arithmetic and complete face-coverage composition.')
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--refresh', action='store_true')
    ap.add_argument('--hub-cap', type=int, choices=(0, 1, 2, 3), default=2)
    args = ap.parse_args()
    assert not (args.verify and args.refresh)
    CAP = args.hub_cap
    TAG += '_cap' + str(CAP) if CAP else '_exact'
    OUT = ROOT / f'results/{TAG}.json'
    saved = json.loads(OUT.read_text()) if args.verify or args.refresh else None
    out = json.loads(json.dumps(run(saved), default=str))
    if args.verify:
        assert out == saved
    else:
        OUT.write_text(json.dumps(out, indent=2) + '\n')
    proof = out['verification']
    print(out['status'], out['certificate']['kind'], proof.get('maximum', proof.get('readable_lower_bound')))
