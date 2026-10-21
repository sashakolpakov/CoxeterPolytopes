#!/usr/bin/env python3
"""Complete-face rank21 scout with normalized A1, I2(5), and H3 crosses.

All bonus variables depend on exact A1 count, H3 presence and exact E7/E8 counts, and capped/exact I2(5) count.
Every common spectator A1/I2(5)/H3/E7/E8 count and every finite forest is accounted for.
The original baseline is unchanged and every added bonus is nonnegative.
"""
import argparse, json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_dimension24_weight_barrier_20261015 as B
B.W.BOUND = 21
B.W.TARGET = F(40)
W = B.W
G = W.G
ROOT = B.ROOT
TAG = 'fresh_sqrt5_dimension21_a4_exceptional_20261016'
OUT = ROOT / f'results/{TAG}.json'
CAP = 2

def hub_context(b, h3, e7, e8, a4):
    return (min(b, CAP) if CAP else b, min(h3,1), e7, e8, min(a4,2))

def hub_contexts():
    return [(b,c,e7,e8,a4) for b,c,e7,e8,a4 in product(range(CAP + 1) if CAP else range(12), range(2), range(4), range(3), range(3)) if 7*e7+8*e8+4*a4<=21]

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

def h3_cross_corners(row, orbits):
    es = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
    out = []
    for a, b in product((0, 2), (1, 3)):
        parts = G.pieces(G.graph(row['white'] + [a, b], es))
        left = next((c for c in parts if a in c))
        right = next((c for c in parts if b in c))
        key = None
        if b not in left:
            lf, rf = (G.family(left), G.family(right))
            comp, root = (right, b) if lf == 'H3' and rf not in ('A1', 'I2(5)', 'H3') else (left, a) if rf == 'H3' and lf not in ('A1', 'I2(5)', 'H3') else (None, None)
            if comp is not None:
                fam = G.family(comp)
                key = 'H3x' + fam + ':' + json.dumps(G.code(comp, (root,)), separators=(',', ':'))
                ids = sorted(comp)
                orbits.setdefault(key, dict(key=key, type=fam, rank=len(comp), nodes=ids, labels=[comp[u].get(v, 2) for u, v in combinations(ids, 2)], marked_root=root))
        out.append(key)
    return out

def a4_corners(row, cross, self_orbits):
    es = dict(zip(combinations(range(row['nodes']),2),row['labels']))
    out=[]
    def rooted(comp,root,prefix,orbits):
        fam=G.family(comp)
        key=prefix+fam+':'+json.dumps(G.code(comp,(root,)),separators=(',',':'))
        ids=sorted(comp)
        orbits.setdefault(key,dict(key=key,type=fam,rank=len(comp),nodes=ids,
            labels=[comp[u].get(v,2) for u,v in combinations(ids,2)],marked_root=root))
        return key
    for a,b in product((0,2),(1,3)):
        parts=G.pieces(G.graph(row['white']+[a,b],es))
        left=next(c for c in parts if a in c);right=next(c for c in parts if b in c)
        keys=[];kind=None
        if b not in left:
            lf,rf=G.family(left),G.family(right)
            if lf==rf=='A4':
                kind='self'
                keys=[rooted(left,a,'A4self',self_orbits),rooted(right,b,'A4self',self_orbits)]
            elif lf=='A4' and rf not in ('A1','I2(5)','H3','A4'):
                kind='cross';keys=[rooted(right,b,'A4x',cross)]
            elif rf=='A4' and lf not in ('A1','I2(5)','H3','A4'):
                kind='cross';keys=[rooted(left,a,'A4x',cross)]
        out.append(dict(kind=kind,keys=keys))
    return out

def build():
    data = W.build()
    oldzi = data['zi']
    orbits = {}
    h3orbits = {}
    a4orbits = {}
    selforbits = {}
    for f in data['faces']:
        f['i25_cross_corners'] = cross_corners(f['row'], orbits)
        f['h3_cross_corners'] = h3_cross_corners(f['row'], h3orbits)
        f['a4_corners'] = a4_corners(f['row'],a4orbits,selforbits)
    original_variables = list(data['variables'])

    def named(k, count):
        return k + '@I25context=' + str(count[0]) + '@H3context=' + str(count[1]) + '@E7context=' + str(count[2]) + '@E8context=' + str(count[3]) + '@A4context=' + str(count[4])
    variables = [(named(k, p), a) for k, a in original_variables for p in hub_contexts()] + [(named(k, p), a) for k, r in sorted(orbits.items()) for a in range(22 - r['rank'] - 2) for p in hub_contexts() if p[0]]
    variables += [(named(k,p),a) for k,r in sorted(h3orbits.items()) for a in range(22-r['rank']-3) for p in hub_contexts() if p[1]]
    variables += [(named(k,p),a) for k,r in sorted(a4orbits.items()) for a in range(22-r['rank']-4) for p in hub_contexts() if p[4]]
    variables += [(named(k,p),a) for k,r in sorted(selforbits.items()) for a in range(22-8) for p in hub_contexts() if p[4]>=2]
    index = {v: i for i, v in enumerate(variables)}
    unique = {}
    cases = 0
    for fi, f in enumerate(data['faces']):
        residual = 23 - f['row']['nodes']
        for outside_a in range(residual + 1):
            for outside_b in range((residual - outside_a) // 2 + 1):
                for outside_h3 in range((residual - outside_a - 2*outside_b) // 3 + 1):
                    for outside_e7 in range((residual - outside_a - 2*outside_b - 3*outside_h3)//7+1):
                        for outside_e8 in range((residual - outside_a - 2*outside_b - 3*outside_h3 - 7*outside_e7)//8+1):
                            for outside_a4 in range((residual - outside_a - 2*outside_b - 3*outside_h3 - 7*outside_e7 - 8*outside_e8)//4+1):
                                if residual - outside_a - 2*outside_b - 3*outside_h3 - 7*outside_e7 - 8*outside_e8 - 4*outside_a4 == 1:
                                    continue
                                co = Counter()
                                baseline = F(0)
                                for c, ak, bk, ck, dk in zip(f['corners'], f['cross_corners'], f['i25_cross_corners'], f['h3_cross_corners'], f['a4_corners']):
                                    a = c['core_forest'].count('A1') + outside_a
                                    b = c['core_forest'].count('I2(5)') + outside_b
                                    h3 = c['core_forest'].count('H3') + outside_h3
                                    e7 = c['core_forest'].count('E7') + outside_e7
                                    e8 = c['core_forest'].count('E8') + outside_e8
                                    a4 = c['core_forest'].count('A4') + outside_a4
                                    baseline += F(c['baseline'])
                                    context = hub_context(b, h3, e7, e8, a4)
                                    if c['orbit'] is not None:
                                        co[index[named(c['orbit'], context), a]] += 1
                                    if ak is not None:
                                        assert a > 0
                                        co[index[named(ak, context), a]] += F(1, a)
                                    if bk is not None:
                                        assert b > 0
                                        co[index[named(bk, context), a]] += F(1, 2*b)
                                    if ck is not None:
                                        assert h3 > 0
                                        co[index[named(ck,context),a]] += F(1,3*h3)
                                    if dk['kind']=='cross':
                                        assert a4>0 and len(dk['keys'])==1
                                        co[index[named(dk['keys'][0],context),a]] += F(1,4*a4)
                                    elif dk['kind']=='self':
                                        assert a4>=2 and len(dk['keys'])==2
                                        for key in dk['keys']:
                                            co[index[named(key,context),a]] += F(1,4*(a4-1))
                                unique.setdefault((tuple(sorted(co.items())), baseline - 1), []).append(dict(face=fi, external_A1=outside_a, external_I25=outside_b, external_H3=outside_h3, external_E7=outside_e7, external_E8=outside_e8, external_A4=outside_a4))
                                cases += 1
    data.update(variables=variables, zi=len(variables), rows=[], rhs=[], meta=[], face_case_count=cases, i25_cross_orbits=list(orbits.values()), h3_cross_orbits=list(h3orbits.values()), a4_cross_orbits=list(a4orbits.values()), a4_self_orbits=list(selforbits.values()), old_variables=oldzi)
    for (co, rhs), witnesses in sorted(unique.items()):
        data['rows'].append({j: -v for j, v in co})
        data['rhs'].append(rhs)
        data['meta'].append(dict(kind='face', witnesses=witnesses))
    graphs = dict(G.models())
    for comp in data['components']:
        comp['old_context_coeffs'] = {p: [Counter({index[named(original_variables[j][0], p), original_variables[j][1]]: v for j, v in co.items()}) for co in comp['coeffs']] for p in hub_contexts()}
        comp['i25_coeffs'] = {p: [Counter() for _ in range(22)] for p in hub_contexts()}
        comp['h3_coeffs'] = {p: [Counter() for _ in range(22)] for p in hub_contexts()}
        if comp['type'] in ('A1', 'I2(5)'):
            continue
        g = graphs[comp['type']]
        for root in g:
            key = 'I25x' + G.family(g) + ':' + json.dumps(G.code(g, (root,)), separators=(',', ':'))
            if key in orbits:
                for a in range(22 - comp['rank'] - 2):
                    for p in hub_contexts():
                        if p[0]:
                            comp['i25_coeffs'][p][a][index[named(key, p), a]] += 1
            if comp['type'] == 'H3':
                continue
            key = 'H3x' + G.family(g) + ':' + json.dumps(G.code(g,(root,)), separators=(',',':'))
            if key in h3orbits:
                for a in range(22-comp['rank']-3):
                    for p in hub_contexts():
                        if p[1]:
                            comp['h3_coeffs'][p][a][index[named(key,p),a]] += 1
    for comp in data['components']:
        comp['a4_coeffs']={p:[Counter() for _ in range(22)] for p in hub_contexts()}
        if comp['type'] in ('A1','I2(5)','H3'):
            continue
        g=graphs[comp['type']]
        for root in g:
            is_self=comp['type']=='A4'
            key=('A4self' if is_self else 'A4x')+G.family(g)+':'+json.dumps(G.code(g,(root,)),separators=(',',':'))
            if key not in (selforbits if is_self else a4orbits):
                continue
            for a in range(22-8 if is_self else 22-comp['rank']-4):
                for context in hub_contexts():
                    if context[4]>=(2 if is_self else 1):
                        comp['a4_coeffs'][context][a][index[named(key,context),a]]+=1
    return data

def coeff(comp, a, context):
    out = Counter(comp['old_context_coeffs'][context][a])
    if context[0]:
        out.update(comp['i25_coeffs'][context][a])
    if context[1]:
        out.update(comp['h3_coeffs'][context][a])
    if context[4]:
        out.update(comp['a4_coeffs'][context][a])
    return out

def maxforests(data, x):
    cs = data['components']
    ia = next(i for i,c in enumerate(cs) if c['type'] == 'A1')
    ib = next(i for i,c in enumerate(cs) if c['type'] == 'I2(5)')
    ie = next(i for i,c in enumerate(cs) if c['type'] == 'H3')
    i7 = next(i for i,c in enumerate(cs) if c['type'] == 'E7')
    i8 = next(i for i,c in enumerate(cs) if c['type'] == 'E8')
    i4 = next(i for i,c in enumerate(cs) if c['type'] == 'A4')
    out = []
    for a in range(22):
        for context in hub_contexts():
            totals = [c['baseline'] + sum(v*x[j] for j,v in coeff(c,a,context).items()) for c in cs]
            dp = {0: (F(0), ())}
            for r in range(1, 22-a):
                for i,c in enumerate(cs):
                    if i in (ia,ib,ie,i7,i8,i4) or c['rank'] > r:
                        continue
                    old = dp.get(r-c['rank'])
                    if old is None:
                        continue
                    candidate = (old[0]+totals[i], tuple(sorted(old[1]+(i,))))
                    if r not in dp or candidate > dp[r]:
                        dp[r] = candidate
            for a4 in range((21-a)//4+1):
                for e8 in range((21-a-4*a4)//8+1):
                    for e7 in range((21-a-4*a4-8*e8)//7+1):
                        for h3 in range((21-a-4*a4-8*e8-7*e7)//3+1):
                            remaining = 21-a-4*a4-8*e8-7*e7-3*h3
                            for b in range(remaining//2+1):
                                if hub_context(b,h3,e7,e8,a4) != context or remaining-2*b not in dp:
                                    continue
                                v,forest = dp[remaining-2*b]
                                out.append(dict(context=[a,b,h3,e7,e8,a4], maximum=v+a*totals[ia]+b*totals[ib]+h3*totals[ie]+e7*totals[i7]+e8*totals[i8]+a4*totals[i4], forest=tuple(sorted(forest+(ia,)*a+(ib,)*b+(ie,)*h3+(i7,)*e7+(i8,)*e8+(i4,)*a4))))
    return out

def forestrow(data, context, forest):
    a,b,h3,e7,e8,a4 = context
    cs = data['components']
    assert sum(cs[i]['rank'] for i in forest) == 21
    assert sum(cs[i]['type']=='A1' for i in forest)==a
    assert sum(cs[i]['type']=='I2(5)' for i in forest)==b
    assert sum(cs[i]['type']=='H3' for i in forest)==h3
    assert sum(cs[i]['type']=='E7' for i in forest)==e7
    assert sum(cs[i]['type']=='E8' for i in forest)==e8
    assert sum(cs[i]['type']=='A4' for i in forest)==a4
    row = Counter({data['zi']: -1})
    rhs = F(0)
    for i in forest:
        row.update(coeff(cs[i], a, hub_context(b,h3,e7,e8,a4)))
        rhs -= cs[i]['baseline']
    return dict(row),rhs

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
        if lp.fun > F(40) + 1e-07:
            for denominator in (1000, 10 ** 6, 10 ** 9, 10 ** 12):
                y = [F(round(float(v)*denominator),denominator) if v < -(1e-7 if m['kind'] == 'face' else 1e-8) else F(0) for v,m in zip(lp.ineqlin.marginals,meta)]
                forest_mass = -sum((y[i] for i, m in enumerate(meta) if m['kind'] == 'forest'))
                if forest_mass > 1:
                    y = [v / forest_mass if m['kind'] == 'forest' else v for v, m in zip(y, meta)]
                # A rounded negligible forest marginal can leave a zero-capacity
                # column. Discard every face term that demands such a column;
                # subsequent uniform shrinking enforces all remaining columns.
                raw_capacity = Counter()
                for r, v, m in zip(rows, y, meta):
                    if m['kind'] == 'forest' and v < 0:
                        for j,c in r.items():
                            if j != data['zi'] and c > 0:
                                raw_capacity[j] -= v*c
                supported = {j for j,v in raw_capacity.items() if v >= F(1,10**7)}
                y = [F(0) if m['kind'] == 'face' and any(c < 0 and j not in supported for j,c in r.items()) else v for r,v,m in zip(rows,y,meta)]
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
                # Numerical dual tolerance can dominate nearly empty columns.
                # Remove affected face terms before the exact global shrink;
                # this weakens the lower bound and never asserts feasibility.
                bad = {j for j,v in demand.items() if v > 0 and capacity[j] < F(999,1000)*v}
                if bad:
                    dropped = [i for i,(r,m,v) in enumerate(zip(rows,meta,y)) if m['kind']=='face' and v and any(j in bad and c < 0 for j,c in r.items())]
                    print('dual tiny-column pruning',len(bad),len(dropped),float(sum(rhs[i]*y[i] for i in dropped)),flush=True)
                    for i in dropped:y[i]=F(0)
                    demand=Counter()
                    for r,v,m in zip(rows,y,meta):
                        if m['kind']=='face':
                            for j,c in r.items():
                                if j != data['zi']:demand[j]+=v*c
                scale = min([F(1)] + [capacity[j] / v for j, v in demand.items() if v > 0])
                y = [v * scale if m['kind'] == 'face' else v for v, m in zip(y, meta)]
                try:
                    lower, slack = B.dual_check(data, rows, rhs, y)
                except AssertionError:
                    print('dual recovery infeasible', denominator, float(scale), flush=True)
                    continue
                print('dual recovery', denominator, float(scale), float(lower), flush=True)
                if lower > F(40):
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
        if maximum <= F(40):
            rounded = [F(-(-v.numerator * 10000 // v.denominator), 10000) for v in x[:-1]]
            upper = max((f['maximum'] for f in maxforests(data, rounded)))
            if upper <= F(40):
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
        assert maximum == F(cert['maximum']) == x[-1] <= F(40)
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
    assert lower == F(cert['lower_bound']) > F(40)
    return dict(lower_bound=str(lower), readable_lower_bound=str(F(lower.numerator * 100 // lower.denominator, 100)), terms=len(rows), exact_weak_duality=True)

def run(saved=None):
    data = build()
    print('I25 model', len(data['faces']), 'faces', data['face_case_count'], 'cases', len(data['rows']), 'unique', data['zi'], 'variables', flush=True)
    cert = solve(data) if saved is None else saved['certificate']
    proof = verify(data, cert)
    manifests = B.sources()
    manifests[str(Path(__file__).relative_to(ROOT))] = sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(status='PASS_EXACT_DIMENSION21_A4_EXCEPTIONAL_MODEL', sources=manifests, dimension=21, target='40', I25_count_context_cap=CAP, I25_context_mode=('exact' if CAP == 0 else 'capped'), H3_context_mode='presence', H3_count_context_cap=1, E7_context_mode='exact', E7_count_context_cap=0, E8_context_mode='exact', E8_count_context_cap=0, A4_context_mode='capped', A4_count_context_cap=2, faces=len(data['faces']), face_cases=data['face_case_count'], unique_face_constraints=len(data['rows']), variables=data['zi'], old_variables=data['old_variables'], connected_orbits=data['orbits'], cross_orbits=data['cross_orbits'], new_cross_orbits=data['i25_cross_orbits'], h3_cross_orbits=data['h3_cross_orbits'], a4_cross_orbits=data['a4_cross_orbits'], a4_self_orbits=data['a4_self_orbits'], variable_order=[dict(orbit=k.rsplit('@I25context=',1)[0], A1_count=a, I25_context=int(k.rsplit('@I25context=',1)[1].split('@')[0]), H3_context=int(k.rsplit('@H3context=',1)[1].split('@')[0]), E7_context=int(k.rsplit('@E7context=',1)[1].split('@')[0]), E8_context=int(k.rsplit('@E8context=',1)[1].split('@')[0]), A4_context=int(k.rsplit('@A4context=',1)[1])) for k, a in data['variables']], certificate=cert, verification=proof, rule='A4-to-other cross weights are v/(4d), excluding destinations A1/I2(5)/H3/A4. Between roots r,s in distinct A4 components, weight (u_r+u_s)/(4(d-1)) for d>=2. These normalize to one bonus per destination/component root, preserving additivity. All bonuses also depend on A4 count capped at two, with every common outside count and exact vertex count included. All bonus variables also depend on exact E7 and E8 counts. Every common external E7/E8 count is included and exact vertex E7/E8 counts are excluded from residual DP then restored. An additional normalized H3-hub cross w/(3c) joins an H3 component to a rooted destination excluding A1/I2(5)/H3, so cross families are disjoint. Every bonus additionally depends on presence of entire H3 components; the common spectator H3 count is enumerated, and vertex forests fix and exclude it from residual DP. Inherited nonnegative q2/H4 baseline; all connected, normalized A1-to-component and normalized I2(5)-to-component bonuses depend on exact A1 count a and hub context min(b,cap), or exact b if cap=0. I2(5) crosses have weight u/(2b), destinations exclude A1/I2(5); A1 crosses have weight v/a. Every common spectator count(a,b) is included; remaining rank0 or>=2 is filled by non-A1/non-I2(5) components. Complete finite-forest separation fixes exact(a,b).', scope='Prototype exact weight-model calculation only. No independent review or global bound claimed; any feasible rule still needs independent arithmetic and complete face-coverage composition.')
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
