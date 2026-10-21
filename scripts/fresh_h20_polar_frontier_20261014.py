#!/usr/bin/env python3
"""Bounded numerical extreme-ray/polar-wall scout on the frozen 98-wall core."""
import argparse, json, math
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_polar_frontier_20261014'
INPUTS = ['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json',
          'results/fresh_h20_e7_bridge_face_20261014_exact.json',
          'results/fresh_h20_bridge_edge_closure_20261014.json']

def sha(p): return sha256((ROOT / p).read_bytes()).hexdigest()

def load(remove=()):
    old, pre, new = [json.loads((ROOT / p).read_text()) for p in INPUTS]
    for d in [old, pre, new]:
        for p, h in d['sources'].items(): assert sha(p) == h, p
    mp.mp.dps = 90
    def rational(a):
        q=F(a);return mp.mpf(q.numerator)/q.denominator
    rad = list(map(mp.sqrt,old['field']['base_square_root_basis']))
    def v0(a): return sum(rational(q) * r for q, r in zip(a, rad))
    d1 = mp.sqrt(v0(old['field']['discriminant']))
    def v1(a): return v0(a[0]) + d1 * v0(a[1])
    d2 = mp.sqrt(v1(pre['field']['new_positive_sqrt_radical']))
    def v2(a): return v1(a[0]) + d2 * v1(a[1])
    d3 = mp.sqrt(v2(new['field']['new_positive_sqrt_radical']))
    def v3(a): return v2(a[0]) + d3 * v2(a[1])
    N = [mp.matrix(v) for v in ([[v1(a) for a in n] for n in old['normals']] +
                 [[v2(a) for a in pre['new_normal']]] +
                 [[v3(a) for a in n] for n in new['new_normals']])]
    Q = mp.matrix([[rational(q) for q in row] for row in old['form']])
    def dot(a,b):return (a.T*Q*b)[0]
    norms=list(map(rational,old['norms']))+[mp.mpf(2)]+[v3(new['new_normal_norm'])]*2
    for n,q in zip(N,norms):assert abs(dot(n,n)-q)<mp.mpf('1e-65')*max(1,q)
    N=[n*mp.sqrt(2/q) for n,q in zip(N,norms)]
    # Positive combination of spanning retained old roots gives a bounded
    # projective section; omit the two extremely boosted newest roots.
    H=sum((n for i,n in enumerate(N[:95]) if i not in remove),mp.zeros(21,1))
    H/=mp.sqrt(-dot(H,H))
    assert max(dot(n,H) for i,n in enumerate(N) if i not in remove)<0
    E=mp.matrix([v2(a) for a in new['strict_escape']['point']])
    assert abs(dot(E,E)-mp.mpf(20181174579467302741)/10**18)<mp.mpf('1e-60')
    U=[]
    for j in range(20):
        v=mp.matrix([int(i==j) for i in range(21)])
        v+=dot(v,H)*H
        for u in U:v-=dot(v,u)*u
        v/=mp.sqrt(dot(v,v));U.append(v)
    assert max(abs(dot(u,v)-int(i==j)) for i,u in enumerate(U) for j,v in enumerate(U))<mp.mpf('1e-60')
    b=np.array([float(-dot(n,H)) for n in N])
    spatial=np.array([[float(dot(n,u)) for u in U] for n in N])
    ee=np.array([float(dot(E,u)) for u in U]);et=float(dot(E,H))
    return spatial,b,ee,et,np.array(list(H),dtype=float),np.array([list(u)for u in U],dtype=float).T

def run(samples, seed, remove=()):
    spatial,b,ee,et,H,U = load(remove)
    # Decode and orthonormalize at90 decimal digits before using float LP.
    # The three final roots are extremely boosted; direct float norms fail.
    retained=[i for i in range(98) if i not in remove]
    assert min(b[retained]) > 0
    A=spatial/b[:,None]
    rng = np.random.default_rng(seed)
    counts, seen, rows = Counter(), set(), []
    directions = [ee, -ee]
    directions += list(rng.normal(size=(samples, 20)))
    directions += [-ee/np.linalg.norm(ee)+scale*rng.normal(size=20)
                   for scale in [.0001,.001,.01,.1,1] for _ in range(64)]
    for j, c in enumerate(directions):
        sol = linprog(c, A_ub=A[retained], b_ub=np.ones(len(retained)), bounds=[(None, None)]*20, method='highs')
        counts['LP_' + str(sol.status)] += 1
        assert sol.success, sol.message
        active = tuple(i for i in retained if abs(A[i]@sol.x-1)<2e-7)
        if active in seen: continue
        seen.add(active)
        q = float(sol.x @ sol.x - 1)
        if q <= 1e-7:
            counts['causal_extreme'] += 1
            continue
        counts['spacelike_extreme'] += 1
        assert len(active) >= 20
        r = (H + U @ sol.x)*math.sqrt(2/q)
        pair = b*(A@sol.x-1)*math.sqrt(2/q)
        # LP-active rows are prescribed zero; no exact cap is accepted here.
        pair[list(active)]=0
        assert max(pair[retained]) < 1e-6
        bad, labels = [], []
        for i, g in enumerate(pair):
            if i in remove:labels.append('removed');continue
            z = -g/2
            if abs(z) < 2e-7: label = 2
            elif z >= 1-2e-7: label = 'parallel_or_dotted'
            elif z < .5-2e-7: label = None
            else:
                m = math.pi/math.acos(z)
                label = round(m) if abs(m-round(m)) < 1e-5 else None
            if label is None: bad.append([i, float(z)])
            labels.append(label)
        counts['polar_pass' if not bad else 'polar_fail'] += 1
        rows.append(dict(active=list(active), normal=r.tolist(), ray_square=q,
                         escape_pairing=float((et+sol.x@ee)*math.sqrt(2/q)), bad_finite_pairs=bad,
                         all_labels=labels))
        if not bad:
            print('NUMERIC_POLAR_PASS', list(active), 'E=', rows[-1]['escape_pairing'], flush=True)
    rows.sort(key=lambda x:(len(x['bad_finite_pairs']),-x['escape_pairing']))
    paths=INPUTS+[f'scripts/{TAG}.py']
    return dict(status='BOUNDED_NUMERICAL_EXTREME_RAY_POLAR_SCOUT',removed_walls=list(remove),
                sources={p:sha(p) for p in paths}, samples=samples, random_seed=seed,
                chart_precision_decimal_digits=90, counts=dict(counts), unique_active_sets=len(seen), candidates=rows,
                scope='Random and escape-biased LP vertices of the projective slice of the frozen98 core after deleting the explicitly listed walls. Radical evaluation and chart formation use90 decimal digits; LP and angle classification remain floating point. Not a full extreme-ray census, exact cap certificate, or finite-volume claim.')

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--samples',type=int,default=3000);p.add_argument('--seed',type=int,default=20261014)
    p.add_argument('--remove',type=int,nargs='*',default=[]);p.add_argument('--verify',action='store_true')
    a=p.parse_args();d=run(a.samples,a.seed,a.remove)
    suffix='_remove_'+'_'.join(map(str,a.remove)) if a.remove else ''
    out=ROOT/f'results/{TAG}{suffix}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'], d['counts'], 'unique',d['unique_active_sets'],flush=True)
