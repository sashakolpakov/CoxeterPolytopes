#!/usr/bin/env python3
"""Exact face identities and a finite-volume counterexample, not a dimension bound.

The all-rank combinatorial proofs are in the companion report. This standalone
Fraction checker verifies product-face counts, coefficient identities, and the
literal Coxeter square pyramid. It imports no repository arithmetic.
"""
import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_finite_volume_simplex_truncation_20261017'


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def partitions(n, low=1):
    if n == 0:
        yield ()
    for a in range(low, n + 1):
        for tail in partitions(n - a, a):
            yield (a,) + tail


def mul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def sub(a, b):
    c = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a): c[i] += x
    for i, x in enumerate(b): c[i] -= x
    while len(c) > 1 and c[-1] == 0: c.pop()
    return c


def scale(a, c):
    return [x * c for x in a]


def link_counts(ranks):
    d = sum(ranks)
    V = math.prod(r + 1 for r in ranks)
    E = F(V * d, 2)
    T = F(V * sum(r * (r - 1) for r in ranks), 6)
    Q = F(V * sum(a * b for a, b in itertools.combinations(ranks, 2)), 4)
    polynomial = [F(1)]
    for r in ranks:
        polynomial = mul(polynomial, [F(math.comb(r + 1, k + 1)) for k in range(r + 1)])
    assert polynomial[0] == V and polynomial[1] == E
    assert (polynomial[2] if d >= 2 else 0) == T + Q
    assert 3 * T + 4 * Q == math.comb(d, 2) * V
    assert all(x.denominator == 1 for x in [E, T, Q])
    demand = 2 * T + Q
    assert demand / V == F(5 * sum(r*r for r in ranks) + 3*d*d - 8*d, 24)
    assert demand / V >= F(d*(d-1), 8)
    return dict(ranks=list(ranks), dimension=d, vertices=V, edges=int(E),
                triangles=int(T), squares=int(Q), two_faces=int(T+Q),
                artificial_demand=int(demand), demand_per_new_vertex=str(demand/V))


def eta(n):
    return F(4*(n-1), n-2) if n % 2 == 0 else F(4*n, n-1)


def budget(n):
    return F(n*(n-6), 8) if n % 2 == 0 else F((n-1)*(n-5), 8)


def all_rank_identities():
    # Polynomial variable m; separately n=2m and n=2m+1.
    records = []
    for parity in [0, 1]:
        n = [F(parity), F(2)]
        d = sub(n, [1])
        lower = scale(mul(d, sub(d, [1])), F(1, 8))
        b = (scale(mul(n, sub(n, [6])), F(1, 8)) if parity == 0
             else scale(mul(sub(n, [1]), sub(n, [5])), F(1, 8)))
        margin = (scale([3*parity+2, 6], F(1, 8)) if parity == 0
                  else scale(sub(n, [1]), F(3, 8)))
        identity = sub(sub(lower, b), margin)
        assert identity == [0]
        records.append(dict(parity=parity, variable='n=2m+parity',
                            zero_coefficients=list(map(str, identity)),
                            artificial_budget_margin_coefficients=list(map(str, margin))))
    # Consistency of the original-face and simple-vertex budgets; no extrapolation.
    for n in range(4, 1001):
        assert (5/eta(n)-1) * math.comb(n, 2) == budget(n)
    return records


def dot(x, y):
    return -x[0]*y[0] + sum(a*b for a, b in zip(x[1:], y[1:]))


def solve(a, b):
    n = len(b)
    m = [[F(x) for x in row] + [F(y)] for row, y in zip(a, b)]
    for i in range(n):
        pivot = next((j for j in range(i, n) if m[j][i]), None)
        if pivot is None: return None
        m[i], m[pivot] = m[pivot], m[i]
        q = m[i][i]
        m[i] = [x/q for x in m[i]]
        for j in range(n):
            if j != i:
                q = m[j][i]
                m[j] = [x-q*y for x, y in zip(m[j], m[i])]
    return [row[-1] for row in m]


def polytope(normals):
    vertices = set()
    for inds in itertools.combinations(range(len(normals)), 3):
        x = solve([normals[i][1:] for i in inds], [normals[i][0] for i in inds])
        if x is None: continue
        v = tuple([F(1)] + x)
        if all(dot(v, r) <= 0 for r in normals): vertices.add(v)
    vertices = sorted(vertices)
    active = [{i for i, r in enumerate(normals) if dot(v, r) == 0} for v in vertices]
    edges = set()
    for i, j in itertools.combinations(range(len(normals)), 2):
        ends = tuple(k for k, a in enumerate(active) if i in a and j in a)
        if len(ends) == 2: edges.add(ends)
    faces = [[k for k, a in enumerate(active) if i in a] for i in range(len(normals))]
    assert all(len(face) >= 3 for face in faces)
    assert len(vertices) - len(edges) + len(faces) == 2
    return vertices, active, sorted(edges), faces


def pyramid():
    normals = [tuple(map(F, r)) for r in
               [(0,0,0,-1), (1,2,0,1), (1,0,2,1), (1,-2,0,1), (1,0,-2,1)]]
    G = [[dot(a,b) for b in normals] for a in normals]
    H = tuple(map(F, (2,0,0,1)))
    assert dot(H,H) == -3 and all(dot(H,r) == -1 for r in normals)
    for i,j in itertools.combinations(range(5),2):
        assert G[i][j] <= 0
        cyc = 4*G[i][j]**2/(G[i][i]*G[j][j])
        assert cyc in [0,1,4] # right, pi/3, parallel
    V, A, E, faces = polytope(normals)
    assert [len(V),len(E),len(faces)] == [5,8,5]
    assert sum(dot(v,v)==0 for v in V) == 1
    assert all(dot(v,v)<=0 for v in V)
    assert all(len(a)==3 for v,a in zip(V,A) if dot(v,v)<0)
    assert sorted(len(f) for f in faces) == [3,3,3,3,4]
    assert all(sum(dot(V[k],V[k])==0 for k in f)==1 for f in faces if len(f)==3)
    # At this finite vertex the two noncompact triangles share the ideal edge.
    v = next(i for i,x in enumerate(V) if x == (1,F(1,2),F(1,2),0))
    nu = next(i for i,x in enumerate(V) if dot(x,x)==0)
    assert A[v] == {0,1,2} and A[nu] == {1,2,3,4}
    assert len(faces[1]) == len(faces[2]) == 3 and len(faces[0]) == 4
    assert tuple(sorted((v,nu))) in E
    cut = tuple(map(F,(1,0,0,2))) # affine truncation z <= t/2
    VV,AA,EE,FF = polytope(normals+[cut])
    assert [len(VV),len(EE),len(FF)] == [8,12,6]
    assert all(dot(x,x)<0 for x in VV)
    assert all(len(a)==3 for a in AA) and all(len(f)==4 for f in FF)
    assert dot(cut,normals[1]) == 1 # truncation is not acute
    link=link_counts((1,1))
    assert [8-5,12-8,6-5] == [link['vertices']-1,link['edges'],link['two_faces']]
    return dict(normals=[[str(x) for x in r] for r in normals],
                gram=[[str(x) for x in row] for row in G],
                strict_center=list(map(str,H)), vertices=[list(map(str,v)) for v in V],
                vertex_norms=[str(dot(v,v)) for v in V], active_sets=[sorted(a) for a in A],
                edges=[list(e) for e in E], faces=faces,
                shared_finite_vertex=v, shared_ideal_vertex=nu,
                original_f_vector=[5,8,5], truncated_f_vector=[8,12,6],
                original_compact_triangles=0, original_noncompact_triangles=4,
                truncated_triangles=0, truncated_squares=6,
                truncation_normal=list(map(str,cut)), positive_truncation_side_product='1')


def run():
    links=[link_counts(p) for d in range(1,17) for p in partitions(d)]
    tables=[]
    for n in [23,24,25,26,163]:
        a=5/eta(n)-1
        assert a>0
        tables.append(dict(dimension=n, eta=str(eta(n)), original_cusp_flag_budget=str(a),
                           finite_vertex_budget=str(budget(n)),
                           original_cusp_budget_per_link_vertex=str(a*F(n-1,2)),
                           artificial_minimum_per_link_vertex=str(F((n-1)*(n-2),8))))
    paths=[f'scripts/{TAG}.py',
           'scripts/fresh_simplex_disjoint_ports_20261017.py',
           'results/fresh_simplex_disjoint_ports_20261017.json',
           'results/fresh_simplex_disjoint_ports_20261017.md']
    return dict(status='PASS_EXACT_FINITE_VOLUME_PORT_TRANSFER_IDENTITIES_AND_SCOPE_COUNTEREXAMPLE',
                scope='Compact simplex faces at finite vertices; exact truncation identities; no new dimension bound.',
                sources={p:sha(p) for p in paths},
                primary_sources={
                    'Balinski1961':'https://msp.org/pjm/1961/11-2/pjm-v11-n2-p03-s.pdf',
                    'Khovanskii1986':'https://www.math.utoronto.ca/askold/1986-Faa-1-english.pdf',
                    'PotyagailoVinberg2005':'https://ems.press/content/serial-article-files/43032'},
                primary_checked_date='2026-10-17',
                all_rank_coefficient_identities=all_rank_identities(),
                product_link_checks=len(links), product_links=links,
                target_budget_table=tables, square_pyramid=pyramid(),
                limitations=['No new reflection-polytope dimension bound.',
                             'No port transitivity for triangles containing ideal vertices.',
                             'Truncation need not be acute or Coxeter.',
                             'Direct simple-at-edges inequality used in weak form.',
                             'Exact bounded partition checks supplement the all-rank proof, not replace it.'])


if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    args=ap.parse_args()
    data=run()
    path=ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text())==data
    else:
        path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(data['status'], 'links=',data['product_link_checks'], 'saved_replay=',args.verify)
