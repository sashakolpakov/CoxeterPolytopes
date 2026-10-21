#!/usr/bin/env python3
"""Exact all-real gate for 308 E7 single-weight/finite-vertex families.

Also covers the 14 corresponding single-weight/source-null-ray families.
No finite Coxeter-label cutoff or numerical acceptance is used.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact_review as K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e7_transverse_20261014'
BASE = 'results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
VERTICES = 'results/fresh_h20_e7_intrinsic_caps_20261014.json'
Z = K.Z


def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()


def radical_sign(a, b, h):
    """Exact sign of a+b*sqrt(h), h rational nonnegative."""
    assert h >= 0
    sa, sb = K.sign(a), K.sign(b)
    if not h or not sb: return sa
    if not sa or sa == sb: return sb
    delta = K.sub(K.multiply(a, a), K.scale(K.multiply(b, b), h))
    return sa*K.sign(delta)


def enc(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, (tuple, list)): return [enc(y) for y in x]
    if isinstance(x, dict): return {str(k): enc(v) for k, v in x.items()}
    return x


def run():
    d, b = [json.loads((ROOT/p).read_text()) for p in [BASE, VERTICES]]
    for doc in (d, b):
        for p, h in doc['sources'].items(): assert digest(p) == h, p
    assert len(d['normals']) == 87 and len(b['normals']) == 85
    assert d['normals'][:85] == b['normals']
    assert d['form'] == b['form'] and d['E7_support'] == b['E7_support']
    assert len(b['face_finite_vertices']) == 44 and len(b['face_ideal_vertices']) == 2
    Q = [[F(x) for x in row] for row in d['form']]
    sparse = [[(j, a) for j, a in enumerate(row) if a] for row in Q]
    cache = {}

    def dot(x, y):
        if y not in cache:
            cache[y] = tuple(tuple(sum(a*y[j][k] for j, a in row) for k in range(8)) for row in sparse)
        out = Z
        for a, bb in zip(x, cache[y]):
            if any(a) and any(bb): out = K.add(out, K.multiply(a, bb))
        return out

    N = list(map(K.vector, d['normals']))
    H = K.vector(d['strict_center'])
    E = K.vector(d['strict_escape']['point'])
    assert d['strict_escape']['point'] == b['strict_escape']['point']
    assert K.sign(dot(H, H)) < 0 and K.sign(dot(E, H)) < 0
    assert dot(E, E) == K.scalar(F(1336028658157091827, 250000000000000000))
    assert all(K.sign(dot(E, n)) < 0 and K.sign(dot(H, n)) < 0 for n in N)
    S = d['E7_support']; assert S == [12, 13, 14, 15, 18, 31, 65]
    G = [[K.rational(dot(N[i], N[j])) for j in S] for i in S]
    assert K.pd(G) and all(G[i][i] == 2 for i in range(7))
    J = K.inverse(G)
    assert all(a > 0 for row in J for a in row)
    degrees = [J[i][i] for i in range(7)]
    assert degrees == list(map(F, ['7/2', '12', '6', '2', '15/2', '3/2', '4']))
    W = [tuple(tuple(sum(row[j]*N[S[j]][a][k] for j in range(7)) for k in range(8)) for a in range(21)) for row in J]
    assert all(dot(W[i], N[S[j]]) == K.scalar(int(i == j)) for i in range(7) for j in range(7))
    assert all(dot(w, w) == K.scalar(degree) for w, degree in zip(W, degrees))
    A = [K.neg(dot(w, E)) for w in W]
    wall_a = [[K.neg(dot(w, n)) for n in N] for w in W]
    assert all(K.sign(a) > 0 for a in A)
    rows = []; vertex_checks = []
    for vi, row in enumerate(b['face_finite_vertices']):
        V = K.vector(row['ambient_point'])
        q = -F(row['norm'])
        assert q > 0 and dot(V, V) == K.scalar(-q)
        assert K.sign(dot(V, H)) < 0 and all(dot(V, w) == Z for w in W)
        beta = K.neg(dot(V, E)); assert K.sign(beta) > 0
        wall_b = [K.neg(dot(V, n)) for n in N]
        assert all(K.sign(x) >= 0 for x in wall_b)
        assert [i for i, x in enumerate(wall_b) if x == Z] == row['ambient_incident']
        incident = row['ambient_incident']; assert len(incident) == 20 and set(S)<=set(incident)
        assert all(wall_a[j][i]==K.scalar(-int(i==S[j])) for j in range(7) for i in incident)
        sq2=(F(0),F(1))+(F(0),)*6
        scales={i:sq2 if i in b['face_facet_walls'][11:13] else K.scalar(1) for i in incident}
        IG = [[K.rational(K.multiply(K.multiply(scales[i],scales[j]),dot(N[i],N[j]))) for j in incident] for i in incident]
        assert K.pd(IG)
        vertex_checks.append(dict(vertex=vi,negative_norm=q,beta=beta,incident=incident))
        for j, degree in enumerate(degrees):
            # r/c=-omega_j+u*V, u²=(degree-2/c²)/q.  This increases
            # with c; any nonzero Coxeter contact with S[j] needs c>=1.
            rec = dict(vertex=vi,E7_position=j,E7_wall=S[j],weight_norm=degree,
                       vertex_negative_norm=q,escape_A=A[j],escape_beta=beta)
            h1 = (degree-2)/q
            if h1 >= 0 and radical_sign(A[j], K.neg(beta), h1) <= 0:
                rec.update(reason='minimum_coefficient_escape_bound',
                           h_at_c1=h1,sign_at_c1=radical_sign(A[j],K.neg(beta),h1))
            else:
                blocker = None
                for i, (a, bb) in enumerate(zip(wall_a[j], wall_b)):
                    if K.sign(a) <= 0: continue
                    if bb == Z:
                        blocker = (i, None); break
                    gap = K.sub(K.multiply(a, beta), K.multiply(A[j], bb))
                    if K.sign(gap) >= 0:
                        blocker = (i, gap); break
                if blocker is not None:
                    i, gap = blocker
                    rec.update(reason='acute_wall_escape_bound',blocking_wall=i,
                               wall_a=wall_a[j][i],wall_b=wall_b[i],cross_product_gap=gap)
                else:
                    # All labels >=4, and parallel/dotted relations, have
                    # c>=sqrt2.  Only c=1 can remain between 1 and sqrt2.
                    h4 = (degree-1)/q
                    sg = radical_sign(A[j], K.neg(beta), h4)
                    assert sg < 0, (vi,j,'large_c_not_excluded')
                    rec.update(reason='Coxeter_coefficient_gap',h_at_sqrt2=h4,sign_at_sqrt2=sg)
                    if h1 <= 0:
                        rec.update(c1_reason='nonpositive_timelike_coefficient_square',h_at_c1=h1)
                    else:
                        positive = [(i, a, bb) for i, (a, bb) in enumerate(zip(wall_a[j],wall_b)) if radical_sign(a,K.neg(bb),h1)>0]
                        assert positive, (vi,j,'c1_not_excluded')
                        i,a,bb=positive[0]
                        rec.update(c1_reason='positive_pairing',blocking_wall=i,wall_a=a,wall_b=bb,h_at_c1=h1,
                                   c1_pairing_sign=radical_sign(a,K.neg(bb),h1))
            rows.append(rec)
    assert len(rows) == 308
    remaining = [(r['vertex'],r['E7_wall']) for r in rows if r['reason']=='Coxeter_coefficient_gap']
    assert remaining == [(40,65),(41,65),(42,65)]
    assert Counter(r['reason'] for r in rows)==dict(minimum_coefficient_escape_bound=174,acute_wall_escape_bound=131,Coxeter_coefficient_gap=3)
    # A null residual gives r=-c*omega+t*nu with c²=2/degree.
    # The two possible coefficient gaps leave only the degree-two node.
    null_rows=[]
    for ni, row in enumerate(b['face_ideal_vertices']):
        nu=K.vector(row['ambient_point'])
        assert dot(nu,nu)==Z and K.sign(dot(nu,H))<0
        assert all(dot(nu,w)==Z for w in W)
        assert all(K.sign(dot(nu,n))<=0 for n in N)
        for j,degree in enumerate(degrees):
            c2=2/degree
            rec=dict(ideal_vertex=ni,E7_position=j,E7_wall=S[j],weight_norm=degree,coefficient_squared=c2)
            if c2<1:rec['reason']='coefficient_between_zero_and_one'
            elif 1<c2<2:rec['reason']='coefficient_between_one_and_sqrt2'
            else:
                assert c2==1 and j==3
                own=85+ni
                assert dot(nu,N[own])==Z and K.neg(dot(W[j],N[own]))==K.scalar(2)
                assert N[own]==K.vadd(K.vscale(W[j],K.scalar(-1)),K.vscale(nu,K.scalar(2)))
                rec.update(reason='positive_pairing_with_retained_same_cusp_filler',blocking_wall=own,pairing=F(2),
                           duplicate_parameter=F(2),duplicate_escape_pairing=dot(E,N[own]))
                assert K.sign(dot(E,N[own]))<0
            null_rows.append(rec)
    assert len(null_rows)==14
    paths=[BASE,VERTICES,f'scripts/{TAG}.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
    return enc(dict(status='PASS_ALL_308_FINITE_RAY_AND_14_NULL_RAY_TRANSVERSE_FAMILIES',
                    sources={p:digest(p) for p in paths},E7_support=S,E7_inverse=J,weights=W,
                    family='r=-c*omega_j+sqrt((c*c*d_j-2)/(-V*V))*V; norm2; positive root; all real c compatible with retained walls',
                    vertex_checks=vertex_checks,counts=dict(Counter(r['reason'] for r in rows)),records=rows,
                    null_counts=dict(Counter(r['reason'] for r in null_rows)),null_records=null_rows,
                    strict_escape_norm=dot(E,E),
                    scope='Every new actual facet in the 308 single-weight/finite-vertex-ray families, Coxeter-compatible with all retained87 walls, has nonpositive pairing with the saved escape. The 14 single-weight/source-null-ray families give no genuinely new compatible facet; duplicate cusp fillers miss the escape. Any finite coherent collection from these families cannot yield finite volume. Mixed E7 weights, residuals along source edges or interiors, and removal of retained walls are outside the statement.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],result['null_counts'])
