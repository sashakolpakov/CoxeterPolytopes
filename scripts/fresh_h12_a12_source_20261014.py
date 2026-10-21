#!/usr/bin/env python3
"""Small complete classical I_(12,1) source with an actual A12 corner.

Standalone rational arithmetic. Enumerates every dominant norm1/2 root
through distance height9, then all91 projective active sets of size12.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h12_a12_source_20261014'


def ip(x, y): return sum(a*b for a, b in zip(x[1:], y[1:]))-x[0]*y[0]


def squares(total, count, maximum):
    if not count:
        if not total: yield ()
        return
    for a in range(min(maximum, isqrt(total)), -1, -1):
        if total-a*a > (count-1)*a*a: continue
        for tail in squares(total-a*a, count-1, a): yield (a,)+tail


def rref(A):
    A = [list(map(F, row)) for row in A]
    pivots = []; k = 0
    for j in range(len(A[0])):
        p = next((i for i in range(k, len(A)) if A[i][j]), None)
        if p is None: continue
        A[k], A[p] = A[p], A[k]
        t = A[k][j]; A[k] = [x/t for x in A[k]]
        for i in range(len(A)):
            if i != k:
                t = A[i][j]; A[i] = [a-t*b for a, b in zip(A[i], A[k])]
        pivots.append(j); k += 1
        if k == len(A): break
    return A, pivots


def solve(A, b):
    n = len(A)
    R, p = rref([list(row)+[v] for row, v in zip(A, b)])
    if p != list(range(n)): return None
    return tuple(row[-1] for row in R)


def determinant(A):
    A = [list(map(F, r)) for r in A]; out = F(1)
    for j in range(len(A)):
        k = next((i for i in range(j, len(A)) if A[i][j]), None)
        if k is None: return F(0)
        if k != j: A[k], A[j] = A[j], A[k]; out = -out
        p = A[j][j]; out *= p
        for i in range(j+1, len(A)):
            t = A[i][j]/p
            A[i] = [a-t*b for a, b in zip(A[i], A[j])]
    return out


def positive_definite(A):
    A = [list(map(F, r)) for r in A]; pivots = []
    for j in range(len(A)):
        p = A[j][j]; assert p > 0
        pivots.append(p)
        for i in range(j+1, len(A)):
            for k in range(i, len(A)):
                A[k][i] = A[i][k] = A[i][k]-A[i][j]*A[j][k]/p
    return pivots


def encode(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [encode(v) for v in x]
    return x


def run():
    N = []
    for i in range(11):
        r = [0]*13; r[i+1] = -1; r[i+2] = 1; N.append(tuple(r))
    N.append((0,)+(0,)*11+(-1,))
    transcripts = []
    levels = sorted((F(t*t, q), t, q) for q in (1, 2) for t in range(1, isqrt(9*q)+1))
    for height, t, q in levels:
        candidates = [(t,)+x for x in squares(t*t+q, 12, isqrt(t*t+q))]
        assert len(set(candidates)) == len(candidates)
        accepted = []; records = []
        for r in candidates:
            assert ip(r, r) == q and gcd(*r) == 1
            assert all(2*a % q == 0 for a in r)
            obstruction = next((i for i, w in enumerate(N) if ip(r, w) > 0), None)
            records.append(dict(root=r, obstruction=obstruction))
            if obstruction is None: accepted.append(r)
        assert all(ip(r, s) <= 0 for r, s in combinations(accepted, 2))
        N += accepted
        transcripts.append(dict(height=height, time=t, norm=q, candidates=records, accepted=accepted))
    assert N[12:] == [(1,)+(1,)*3+(0,)*9, (3,)+(1,)*11+(0,)]
    assert len(N) == 14 and sum(len(r['candidates']) for r in transcripts) == 22
    G = [[ip(r, s) for s in N] for r in N]; counts = Counter()
    for i, j in combinations(range(14), 2):
        g = G[i][j]; assert g <= 0
        if not g: label = 2
        else:
            cosine2 = F(g*g, G[i][i]*G[j][j])
            assert cosine2 in (F(1, 4), F(1, 2))
            label = 3 if cosine2 == F(1, 4) else 4
        counts[str(label)] += 1
    H = (100,)+tuple(range(12, 0, -1))
    assert ip(H, H) < 0 and all(ip(H, r) < 0 for r in N)
    facets = []
    for r in N:
        h = F(ip(H, r), ip(r, r)); p = tuple(F(a)-h*b for a, b in zip(H, r))
        assert ip(p, p) < 0 and ip(p, H) < 0
        assert [i for i, s in enumerate(N) if ip(p, s) == 0] == [len(facets)]
        assert all(ip(p, s) <= 0 for s in N); facets.append(p)

    # x0=1 is bounded: x1>=...>=x12>=0 and x1+x2+x3<=1.
    # Thus every extreme point has a rank12 active subset of the fourteen rows.
    vertices = {}; active_checks = 0
    for J in combinations(range(14), 12):
        active_checks += 1
        x = solve([N[i][1:] for i in J], [N[i][0] for i in J])
        if x is None: continue
        p = (F(1),)+x
        if any(ip(p, r) > 0 for r in N): continue
        active = [i for i, r in enumerate(N) if ip(p, r) == 0]
        q = ip(p, p); assert q <= 0
        if q < 0:
            assert len(active) == 12
            positive_definite([[G[i][j] for j in active] for i in active])
        vertices[p] = dict(point=p, norm=q, incident=active)
    assert active_checks == 91 and len(vertices) == 37
    assert sum(row['norm'] < 0 for row in vertices.values()) == 35
    cusps = []
    for point, row in vertices.items():
        if row['norm']: continue
        J = row['incident']; unseen = set(J); blocks = []; quotient_rank = 0
        while unseen:
            todo = [min(unseen)]; component = set()
            while todo:
                i = todo.pop()
                if i in component: continue
                component.add(i)
                todo.extend(j for j in J if j != i and G[i][j] and j not in component)
            unseen -= component; K = sorted(component)
            B = [[G[i][j] for j in K] for i in K]
            RR, piv = rref(B); assert len(piv) == len(K)-1
            free = next(i for i in range(len(K)) if i not in piv)
            weights = [F(0)]*len(K); weights[free] = 1
            for rr, i in zip(RR, piv): weights[i] = -rr[free]
            if all(x < 0 for x in weights): weights = [-x for x in weights]
            assert all(x > 0 for x in weights)
            minor = [[B[i][j] for j in range(len(K)) if j != free] for i in range(len(K)) if i != free]
            positive_definite(minor)
            total = tuple(sum(weights[j]*N[k][i] for j, k in enumerate(K)) for i in range(13))
            mu = total[0]; assert mu > 0 and total == tuple(mu*x for x in point)
            quotient_rank += len(K)-1
            blocks.append(dict(incident=K, positive_kernel=weights, null_multiple=mu))
        assert quotient_rank == 11
        cusps.append(dict(point=point, incident=J, blocks=blocks, euclidean_quotient_rank=quotient_rank))
    assert len(cusps) == 2

    chain = list(range(11))+[13]
    A = [[G[i][j] for j in chain] for i in chain]
    assert A == [[2 if i == j else -int(abs(i-j) == 1) for j in range(12)] for i in range(12)]
    v = (11,)+(3,)*12
    assert ip(v, v) == -13 and all(ip(v, N[i]) == 0 for i in chain)
    assert ip(v, N[12]) == -2 and ip(v, N[11]) == -3
    assert (F(1),)+(F(3, 11),)*12 in vertices
    height = F(v[0]*v[0], -ip(v, v))-1
    assert height == F(108, 13) and height < 9
    weights = []
    for k in range(1, 13):
        coeff = [F(min(i, k)*(13-max(i, k)), 13) for i in range(1, 13)]
        w = tuple(sum(c*N[j][i] for c, j in zip(coeff, chain)) for i in range(13))
        assert all(ip(w, N[j]) == int(n == k-1) for n, j in enumerate(chain))
        weights.append(w)
    T = tuple(F(x, 13) for x in v)
    assert ip(T, T) == -F(1, 13)
    cap_r = tuple(-x+2*t for x, t in zip(weights[2], T))
    cap_b = tuple(-x+3*t for x, t in zip(weights[10], T))
    assert cap_r == N[12] and cap_b == N[11] and ip(cap_r, cap_b) == 0

    # The glue statement is witnessed in actual integer I_(12,1) coordinates.
    glue = tuple(x-5*t for x, t in zip(weights[0], T))
    assert all(x.denominator == 1 for x in glue) and ip(glue, glue) == -1
    basis = [N[i] for i in chain[:11]]+[v, glue]
    basis_det = determinant(list(zip(*basis))); assert abs(basis_det) == 1
    assert ip(glue, v) == 5 and ip(v, v) == -13 and gcd(5, 13) == 1
    assert tuple(13*g-13*l+5*x for g, l, x in zip(glue, weights[0], v)) == (0,)*13
    assert determinant(A) == 13

    path = f'scripts/{TAG}.py'
    return encode(dict(status='PASS_COMPLETE_I12_1_SOURCE_WITH_ACTUAL_A12_CORNER',
        sources={path: sha256((ROOT/path).read_bytes()).hexdigest()}, form_diagonal=[-1]+[1]*12,
        root_norm_alphabet=[1, 2], complete_vinberg_height=9, dominant_candidate_count=22,
        vinberg_levels=transcripts, normals=N, gram=G, pair_labels=dict(counts),
        strict_center=H, strict_facet_points=facets, projective_active_sets=91,
        projective_vertices=list(vertices.values()), finite_vertex_count=35, ideal_vertex_count=2,
        bounded_cusp_links=cusps, corner=dict(incident=chain, point=v, norm=-13, distance_height=height, type='A12'),
        normalized_tangent_time=T, fundamental_weights=weights,
        two_nonincident_facets=dict(norm2_root=N[12], norm1_root=N[11], formulas=['-omega3+2T','-omega11+3T'], mutual_pair=0),
        unimodular_glue=dict(root_basis=chain, glue=glue, glue_formula='omega1-5v/13', glue_norm=-1,
                            basis=basis, integral_basis_determinant=basis_det, perpendicular_exactly='A12'),
        scope='Complete classical14-wall finite-volume H12 source,35finite vertices and2bounded cusp links. A12 corner certified within complete distance height9. Supplies two intrinsic H12 face-completion normals; no ambient H20 lift or compatibility claim.'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    d = run(); path = ROOT/f'results/{TAG}.json'
    if args.verify: assert d == json.loads(path.read_text())
    else: path.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['dominant_candidate_count'], d['finite_vertex_count'], d['ideal_vertex_count'])
