#!/usr/bin/env python3
"""Separate nested-quadratic replay of the frozen changed-support checkpoint.

No producer imports; signs use exact squaring in Q(sqrt2)(sqrt79), rather
than the producer's fixed-point radical intervals. Old68 Coxeter relations
and the inherited cusp-link completeness remain certified dependencies.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_contact_closure_20261014_replacement'
ZERO = (F(0),) * 4


def add(a, b): return tuple(x+y for x, y in zip(a, b))
def neg(a): return tuple(-x for x in a)
def scale(a, b): return tuple(x*b for x in a)
def scalar(a): return (F(a), F(0), F(0), F(0))
def qmul(a, b): return (a[0]*b[0]+2*a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def qadd(a, b): return (a[0]+b[0], a[1]+b[1])


def mul(a, b):
    u, v, x, y = a[:2], a[2:], b[:2], b[2:]
    vy = qmul(v, y)
    return qadd(qmul(u, x), (79*vy[0], 79*vy[1])) + qadd(qmul(u, y), qmul(v, x))


def rs(x): return (x > 0)-(x < 0)


def qs(x):
    a, b = x
    if not b: return rs(a)
    if not a or rs(a) == rs(b): return rs(b)
    d = a*a-2*b*b
    assert d
    return rs(a)*rs(d)


def sign(x):
    u, v = x[:2], x[2:]
    su, sv = qs(u), qs(v)
    if not sv: return su
    if not su or su == sv: return sv
    vv = qmul(v, v)
    d = qadd(qmul(u, u), (-79*vv[0], -79*vv[1]))
    sd = qs(d)
    assert sd  # 79 is not a square in Q(sqrt2).
    return su*sd


def pd(A):
    A = [row[:] for row in A]
    for j in range(len(A)):
        pivot = A[j][j]
        if pivot <= 0: return False
        for i in range(j+1, len(A)):
            for k in range(i, len(A)):
                A[k][i] = A[i][k] = A[i][k]-A[i][j]*A[j][k]/pivot
    return True


def rank(A):
    A = [row[:] for row in A]
    r = 0
    for j in range(len(A[0])):
        k = next((k for k in range(r, len(A)) if A[k][j]), None)
        if k is None: continue
        A[r], A[k] = A[k], A[r]
        p = A[r][j]
        A[r] = [x/p for x in A[r]]
        for k in range(r+1, len(A)):
            p = A[k][j]
            A[k] = [x-p*y for x, y in zip(A[k], A[r])]
        r += 1
        if r == len(A): break
    return r


def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p): return json.loads((ROOT/p).read_text())


def run():
    d = load(f'results/{TAG}.json')
    basepath = 'results/fresh_h20_contact_closure_20261014.json'
    base = load(basepath)
    for doc in (d, base):
        for path, h in doc['sources'].items(): assert digest(path) == h, path
    Q = [[F(x) for x in row] for row in d['form']]
    assert Q == [[F(x['a']) for x in row] for row in base['form']]
    assert all(F(x['b']) == 0 for row in base['form'] for x in row)
    nonzero = [[(j, x) for j, x in enumerate(row) if x] for row in Q]
    vec = lambda row: [tuple(map(F, x)) for x in row]
    N = list(map(vec, d['normals']))
    old = [[(F(x['a']), F(x['b']), F(0), F(0)) for x in row] for row in base['normals']]
    cache = {}

    def cov(v):
        key = tuple(v)
        if key not in cache:
            cache[key] = [tuple(sum((x*v[j][k] for j, x in row), F(0)) for k in range(4)) for row in nonzero]
        return cache[key]

    def ip(v, w):
        out = ZERO
        for x, y in zip(v, cov(w)):
            if x != ZERO and y != ZERO: out = add(out, mul(x, y))
        return out

    keep = d['retained_predecessor_walls']
    assert d['removed_predecessor_walls'] == list(range(61, 76, 2))
    assert keep == [i for i in range(76) if i not in range(61, 76, 2)]
    assert N[:68] == [old[i] for i in keep]
    recipes = d['new_recipes']
    assert [r['candidate_index'] for r in recipes] == [1, 3, 5, 7]
    for i, r in enumerate(recipes):
        u = list(map(F, r['rational_offset']))
        w = list(map(F, r['rational_direction']))
        ta, tb = F(r['parameter_a']), F(r['parameter_b'])
        sg = 1 if i < 2 else -1
        assert (ta, tb) == (F(sg*107, 54), F(sg*11, 54))
        v = [(a+b*ta, F(0), b*tb, F(0)) for a, b in zip(u, w)]
        assert v == N[68+i]
        U, W = list(map(scalar, u)), list(map(scalar, w))
        abc = [ip(W, W), scale(ip(U, W), 2), add(ip(U, U), scalar(-2))]
        assert abc == list(map(scalar, [F(54, 242), F(-214*sg, 242), F(35, 242)]))
        assert abc == list(map(scalar, map(F, r['norm_polynomial'])))
        J = r['support_prior_indices']+[r['anchor']]
        assert len(J) == len(set(J)) == 20
        rows = [cov(old[j]) for j in J]
        assert all(x[1:] == (0, 0, 0) for row in rows for x in row)
        assert rank([[x[0] for x in row] for row in rows]) == 20
        for j in J:
            assert ip(v, old[j]) == scalar(-int(j in r['profile'] or j == r['anchor']))
        flipped = v[:]
        flipped[19] = neg(flipped[19])
        assert flipped == N[72+i]

    G = [[ZERO]*76 for _ in N]
    for i in range(76):
        for j in range(i+1): G[i][j] = G[j][i] = ip(N[i], N[j])
    norms = list(map(F, d['norms']))
    assert all(G[i][i] == scalar(norms[i]) and norms[i] > 0 for i in range(76))
    assert norms[68:] == [2]*8
    counts = Counter()
    for i in range(68, 76):
        for j in range(i):
            g = G[i][j]
            assert sign(g) <= 0
            if g == ZERO: counts['right'] += 1
            elif mul(g, g) == scalar(norms[i]*norms[j]/4): counts['pi3'] += 1
            else:
                assert sign(add(mul(g, g), scalar(-norms[i]*norms[j]))) > 0
                counts['dotted'] += 1
    assert dict(counts) == d['new_pair_counts'] == {'right': 144, 'pi3': 16, 'dotted': 412}
    assert all(sign(G[i][j]) <= 0 for i, j in combinations(range(76), 2))
    H = list(map(scalar, map(F, d['strict_center'])))
    assert sign(ip(H, H)) < 0 and all(sign(ip(H, r)) < 0 for r in N)
    # This strict timelike center and the acute Gram give all actual facets
    # by orthogonal projection; no stored facet-point assertion is needed.

    def gram(J):
        assert all(G[i][j][1:] == (0, 0, 0) for i in J for j in J)
        return [[G[i][j][0] for j in J] for i in J]

    vertices = d['finite_vertices']
    assert len(vertices) == 16
    for v in vertices:
        J = v['incident']; p = vec(v['point'])
        assert len(J) == len(set(J)) == 20 and pd(gram(J))
        assert ip(p, p) == tuple(map(F, v['norm'])) and sign(ip(p, p)) < 0
        assert sign(ip(p, H)) < 0
        for j, r in enumerate(N): assert ip(p, r) == ZERO if j in J else sign(ip(p, r)) < 0
    assert len({tuple(sorted(v['incident'])) for v in vertices}) == 16
    clique_tests = 0
    for row in d['new_cap_vertex_census']:
        k = row['new_wall']
        finite = [j for j in range(76) if j != k and sign(add(mul(G[k][j], G[k][j]), scalar(-norms[k]*norms[j]))) < 0]
        assert finite == row['finite_neighbors'] and len(finite) == 20
        accepted = []
        for J0 in combinations(finite, 19):
            clique_tests += 1
            J = list(J0)+[k]
            if any(sign(add(mul(G[i][j], G[i][j]), scalar(-norms[i]*norms[j]))) > 0 for i, j in combinations(J, 2)): continue
            if pd(gram(J)): accepted.append(set(J))
        assert len(accepted) == 1 and accepted[0] == set(vertices[row['actual_vertex']]['incident'])
    assert len(d['new_cap_vertex_census']) == 8 and clique_tests == 160
    edges = []
    for i, j in combinations(range(16), 2):
        S = sorted(set(vertices[i]['incident']) & set(vertices[j]['incident']))
        if len(S) == 19:
            assert pd(gram(S)); edges.append(dict(vertices=[i, j], support=S))
    assert edges == d['finite_edges'] and len(edges) == 8

    E = list(map(scalar, map(F, d['strict_escape']['point'])))
    assert ip(E, E) == scalar(F(269344645863312077, 40000000000000000))
    assert ip(E, E) == tuple(map(F, d['strict_escape']['norm']))
    assert sign(ip(E, H)) < 0 and all(sign(ip(E, r)) < 0 for r in N)
    oldE = [(F(x['a']), F(x['b']), F(0), F(0)) for x in base['strict_escape']['point']]
    for row in d['previous_escape_sign_images_cut']:
        e = oldE[:]
        for j, s in zip((19, 20), row['signs']): e[j] = scale(e[j], s)
        assert [i for i in range(68, 76) if sign(ip(N[i], e)) > 0] == row['owners']
        assert row['owners']
    sourcepath = 'results/fresh_dimension_construction_routes_20261014.json'
    source = load(sourcepath)
    assert len(source['repaired_cusps']) == len(d['preserved_three_cusp_links']) == 3
    for r, saved in zip(source['repaired_cusps'], d['preserved_three_cusp_links']):
        delta = [(F(x['a']), F(x['b']), F(0), F(0)) for x in r['delta']]
        assert ip(delta, delta) == ZERO
        vals = [ip(n, delta) for n in N[68:]]
        assert vals == list(map(tuple, (map(F, a) for a in saved['new_pairings'])))
        assert all(sign(a) < 0 for a in vals)

    paths = [f'scripts/{TAG}.py', f'results/{TAG}.json', basepath, sourcepath, f'scripts/{TAG}_review.py']
    return dict(status='PASS_SEPARATE_NESTED_QUADRATIC_REPLACEMENT_GEOMETRY_REPLAY',
                sources={p: digest(p) for p in paths},
                new_recipe_ranks=[20]*4, new_pair_counts=dict(counts),
                vertices_checked=16, vertex_inequalities_checked=16*76,
                per_new_cap_subsets_checked=160, finite_edges_checked=8,
                strict_center_and_all_facets=True, preserved_cusp_rays=3,
                escape_norm='269344645863312077/40000000000000000',
                scope='Separate Fraction nested-quadratic arithmetic with exact squaring signs and no producer imports. Reconstructs all8 new normals,572 new pairs,16 specified vertices, each new-cap finite-vertex census,8 edges, old escape cuts and new strict escape. Old68 Coxeter validity, completeness of three inherited cusp links and conditional nonquasiarithmetic remain dependencies. No whole chamber vertex or ideal census.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    d = run(); path = ROOT/f'results/{TAG}_review.json'
    if args.verify: assert d == json.loads(path.read_text())
    else: path.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['new_pair_counts'], d['finite_edges_checked'])
