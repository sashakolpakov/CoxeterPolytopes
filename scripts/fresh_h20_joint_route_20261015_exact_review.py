#!/usr/bin/env python3
"""Independent Fraction biquadratic replay of the43-wall cap and its finite edge.

No producer, helper, symbolic algebra, floating point or project imports.
Signs use rational square-root isolating intervals; inverses use Galois norms.
Finite20-supports use full matrix LDL, not the producer's forest algorithm.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_joint_route_20261015_exact_review'
PRIMARY = 'results/fresh_h20_joint_route_20261015_exact.json'
SOURCE = 'results/fresh_h20_a1a1_h18_coupling_20261015.json'
ZERO = (F(0),) * 4
ONE = (F(1), F(0), F(0), F(0))


def dec(x):
    if isinstance(x, (int, str, F)):
        return (F(x), F(0), F(0), F(0))
    if isinstance(x, dict):
        assert set(x) == {'a', 'b'}
        return (F(x['a']), F(x['b']), F(0), F(0))
    assert len(x) in (2, 4)
    return tuple(map(F, x)) + (F(0),) * (4-len(x))


def add(x, y):
    return tuple(a+b for a, b in zip(x, y))


def neg(x):
    return tuple(-a for a in x)


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    z = [F(0)] * 4
    for i, a in enumerate(x):
        if not a:
            continue
        for j, b in enumerate(y):
            if b:
                z[i ^ j] += a*b*(2 if i & j & 1 else 1)*(3 if i & j & 2 else 1)
    return tuple(z)


@lru_cache(None)
def inv(x):
    assert x != ZERO
    c = ONE
    for sigma in (1, 2, 3):
        c = mul(c, tuple(a*(-1 if (i & sigma).bit_count() % 2 else 1) for i, a in enumerate(x)))
    n = mul(x, c)
    assert n[1:] == (0, 0, 0) and n[0]
    return tuple(a/n[0] for a in c)


@lru_cache(None)
def radical_interval(d, bits):
    scale = 1 << bits
    q = isqrt(d * scale * scale)
    assert q*q <= d*scale*scale < (q+1)*(q+1)
    return F(q, scale), F(q+1, scale)


SIGN_PRECISIONS = Counter()


@lru_cache(None)
def sign(x):
    if x == ZERO:
        return 0
    for bits in (64, 128, 256, 512):
        lo = hi = x[0]
        for a, d in zip(x[1:], (2, 3, 6)):
            l, u = radical_interval(d, bits)
            lo += a * (l if a >= 0 else u)
            hi += a * (u if a >= 0 else l)
        if lo > 0 or hi < 0:
            SIGN_PRECISIONS[bits] += 1
            return 1 if lo > 0 else -1
    raise AssertionError(('unresolved exact sign', x))


def dot(x, y):
    z = ZERO
    for a, b in zip(x, y):
        z = add(z, mul(a, b))
    return z


def mv(M, x):
    return [dot(row, x) for row in M]


def inertia(M, require_positive=False):
    """Symmetric congruence; a zero diagonal is repaired by a literal basis sum."""
    a = [row[:] for row in M]
    counts = Counter()
    while a:
        n = len(a)
        j = next((i for i in range(n) if a[i][i] != ZERO), None)
        if j is None:
            pair = next(((i, j) for i, j in combinations(range(n), 2) if a[i][j] != ZERO), None)
            if pair is None:
                counts[0] += n
                break
            i, j = pair
            # Replace basis vector i by e_i+e_j; this is invertible.
            for k in range(n):
                a[i][k] = add(a[i][k], a[j][k])
            for k in range(n):
                a[k][i] = add(a[k][i], a[k][j])
            j = i
        if j:
            a[0], a[j] = a[j], a[0]
            for row in a:
                row[0], row[j] = row[j], row[0]
        s = sign(a[0][0])
        if require_positive and s <= 0:
            return False
        counts[s] += 1
        q = inv(a[0][0])
        a = [[sub(a[i][j], mul(mul(a[i][0], q), a[0][j])) for j in range(1, n)] for i in range(1, n)]
    if require_positive:
        return counts[0] == 0 and counts[-1] == 0
    return counts[1], counts[-1], counts[0]


def rank(M):
    a = [row[:] for row in M]
    k = 0
    for j in range(len(a[0])):
        i = next((i for i in range(k, len(a)) if a[i][j] != ZERO), None)
        if i is None:
            continue
        a[k], a[i] = a[i], a[k]
        q = inv(a[k][j])
        a[k] = [mul(v, q) for v in a[k]]
        for i in range(k+1, len(a)):
            q = a[i][j]
            if q != ZERO:
                a[i] = [sub(v, mul(q, w)) for v, w in zip(a[i], a[k])]
        k += 1
        if k == len(a):
            break
    return k


def run():
    def sha(p):
        return sha256((ROOT/p).read_bytes()).hexdigest()
    seen = set()
    def read(p):
        d = json.loads((ROOT/p).read_text())
        if p not in seen:
            seen.add(p)
            for q, h in d.get('sources', {}).items():
                assert sha(q) == h, q
                if q.endswith('.json'):
                    read(q)
        return d
    d, old = map(read, (PRIMARY, SOURCE))
    assert d['status'] == 'PASS_EXACT_NONINDUCING_H20_CAP_AND_FINITE_EDGE_WITH_RENEWED_ESCAPE'
    assert d['field_basis'] == ['1', 'sqrt2', 'sqrt3', 'sqrt6']
    Q = [[dec(a) for a in row] for row in d['form']]
    N = [[dec(a) for a in row] for row in d['normals']]
    Q0 = [[dec(a) for a in row] for row in old['source_form']]
    N0 = [[dec(a) for a in row] for row in old['source_normals']]
    assert len(Q) == 21 and len(Q0) == 19 and len(N0) == 40 and len(N) == 43
    assert all(Q[i][j] == Q[j][i] for i in range(21) for j in range(21))
    expected = [row+[ZERO, ZERO] for row in Q0] + [[ZERO]*19+[dec(2), ZERO], [ZERO]*20+[dec(2)]]
    assert Q == expected and N[:40] == [row+[ZERO, ZERO] for row in N0]
    assert N[40] == [ZERO]*19+[ONE, ZERO] and N[41] == [ZERO]*20+[ONE]
    assert inertia(Q) == (20, 1, 0) and inertia(Q0) == (18, 1, 0)
    cov = [mv(Q, n) for n in N]
    G = [[dot(a, c) for c in cov] for a in N]
    assert G == [[dec(a) for a in row] for row in d['gram']]
    norms = [G[i][i] for i in range(43)]
    assert all(sign(a) > 0 for a in norms) and norms[42] == dec(2)
    recipe = d['cap_recipe']
    rights = [1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,21,36]
    assert recipe['source_rights'] == rights and recipe['source_anchors'] == [12, 18]
    Y = list(map(dec, recipe['source_projection']))
    C0 = list(map(dec, old['source_center']))
    old_cov = [mv(Q0, n) for n in N0]
    assert rank([old_cov[i] for i in rights+[12, 18]]) == 19
    assert [dot(Y, old_cov[i]) for i in rights+[12, 18]] == [ZERO]*17+[dec(-1), dec(-1)]
    assert dot(Y, mv(Q0, Y)) == dec('-3/2') and sign(dot(Y, mv(Q0, C0))) < 0
    assert N[42] == Y+[(F(0), F(0), F(-1,2), F(0)), dec(-1)]
    labels = {}
    totals, new = Counter(), Counter()
    for i, j in combinations(range(43), 2):
        g = G[i][j]
        assert sign(g) <= 0
        t = mul(dec(4), mul(mul(g, g), inv(mul(norms[i], norms[j]))))
        if g == ZERO:
            label = 'right'
        elif t in (dec(1), dec(2), dec(3), dec(4)):
            label = {1:'pi3', 2:'pi4', 3:'pi6', 4:'parallel'}[int(t[0])]
        else:
            assert sign(sub(t, dec(4))) > 0
            label = 'dotted'
        labels[i,j] = label
        totals[label] += 1
        if j == 42:
            new[label] += 1
    assert dict(totals) == d['all_pair_counts'] and dict(new) == d['new_pair_counts']
    C = list(map(dec, d['strict_center']))
    cp = [dot(C, c) for c in cov]
    assert cp == list(map(dec, d['center_pairings']))
    assert dot(C, mv(Q, C)) == dec(d['center_norm']) and sign(dot(C, mv(Q, C))) < 0
    assert all(sign(a) < 0 for a in cp)
    for i, point in enumerate(d['facet_points']):
        v = list(map(dec, point))
        assert sign(dot(v, mv(Q, v))) < 0 and sign(dot(v, mv(Q, C))) < 0
        assert [sign(dot(v, c)) for c in cov] == [-1]*i+[0]+[-1]*(42-i)
    finite = {'right', 'pi3', 'pi4', 'pi6'}
    neighbors = [i for i in range(42) if labels[i,42] in finite]
    assert neighbors == d['complete_new_cap_vertex_census']['finite_neighbors']
    accepted, finite_pair_cliques = [], 0
    for subset in combinations(neighbors, 19):
        ids = list(subset)+[42]
        if not all(labels[min(i,j),max(i,j)] in finite for i,j in combinations(ids,2)):
            continue
        finite_pair_cliques += 1
        if inertia([[G[i][j] for j in ids] for i in ids], require_positive=True):
            accepted.append(ids)
    assert len(list(combinations(neighbors,19))) == 1540
    assert accepted == d['complete_new_cap_vertex_census']['positive_supports'] and len(accepted) == 1
    endpoints = d['finite_edge']['endpoints']
    assert len(endpoints) == 2
    for row in endpoints:
        v = list(map(dec, row['point']));S = row['support']
        assert len(S) == 20 and inertia([[G[i][j] for j in S] for i in S], True)
        assert dot(v,mv(Q,C)) == dec(-1)
        assert dot(v,mv(Q,v)) == dec(row['norm']) and sign(dec(row['norm'])) < 0
        p = [dot(v,c) for c in cov]
        assert p == list(map(dec,row['pairings']))
        assert [sign(a) for a in p] == [0 if i in S else -1 for i in range(43)]
    edge = sorted(set(endpoints[0]['support']) & set(endpoints[1]['support']))
    assert edge == d['finite_edge']['shared_support'] and len(edge) == 19
    assert inertia([[G[i][j] for j in edge] for i in edge], True)
    assert set(endpoints[0]['support'])-set(edge) == {41}
    assert set(endpoints[1]['support'])-set(edge) == {42}
    Eold = list(map(dec,old['strict_escape']))
    assert dot(Eold,cov[42]) == dec(d['old_escape_cut']) and sign(dec(d['old_escape_cut'])) > 0
    E = list(map(dec,d['strict_escape']))
    assert E == Eold[:19]+[dec('-3/4'),dec('-1/1000')]
    ep = [dot(E,c) for c in cov]
    assert ep == list(map(dec,d['escape_pairings'])) and all(sign(a)<0 for a in ep)
    assert dot(E,mv(Q,E)) == dec(d['escape_norm']) and sign(dec(d['escape_norm'])) > 0
    assert dot(E,mv(Q,C)) == dec(d['escape_center_pairing']) and sign(dec(d['escape_center_pairing'])) < 0
    face = list(map(dec,d['full_H18_source_face']['point']))
    assert face == C0+[ZERO,ZERO]
    assert [sign(dot(face,c)) for c in cov] == [-1]*40+[0,0,-1]
    cycle = mul(dec(4),mul(mul(G[0][36],G[0][36]),inv(mul(norms[0],norms[36]))))
    minor = sub(mul(norms[38],norms[39]),mul(G[38][39],G[38][39]))
    assert cycle == dec([12,8]) == dec(d['conditional_nonquasi']['cartan_two_cycle'])
    assert minor == dec(-32) == dec(d['conditional_nonquasi']['determinant'])
    paths = [PRIMARY,SOURCE,f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_EXACT_H20_CAP_EDGE_ESCAPE_REVIEW',
                sources={p:sha(p) for p in paths}, checked_source_manifests=len(seen),
                arithmetic='Standalone Fraction basis1,sqrt2,sqrt3,sqrt6; Galois-norm inversion, rational square-root sign intervals, full matrix LDL. No producer/helper imports.',
                form_signature=[20,1,0], facets=43, pair_count=903,
                all_pair_counts=dict(totals), new_pair_counts=dict(new),
                actual_facet_witnesses=43, cap_recipe_unique_equations=19,
                candidate_new_vertex_supports=1540, finite_pair_cliques=finite_pair_cliques,
                positive_new_vertex_supports=accepted, exact_finite_edge_endpoints=2,
                strict_escape=True, old_escape_cut=True,
                field_sign_precision_histogram={str(k):v for k,v in SIGN_PRECISIONS.items()},
                scope='Complete new-cap finite-vertex census, one exact finite edge, all43 actual facets and903 pairs, unique cap recipe, and renewed strict escape independently replayed. No full ambient vertex/edge census. The conditional non-quasi-arithmetic numbers are checked; the external Vinberg criterion and classification are inherited theorem dependencies, not independently established here.')


if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a = p.parse_args()
    data = run();out = ROOT/f'results/{TAG}.json'
    if a.verify:
        assert data == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['finite_pair_cliques'],data['field_sign_precision_histogram'])
