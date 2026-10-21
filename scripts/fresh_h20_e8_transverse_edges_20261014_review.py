#!/usr/bin/env python3
"""Independent interval replay of the 280 transverse-family certificates.

Imports only the earlier independent sparse Fraction arithmetic checker,
not the producer or its nested-quadratic sign implementation.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import fresh_h20_e8_cusp_completion_20261014_review as K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e8_transverse_edges_20261014_review'
PRIMARY = 'results/fresh_h20_e8_transverse_edges_20261014.json'
BASE = 'results/fresh_h20_e8_cusp_completion_20261014.json'
FACE = 'results/fresh_h20_contact_closure_20261014_h12_replacement.json'


def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def vec(row): return [tuple(map(F, x)) for x in row]
def sub(a, b): return K.add(a, K.scale(b, -1))


def interval(x):
    return (sum(a*(b[0] if a >= 0 else b[1]) for a, b in zip(x, K.BOUNDS)),
            sum(a*(b[1] if a >= 0 else b[0]) for a, b in zip(x, K.BOUNDS)))


def sqrt_interval(h):
    assert h >= 0
    n = isqrt(h.numerator*K.DEN*K.DEN//h.denominator)
    lo = F(n, K.DEN)
    return lo, lo if lo*lo == h else F(n+1, K.DEN)


def radical_interval(a, b, h):
    al, ah = interval(a); bl, bh = interval(b); hl, hh = sqrt_interval(h)
    p = [bl*hl, bl*hh, bh*hl, bh*hh]
    return al+min(p), ah+max(p)


def run():
    doc, d, f = [json.loads((ROOT/p).read_text()) for p in [PRIMARY, BASE, FACE]]
    for data in [doc, d, f]:
        for p, h in data['sources'].items(): assert sha(p) == h, p
    Q = [[F(x) for x in row] for row in d['form']]; cache = {}

    def dot(x, y):
        key = tuple(y)
        if key not in cache:
            cache[key] = [tuple(sum(Q[i][j]*y[j][k] for j in range(21)) for k in range(4)) for i in range(21)]
        g = K.Z
        for a, b in zip(x, cache[key]): g = K.add(g, K.mul(a, b))
        return g

    N = list(map(vec, d['normals'])); E = vec(d['strict_escape']['point'])
    V = [vec(row['point']) for row in f['H12_face_finite_vertices']]
    J = d['E8_support']; B = [[dot(N[i], N[j]) for j in J] for i in J]
    assert all(x[1:] == (0, 0, 0) for row in B for x in row)
    I = K.matrix_inverse([[x[0] for x in row] for row in B])
    weights = [[tuple(sum(I[i][j]*N[J[j]][a][k] for j in range(8)) for k in range(4)) for a in range(21)] for i in range(8)]
    assert K.sign(dot(E, E)) > 0
    counts = Counter(); seen = set()
    for row in doc['records']:
        vi, j = row['vertex'], row['E8_position']; assert (vi, j) not in seen; seen.add((vi, j))
        w, v = weights[j], V[vi]
        degree = I[j][j]; q = K.scale(dot(v, v), -1)
        assert q[1:] == (0, 0, 0) and q[0] > 0
        q = q[0]; A = K.scale(dot(w, E), -1); beta = K.scale(dot(v, E), -1)
        assert K.sign(A) > 0 and K.sign(beta) > 0 and dot(w, v) == K.Z
        assert tuple(map(F, row['escape_A'])) == A and tuple(map(F, row['escape_beta'])) == beta
        assert F(row['weight_norm']) == degree and F(row['vertex_negative_norm']) == q
        h1 = (degree-2)/q; reason = row['reason']; counts[reason] += 1
        if reason == 'minimum_coefficient_escape_bound':
            assert radical_interval(A, K.scale(beta, -1), h1)[1] < 0
            continue
        i = row['blocking_wall']; a = K.scale(dot(w, N[i]), -1); b = K.scale(dot(v, N[i]), -1)
        assert tuple(map(F, row['wall_a'])) == a and tuple(map(F, row['wall_b'])) == b
        if reason == 'acute_wall_escape_bound':
            assert K.sign(a) > 0 and K.sign(b) >= 0
            if b != K.Z:
                gap = sub(K.mul(a, beta), K.mul(A, b))
                assert tuple(map(F, row['cross_product_gap'])) == gap and K.sign(gap) >= 0
            continue
        assert reason == 'Coxeter_coefficient_gap'
        assert radical_interval(A, K.scale(beta, -1), (degree-1)/q)[1] < 0
        lo, hi = radical_interval(a, K.scale(b, -1), h1)
        if row['c1_reason'] == 'nonacute': assert lo > 0
        else:
            assert row['c1_reason'] == 'forbidden_angle_between_pi3_and_pi2'
            assert hi < 0 and lo*lo < F(d['norms'][i])/2
    assert seen == {(i, j) for i in range(35) for j in range(8)}
    assert dict(counts) == doc['counts']
    paths = [PRIMARY, BASE, FACE, f'scripts/{TAG}.py', 'scripts/fresh_h20_e8_cusp_completion_20261014_review.py']
    return dict(status='PASS_INDEPENDENT_INTERVAL_TRANSVERSE_FAMILY_REVIEW', sources={p: sha(p) for p in paths}, records=280, counts=dict(counts), scope='Reconstructs E8 inverse, weights, vertex norms, all certificate coefficients and inequalities with independent sparse field multiplication and rational radical intervals. Inherited79 chamber and35 vertex completeness remain certified dependencies. This does not review unrestricted caps or finite-volume completion.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    d = run(); p = ROOT/f'results/{TAG}.json'
    if args.verify: assert d == json.loads(p.read_text())
    else: p.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], d['counts'])
