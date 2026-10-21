#!/usr/bin/env python3
"""Separate sparse-Fraction/interval replay of eleven all-y dual bounds."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e8_cusp_completion_20261014_review as K

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_e8_completion_gate_20261014_review'
PRIMARY = 'results/fresh_h20_e8_completion_gate_20261014.json'
BASE = 'results/fresh_h20_e8_cusp_completion_20261014.json'
FACE = 'results/fresh_h20_contact_closure_20261014_h12_replacement.json'


def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def vec(r): return [tuple(map(F, x)) for x in r]
def fld(x): return tuple(map(F, x))
def sub(a, b): return K.add(a, K.scale(b, -1))


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

    def linear(vs, cs):
        out = [K.Z]*21
        for v, c in zip(vs, cs): out = [K.add(a, K.mul(b, c)) for a, b in zip(out, v)]
        return out

    N = list(map(vec, d['normals'])); E = vec(d['strict_escape']['point']); J = d['E8_support']
    B = [[dot(N[i], N[j]) for j in J] for i in J]
    assert all(x[1:] == (0, 0, 0) for row in B for x in row)
    B = [[x[0] for x in row] for row in B]; K.pd(B); I = K.matrix_inverse(B)
    W = [linear([N[j] for j in J], list(map(K.scalar, row))) for row in I]
    assert W == list(map(vec, doc['E8_fundamental_weights']))

    def project(x): return [sub(a, b) for a, b in zip(x, linear(W, [dot(x, N[j]) for j in J]))]

    HP = project(vec(d['strict_center'])); assert K.sign(dot(HP, HP)) < 0
    EH = project(E); assert K.sign(dot(EH, EH)) < 0 and K.sign(dot(EH, HP)) < 0
    domains = []
    for row in doc['dual_certificates']:
        X = E[:]
        for i, value in row['nonnegative_wall_multipliers']:
            lam = F(value); assert lam > 0 and 0 <= i < 79
            X = [sub(a, K.scale(b, lam)) for a, b in zip(X, N[i])]
        assert X == vec(row['adjusted_escape'])
        Y = project(X); assert Y == vec(row['future_timelike_projection'])
        kappa = K.scale(dot(Y, Y), -1); assert kappa == fld(row['kappa'])
        assert K.sign(kappa) > 0 and K.sign(dot(Y, HP)) < 0
        a = [K.scale(dot(w, X), -1) for w in W]; assert a == list(map(fld, row['linear_coefficients']))
        M = [[sub(K.scale(kappa, I[i][j]), K.mul(a[i], a[j])) for j in range(8)] for i in range(8)]
        assert M == [list(map(fld, rr)) for rr in row['quadratic_matrix']]
        S = row['support']; kind = row['kind']
        if kind == 'single_negative_linear':
            assert len(S) == 1 and K.sign(a[S[0]]) < 0
            margin = K.scale(a[S[0]], -1)
        elif kind == 'single_threshold':
            assert len(S) == 1; i = S[0]; lower = F(row['minimum_squared_magnitude'])
            assert lower in [1, 2] and K.sign(M[i][i]) > 0
            margin = sub(K.scale(M[i][i], lower), K.scale(kappa, 2))
        else:
            assert len(S) == 2
            if kind == 'independent_magnitudes_at_least_one':
                assert all(K.sign(M[i][j]) > 0 for i in S for j in S)
            else: assert kind == 'equal_magnitudes_at_least_one' and S == [0, 7]
            total = K.Z
            for i in S:
                for j in S: total = K.add(total, M[i][j])
            assert K.sign(total) > 0
            margin = sub(total, K.scale(kappa, 2))
        assert margin == fld(row['strict_algebraic_margin']) and K.sign(margin) > 0
        domains.append(dict(kind=kind, support=S, minimum_squared_magnitude=row['minimum_squared_magnitude']))
    assert len(domains) == 11
    w0 = [K.scale(x, -1) for x in W[0]]; C = dot(w0, E)
    assert dot(w0, w0) == K.scalar(2) and K.sign(C) > 0
    for k, row in enumerate(doc['minimal_weight_pi3_null_ray_bounds']):
        nu = vec(row['point']); assert nu == vec(d['cusp_rays'][k]) == vec(f['H12_face_ideal_vertices'][k]['point'])
        assert dot(nu, nu) == K.Z and K.sign(dot(nu, HP)) < 0
        slope = dot(nu, E); assert slope == fld(row['escape_slope']) and K.sign(slope) < 0
        rec = row['witness']; n = N[rec['wall']]
        denominator = K.scale(dot(nu, n), -1); offset = dot(w0, n)
        assert denominator == fld(rec['positive_denominator']) and offset == fld(rec['positive_offset'])
        assert K.sign(denominator) > 0 and K.sign(offset) > 0
        upper = K.add(K.mul(denominator, C), K.mul(offset, slope))
        assert upper == fld(rec['strict_upper_numerator']) and K.sign(upper) < 0
    assert len(doc['minimal_weight_pi3_null_ray_bounds']) == 2
    paths = [PRIMARY, BASE, FACE, f'scripts/{TAG}.py', 'scripts/fresh_h20_e8_cusp_completion_20261014_review.py']
    return dict(status='PASS_SEPARATE_INTERVAL_ALL_Y_DUAL_BOUND_REVIEW', sources={p: sha(p) for p in paths}, dual_domains=domains, null_ray_bounds=2, scope='Reconstructs all11 nonnegative dual recipes and all field inequalities using sparse radical multiplication and160bit rational intervals; no optimizer or producer imports. Source-cone completeness and embedding remain certified inherited dependencies, not a new enumeration. The bounds concern exactly the stated contact domains.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    d = run(); p = ROOT/f'results/{TAG}.json'
    if args.verify: assert d == json.loads(p.read_text())
    else: p.write_text(json.dumps(d, indent=2)+'\n')
    print(d['status'], len(d['dual_domains']))
