#!/usr/bin/env python3
"""Separate SymPy verification of the saved single-factor geometric data.

No producer/helper arithmetic is imported. Continuous C2/A2 coverage is
reviewed in prose; this checks every saved source face and sum-plane ridge.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_single_factor_links_20261013'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'


def sha(p):
    return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    packet = json.loads((ROOT/'results'/f'{TAG}.json').read_text())
    assert all(sha(p) == h for p, h in packet['sources'].items())
    raw = json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q = -s.Matrix(raw['basis_gram']).applyfunc(s.Rational)
    W = s.Matrix(raw['wall_coordinates_as_columns_transposed']).applyfunc(s.Rational)
    C = W*Q
    G = C*W.T
    def vec(x):
        return s.Matrix(list(map(s.Rational, x)))
    def pair(a, b):
        return (a.T*Q*b)[0]
    faces = {}
    def traverse(x):
        if isinstance(x, dict):
            if {'walls', 'point', 'norm', 'future_pair', 'LDL'} <= set(x):
                key = tuple(x['walls'])
                if key in faces:
                    assert faces[key] == x
                else:
                    p = vec(x['point'])
                    values = C*p
                    assert pair(p, p) == s.Rational(x['norm']) < 0
                    assert all(v <= 0 for v in values)
                    assert [i for i, v in enumerate(values) if v == 0] == list(key)
                    M = G.extract(key, key)
                    assert all(M[:i, :i].det() > 0 for i in range(1, len(key)+1))
                    assert max(v for v in values if v) == s.Rational(x['maximum_other_pairing'])
                    faces[key] = x
            for value in x.values():
                traverse(value)
        elif isinstance(x, list):
            for value in x:
                traverse(value)
    traverse(packet)
    t = s.Symbol('t', real=True)
    for row in packet['sum_planes']:
        a, D = vec(row['normal_constant']), vec(row['normal_shift'])
        k = len(row['selected_leaves'])
        assert pair(a, a) == 2*k and pair(D, D) == pair(a, D) == 0
        e = a+t*D
        X, Y = vec(row['ridge_point_at_minus_one']), vec(row['ridge_point_at_zero'])
        assert pair(X, X) < 0 and pair(Y, Y) < 0 and pair(X, Y) < 0
        point = -t*X+(1+t)*Y
        assert s.expand(pair(point, e)) == 0
        i, j = row['ridge_walls']
        assert set(n for n, v in enumerate(C*X) if not v) & set(n for n, v in enumerate(C*Y) if not v) == {i, j}
        assert all(v <= 0 for v in C*X) and all(v <= 0 for v in C*Y)
        ep = C*e
        p = s.expand(G[i, j]-ep[i]*ep[j]/(2*k))
        ni = s.expand(G[i, i]-ep[i]**2/(2*k))
        nj = s.expand(G[j, j]-ep[j]**2/(2*k))
        assert s.simplify(p-(1+t)/k) == 0
        assert ni == 2-2/s.Integer(k)
        assert s.simplify(nj-(2-(1+t)**2/(2*k))) == 0
        assert s.simplify(ni*nj-p*p-(4-4/s.Integer(k)-(1+t)**2/k)) == 0
        assert 4-5/s.Integer(k) > 0
        # This proves positive projected norm/determinant and an obtuse
        # actual ridge for every -1<t<0, using convexity of the future cone.
    for row in packet['full_C2_spans']:
        a, b, c = row['walls']
        j, k = row['auxiliary_pair']
        assert G.extract([a, b, c], [a, b, c]) == s.Matrix([[2, 0, -2], [0, 2, -2], [-2, -2, 4]])
        assert [G[j, i] for i in [a, b, c, k]] == [-1, 0, 0, 0]
        assert [G[k, i] for i in [a, b, c]] == [0, -1, 0]
        D = vec(row['null_sum'])
        for r in row['triangle_face_points']:
            p = vec(r['normalized_point'])
            assert pair(p, D) == -1 and pair(p, W[r['only_nonzero_factor_wall'], :].T) == -1
        for r in row['same_sign_faces']:
            X, Y = vec(r['normalized_left']), vec(r['normalized_right'])
            assert pair(X, D) == pair(Y, D) == -1
            assert pair(X, W[r['large_coefficient_wall'], :].T) == -1
            assert pair(Y, W[c, :].T) == -1
    for row in packet['full_A2_spans']:
        ids = row['factor']
        assert G.extract(ids, ids) == s.Matrix([[4 if i == j else -2 for j in range(3)] for i in range(3)])
        a, b = row['ridge_walls']
        negative = row['negative_wall']
        assert G[a, b] == -1 and all(G[a, i] == 0 for i in ids)
        assert G[b, negative] == -2 and all(G[b, i] == 0 for i in row['positive_walls'])
    assert len(faces) == packet['counts']['exact_ambient_faces']
    sources = [RAW, f'results/{TAG}.json', f'scripts/{TAG}.py', f'scripts/{TAG}_review.py']
    return {'status': 'PASS_INDEPENDENT_SYMPY_SINGLE_FACTOR_GEOMETRIC_DATA',
            'sources': {p: sha(p) for p in sources},
            'exact_source_faces': len(faces), 'all_real_sum_plane_ridges': len(packet['sum_planes']),
            'C2_full_span_auxiliary_configurations': len(packet['full_C2_spans']),
            'A2_full_span_orientation_configurations': len(packet['full_A2_spans']),
            'scope': 'All saved ambient face points checked against all210 walls by separate SymPy arithmetic; all79 real-parameter sum-plane ridge identities verified symbolically; all rank-two auxiliary configurations checked. Full parameter coverage and rank-two inequalities additionally require the mathematical report review.'}


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    p = ROOT/'results'/f'{TAG}_review.json'
    if args.verify:
        assert json.loads(p.read_text()) == out
    else:
        p.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['exact_source_faces'])
