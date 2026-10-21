#!/usr/bin/env python3
"""Independent physical identity/sign/rank replay for the S75 escape dual.

No producer or project arithmetic imports. Squarefree-radical multiplication,
rational interval signs and one exact finite-field rank witness suffice.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_s75_completion_20261014_gate_review'
SOURCE = 'results/fresh_h20_s75_completion_20261014_gate.json'
PRIMARY = 'scripts/fresh_h20_s75_completion_20261014_gate.py'
BASE = 'results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
CURRENT = 'results/fresh_h20_bridge_edge_closure_20261014.json'
PROFILE = 'results/fresh_h20_s75_projection_gate_20261014.json'
RAD = [1, 2, 79, 158, 11, 22, 869, 1738, 3, 6, 237, 474, 33, 66, 2607, 5214]
INDEX = {r: i for i, r in enumerate(RAD)}
Z = (F(0),) * 16

def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()

def mul(a, b):
    out = [F(0)] * 16
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    g = gcd(RAD[i], RAD[j])
                    out[INDEX[RAD[i] * RAD[j] // (g*g)]] += x*y*g
    return tuple(out)

def interval(a, bits=192):
    den = 1 << bits
    lo = hi = F(0)
    for r, q in zip(RAD, a):
        n = isqrt(r*den*den)
        l, h = F(n, den), F(n + (n*n != r*den*den), den)
        lo += q * (l if q >= 0 else h)
        hi += q * (h if q >= 0 else l)
    return lo, hi

def mod_rank(matrix, p):
    a = [row[:] for row in matrix]
    row = 0
    pivots = []
    for j in range(len(a[0])):
        k = next((k for k in range(row, len(a)) if a[k][j] % p), None)
        if k is None:
            continue
        a[row], a[k] = a[k], a[row]
        inv = pow(a[row][j], -1, p)
        a[row] = [(x*inv) % p for x in a[row]]
        for k in range(row+1, len(a)):
            q = a[k][j]
            a[k] = [(x-q*y) % p for x, y in zip(a[k], a[row])]
        pivots.append(j)
        row += 1
    return row, pivots

def encode(v):
    if isinstance(v, F):
        return str(v)
    if isinstance(v, (tuple, list)):
        return [encode(x) for x in v]
    if isinstance(v, dict):
        return {str(k): encode(x) for k, x in v.items()}
    return v

def run():
    d, old, cur, profile = [json.loads((ROOT/p).read_text()) for p in [SOURCE, BASE, CURRENT, PROFILE]]
    for doc in [d, old, cur, profile]:
        for path, value in doc['sources'].items():
            assert sha(path) == value, path
    assert d['field_square_root_basis'] == old['field']['base_square_root_basis'] == RAD
    rows = d['identity_coefficients']
    ids = [r['wall'] for r in rows]
    assert len(ids) == len(set(ids)) == 21
    coefficients = [tuple(map(F, r['coefficient'])) for r in rows]
    vectors = []
    for i in ids:
        v = []
        for a, b in old['normals'][i]:
            assert all(F(x) == 0 for x in b)
            v.append(tuple(map(F, a)))
        vectors.append(v)
    E = []
    for ab, cd in cur['strict_escape']['point']:
        a, b = ab
        c, e = cd
        assert all(F(x) == 0 for x in b+c+e+a[1:])
        E.append(F(a[0]))
    assert E == list(map(F, d['escape_point']))
    for k in range(21):
        value = [F(0)] * 16
        for c, v in zip(coefficients, vectors):
            value = [x+y for x, y in zip(value, mul(c, v[k]))]
        assert value == [E[k]] + [F(0)] * 15
    enclosures = [interval(c) for c in coefficients]
    signs = [1 if lo > 0 else -1 if hi < 0 else 0 for lo, hi in enclosures]
    assert all(signs) and signs == [r['sign'] for r in rows]
    negative = [i for i, s in zip(ids, signs) if s < 0]
    positive = [i for i, s in zip(ids, signs) if s > 0]
    assert negative == [8, 15] and len(positive) == 19
    assert all(signs[ids.index(i)] > 0 for i in d['positive_wall_indices'])
    # A homomorphism into a finite field with nonzero determinant proves
    # the physical 21 columns independent over the characteristic-zero field.
    denominators = [x.denominator for v in vectors for a in v for x in a]
    chosen = None
    for p in range(101, 10000, 2):
        if any(p % k == 0 for k in range(2, isqrt(p)+1)):
            continue
        if any(q % p == 0 for q in denominators):
            continue
        roots = []
        for q in [2, 79, 11, 3]:
            root = next((x for x in range(p) if x*x % p == q % p), None)
            if root is None:
                break
            roots.append(root)
        if len(roots) != 4:
            continue
        values = []
        for i in range(16):
            x = 1
            for j, root in enumerate(roots):
                if i & (1 << j):
                    x = x*root % p
            values.append(x)
        assert all(x*x % p == q % p for x, q in zip(values, RAD))
        for i in range(16):
            for j in range(16):
                g = gcd(RAD[i], RAD[j])
                k = INDEX[RAD[i]*RAD[j]//(g*g)]
                assert values[i]*values[j] % p == g*values[k] % p
        def ev(a):
            return sum(x.numerator*pow(x.denominator, -1, p)*y for x, y in zip(a, values)) % p
        matrix = [[ev(v[k]) for v in vectors] for k in range(21)]
        rank, pivots = mod_rank(matrix, p)
        if rank == 21:
            chosen = dict(prime=p, generator_square_roots=roots, basis_images=values,
                          matrix=matrix, rank=rank, pivot_columns=pivots)
            break
    assert chosen is not None
    G = [list(map(F, row)) for row in d['support_gram']]
    J = [list(map(F, row)) for row in d['support_inverse']]
    assert all(sum(G[i][k]*J[k][j] for k in range(18)) == int(i == j)
               for i in range(18) for j in range(18))
    assert all(q >= 0 for row in J for q in row)
    forced = [d['support'][i] for i in range(18) if J[i][i] >= 2]
    assert forced == d['mandatory_right_walls'] == profile['forced_right_walls']
    assert [J[d['support'].index(i)][d['support'].index(i)] for i in negative] == [F(24, 11), F(4)]
    Q = [list(map(F, row)) for row in old['form']]
    qE = sum(E[i]*Q[i][j]*E[j] for i in range(21) for j in range(21))
    assert qE == F(d['escape_norm']) == F(20181174579467302741, 10**18) > 0
    paths = [SOURCE, PRIMARY, BASE, CURRENT, PROFILE, f'scripts/{TAG}.py']
    return encode(dict(status='PASS_INDEPENDENT_S75_ESCAPE_DUAL', sources={p: sha(p) for p in paths},
        identity_coordinates_verified=21, negative_coefficient_walls=negative,
        positive_coefficient_walls=positive, coefficient_enclosures_192_bits=enclosures,
        physical_rank_witness=chosen, escape_norm=qE, mandatory_right_walls=forced,
        stronger_corollary='Any nonzero normal acute to the 19 positive-coefficient walls and right to walls8,15 has strictly negative escape pairing. Equality would force orthogonality to all21 independent roots.',
        scope='Independent coefficient identity, sign, modular physical-rank and Schur-gate checks; no replay of the old98 pair or vertex census.'))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    result = run()
    target = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], 'rank=', result['physical_rank_witness']['rank'])
