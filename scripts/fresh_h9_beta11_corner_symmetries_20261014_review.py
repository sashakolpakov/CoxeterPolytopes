#!/usr/bin/env python3
"""Separate exact verification of the full-chamber infinite-order symmetry."""
import argparse
import json
from hashlib import sha256
from pathlib import Path
from fresh_h9_h8_cap_height_20261014_norm_pool_review import (
    K, dot, mm, transpose, determinant)

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_beta11_corner_symmetries_20261014_review'
BASE = 'results/fresh_h9_beta11_corner_symmetries_20261014.json'
PREFIX = 'results/fresh_h9_beta11_vinberg_h300_20261014.json'
GEOMETRY = 'results/fresh_h9_beta11_geometry_h300_20261014.json'
HELPER = 'scripts/fresh_h9_h8_cap_height_20261014_norm_pool_review.py'

def sha(p): return sha256((ROOT / p).read_bytes()).hexdigest()
def read(p):
    d = json.loads((ROOT / p).read_text())
    for f, h in d['sources'].items(): assert sha(f) == h, f
    return d
def dec(x): return K(x['a'], x['b'])
def vec(x): return list(map(dec, x))
def mat(x): return list(map(vec, x))
def mv(a, v): return [dot(r, v) for r in a]
def identity(n): return [[K(i == j) for j in range(n)] for i in range(n)]
def enc(x):
    if isinstance(x, K): return dict(a=str(x.a), b=str(x.b))
    if isinstance(x, (list, tuple)): return [enc(z) for z in x]
    if isinstance(x, dict): return {k: enc(v) for k, v in x.items()}
    return x

def run():
    d, p, g = map(read, (BASE, PREFIX, GEOMETRY))
    G, N, H = mat(p['form']), mat(p['roots']), vec(p['basepoint'])
    T, TI = mat(d['matrix']), mat(d['inverse_matrix'])
    pair = lambda x, y: dot(x, mv(G, y))
    assert len(G) == 10 and len(N) == 35 and p['height_bound'] == '300'
    assert mm(T, TI) == mm(TI, T) == identity(10)
    assert mm(transpose(T), mm(G, T)) == G
    assert all(z.integral() for row in T + TI for z in row)
    assert determinant(T) == 1
    trace = sum((T[i][i] for i in range(10)), K())
    assert trace == K(4, 8) == dec(d['trace']) and (trace - 10).sign() > 0

    # Check both signatures directly from the literal orthogonal block form.
    phi = K(0, 1)
    assert G[0][0] == -2 * phi and G[0][9] == G[9][0] == -1 and G[9][9] == 2
    assert all(not G[i][j] for i in (0, 9) for j in range(1, 9))
    E8 = [row[1:9] for row in G[1:9]]
    minors = [determinant([row[:k] for row in E8[:k]]) for k in range(1, 9)]
    assert all(z.sign() > 0 and z.conjugate().sign() > 0 for z in minors)
    binary_det = G[0][0] * G[9][9] - G[0][9] * G[9][0]
    assert binary_det.sign() < 0 and binary_det.conjugate().sign() > 0
    assert G[0][0].conjugate().sign() > 0

    sids = list(range(9))
    images = [7, 4, 6, 5, 20, 29, 30, 11, 2]
    target = g['finite_vertices'][735]
    v = vec(target['point'])
    assert sorted(images) == target['incident']
    assert d['source_corner']['incident'] == sids
    assert d['target_corner']['ordered_images'] == images
    assert H == vec(d['source_corner']['point']) and v == vec(d['target_corner']['point'])
    assert H == [K(2)] + [K()] * 8 + [K(1)]
    assert pair(H, H).sign() < 0 and pair(v, v).sign() < 0 and pair(v, H).sign() < 0
    hp, vp = [pair(H, n) for n in N], [pair(v, n) for n in N]
    assert all(z.sign() <= 0 for z in hp + vp)
    assert [i for i, z in enumerate(hp) if not z] == sids
    assert [i for i, z in enumerate(vp) if not z] == sorted(images)
    B = [[pair(N[i], N[j]) for j in sids] for i in sids]
    assert B == [[pair(N[i], N[j]) for j in images] for i in images]
    assert all(B[i][j].sign() <= 0 for i in range(9) for j in range(i))
    spherical_minors = [determinant([row[:k] for row in B[:k]]) for k in range(1, 10)]
    assert all(z.sign() > 0 for z in spherical_minors)
    J = mat(d['inverse_gram'])
    assert mm(B, J) == mm(J, B) == identity(9)
    assert all(z.sign() >= 0 for row in J for z in row)

    a = pair(H, v) / pair(v, v)
    b = mv(J, [hp[i] for i in images])
    assert a.sign() > 0 and all(z.sign() <= 0 for z in b)
    assert [a * v[k] + sum((b[j] * N[images[j]][k] for j in range(9)), K())
            for k in range(10)] == H
    assert a == dec(d['preservation_decomposition']['a'])
    assert b == vec(d['preservation_decomposition']['coefficients'])
    scale = K(8, 14)
    assert scale.sign() > 0 and mv(T, H) == [scale * z for z in v]
    assert all(mv(T, N[i]) == N[j] for i, j in enumerate(images))
    assert all(pair(N[i], N[j]).sign() <= 0 for i in range(35) for j in range(i))
    C = vec(p['strict_center'])
    assert pair(C, C).sign() < 0 and pair(C, H).sign() < 0
    assert all(pair(C, n).sign() < 0 for n in N)
    return enc(dict(
        status='PASS_INDEPENDENT_INFINITE_FULL_CHAMBER_STABILIZER',
        sources={f: sha(f) for f in (BASE, PREFIX, GEOMETRY, HELPER, f'scripts/{TAG}.py')},
        arithmetic='Separate Fraction field and determinant implementation; no primary symmetry or Vinberg arithmetic imports.',
        determinant=K(1), trace=trace, trace_minus_dimension=trace-10,
        signatures=[[9, 1], [10, 0]], source_incident=sids, target_incident=sorted(images),
        outward_root_images=images, positive_time_scale=scale,
        spherical_leading_minors=spherical_minors, preservation_a=a, preservation_b=b,
        checks=dict(isometry_entries=100, integral_matrix_and_inverse_entries=200,
                    corner_wall_signs=70, outward_root_images=9),
        conceptual_proof_review='PASS: acute finite-vertex preservation, complete zero sector, local finiteness, oriented integral sector map, infinite order and infinite chamber volume.',
        conclusion='The full integral reflection group of this literal discriminant-11 lattice has infinite covolume. No finite-volume Coxeter polytope can have all wall reflections preserving this lattice.',
        scope='Independent witness verification and mathematical proof review. Root-prefix completeness and the norm alphabet are separately reviewed dependencies. No complete corner or symmetry census is repeated. Arbitrary geometric walls remain outside this exclusion.'))

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify', action='store_true'); args = ap.parse_args()
    result = run(); out = ROOT / f'results/{TAG}.json'
    if args.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])
