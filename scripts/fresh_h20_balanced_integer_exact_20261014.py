#!/usr/bin/env python3
"""Exact810-wall extension of the balanced bridge chamber; still infinite."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from flint import fmpq as Q, fmpq_mat as M
from fresh_h20_balanced_20261014 import exact_data

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_balanced_integer_exact_20261014'
SOURCES = ['results/fresh_h20_balanced_integer_caps_20261014.json',
           'results/fresh_h20_balanced_integer_joint_20261014_motif.json',
           'results/fresh_h20_balanced_20261014_bridge_exact.json',
           'scripts/fresh_h20_balanced_20261014.py', f'scripts/{TAG}.py']


def run():
    scout, joint, motif = [json.loads((ROOT/p).read_text()) for p in SOURCES[:3]]
    for doc in (scout, joint, motif):
        for p, h in doc['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    E = exact_data()
    B = M(21, 21)
    for i, row in enumerate(E['spatial_Gram']):
        for j, x in enumerate(row):
            B[i, j] = Q(x.numerator, x.denominator)
    B[19, 20] = B[20, 19] = -1
    inv = M([[B[i, j] for j in range(19)] for i in range(19)]).inv()
    def rat(x):
        return Q(x.numerator, x.denominator)
    seed = [[rat(x) for x in r]+[Q(1), Q(0)] for r in E['seed_vectors']]
    vertices, local_vertices = [], []
    for a in range(7):
        for b in range(7):
            for c in range(8):
                local_vertices.append((a, b, c))
                vertices.append([sum(rat(E['factor_vertices'][k][v][j]) for k, v in enumerate((a, b, c))) for j in range(19)])
    def dot(c, d):
        return (M([c+[Q(0), Q(0)]])*B*M([d+[Q(0), Q(0)]]).transpose())[0, 0]
    def cap(c, radius2):
        return c+[(dot(c, c)-radius2)/2, Q(1)]
    roots = seed+[[Q(0)]*19+[Q(-2), Q(1)]]+[cap(v, Q(28)) for v in vertices]
    roots += [cap(list(map(Q, r['center'])), Q(r['radius'])**2) for r in motif['new_caps']]
    records = []
    for ti, pi in joint['best']['selected']:
        base = scout['recipes'][pi]['gaps']
        gaps = []
        for g, shift in zip(base, local_vertices[ti]):
            gaps += [g[(j-shift) % len(g)] for j in range(len(g))]
        rhs = M([[1-gaps[j]] for j in range(22) if j not in (0, 7, 14)])
        coeff = inv*rhs
        center = [coeff[j, 0] for j in range(19)]
        root = cap(center, Q(4))
        assert [1-dot(center, s[:19]) for s in seed] == gaps
        x = [Q(2, 5)*v for v in vertices[ti]]
        power = dot(x, x)+dot(center, center)-2*dot(x, center)-4
        assert power < 0
        records.append(dict(probe=ti, base_recipe=pi, seed_gaps=gaps,
                            center=list(map(str, center)), probe_power=str(power)))
        roots.append(root)
    assert len(records) == 392 and len(roots) == 810
    assert len(set(tuple(r) for r in roots)) == 810
    W = M(roots)
    G = W*B*W.transpose()
    counts = Counter()
    for i in range(810):
        assert G[i, i] > 0
        for j in range(i):
            g, den2 = G[i, j], G[i, i]*G[j, j]
            if not g:
                typ = 'right'
            elif g < 0 and 4*g*g == den2:
                typ = 'pi3'
            else:
                assert g < 0 and g*g >= den2, (i, j, g, den2)
                typ = 'parallel' if g*g == den2 else 'dotted'
            counts[typ] += 1
    H = M([[Q(0)]*19+[Q(200), Q(1)]])
    HW = H*B*W.transpose()
    assert all(HW[0, i] < 0 for i in range(810))
    for i in range(810):
        assert -400-HW[0, i]**2/G[i, i] < 0
        assert all(HW[0, j]-HW[0, i]*G[i, j]/G[i, i] < 0 for j in range(810) if i != j)
    P = M([v+[Q(112), Q(1)] for v in vertices])
    PW = P*B*W.transpose()
    for i in range(392):
        assert sum(PW[i, j] == 0 for j in range(810)) == 20
        assert all(PW[i, j] <= 0 for j in range(810))
    x = [Q(3, 5)*v for v in vertices[0]]
    x2 = dot(x, x)
    escape = x+[(x2-1)/2, Q(1)]
    X = M([escape])
    assert (X*B*X.transpose())[0, 0] == 1
    assert (X*B*H.transpose())[0, 0] < 0
    assert all(v < 0 for v in list(X*B*W.transpose()))
    return dict(status='PASS_EXACT_810_FACETS_AND_ALL392_PROBE_CUTS',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
                wall_count=810, pair_count=810*809//2, pair_types=dict(counts),
                actual_facet_witnesses=810, preserved_finite_vertices=392,
                covered_probe_multiple='2/5', covered_probe_count=392,
                new_caps=records, walls=[list(map(str, r)) for r in roots],
                strict_escape=dict(spatial_multiple='3/5', old_vertex=0,
                                   point=list(map(str, escape)), norm=1),
                scope='Exact coherent extension retaining all418 previous facets. Every .4 old-vertex probe is cut, but a .6 probe gives a strict infinite-volume escape. No complete new vertex/cusp enumeration or finite-volume claim.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert json.loads(path.read_text()) == out
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['pair_types'], flush=True)
