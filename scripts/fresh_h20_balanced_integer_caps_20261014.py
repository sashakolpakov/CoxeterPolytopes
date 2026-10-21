#!/usr/bin/env python3
"""Fixed radius-two caps with integer seed gaps, proposed numerically.

The product composition search exactly describes seed contacts; floating
pair and coverage checks only propose candidates for separate exact review.
"""
import json
from itertools import product
from pathlib import Path
import numpy as np
from fresh_h20_balanced_20261014 import ROOT, data, sha

TAG = 'fresh_h20_balanced_integer_caps_20261014'


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for i in range(total+1):
            for rest in compositions(total-i, length-1):
                yield (i,)+rest


def rational_cosine_allowed(c):
    return (abs(c) < 1e-8) | (abs(c-.5) < 1e-8) | (c >= 1-1e-8)


def run():
    d = data()
    vv = np.array([v['point'] for v in d['vertices']])
    x = .4*vv[0]
    prepared = []
    for n, lo in ((6, 0), (6, 6), (7, 12)):
        inverse = np.linalg.inv(d['spatial_Gram'][lo:lo+n, lo:lo+n])
        gaps = np.array(list(compositions(n+1, n+1)), dtype=int)
        scalar = 1-gaps[:, 1:]
        coeff = scalar@inverse
        xx = .4*np.ones(n)@inverse
        dist = np.einsum('ij,ij->i', scalar-.4, coeff-xx)
        keep = np.where(dist < 4+1e-9)[0]
        prepared.append(dict(gaps=gaps[keep], coeff=coeff[keep], dist=dist[keep],
                             full_count=len(gaps)))
    recipes, coefficient_rows = [], []
    for i, j, k in product(*[range(len(p['dist'])) for p in prepared]):
        idx = (i, j, k)
        power = sum(p['dist'][a] for p, a in zip(prepared, idx))-4
        if power >= -1e-8:
            continue
        coefficient_rows.append(np.concatenate([p['coeff'][a] for p, a in zip(prepared, idx)]))
        recipes.append(dict(gaps=[p['gaps'][a].tolist() for p, a in zip(prepared, idx)],
                            probe_power=float(power)))
    coeff = np.array(coefficient_rows)
    C = coeff@d['cholesky']
    norm = np.sum(C*C, axis=1)
    power0 = np.sum((C-x)**2, axis=1)-4
    assert np.max(abs(power0-np.array([r['probe_power'] for r in recipes]))) < 1e-8
    central_cos = (norm-8)/8
    central_ok = rational_cosine_allowed(central_cos)
    # Rational distance andR2 force right or dotted contact withRsqrt28:
    # no nonzero rational-angle cosine has square of the form rational_square/7.
    vertex_dist = norm[:, None]+196-2*C@vv.T
    vertex_ok = ((abs(vertex_dist-32) < 1e-8) | (vertex_dist > 32+8*np.sqrt(7)-1e-8)).all(axis=1)
    core_ok = central_ok & vertex_ok
    motif = json.loads((ROOT/'results/fresh_h20_balanced_20261014_bridge_exact.json').read_text())
    old_c = np.array([[float(__import__('fractions').Fraction(v)) for v in r['center']]
                      for r in motif['new_caps']])@d['cholesky']
    old_r = np.array([float(__import__('fractions').Fraction(r['radius'])) for r in motif['new_caps']])
    pair_cos = (norm[:, None]+np.sum(old_c*old_c, axis=1)-2*C@old_c.T-4-old_r**2)/(4*old_r)
    motif_gate = rational_cosine_allowed(pair_cos)
    accepted = np.where(core_ok & motif_gate.all(axis=1))[0]
    partial = {str(mask): np.where(core_ok & np.all(motif_gate[:, [j for j in range(3) if mask >> j & 1]], axis=1))[0].tolist()
               for mask in range(8)}
    for i, row in enumerate(recipes):
        row.update(center_square=float(norm[i]), central_cosine=float(central_cos[i]),
                   central_ok=bool(central_ok[i]), full_vertex_layer_ok=bool(vertex_ok[i]),
                   bridge_pair_cosines=pair_cos[i].tolist())
    arrays = f'data/enumerations/{TAG}.npz'
    np.savez_compressed(ROOT/arrays, centers=C, coefficients=coeff,
                        vertex_distances=vertex_dist, motif_pair_cosines=pair_cos,
                        core_ok=core_ok, all418_ok=core_ok & motif_gate.all(axis=1))
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_h20_balanced_20261014.py',
             'results/fresh_h20_balanced_20261014_bridge_exact.json', arrays]
    out = dict(status='NUMERICAL_RADIUS_TWO_INTEGER_GAP_CUTTERS',
               sources={p: sha(p) for p in paths}, radius=2, probe_multiple='2/5',
               factor_composition_counts=[p['full_count'] for p in prepared],
               probe_cutters=len(C), core_compatible=int(sum(core_ok)),
               all418_compatible=len(accepted), all418_indices=accepted.tolist(),
               motif_subset_compatible_counts={k: len(v) for k, v in partial.items()},
               motif_subset_compatible_indices=partial, recipes=recipes,
               power_range=[float(power0.min()), float(power0.max())],
               center_square_range=[float(norm.min()), float(norm.max())],
               scope='All nonnegative integer-gap compositions summing7,7,8 at radius2, numerically restricted to strict cutters ofone .4vertex probe. Central and392 Rsqrt28 vertex-layer gates; optional existing3-cap motif gates. This does not test mutual compatibility of proposed caps or prove finite volume.')
    (ROOT/f'results/{TAG}.json').write_text(json.dumps(out, indent=2)+'\n')
    print({k: v for k, v in out.items() if k not in ('sources', 'recipes', 'scope', 'motif_subset_compatible_indices', 'all418_indices')}, flush=True)


if __name__ == '__main__':
    run()
