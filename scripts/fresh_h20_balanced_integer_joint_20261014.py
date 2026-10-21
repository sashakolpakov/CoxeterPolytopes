#!/usr/bin/env python3
"""Joint integer-cap covering scout, with exact rational distance-grid gates.

The candidate geometry and greedy search are numerical; selected proposals
must be independently exactified before being called coherent chambers.
"""
import argparse
import json
import time
from fractions import Fraction
from itertools import product
import numpy as np
from fresh_h20_balanced_20261014 import ROOT, data, sha

TAG = 'fresh_h20_balanced_integer_joint_20261014'
SOURCE = 'results/fresh_h20_balanced_integer_caps_20261014.json'


def prepare():
    d = data()
    src = json.loads((ROOT/SOURCE).read_text())
    recipes = src['recipes']
    base = [np.array([r['gaps'][f] for r in recipes]) for f in range(3)]
    C = []
    for local in product(range(7), range(7), range(8)):
        parts = []
        for f, (n, lo) in enumerate(((6, 0), (6, 6), (7, 12))):
            gaps = np.roll(base[f], local[f], axis=1)
            parts.append((1-gaps[:, 1:])@np.linalg.inv(d['spatial_Gram'][lo:lo+n, lo:lo+n]))
        C.append(np.concatenate(parts, axis=1)@d['cholesky'])
    C = np.array(C)
    X = .4*np.array([v['point'] for v in d['vertices']])
    assert C.shape == (392, 555, 19)
    assert np.max(np.sum((C-X[:, None, :])**2, axis=2)-4) < -1e-5
    return d, C, X


def run(restarts, keep_motif):
    d, C, X = prepare()
    flat = C.reshape(-1, 19)
    norms = np.sum(flat*flat, axis=1)
    allowed = np.ones((392, 555), dtype=bool)
    if keep_motif:
        motif = json.loads((ROOT/'results/fresh_h20_balanced_20261014_bridge_exact.json').read_text())
        for r in motif['new_caps']:
            c = np.array([float(Fraction(x)) for x in r['center']])@d['cholesky']
            rad = float(Fraction(r['radius']))
            cos = (norms+c@c-2*flat@c-4-rad*rad)/(4*rad)
            gate = (abs(cos) < 1e-7) | (abs(cos-.5) < 1e-7) | (cos >= 1-1e-7)
            allowed &= gate.reshape(392, 555)
    initial_sizes = allowed.sum(axis=1)
    assert initial_sizes.min() > 0
    pair_examples = []
    for ti in (1, 8, 56):
        dist = np.sum(C[0]*C[0], axis=1)[:, None]+np.sum(C[ti]*C[ti], axis=1)-2*C[0]@C[ti].T
        grid = np.rint(28*dist).astype(int)
        assert np.max(abs(28*dist-grid)) < 1e-7
        compatible = (grid == 0) | (grid == 224) | (grid == 336) | (grid >= 448)
        subset = compatible & allowed[0][:, None] & allowed[ti][None, :]
        pair_examples.append(dict(probes=[0, ti], compatible_pairs=int(subset.sum())))
    best = None
    attempts = []
    rng = np.random.default_rng(20261014)
    began = time.time()
    for attempt in range(restarts):
        live = allowed.copy()
        assigned = np.zeros(392, dtype=bool)
        selection = []
        while not assigned.all():
            sizes = live.sum(axis=1)
            unassigned = np.where(~assigned)[0]
            least = sizes[unassigned].min()
            if least == 0:
                break
            choices = unassigned[sizes[unassigned] == least]
            ti = int(choices[0] if attempt == 0 else rng.choice(choices))
            choices = np.where(live[ti])[0]
            if attempt == 0:
                distances = np.sum((C[ti, choices]-X[ti])**2, axis=1)
                pi = int(choices[np.argmin(distances)])
            else:
                pi = int(rng.choice(choices))
            v = C[ti, pi]
            dist = norms+v@v-2*flat@v
            grid = np.rint(28*dist).astype(int)
            assert np.max(abs(28*dist-grid)) < 1e-7
            compatible = (grid == 0) | (grid == 224) | (grid == 336) | (grid >= 448)
            live &= compatible.reshape(392, 555)
            assigned[ti] = True
            selection.append([ti, pi])
            # A single cap may already cover multiple required probes.
            assigned |= np.sum((X-v)**2, axis=1) < 4-1e-8
        row = dict(attempt=attempt, covered_probes=int(assigned.sum()),
                   chosen_caps=len(selection), zero_domains=np.where(~assigned & (live.sum(axis=1) == 0))[0].tolist(),
                   selected=selection)
        attempts.append(row)
        if best is None or row['covered_probes'] > best['covered_probes']:
            best = row
        print('attempt', attempt, 'covered', row['covered_probes'], 'caps', len(selection),
              'elapsed', round(time.time()-began, 2), flush=True)
        if assigned.all():
            break
    paths = [SOURCE, 'scripts/fresh_h20_balanced_20261014.py', f'scripts/{TAG}.py']
    if keep_motif:
        paths.append('results/fresh_h20_balanced_20261014_bridge_exact.json')
    out = dict(status='BOUNDED_GREEDY_INTEGER_CAP_JOINT_COVER',
               sources={p: sha(p) for p in paths}, keep_old_three_cap_motif=keep_motif,
               candidate_occurrences=len(flat), initial_domain_range=[int(initial_sizes.min()), int(initial_sizes.max())],
               nearest_probe_pair_screen=pair_examples, best=best, attempts=attempts,
               scope='Finite555x392 integer-gap cap occurrences. Rational squared-distance grid checked numerically before integer compatibility gates. Greedy failure is not infeasibility. Probe coverage is not finite volume. Exact geometry replays remain required.')
    suffix = 'motif' if keep_motif else 'core'
    (ROOT/f'results/{TAG}_{suffix}.json').write_text(json.dumps(out, indent=2)+'\n')
    print('BEST', best['covered_probes'], 'of392; initial domains', out['initial_domain_range'], flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--restarts', type=int, default=12)
    ap.add_argument('--keep-motif', action='store_true')
    a = ap.parse_args()
    run(a.restarts, a.keep_motif)
