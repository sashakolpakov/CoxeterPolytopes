#!/usr/bin/env python3
"""Finite-angle endpoint/central anchors for the balanced intrinsic cusp.

Numerical proposal bank only. No lattice condition and no volume inference.
"""
import argparse
import json
import math
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h20_balanced_20261014 import ROOT, TAG, data, sha

PROFILE = 'results/fresh_h20_balanced_endpoint_20261014.json'


def allowed(x, cosines, tol=1e-8):
    x = np.asarray(x)
    out = x >= 1-tol
    for c in cosines:
        out |= abs(x-c) < tol
    return out


def run(mode):
    scales = (1., 1., (1+math.sqrt(5))/2 if mode == 'golden' else 1.)
    d = data(scales)
    p = json.loads((ROOT/PROFILE).read_text())
    labels = list(range(2, 17))+[20, 24, 30, 60]
    cosines = [0. if m == 2 else math.cos(math.pi/m) for m in labels]+[1.]
    anchors = list(zip(labels+['parallel'], cosines))
    profiles = p['types'][0]['profiles']
    z0 = d['vertices'][0]['point']
    ee = d['vectors'][p['types'][0]['edge_support']]
    A = np.array([7, 7, 8])*scales
    prepared = []
    for prof in profiles:
        c = np.array(prof['cosine_vector_approx'])
        y = np.linalg.solve(ee, c)
        pp, zy, zn = float(y@y), float(z0@y), float(z0@z0)
        B = np.array([r['subtract_approx'] for r in prof['omitted_contact_ratios']])
        trial = []
        for f in range(3):
            for label, cc in anchors:
                if B[f]+cc > 1e-10:
                    trial.append((float(A[f]/(B[f]+cc)), ['omitted', f, label]))
        for label, cc in anchors:
            aa, bb, dd = 1-pp, 2*zy+4*cc, zn-4
            if abs(aa) < 1e-9:
                if bb > 1e-9:
                    trial.append((dd/bb, ['central', label]))
            else:
                rr = 2*dd/(math.sqrt(bb*bb+4*aa*dd)+bb)
                trial.append((rr, ['central', label]))
        radii = []
        for R, anchor in sorted(trial):
            if R <= 0 or not allowed(A/R-B, cosines).all():
                continue
            central = (zn-2*R*zy+(pp-1)*R*R-4)/(4*R)
            if not bool(allowed(central, cosines)):
                continue
            if any(abs(R-r[0]) < 1e-8 for r in radii):
                continue
            radii.append((R, anchor))
        prepared.append(radii)
    centers, radii, recipes = [], [], []
    seen = set()
    for ti, tr in enumerate(p['literal_edge_transports']):
        assert tr['omitted'] == d['vertices'][ti]['omitted']
        v = d['vectors'][tr['edge_support_in_profile_order']]
        z = d['vertices'][ti]['point']
        yy = np.linalg.solve(v, np.array([q['cosine_vector_approx'] for q in profiles]).T).T
        for prof, y, rr in zip(profiles, yy, prepared):
            for R, anchor in rr:
                w = z-R*y
                assert allowed((d['offsets']-d['vectors']@w)/R, cosines, 1e-7).all()
                key = tuple(np.rint(np.append(w, R)*1e8).astype(np.int64))
                if key in seen:
                    continue
                seen.add(key)
                centers.append(w)
                radii.append(R)
                recipes.append(dict(transport=ti, profile=prof['id'], radius=R,
                                    anchor=anchor, endpoint=prof['endpoint']))
    C, R = np.array(centers), np.array(radii)
    V = np.array([v['point'] for v in d['vertices']])
    assert len(R)
    X = np.concatenate([np.zeros((1, 19))]+[V*s for s in [.125, .25, .5, .75, 1.]])
    minima, nearest = [], []
    for x in X:
        power = np.sum((C-x)**2, axis=1)-R*R
        i = int(np.argmin(power))
        minima.append(min(float(power[i]), float(x@x-4)))
        nearest.append(i)
    worst = int(np.argmax(minima))
    arrays = f'data/enumerations/{TAG}_{mode}_endpoint_bank.npz'
    np.savez_compressed(ROOT/arrays, centers=C, radii=R, probes=X,
                        minimum_powers=np.array(minima), nearest=np.array(nearest))
    paths = [PROFILE, f'scripts/{TAG}.py', str(Path(__file__).resolve().relative_to(ROOT)), arrays]
    out = dict(status='NUMERICAL_BALANCED_ENDPOINT_AND_CENTRAL_ANCHOR_BANK',
               sources={s: sha(s) for s in paths}, mode=mode, scales=list(scales),
               finite_labels=labels, profile_radius_counts=list(map(len, prepared)),
               bank_size=len(R), endpoint_types=dict(Counter(r['endpoint'] for r in recipes)),
               radius_range=[float(R.min()), float(R.max())], probe_count=len(X),
               uncovered_probe_count=sum(v > 1e-7 for v in minima),
               worst_probe=dict(index=worst, point=X[worst].tolist(), minimum_power=minima[worst],
                                closest_endpoint_cap=nearest[worst]), recipes=recipes,
               scope='Finite numerical anchors at displayed omitted or central contacts, plus parallel boundaries. Every cap passes seed and central pair gates, but caps are not yet mutually checked. Probe coverage is not a finite-volume certificate; PSD endpoints have no bounded-link certificate.')
    target = ROOT/f'results/{TAG}_{mode}_endpoint_scout.json'
    target.write_text(json.dumps(out, indent=2)+'\n')
    print({k: v for k, v in out.items() if k not in ('sources', 'recipes', 'scope')}, flush=True)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['unit', 'golden'], default='unit')
    run(ap.parse_args().mode)
