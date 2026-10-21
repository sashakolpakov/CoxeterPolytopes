#!/usr/bin/env python3
"""Numerical contact search with a variable Coxeter central-cap radius.

Searches identical scaled-vertex caps, not arbitrary cap pairs or completions.
The output is a bounded numerical scout and cannot certify nonexistence.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from flint import arb, ctx
from fresh_h20_balanced_20261014 import ROOT, data, sha

TAG = 'fresh_h20_balanced_central_radius_20261014'


def review_near_hit(row):
    ctx.prec = 192
    def cosine(m):
        return arb(1) if m == 'parallel' else (arb.pi()/m).cos()
    ck = cosine(row['central_seed_label'])
    cn = cosine(row['bridge_seed_label'])
    cm = cosine(row['bridge_central_label'])
    r0 = 1/ck
    aa = 196-r0*r0
    bb = 2*r0*r0+2*r0*cm/cn
    cc = 1/(cn*cn)+r0*r0+2*r0*cm/cn
    tt = (bb+(bb*bb+4*aa*cc).sqrt())/(2*aa)
    paircos = tt*tt*cn*cn*row['vertex_separation_squared']/2-1
    error = paircos-cosine(row['bridge_pair_label'])
    assert not error.contains(0)
    return dict(**row, certified_difference_ball=str(error),
                verdict='REJECTED_NONZERO_ARB_INTERVAL')


def run(limit):
    d = data()
    vv = np.array([v['point'] for v in d['vertices']])
    distances = sorted(set(map(int, np.rint(np.sum((vv-vv[0])**2, axis=1))))-{0})
    labels = list(range(3, limit+1))+['parallel']
    cosines = np.array([math.cos(math.pi/m) if isinstance(m, int) else 1. for m in labels])
    mlabs = [2]+labels
    mc = np.r_[0., cosines]
    nn, mm = np.meshgrid(cosines, mc, indexing='ij')
    records, closest = [], None
    examined = 0
    for ki, ck in enumerate(cosines):
        r0 = 1/ck
        aa = 196-r0*r0
        bb = 2*r0*r0+2*r0*mm/nn
        cc = 1/(nn*nn)+r0*r0+2*r0*mm/nn
        tt = (bb+np.sqrt(bb*bb+4*aa*cc))/(2*aa)
        gamma = tt/(1+tt)
        rr = 1/(nn*(1+tt))
        omitted = [(1+6*gamma)/rr, (1+7*gamma)/rr]
        # Retained contacts and central contacts hold by construction.
        # For this bounded bank both omitted-factor contacts are dotted.
        assert all(np.min(o) > 1 for o in omitted)
        lam = tt*tt*nn*nn
        examined += lam.size
        for ds in distances:
            paircos = lam*ds/2-1
            ii, jj = np.where((paircos >= -1e-12)&(paircos < 1-1e-12))
            if not len(ii):
                continue
            c = paircos[ii, jj]
            real_label = np.pi/np.arccos(np.clip(c, -1, 1))
            rounded = np.maximum(np.rint(real_label).astype(int), 2)
            target = np.cos(np.pi/rounded)
            error = np.abs(c-target)
            ix = int(np.argmin(error))
            def rec(p):
                ni, mi = int(ii[p]), int(jj[p])
                return dict(central_seed_label=labels[ki], bridge_seed_label=labels[ni],
                            bridge_central_label=mlabs[mi], bridge_pair_label=int(rounded[p]),
                            vertex_separation_squared=int(ds), central_radius=float(r0),
                            gamma=float(gamma[ni, mi]), radius=float(rr[ni, mi]),
                            pair_cosine=float(c[p]), cosine_error=float(error[p]))
            if closest is None or error[ix] < closest['cosine_error']:
                closest = rec(ix)
            records += [rec(int(p)) for p in np.where(error < 1e-10)[0]]
    sources = [f'scripts/{TAG}.py', 'scripts/fresh_h20_balanced_20261014.py']
    out = dict(status='BOUNDED_NUMERICAL_IDENTICAL_VERTEX_BRIDGE_CONTACT_SCOUT',
               sources={p: sha(p) for p in sources}, finite_label_bound=limit,
               cap_profiles=examined, vertex_distance_alphabet=distances,
               near_hit_reviews=[review_near_hit(r) for r in records],
               surviving_near_hits=0, closest_approximate_contact=closest,
               scope='Unit A6+A6+A7 cusp. R0=1/cos(pi/k); centers gamma*z_v; retained seed label n; central label m; all k,n=3..bound or parallel, m=2..bound or parallel. Searches pairs sharing one profile, with any inferred integer contact label. Near hits below1e-10 are individually rejected by192-bit Arb intervals. The broad float sweep is not an exact exclusion; arbitrary-profile pairs and finite volume are not checked.')
    path = ROOT/f'results/{TAG}.json'
    path.write_text(json.dumps(out, indent=2)+'\n')
    print({k: v for k, v in out.items() if k not in ('sources', 'scope', 'near_hit_reviews')}, flush=True)
    print('numerical_hits', len(records), flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=100)
    run(ap.parse_args().limit)
