#!/usr/bin/env python3
"""Complete sampled cusp orbits and test their height-one H22 caps.

Integer orbit construction and compatibility checks are exact. The final
volume screen is numerical and is separately labelled.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fresh_h22_borcherds_extension_20261013 import prepare

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'fresh_h22_cusp_orbits_20261013'
SAMPLE = 'data/enumerations/fresh_h22_borcherds_extension_20261013_scout.json'
SYMS = 'data/enumerations/fresh_h20_global_interior_sections_20261013.json'
OLD = 'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz'
ESCAPE = 'results/fresh_h22_symmetry_orbit_verify_20261013.json'


def read(p):
    return json.loads((ROOT/p).read_text())


def run():
    Q, W, G, h, X, A, b, Qn, R = prepare()
    generators = [np.array(x['matrix'], dtype=np.int64) for x in read(SYMS)['ambient_symmetry_generators']]
    U = np.load(ROOT/OLD)['twice_old_projections']
    candidates = [tuple(x['vector']) for x in read(SAMPLE)['old_vertices'] if Fraction(x['norm']) == 0]
    # Prefer the orbit that cuts the old exact escape.
    saved = read(ESCAPE)['strict_escape']
    wnum = np.array(list(map(int, saved['numerator'])), dtype=np.int64)
    candidates.sort(key=lambda v: int(np.array(v)@Q@wnum[:22]), reverse=True)
    done = set(); accepted = []; heights = []; records = []
    for seed in candidates:
        if seed in done:
            continue
        v = np.array(seed, dtype=np.int64)
        assert int(v@Q@v) == 0 and np.max(W@Q@v) == 0 and int(h@Q@v) < 0
        orbit = [v]; seen = {seed}; cursor = 0
        while cursor < len(orbit):
            s = orbit[cursor]; cursor += 1
            for g in generators:
                t = g@s; key = tuple(t)
                if key not in seen:
                    seen.add(key); orbit.append(t)
            assert len(orbit) <= 241920
        done.update(seen)
        D = np.array(orbit, dtype=np.int64)
        assert np.all(np.einsum('ij,jk,ik->i', D, Q, D) == 0)
        hp = int(h@Q@v); pp = U@(Q@v)
        oldpairs = W@Q@v
        # These strict bounds give dotted pairs for both cap height signs,
        # uniformly in the old center-cap parameter m>=3.
        compatible = int(max(pp)) <= -6 and hp <= -6
        height = 1 if compatible else 0.25
        assert compatible or (int(max(pp)) <= -2 and hp <= -2)
        rec = {'representative': list(map(int, v)), 'size': len(D), 'center_pairing': hp,
               'maximum_twice_timelike_orbit_pairing': int(max(pp)), 'height_one_compatible': compatible,
               'height':str(Fraction(height)),
               'old_pair_counts': dict(Counter(map(str, oldpairs))),
               'old_incident_walls': np.flatnonzero(oldpairs == 0).tolist(),
               'old_escape_signed_cut_numerator': int(v@Q@wnum[:22]+abs(wnum[-1]))}
        rec['array_offset'] = len(accepted)
        accepted.extend(orbit); heights.extend([height]*len(orbit))
        records.append(rec)
        print('cusp orbit', len(records), len(D), hp, int(max(pp)), compatible, 'total', len(accepted), flush=True)
    D = np.array(accepted, dtype=np.int64)
    array = ROOT/'data/enumerations'/f'{PREFIX}.npz'
    heights = np.array(heights)
    np.savez_compressed(array, old_null_projections=D, four_heights=(4*heights).astype(np.int64))
    cuts = D@Q@wnum[:22]+heights*abs(wnum[-1])
    print('cusp directions', len(D), 'cutters of exact escape', sum(cuts>0), flush=True)
    out = {'status':'CUSP_ORBITS_WITH_NUMERICAL_VOLUME_SCREEN',
           'sources':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (SAMPLE,SYMS,OLD,ESCAPE)},
           'orbits':records,'accepted_null_directions':len(D),'signed_cusp_caps':2*len(D),
           'total_wall_count':121172+2*len(D),'old_escape_cutters':int(sum(cuts>0)),
           'array':str(array.relative_to(ROOT)),'array_sha256':hashlib.sha256(array.read_bytes()).hexdigest(),
           'runs':[],'scope':'Complete symmetry orbits of sampled old cusps, not a complete old cusp census. Numerical volume screens require separate exactification. Cusp caps have exact integral old pairings; all pair and bounded-link proofs require standalone replay.'}
    output = ROOT/'data/enumerations'/f'{PREFIX}.json'
    output.write_text(json.dumps(out,indent=2)+'\n')
    # Signed caps (delta,+/-1), along with the old complete timelike orbit.
    B = U@Q@X/2; bb = -U@Q@h/(2*np.sqrt(14))
    C = D@Q@X; cc = -D@Q@h/np.sqrt(14)
    # Constraint generation avoids passing millions of dense rows to HiGHS.
    P = np.vstack([B,C]); limits = np.r_[bb,cc]
    ht = np.r_[np.full(len(B),np.sqrt(2.5)),heights]
    rng = np.random.default_rng(221310)
    runs = []
    for m in (3,4,5,6,8,12):
        height = np.sqrt(14+1/(2*np.cos(np.pi/m)**2))
        AA = np.vstack([np.column_stack([A,np.zeros(len(A))]),[[0]*21+[height],[0]*21+[-height]]]); rhs = np.r_[b,np.sqrt(14),np.sqrt(14)]
        active = set()
        escape = None; maximum = -1
        for k in range(30):
            objective = rng.normal(size=22)
            for iteration in range(150):
                chosen = sorted(active)
                extra = np.array([np.r_[P[i],s*ht[i]]for i,s in chosen]).reshape((-1,22))
                p = linprog(objective,A_ub=np.vstack([AA,extra]),b_ub=np.r_[rhs,[limits[i]for i,s in chosen]],bounds=[(None,None)]*22,method='highs')
                assert p.success,p.message
                violations = P@p.x[:21]+ht*abs(p.x[-1])-limits
                if max(violations)<2e-7:
                    break
                ids = np.flatnonzero(violations>2e-7)
                ids = ids[np.argsort(violations[ids])[-100:]]
                sign = 1 if p.x[-1]>=0 else -1
                before = len(active);active.update((int(i),sign)for i in ids)
                assert len(active)>before
            else:
                raise RuntimeError('bounded constraint-generation limit reached')
            if not p.success:
                continue
            q = float(p.x@p.x-1); maximum = max(maximum,q)
            if q > 1e-7:
                point = np.r_[h/np.sqrt(14)+X@p.x[:21],p.x[-1]]
                escape = {'point':point.tolist(),'norm':q,
                          'active_old':np.flatnonzero(abs(A@p.x[:21]-b)<1e-6).tolist()}
                break
        runs.append({'m':m,'largest_sample_square':maximum,'escape':escape,'active_constraints':len(active)})
        out['runs']=runs;output.write_text(json.dumps(out,indent=2)+'\n')
        print('m',m,'numerical escape',None if escape is None else escape['norm'],flush=True)
    out['completed']=True;output.write_text(json.dumps(out,indent=2)+'\n')


if __name__ == '__main__':
    run()
