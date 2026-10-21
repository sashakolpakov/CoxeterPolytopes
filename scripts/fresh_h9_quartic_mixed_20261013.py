#!/usr/bin/env python3
"""Complete both-F4-block census, with exact unit-balanced coefficient caps."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path
import time

from fresh_h9_glued_host_20261013 import Host, geometry
from fresh_h9_glued_host_20261013_arithmetic import *

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'results/fresh_h9_glued_host_20261013_h9_h100.json'
PREVIOUS = ROOT/'results/fresh_h9_glued_host_20261013_symmetry.json'
UNITS = (T-1, T**3-3*T-1, -T**3+3*T-1)
assert all(abs(u.norm()) == 1 and u.inverse().is_integral for u in UNITS)


@lru_cache(None)
def balanced_alphabet(caps):
    """Unit choice is heuristic; the resulting finite enumeration is exact.

    x=u*y with u an integral unit is a bijection of O_K. We enumerate all y
    in the exact trace ball implied by y_e^2 <= cap_e/u_e^2, and finally check
    every original embedding inequality. No rounded numerical bound is used.
    """
    if any(c.sign(e) < 0 for e, c in enumerate(caps)):
        return ()
    if any(c.sign(e) == 0 for e, c in enumerate(caps)):
        return (Z,)
    import numpy as np
    values = np.array([c.real(e) for e, c in enumerate(caps)])
    logs = np.array([[np.log(abs(u.real(e))) for u in UNITS] for e in range(4)])
    target = (np.log(values)-np.log(values).mean())/2
    estimate = np.linalg.lstsq(logs, target, rcond=None)[0]
    rounded = [int(round(x)) for x in estimate]
    exponents = min((tuple(k+d for k, d in zip(rounded, delta))
                     for delta in product((-1,0,1), repeat=3)),
                    key=lambda ks: float(sum(values/np.exp(2*logs@np.array(ks)))))
    unit = O
    for u, k in zip(UNITS, exponents):
        unit *= u**k
    if unit < 0:
        unit = -unit
    assert unit.is_integral and unit.inverse().is_integral
    transformed = tuple(c/(unit*unit) for c in caps)
    bound = sum(ceil_embedding(c, e) for e, c in enumerate(transformed))
    found = []
    for y, ysquare in trace_ball(bound):
        if any((c-ysquare).sign(e) < 0 for e, c in enumerate(transformed)):
            continue
        x = unit*y
        assert x.is_integral and x >= 0
        assert all((c-x*x).sign(e) >= 0 for e, c in enumerate(caps))
        found.append(x)
    return tuple(sorted(found))


@lru_cache(None)
def block_census(energy):
    alphabets = [balanced_alphabet((energy/FI[j][j], 2*energy, 2*energy, 2*energy))
                 for j in range(4)]
    found, nodes = [], 0

    def walk(prefix):
        nonlocal nodes
        nodes += 1
        k = len(prefix)
        if k:
            lower = pair(prefix, SCHUR[k-1], prefix)
            if not le_all(lower, energy):
                return
            primary_lower = pair(prefix, tuple(row[:k] for row in FI[:k]), prefix)
            if primary_lower > energy:
                return
        if k == 4:
            vector = mv(FI, prefix)
            assert all(z.is_integral for z in vector)
            found.append((lower, vector))
            return
        for value in alphabets[k]:
            walk((*prefix, value))

    walk(())
    return tuple(found), {'coordinate_alphabets': list(map(len, alphabets)),
                          'block_nodes': nodes, 'block_vectors': len(found)}


def coxeter(root, q, covectors, norms):
    for cov, nq in zip(covectors, norms):
        p = dot(root, cov)
        if p > 0:
            return False
        ratio = 4*p*p/(q*nq)
        if ratio < 4 and ratio not in ALLOWED:
            return False
    return True


def run(height, output):
    started = time.monotonic()
    old_data = json.loads(SOURCE.read_text())
    previous_data = json.loads(PREVIOUS.read_text())
    g = Host(2).g
    old = [tuple(map(K.from_json, r)) for r in old_data['all_walls']]
    previous = [tuple(map(K.from_json, r)) for r in previous_data['all_walls']]
    target = tuple(map(K.from_json, previous_data['geometry']['vector']))
    oc = [mv(g,r) for r in old]; oq = [dot(r,c) for r,c in zip(old,oc)]
    pc = [mv(g,r) for r in previous]; pq = [dot(r,c) for r,c in zip(previous,pc)]
    old_keys = {ray_key(r) for r in old}
    previous_keys = {ray_key(r) for r in previous}
    controls = sorted((ALPHA*x*x/q, q, x) for q in NORMS
                      for x in balanced_alphabet((height*q/ALPHA, -q/ALPHA, -q/ALPHA, -q/ALPHA))
                      if x > 0)
    out = {'status':'RUNNING','completed':False,'source':str(SOURCE.relative_to(ROOT)),
           'previous_source':str(PREVIOUS.relative_to(ROOT)), 'height_bound':height,
           'norm_pool':enc(NORMS),'control_count':len(controls), 'levels':[],
           'compatible_dictionary':[]}
    print(json.dumps({'controls':len(controls),'height':height}), flush=True)
    dictionary = {}
    for h,q,x in controls:
        energy = q+ALPHA*x*x
        blocks, stats = block_census(energy)
        scalars = {z*z:z for z in balanced_alphabet((energy,)*4)}
        counts = Counter()
        candidates = []
        for (e1,v1),(e2,v2) in product(blocks, repeat=2):
            rem = energy-e1-e2
            if rem not in scalars:
                continue
            r = (x,*v1,*v2,scalars[rem])
            assert pair(r,g,r) == q
            if all((z/T).is_integral for z in r):
                counts['imprimitive'] += 1
                continue
            assert all((2*z/q).is_integral for z in mv(g,r))
            candidates.append(r)
            mixed = any(v1) and any(v2)
            counts['primitive_candidates'] += 1
            counts['mixed_candidates'] += mixed
            if not coxeter(r,q,oc,oq):
                continue
            key = ray_key(r)
            if key in old_keys:
                continue
            already_previous = key in previous_keys
            good_previous = already_previous or coxeter(r,q,pc,pq)
            cuts = pair(r,g,target) > 0
            counts['old14_compatible'] += 1
            counts['old22_compatible'] += good_previous
            counts['cuts_old22_escape'] += cuts
            counts['mixed_old14_compatible'] += mixed
            dictionary[key] = {'height':h.to_json(),'norm':q.to_json(),'vector':enc(r),
                               'both_blocks_active':mixed,'compatible_with_old22':good_previous,
                               'already_in_old22':already_previous,
                               'old22_escape_pairing':pair(r,g,target).to_json()}
        out['levels'].append({'height':h.to_json(),'norm':q.to_json(),'time':x.to_json(),
                              **stats,**dict(counts),'candidates':menc(candidates)})
        out['compatible_dictionary'] = list(dictionary.values())
        output.write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({'height':round(h.real(),6),'blocks':len(blocks),**dict(counts),
                          'elapsed':round(time.monotonic()-started,2)}),flush=True)
    rows = sorted(dictionary.values(),key=lambda row:(K.from_json(row['height']),tuple(map(K.from_json,row['vector']))))
    walls, cov, norms = list(old), list(oc), list(oq)
    chosen, rejected = [], []
    for i,row in enumerate(rows):
        r,q = tuple(map(K.from_json,row['vector'])),K.from_json(row['norm'])
        if coxeter(r,q,cov,norms):
            chosen.append(i);walls.append(r);cov.append(mv(g,r));norms.append(q)
        else:
            rejected.append(i)
    inside = tuple(map(K.from_json,old_data['strict_interior']))
    while any(dot(inside,c)>=0 for c in cov):
        inside=(inside[0],*(z/2 for z in inside[1:]))
    assert pair(inside,g,inside)<0
    assert all(coxeter(r,norms[i],cov[:i],norms[:i]) for i,r in enumerate(walls))
    geom=geometry(walls,g,inside,9)
    out.update(status='COMPLETE_BOUNDED_TWO_F4_CENSUS',completed=True,
               compatible_dictionary=rows,selected_indices=chosen,rejected_indices=rejected,
               all_walls=menc(walls),wall_count=len(walls),strict_interior=enc(inside),geometry=geom,
               seconds=round(time.monotonic()-started,3),
               scope='Complete dominant primitive-root enumeration in both F4 blocks and scalar through the stated bound, using all five necessary norm classes. Dictionary retains every root compatible with the original 14-wall prefix, regardless of target escape. Height-ordered selection checks all Coxeter pairs. Full-host reflectivity remains unresolved unless separately proved.')
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'dictionary':len(rows),'selected':len(chosen),
                      'walls':len(walls),'geometry':geom['status'],'seconds':out['seconds']}),flush=True)


if __name__ == '__main__':
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('--height',type=int,default=300)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    run(a.height,a.output)
