#!/usr/bin/env python3
"""Numerical positive-construction scout using enlarged ambient mirror orbits.

The original chamber is replaced: only its valid projected cusp walls are
fixed. This is not an exact or exhaustive mirror-bank exclusion.
"""
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_quartic_orbit_completion_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
LOCAL='results/fresh_h20_rank3_omitted_sections_20261014.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'

def sha(p):
    return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def mpr(q):
    q=Fraction(q)
    return mp.mpf(q.numerator)/q.denominator

def run():
    mp.mp.dps=90
    raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
    Q=-np.array(raw['basis_gram'],dtype=np.int64)
    W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
    G=W@Q@W.T; d=np.diag(G)
    h=W.sum(axis=0);h=h//np.gcd.reduce(h);hq=Q@h
    provenance={tuple(r):(i,-1) for i,r in enumerate(W)}
    for generator,r in enumerate(W):
        assert np.all((2*G[generator])%d[generator]==0)
        reflected=W-(2*G[generator]//d[generator])[:,None]*r
        for original,v in enumerate(reflected):
            if v@hq>0:
                v=-v
            provenance.setdefault(tuple(v),(original,generator))
    bank=np.array(sorted(provenance),dtype=np.int64)
    norms=np.einsum('ij,ij->i',bank@Q,bank)
    assert len(bank)==29043 and set(norms)=={2,4}
    recipes=np.array([provenance[tuple(v)] for v in bank],dtype=np.int64)
    bankfile=f'data/enumerations/{TAG}_bank.npz'
    np.savez_compressed(ROOT/bankfile,roots=bank,source_reflection=recipes)
    walk=json.loads((ROOT/WALK).read_text())['corner_orbits']
    local=json.loads((ROOT/LOCAL).read_text())
    # cos²(pi/m) in a quartic field has degree dividing four. This is a
    # necessary superset of the actual field's cyclotomic subfields.
    labels=np.array([2,3,4,5,6,8,10,12,15,16,20,24,30])
    targets=np.cos(np.pi/labels)**2;targets[0]=0
    params=['1/10','1/4','1/2','3/4','9/10']
    output=[]; counts=Counter()
    for index,row in enumerate(local['directions']):
        if len(row['minimal_polynomial'])!=5:
            continue
        spec=local['systems'][row['system']]
        poly=list(map(mpr,row['minimal_polynomial']))
        ab=list(map(mpr,row['root_interval']))
        alpha=mp.findroot(lambda z:mp.polyval(poly,z),sum(ab)/2,
                          tol=mp.mpf('1e-75'),maxsteps=100)
        assert ab[0]<=alpha<=ab[1]
        ev=lambda cs:mp.polyval(list(map(mpr,cs)),alpha)
        e0=[ev(p) for p in row['normal_constant']]
        shift=list(map(mpr,row['null_shift']))
        lower=ev(row['parameter_lower_magnitude'])
        seed=[i for i in walk[spec['cusp']]['active_walls'] if i!=spec['omitted_wall']]
        delta=np.array(list(map(float,walk[spec['cusp']]['primitive_ray'])))
        cusp_pair=bank@(Q@delta)
        orientation=np.where(cusp_pair>0,-1,1)
        oriented=bank*orientation[:,None]
        caps=cusp_pair!=0
        pair=oriented@Q@W[seed].T
        seed_norms=d[seed].astype(float)
        for uu in params:
            u=mpr(uu)
            en=[a-u*lower*b for a,b in zip(e0,shift)]
            e=np.array(list(map(float,en)))
            q=float(ev(row['normal_square']))
            assert abs(e@Q@e-q)<1e-8*max(1,abs(q))
            se=W[seed]@Q@e;be=oriented@Q@e
            sn=seed_norms-se*se/q
            bn=norms-be*be/q
            assert np.min(sn)>1e-8
            pp=pair-be[:,None]*se[None,:]/q
            indices=np.flatnonzero(caps & (bn>1e-8) & (np.max(pp,axis=1)<1e-8))
            compatible=[]; borderline=[]
            for ix in indices:
                c=pp[ix]**2/(bn[ix]*sn)
                finite=c<1-1e-8
                diffs=np.min(np.abs(c[finite,None]-targets[None,:]),axis=1)
                if len(diffs)==0 or np.max(diffs)<1e-8:
                    compatible.append(int(ix))
                elif np.max(diffs)<1e-5:
                    borderline.append(int(ix))
            counts['tested_sections']+=1
            counts['seed_acute_cap_proposals']+=len(indices)
            counts['Coxeter_compatible_cap_proposals']+=len(compatible)
            counts['near_angle_remainders']+=len(borderline)
            output.append({'local_record':index,'parameter_u':uu,'cusp':spec['cusp'],
                           'seed_walls':seed,'positive_acute_caps':len(indices),
                           'compatible_cap_indices':compatible,
                           'near_angle_cap_indices':borderline})
        counts['quartic_direction_records']+=1
        print('direction',index,'done; cap proposals',counts['Coxeter_compatible_cap_proposals'],flush=True)
    sources=[RAW,LOCAL,WALK,bankfile,f'scripts/{TAG}.py']
    out={'status':'COMPLETED_NUMERICAL_ENLARGED_MIRROR_CONSTRUCTION_SCOUT',
         'sources':{p:sha(p) for p in sources},'bank_size':len(bank),
         'bank_root_norm_counts':dict(Counter(map(int,norms))),
         'tested_parameters':params,'angle_label_superset':list(map(int,labels)),
         'counts':dict(counts),'sections':output,
         'scope':'Positive construction attempt retaining the known Coxeter cusp link while allowing other facets from original roots and one ambient simple reflection. All44 quartic local records at five interior parameters. Numerical compatibility only; no exhaustive or exact exclusion of these normals, other parameters, or the full ambient reflection orbit.'}
    (ROOT/'results'/f'{TAG}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status']);print(out['counts'])

if __name__=='__main__':
    run()
