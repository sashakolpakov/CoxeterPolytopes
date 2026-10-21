#!/usr/bin/env python3
"""Numerical fixed-chart coverage of the exact dominant-axis cap.
Qhull counts and Sobol angular integrals are diagnostics, not exact cone proofs.
"""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scipy.stats import qmc,norm
from fresh_h20_triangle_adjacent_h3_20261016 import geometry,close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_dominant142_cap_20261016_residual_scout'
def load(p):return json.loads((ROOT/p).read_text())
def mid(x):return sum(float(F(t))for t in x)/2

def run():
    G=geometry();P='results/fresh_h20_dominant142_cap_20261016.json';p=load(P)
    assert p['status']=='PASS_EXACT_154_WALL_DOMINANT_AXIS_PI5_CAP'
    cutoff=mid(p['definition']['axis_cutoff_bounds']);oldsrc='results/fresh_h20_pi8_triangle_face_20261016_residual_scout.json';old=load(oldsrc)['new153']
    A,b=G['A'],G['b'];AA=np.r_[A,[[0.,1.,0.,0.]]];bb=np.r_[b,cutoff];labels=G['labels']+[153]
    vv=[]
    for z in HalfspaceIntersection(np.c_[AA,-bb],np.full(4,1e-5)).intersections:
        if any(np.linalg.norm(z-w['z'])<1e-7 for w in vv):continue
        act=[labels[i]for i in range(len(bb))if abs(AA[i]@z-bb[i])<1e-7]
        vv.append(dict(z=z.tolist(),squared_radius=float(z@z),active=sorted(set(act))))
    vv.sort(key=lambda v:(v['active'],v['z']))
    counts=dict(total=len(vv),inside=sum(v['squared_radius']<1-1e-6 for v in vv),near_sphere=sum(abs(v['squared_radius']-1)<=1e-6 for v in vv),outside=sum(v['squared_radius']>1+1e-6 for v in vv))
    assert counts==dict(total=168,inside=18,near_sphere=19,outside=131)
    cut=[v for v in old['vertices']if v['z'][1]>cutoff+1e-8];new=[v for v in vv if 153 in v['active']]
    assert len(cut)==7 and all(v['squared_radius']>1+1e-6 for v in cut)
    assert len(new)==12 and sum(v['squared_radius']<1-1e-6 for v in new)==2
    # Normalize non-coordinate inequalities; radial support of P in its positive orthant.
    rows=np.where(b>1e-9)[0];An=A[rows]/b[rows,None]
    U=qmc.Sobol(4,scramble=True,seed=20261016).random_base2(20);Z=np.abs(norm.ppf(np.clip(U,1e-12,1-1e-12)));Z/=np.linalg.norm(Z,axis=1)[:,None]
    assert np.all(np.isfinite(Z))
    oldrad=[];newrad=[]
    for j in range(0,len(Z),8192):
        V=Z[j:j+8192];den=np.max(V@An.T,axis=1);assert min(den)>0
        r=1/den;oldrad.extend(r.tolist());newrad.extend(np.minimum(r,cutoff/V[:,1]).tolist())
    oldrad=np.array(oldrad);newrad=np.array(newrad);scores=[]
    for power in [18,20]:
        n=1<<power;v0=np.pi**2/32*np.mean(np.maximum(oldrad[:n]**4-1,0));v1=np.pi**2/32*np.mean(np.maximum(newrad[:n]**4-1,0))
        scores.append(dict(power=power,old_exterior_angular_fraction=float(np.mean(oldrad[:n]>1)),new_exterior_angular_fraction=float(np.mean(newrad[:n]>1)),old_exterior_euclidean_volume=float(v0),new_exterior_euclidean_volume=float(v1),estimated_relative_volume_reduction=float(1-v1/v0)))
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_triangle_adjacent_h3_20261016.py',P,oldsrc,'results/fresh_h20_triangle_adjacent_h3_20261016.json']
    return dict(status='PASS_BOUNDED_NUMERICAL_DOMINANT_AXIS_REGION_SCORE',sources={x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest()for x in src},axis_cutoff=cutoff,
      old153=dict(counts=old['counts'],maximum_squared_radius=max(v['squared_radius']for v in old['vertices'])),
      new154=dict(counts=counts,maximum_squared_radius=max(v['squared_radius']for v in vv),vertices=vv),
      old_exterior_vertices_cut=cut,new_facet_vertices=new,common_direction_scores=scores,
      method='Positive-orthant sphere directions from fixed-seed scrambled Sobol inverse-normal samples. Radial exterior4-volume=(pi^2/32)*mean(max(radius^4-1,0)); same directions and chart for both chambers. Qhull and 1e-6 sphere band.',
      scope='Numerical diagnostics only: no exact ideal/exterior census, volume value or finite-volume assertion. Primary exact gates and surviving E6 are separate.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['new154']['counts'],d['common_direction_scores'][-1],flush=True)
