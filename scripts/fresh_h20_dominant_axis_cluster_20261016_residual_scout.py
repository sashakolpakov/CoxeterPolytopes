#!/usr/bin/env python3
"""Numerical whole-region score for the exact compatible two-axis cluster."""
import argparse,hashlib,json
from pathlib import Path
from fractions import Fraction as F
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scipy.stats import qmc,norm
from fresh_h20_triangle_adjacent_h3_20261016 import geometry,close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_dominant_axis_cluster_20261016_residual_scout'
def load(p):return json.loads((ROOT/p).read_text())
def mid(x):return sum(float(F(t))for t in x)/2

def run():
    G=geometry();P='results/fresh_h20_dominant_axis_cluster_20261016.json';D='results/fresh_h20_dominant142_cap_20261016.json'
    p=load(P);d=load(D);assert p['status']=='PASS_EXACT_156_WALL_COORDINATED_DOMINANT_AXIS_CLUSTER'
    a=mid(d['definition']['axis_cutoff_bounds']);b=mid(p['definition']['axis_cutoff_bounds']);oldsrc='results/fresh_h20_dominant142_cap_20261016_residual_scout.json';old=load(oldsrc)['new154']
    A,bb0=G['A'],G['b'];AA=np.r_[A,np.eye(4)[[1,0,2]]];bb=np.r_[bb0,a,b,b];labels=G['labels']+[153,154,155]
    vv=[]
    for z in HalfspaceIntersection(np.c_[AA,-bb],np.full(4,1e-5)).intersections:
        if any(np.linalg.norm(z-w['z'])<1e-7 for w in vv):continue
        act=[labels[i]for i in range(len(bb))if abs(AA[i]@z-bb[i])<1e-7]
        vv.append(dict(z=z.tolist(),squared_radius=float(z@z),active=sorted(set(act))))
    vv.sort(key=lambda v:(v['active'],v['z']))
    counts=dict(total=len(vv),inside=sum(v['squared_radius']<1-1e-6 for v in vv),near_sphere=sum(abs(v['squared_radius']-1)<=1e-6 for v in vv),outside=sum(v['squared_radius']>1+1e-6 for v in vv))
    assert counts==dict(total=176,inside=22,near_sphere=19,outside=135)
    cut=[v for v in old['vertices']if max(v['z'][0],v['z'][2])>b+1e-8];new=[v for v in vv if set(v['active'])&{154,155}]
    assert len(cut)==8 and all(v['squared_radius']>1+1e-6 for v in cut)
    assert len(new)==16 and sum(v['squared_radius']<1-1e-6 for v in new)==4
    rows=np.where(bb0>1e-9)[0];An=A[rows]/bb0[rows,None]
    U=qmc.Sobol(4,scramble=True,seed=20261016).random_base2(20);Z=np.abs(norm.ppf(np.clip(U,1e-12,1-1e-12)));Z/=np.linalg.norm(Z,axis=1)[:,None];assert np.all(np.isfinite(Z))
    rr=[[],[],[]]
    for j in range(0,len(Z),8192):
        V=Z[j:j+8192];den=np.max(V@An.T,axis=1);assert min(den)>0
        r0=1/den;r1=np.minimum(r0,a/V[:,1]);r2=np.minimum(r1,np.minimum(b/V[:,0],b/V[:,2]))
        for bank,x in zip(rr,[r0,r1,r2]):bank.extend(x.tolist())
    rr=list(map(np.array,rr));scores=[]
    for power in [18,20]:
        n=1<<power;row=dict(power=power)
        for name,r in zip(['old153','old154','new156'],rr):
            row[name]=dict(exterior_angular_fraction=float(np.mean(r[:n]>1)),exterior_euclidean_volume=float(np.pi**2/32*np.mean(np.maximum(r[:n]**4-1,0))))
        row['relative_volume_reduction_from153']=1-row['new156']['exterior_euclidean_volume']/row['old153']['exterior_euclidean_volume'];scores.append(row)
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_triangle_adjacent_h3_20261016.py',P,D,oldsrc,'results/fresh_h20_triangle_adjacent_h3_20261016.json']
    return dict(status='PASS_BOUNDED_NUMERICAL_COORDINATED_AXIS_REGION_SCORE',sources={x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest()for x in src},axis_cutoffs=[b,a,b],
      old154=dict(counts=old['counts'],maximum_squared_radius=old['maximum_squared_radius']),
      new156=dict(counts=counts,maximum_squared_radius=max(v['squared_radius']for v in vv),vertices=vv),
      old_exterior_vertices_cut=cut,new_facet_vertices=new,common_direction_scores=scores,
      method='Same fixed-seed scrambled2^20 Sobol inverse-normal sphere directions and fixed C chart at153,154,156. Exterior4-volume=(pi^2/32)*mean(max(radius^4-1,0)); Qhull sphere band1e-6.',
      scope='Numerical whole-region diagnostics only; no exact exterior/ideal census or rigorous integration error bounds. Exact new wall/corner checks and surviving escape are separate.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['new156']['counts'],d['common_direction_scores'][-1],flush=True)
