#!/usr/bin/env python3
"""Bounded numerical residual-face comparison151 versus optional152.

Floating-point Qhull vertices and near-sphere tolerance are diagnostics,
not exact exterior/ideal counts or a volume certificate.
"""
import json,numpy as np
from fractions import Fraction as F
from scipy.optimize import linprog
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(tag):return json.loads((ROOT/'results'/f'{tag}.json').read_text())
def mid(x):return sum(float(F(t))for t in x)/2

import argparse,hashlib
from scipy.spatial import HalfspaceIntersection
TAG='fresh_h20_post_cusp_escape_cap_20261016_residual_scout'
def run():
    D=load('fresh_h20_residual_h4_corner_20261016');M=load('fresh_h20_mixed_isolated_cap_20261016');P=load('fresh_h20_h4_parallel_partner_20261016');T=load('fresh_h20_second_cusp_completion_20261016')
    R=np.array([[float(F(x))for x in row]for row in D['definition']['four_recipes']]);q=np.array([float(F(x))for x in D['definition']['four_norms']]);C=np.array([float(F(x))for x in D['new_finite_corner']['chart']]);k=-float(F(D['new_finite_corner']['norm']));g=np.diag([-1,1,1,.5,.5])
    cp=np.array([mid(x)for x in D['new_finite_corner']['all145_pair_bounds']]);pp=np.zeros((145,4))
    for row in M['pairs']:pp[row['wall'],0]=mid(row['pair_bounds'])
    for row in D['pairs']:
     a,b=row['walls'];pp[a,b-141]=mid(row['pair_bounds'])
    for i in range(4):pp[141+i,i]=q[i]
    N=np.array([sum((pp[j,i]/q[i])*R[i]for i in range(4))-cp[j]*C/k for j in range(145)])
    assert np.max(np.abs(N[141:]-R))<1e-8
    nu=C/np.sqrt(k)-R[3]/np.sqrt(q[3]);z=mid(P['definition']['z_bounds']);r=-np.sqrt(2/q[0])*R[0]+np.sqrt(2)*z*nu
    N=np.r_[N,[r],[3*nu],[3*nu],[3*nu],[-R[1]+5*nu],[-R[2]+6*nu]]
    norms=np.r_[np.full(140,2.),1.46,q,[2,2,2,2,q[1],q[2]]];assert len(N)==len(norms)==151
    E5=np.array([1,-.01,-.98,-.3,-.04]);E4=E5.copy();E4[-1]=-.05
    S16=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,61,87]

    def census(normals):
        B=-(R@g@normals.T).T*np.sqrt(k/q);bb=-normals@g@C
        active=[i for i in range(len(normals))if i not in S16 and i not in range(141,145)]
        assert min(bb[active])>0
        half=np.r_[np.c_[B[active]/bb[active,None],-np.ones(len(active))],np.c_[-np.eye(4),np.zeros(4)]]
        inner=np.full(4,1e-5);assert max(half[:,:4]@inner+half[:,4])<0
        vv=[]
        for zz in HalfspaceIntersection(half,inner).intersections:
            act=[i for i in active if abs(B[i]@zz-bb[i])<1e-8]+[141+i for i,x in enumerate(zz)if abs(x)<1e-7]
            vv.append(dict(z=zz.tolist(),squared_radius=float(zz@zz),active=sorted(set(act))))
        vv.sort(key=lambda x:(x['active'],x['z']))
        counts=dict(total=len(vv),inside=sum(x['squared_radius']<1-1e-6 for x in vv),near_sphere=sum(abs(x['squared_radius']-1)<=1e-6 for x in vv),outside=sum(x['squared_radius']>1+1e-6 for x in vv))
        return vv,counts
    old,oc=census(N)
    recipe=[F(60),F(-3,2),F(-59),F(-15),F(-3)]
    rr=np.array(recipe,float);new,nc=census(np.r_[N,[rr]])
    scores=[]
    for v in old:
        point=C-(np.array(v['z'])*np.sqrt(k/q))@R
        pair=float(point@g@rr)
        scores.append(dict(active=v['active'],squared_radius=v['squared_radius'],new_cap_pair=pair,cut=pair>1e-10))
    nf=[v for v in new if 151 in v['active']]
    assert oc==dict(total=150,inside=13,near_sphere=19,outside=118)
    assert nc==dict(total=156,inside=13,near_sphere=19,outside=124)
    assert sum(v['cut']for v in scores)==6 and all(v['squared_radius']>1+1e-6 for v in scores if v['cut'])
    assert len(nf)==12 and all(v['squared_radius']>1+1e-6 for v in nf)
    src=[f'scripts/{TAG}.py']+[f'results/{x}.json'for x in ['fresh_h20_residual_h4_corner_20261016','fresh_h20_mixed_isolated_cap_20261016','fresh_h20_h4_parallel_partner_20261016','fresh_h20_second_cusp_completion_20261016','fresh_h20_post_cusp_escape_cap_20261016']]
    return dict(status='PASS_BOUNDED_NUMERICAL_RESIDUAL_FACE_COMPARISON',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},
       method='Projected all151 inequalities into orthogonal residualC+4caps;scipyQhullin4D;near-sphere threshold1e-6;boundmidpointsonly',
       old151=dict(counts=oc,vertices=old),optional152=dict(counts=nc,vertices=new),old_vertex_cut_scores=scores,
       new_facet_vertices=nf,cut_exterior_vertices=6,new_exterior_vertices=12,
       scope='Numericaldiagnostic,notexactconecompleteness/idealclassification/volumeproof. Optional152cutsE5regionbutincreasesexteriorvertexcount;151remainsdefaultbaseforcoordinatedsearch.')
def close(a,b):
    if isinstance(a,float):return isinstance(b,(int,float))and abs(a-b)<=1e-7*(1+abs(a))
    if isinstance(a,dict):return a.keys()==b.keys()and all(close(v,b[k])for k,v in a.items())
    if isinstance(a,list):return len(a)==len(b)and all(close(x,y)for x,y in zip(a,b))
    return a==b
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();data=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert close(data,json.loads(p.read_text()))
    else:p.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['old151']['counts'],data['optional152']['counts'],flush=True)
