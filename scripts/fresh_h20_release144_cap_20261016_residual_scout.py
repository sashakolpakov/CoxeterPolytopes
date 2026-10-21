#!/usr/bin/env python3
"""Bounded numerical comparison of three release144 recipes; exact wall separate."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.spatial import HalfspaceIntersection
from scipy.stats import qmc,norm
from fresh_h20_triangle_adjacent_h3_20261016 import geometry,mid,close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_release144_cap_20261016_residual_scout'
P='results/fresh_h20_dominant142_cap_20261016.json'
Q='results/fresh_h20_dominant_axis_cluster_20261016.json'
EXACT='results/fresh_h20_release144_cap_20261016.json'
def load(p):return json.loads((ROOT/p).read_text())
def geometry156():
    G=geometry();N,g,q,R,C,k,H=map(G.get,['N','g','q','R','C','k','H']);p,z=load(P),load(Q)
    s=mid(p['definition']['s_bounds']);c=mid(p['definition']['c_bounds']);t=mid(z['definition']['s_star_bounds']);ct=mid(z['definition']['c_star_bounds'])
    e0=C/np.sqrt(k);ee=R/np.sqrt(q)[:,None]
    N=np.r_[N,[np.sqrt(2)*(-c*ee[1]+s*e0),np.sqrt(2)*(-ct*ee[0]+t*e0),np.sqrt(2)*(-ct*ee[2]+t*e0)]]
    G.update(N=N,norms=np.r_[G['norms'],2,2,2],A=np.r_[G['A'],np.eye(4)[[1,0,2]]],b=np.r_[G['b'],s/c,t/ct,t/ct],labels=G['labels']+[153,154,155],e0=e0,ee=ee,boosts=np.array([t,s,t]))
    return G
def gates(G,r):
    co=-(G['N']@G['g']@r)/np.sqrt(2*G['norms']);bad=[]
    for j,z in enumerate(co):
        if abs(z)<1e-6 or z>=1-1e-7:continue
        if z<0:bad.append(dict(wall=j,kind='obtuse',cosine=float(z)));continue
        m=np.pi/np.arccos(z)
        if m<2.99999 or abs(m-round(m))>1e-4:bad.append(dict(wall=j,kind='nonintegral_label',cosine=float(z),label=float(m)))
    return bad
def constraint(G,r):
    den=-G['C']@G['g']@r;assert den>0
    return -(G['R']@G['g']@r)*np.sqrt(G['k']/G['q'])/den
def census(A,b,ids):
    vv=[]
    for z in HalfspaceIntersection(np.c_[A,-b],np.full(4,1e-5)).intersections:
        if any(np.linalg.norm(z-np.array(v['z']))<1e-7 for v in vv):continue
        vv.append(dict(z=z.tolist(),squared_radius=float(z@z),active=sorted(set(ids[i]for i in range(len(b))if abs(A[i]@z-b[i])<1e-7))))
    vv.sort(key=lambda v:(v['active'],v['z']))
    return dict(counts=dict(total=len(vv),inside=sum(v['squared_radius']<1-1e-6 for v in vv),near_sphere=sum(abs(v['squared_radius']-1)<=1e-6 for v in vv),outside=sum(v['squared_radius']>1+1e-6 for v in vv)),maximum_squared_radius=max(v['squared_radius']for v in vv),vertices=vv)
def directions(power=20):
    Z=qmc.Sobol(4,scramble=True,seed=20261016).random_base2(power);Z=np.abs(norm.ppf(np.clip(Z,1e-12,1-1e-12)));Z/=np.linalg.norm(Z,axis=1)[:,None];assert np.all(np.isfinite(Z));return Z
def radial(A,b,Z):
    ix=np.where(b>1e-9)[0];An=A[ix]/b[ix,None];rr=[]
    for j in range(0,len(Z),8192):
        den=np.max(Z[j:j+8192]@An.T,axis=1);assert min(den)>0;rr.extend(1/den)
    return np.array(rr)
def score(rr):return dict(exterior_angular_fraction=float(np.mean(rr>1)),exterior_euclidean_volume=float(np.pi**2/32*np.mean(np.maximum(rr**4-1,0))))
def recipe(G,m):
    N,g=G['N'],G['g'];U=N[[95,142,145,151]];inv=np.linalg.inv(U@g@U.T)
    f=np.array([0,-np.sqrt(2*G['norms'][142])*np.cos(np.pi/m),-1,0]);p=(inv@f)@U
    seed=np.array([0.,0,0,0,1]);k=seed-(inv@(U@g@seed))@U;assert k@g@k>0
    lam=np.sqrt((2-p@g@p)/(k@g@k));E6=np.array([1,-.01,-.98,-.28,-.05]);success=[]
    for sign in [-1,1]:
        r=p+sign*lam*k
        if not gates(G,r)and r@g@G['H']<0 and r@g@E6>0:success.append(r)
    assert len(success)==1;return success[0]
def run():
    G=geometry156();doc=load(EXACT);assert doc['status']=='PASS_EXACT_157_WALL_RELEASE144_E6_CUTTER'
    A,b,ids=G['A'],G['b'],G['labels'];Z=directions();rr=radial(A,b,Z);rows=[]
    E6=np.array([1,-.01,-.98,-.28,-.05]);E7=G['C']-(np.array([.7,.4,.7,.4])*np.sqrt(G['k']/G['q']))@G['R']
    for m in [3,4,5]:
        r=recipe(G,m);ar=constraint(G,r);vv=census(np.r_[A,[ar]],np.r_[b,1],ids+[156]);rnew=np.minimum(rr,1/np.maximum(Z@ar,1e-30))
        rows.append(dict(labels=[2,m,3,2],contacts=[95,142,145,151],chart=r.tolist(),norm=float(r@G['g']@r),failed_gates=gates(G,r),E6_pair=float(r@G['g']@E6),E7_pair=float(r@G['g']@E7),residual=vv,score=score(rnew)))
    assert rows[0]['residual']['counts']==dict(total=181,inside=24,near_sphere=19,outside=138)
    assert rows[0]['score']['exterior_euclidean_volume']<min(x['score']['exterior_euclidean_volume']for x in rows[1:])
    face=[v for v in rows[0]['residual']['vertices']if {151,156}<=set(v['active'])];assert len(face)==5
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_triangle_adjacent_h3_20261016.py',P,Q,EXACT,'results/fresh_h20_triangle_adjacent_h3_20261016.json']
    return dict(status='PASS_BOUNDED_RELEASE144_THREE_RECIPE_REGION_COMPARISON',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},old156=dict(residual=census(A,b,ids),score=score(rr)),candidates=rows,selected_candidate=0,selected_new_H2_face=dict(support='S16+{151,156}',vertices=face),method='Fixed C chart; all retained residual inequalities. Same fixed-seed scrambled2^20 Sobol directions; numerical Qhull sphere band1e-6.',scope='Three prescribed recipes only. No all-cap completeness, exact exterior/ideal census or rigorous quadrature error bound. Only m3 is promoted by the separate exact certificate.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert close(d,json.loads(p.read_text()))
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['candidates'][0]['residual']['counts'],flush=True)
