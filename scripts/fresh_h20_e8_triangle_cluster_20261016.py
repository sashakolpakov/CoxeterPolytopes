#!/usr/bin/env python3
"""Bounded numerical E8-facet triangle construction, full ambient old gates."""
import argparse,hashlib,json
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection
import fresh_h20_spherical_attachment_profiles_20261016_scout as S
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_triangle_cluster_20261016'
BASE='results/fresh_h20_spherical_attachment_profiles_20261016_scout.json'
def load(p):return json.loads((ROOT/p).read_text())
def census(G,Y):
 A=np.r_[G['A'],Y[:,1:]];b=np.r_[G['b'],Y[:,0]];nn=np.linalg.norm(A,axis=1)
 lp=linprog([0,0,0,-1],A_ub=np.c_[A/nn[:,None],np.ones(len(b))],b_ub=b/nn,bounds=[(None,None)]*3+[(0,None)],method='highs')
 if not lp.success or lp.x[3]<1e-9:return dict(empty=True)
 v=HalfspaceIntersection(np.c_[A,-b],lp.x[:3]).intersections;r=np.sum(v*v,axis=1)
 return dict(vertices=len(v),finite=int(sum(r<1-1e-6)),exterior=int(sum(r>1+1e-6)),near_sphere=int(sum(abs(r-1)<=1e-6)),maximum_squared_radius=float(max(r)),center=lp.x[:3].tolist())
def run():
 D,G=S.geometry();B=load(BASE);m=G['m'];notT=[j for j in range(158)if j not in G['T']]
 # Norm2 ambient B3 on the three partners; their A10 projection is removed.
 M=np.array([[12/11,-1/11,-np.sqrt(2)],[-1/11,12/11,-1],[-np.sqrt(2),-1,2.]])
 a=-np.linalg.solve(M,np.ones(3));hn=a@M@a;assert hn<0 and min(a)>0
 t=1/np.sqrt(-hn);P=M+t*t;ev,V=np.linalg.eigh(P);assert abs(ev[0])<1e-12 and ev[1]>0
 spatial=V[:,1:]*np.sqrt(ev[1:])[None,:]
 profileids=[38,42,56];O=np.array([S.offsets(G,D['census']['profiles'][i])for i in profileids])
 assert np.max(abs(O[:,G['T']]))<3
 theta=np.arange(1440)*2*np.pi/1440;rot=np.array([[[np.cos(x),-np.sin(x)],[np.sin(x),np.cos(x)]]for x in theta])
 mults=[1,1.25,1.5,2,3,4,8,16];out=[];best=[];stats=Counter()
 for cap in B['E8_all_other_dotted_caps']:
  r0=np.array(cap['residual4']);direction=r0[1:]/np.linalg.norm(r0[1:]);h=1/30
  off=S.offsets(G,D['census']['profiles'][cap['profile']])
  for factor in mults:
   k=cap['k']*factor;r=np.sqrt(h)/2*np.r_[k-1/k,(k+1/k)*direction]
   assert S.gates(G,r[None,:],off)[0][0]
   H=np.array([1.,0,0,0]);Hf=H-(H@m@r)/h*r;ht=Hf/np.sqrt(-Hf@m@Hf)
   Z=null_space(np.array([r@m,ht@m])).T;ee,VV=np.linalg.eigh(Z@m@Z.T);assert min(ee)>0
   frame=np.r_[[ht],(VV/np.sqrt(ee)[None,:]).T@Z]
   ps=np.einsum('ij,tjk->tik',spatial,rot)
   W=t*frame[0]+np.einsum('tij,jk->tik',ps,frame[1:])
   co=-(np.einsum('tij,jk,lk->til',W,m,G['Y'])+O[None,:,:])/np.sqrt(2*G['q'])[None,None,:]
   co[:,:,G['T']]=0
   ok=(abs(co)<2e-6)|(co>=1-2e-6)
   for lab in [3,4,5,6,7,8,9,10,12,16,20,24]:ok|=abs(co-np.cos(np.pi/lab))<2e-6
   good=np.where(np.all(ok,axis=(1,2)))[0];stats['rotations']+=len(theta);stats['accepted_rotations']+=len(good)
   # Record the best all-dotted margin even when no Coxeter cluster survives.
   margins=np.min(co[:,:,notT],axis=(1,2));ii=int(np.argmax(margins))
   best.append(dict(profile=cap['profile'],direction=cap['direction'],k_factor=factor,theta_index=ii,minimum_nonsupport_cosine=float(margins[ii]),worst_pair=list(map(int,np.unravel_index(np.argmin(co[ii,:,notT].T),co[ii,:,notT].T.shape)))))
   if len(good):
    # Keep every successful grid index, but only one literal representative.
    ii=int(good[np.argmax(margins[good])])
    vals=W[ii]@m@W[ii].T;assert np.max(abs(vals-M))<1e-5
    assert np.max(abs(W[ii]@m@r))<1e-5
    out.append(dict(E8_profile=cap['profile'],direction=cap['direction'],k_factor=factor,k=float(k),accepted_theta_indices=good.tolist(),theta_index=ii,theta=float(theta[ii]),partner_profiles=profileids,residual4=np.r_[[r],W[ii]].tolist(),minimum_nonsupport_cosine=float(margins[ii])))
   print('FAMILY',cap['profile'],cap['direction'],factor,'accepted',len(good),'best',float(max(margins)),flush=True)
 best.sort(key=lambda x:-x['minimum_nonsupport_cosine'])
 # Evaluate a bounded representative per cap/direction/boost family.
 geometry=[]
 for i,rec in enumerate(out):geometry.append(dict(record=i,face=census(G,np.array(rec['residual4']))))
 src=[f'scripts/{TAG}.py',BASE,S.INPUT,f'scripts/{S.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_E8_TRIANGLE_CLUSTER_SEARCH',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},partner_projected_gram=M.tolist(),partner_full_gram=[[2,0,-float(np.sqrt(2))],[0,2,-1],[-float(np.sqrt(2)),-1,2]],triangle_center_coefficients=a.tolist(),angle_count=len(theta),boost_multipliers=mults,statistics=dict(stats),accepted=out,best_margins=best,geometry=geometry,scope='Numerical bounded 36 initial E8 rays,8 boosts,1440 rotations; all158 native ambient gates. No exhaustive parameter or exact-wall claim. Target cap facet has support A10+E8 and triangular partner types E8+A12,E8+B12,E8+A12.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,load(str(out.relative_to(ROOT))))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],flush=True)
