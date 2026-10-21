#!/usr/bin/env python3
"""Bounded discrete source-corner contacts, allowing interior causal residuals."""
import json,math,time
from itertools import combinations,product
from collections import Counter
import numpy as np
import fresh_h20_e7_causal_residuals_20261014_scout as B
def run():
 start=time.time();d,f,S,J,N,V,Q,VG,VP,VE,SP,SE,edges=B.data()
 A=f['face_facet_walls'][:13];C=N[A]@Q@N[A].T;C[np.abs(C)<1e-8]=0;JI=np.linalg.inv(C)
 weights=JI@N[A];WP=weights@Q@N.T;E=np.array([sum(float(B.F(x))*math.sqrt(k)for x,k in zip(v,[1,2,79,158,11,22,869,1738]))for v in d['strict_escape']['point']]);WE=weights@Q@E
 T=np.array([sum(float(B.F(x))*math.sqrt(k)for x,k in zip(v,[1,2,79,158,11,22,869,1738]))for v in f['face_basis'][-1]])/math.sqrt(2)
 assert abs(T@Q@T+2)<1e-8;TP=T@Q@N.T;TE=T@Q@E
 values=[1,math.sqrt(2),math.sqrt(3),2,3,4];recipes=[];rows=[]
 for size in [1,2,3]:
  for support in combinations(range(13),size):
   for choices in product(range(len(values)),repeat=size):
    row=np.zeros(13);row[list(support)]=[values[i]for i in choices];rows.append(row);recipes.append([list(support),list(choices)])
 D=np.array(rows);spaceq=np.einsum('ij,ij->i',D@JI,D);pairs=-D@WP;target=-D@WE
 profiles=[]
 for size in [1,2,3]:
  for support in combinations(range(7),size):
   c=np.zeros(7);c[list(support)]=1;profiles.append((list(support),[3]*size,c))
 for label in [4,5,6]:
  c=np.zeros(7);c[6]=2*np.cos(np.pi/label);profiles.append(([6],[label],c))
 norms=np.array(d['norms'],float);den=np.sqrt(2*norms);cosines=np.cos(np.pi/np.arange(2,65));counts=Counter();out=[];seen=set();best=[];backtracks=[]
 for pi,(support,labels,c)in enumerate(profiles):
  wc=c@J;mass=wc@c-2
  if mass<1e-8:continue
  t=np.sqrt((spaceq+mass)/2);cut=target+t*TE-wc@SE;ids=np.flatnonzero(cut>1e-7);counts['escape_cutters']+=len(ids)
  pair=pairs[ids]+t[ids,None]*TP[None,:]-(wc@SP)[None,:];acute=np.max(pair,axis=1)<1e-7;ids=ids[acute];pair=pair[acute];counts['acute_cutters']+=len(ids)
  for idx,p in zip(ids,pair):
   cos=-p/den;finite=(cos>1e-7)&(cos<1-1e-7);err=np.min(abs(cos[finite,None]-cosines[None,:]),axis=1);bad=np.flatnonzero(finite)[err>1e-7]
   if len(bad):
    if set(bad).issubset({61,75}):
     r=-wc@N[S]-D[idx]@weights+t[idx]*T
     backtracks.append(dict(E7_support=support,E7_labels=labels,source_recipe=recipes[idx],time_coefficient=float(t[idx]),normal=r.tolist(),pairings=p.tolist(),escape_pairing=float(cut[idx]),removed_walls=bad.tolist()))
    if len(best)<30 or len(bad)<best[-1]['bad_count']:
     best.append(dict(E7_support=support,E7_labels=labels,source_recipe=recipes[idx],bad_count=len(bad),bad_walls=bad.tolist(),bad_cosines=cos[bad].tolist(),escape_pairing=float(cut[idx])));best.sort(key=lambda x:x['bad_count']);best=best[:30]
    continue
   r=-wc@N[S]-D[idx]@weights+t[idx]*T;key=tuple(np.round(r,8))
   if key in seen:continue
   seen.add(key);assert abs(r@Q@r-2)<1e-4
   out.append(dict(E7_support=support,E7_labels=labels,source_recipe=recipes[idx],time_coefficient=float(t[idx]),normal=r.tolist(),pairings=p.tolist(),escape_pairing=float(cut[idx])))
  print(pi,support,labels,'acute',len(ids),'candidates',len(out),flush=True)
 paths=[B.BASE,B.FACE,f'scripts/{B.TAG}_scout.py',f'scripts/{B.TAG}_interior_scout.py']
 bp=np.array([r['normal']for r in backtracks]);pair=bp@Q@bp.T;compatible=(pair<-2+1e-7)|(np.min(abs(pair[:,:,None]/2+cosines[None,None,:]),axis=2)<1e-7)
 np.fill_diagonal(compatible,False)
 finite_edges=[[int(i),int(j),float(pair[i,j])]for i,j in zip(*np.where(np.triu(compatible&(pair>-2+1e-7),1)))]
 result=dict(status='BOUNDED_NUMERICAL_SOURCE_CORNER_CONTACT_INTERIOR_SCOUT',sources={x:B.sha(x)for x in paths},source_contact_values=values,source_support_sizes=[1,2,3],source_profile_count=len(D),E7_profile_count=len(profiles),counts=dict(counts),candidates=out,near_misses=best,delete61_75_candidates=backtracks,backtrack_compatible_pair_count=int(compatible.sum()//2),backtrack_finite_pairs=finite_edges,elapsed_seconds=time.time()-start,scope='Listed sparse finite/parallel/dotted source-corner coefficients and E7 profiles only. Norm fixes positive timelike coefficient. All reported acceptances require separate exactification; zero numerical hits is no continuous/all-label theorem.')
 (B.ROOT/f'results/{B.TAG}_interior_scout.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
 d=run();print(d['status'],len(d['candidates']))
