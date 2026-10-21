#!/usr/bin/env python3
"""Finite-event scout coordinated by the exact171 open-edge incidence ledger.

This is a bounded numerical scout, not all-label or continuous exhaustion.
"""
import argparse,hashlib,itertools,json,math
from collections import Counter
import numpy as np
import fresh_h20_replaced_old_edge_profiles_20261017 as O
import fresh_h20_replacement_boundary_20261017 as L
ROOT=O.ROOT;TAG='fresh_h20_replacement_outer_patch_20261017';INPUT=f'results/{L.TAG}.json'
FOCUS=[[95,157],[145,157],[145,152],[95,171],[152,172],[171,172],[161,171]]

def run():
 D=O.P.P.load(INPUT);GG=D['geometry'];Y=np.array(GG['residual_normals']);CT=np.array(GG['T_products']);J=np.array(GG['T_inverse']);q=np.array(GG['native_norms']);m=np.array(GG['residual_form']);T=D['support'];G=dict(T=T,m=m);base=D['stages'][-1];ret=base['retained_ids'];sides=base['actual_facet_ids'];V=np.array([[1.]+v['point']for v in base['poly_vertices']if v['kind']=='finite']);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
 old=O.P.P.load('results/fresh_h20_replaced_old_edge_profiles_20261017.json');cache={tuple(z['edge']):z for z in old['endpoint_censuses']};full=Y@m@Y.T+CT@J@CT.T;censuses=[]
 for edge in FOCUS:
  if tuple(edge)in cache:cen=cache[tuple(edge)]
  else:W=T+edge;cen=O.census(W,full[np.ix_(W,W)]);cen['edge']=edge
  censuses.append(cen);print('CENSUS',edge,len(cen['profiles']),flush=True)
 stats=Counter();candidates=[];trials=[];near=[];openpairs=[z['pair']for z in base['open_edges']]
 for ci,cen in enumerate(censuses):
  rows=list(cen['profiles'])
  for j in cen['parametric_singletons']:
   rows.extend(dict(id=f'single{j}_{lab}',contacts=[[j,lab]])for lab in range(3,13))
  for row in rows:
   for sg,p,a,b,dt,offset in O.conic(cen,row,G,Y,CT,J,q):
    aa=a@m@Y.T;bb=b@m@Y.T;cc=p@m@Y.T+offset;hh=2-dt@J@dt
    for anchor in sides:
     if anchor in cen['edge']:continue
     for lab in range(2,7):
      stats['equations']+=1;target=-math.sqrt(2*q[anchor])*(0 if lab==2 else math.cos(math.pi/lab))
      for k in O.B.roots(aa[anchor],cc[anchor]-target,bb[anchor]):
       stats['real_branches']+=1;r=p+k*a+b/k;ps=r@m@Y.T+offset;co=-ps/np.sqrt(2*q)
       if max(r@m@V.T)>=-1e-8:stats['old_corner_violation']+=1;continue
       if hh*(Y[anchor]@m@Y[anchor])-(r@m@Y[anchor])**2<=1e-8:stats['anchor_misses_target_facet']+=1;continue
       labs=[O.R.classify(float(x))for x in co];bad=[j for j in ret if labs[j]in('obtuse','angle_gap')];trials.append([ci,row['id'],sg,anchor,lab,float(k),bad]);stats['actual_anchor_events']+=1
       rec=dict(census=ci,edge=cen['edge'],profile=row['id'],contacts=row['contacts'],branch=sg,anchor=anchor,label=lab,k=float(k),residual4=r.tolist(),T_pairings=dt.tolist(),failed_gates=[dict(wall=j,kind=labs[j],cosine=float(co[j]))for j in bad])
       if bad:
        if len(bad)<=2:near.append(rec)
        continue
       if any(np.linalg.norm(r-z['residual4'])<1e-6 and np.linalg.norm(dt-z['T_pairings'])<1e-6 for z in candidates):stats['duplicate']+=1;continue
       closed=[];corners=[]
       for edge in openpairs:
        W=T+edge;g=np.zeros((20,20));g[:19,:19]=full[np.ix_(W,W)];g[19,:19]=g[:19,19]=ps[W];g[19,19]=2
        if min(np.linalg.eigvalsh(g))<=1e-8:continue
        mat=Y[edge+[0]][:0] # empty placeholder avoids accidental native normalization.
        ns=np.r_[Y[edge],[r]];xx=np.linalg.solve(ns[:,1:],ns[:,0]);rr=float(xx@xx)
        if rr>=1-1e-7:continue
        if max(Y[ret,1:]@xx-Y[ret,0])>1e-6:continue
        closed.append(edge);corners.append(dict(active=edge+[173],point=xx.tolist(),radius_squared=rr))
       rec.update(id=len(candidates),all_retained_labels=[labs[j]for j in ret],closes_old_open_edges=closed,new_endpoint_corners=corners)
       candidates.append(rec);stats['accepted_unique']+=1
  print('SOLVED',ci,'candidates',len(candidates),dict(stats),flush=True)
 # Exact-angle targets define only numerical discovery here. Pairwise
 # compatibility and actual finite new-neighbor incidence are both recorded.
 pairs=[]
 for a,b in itertools.combinations(candidates,2):
  ra,rb=np.array(a['residual4']),np.array(b['residual4']);da,db=np.array(a['T_pairings']),np.array(b['T_pairings']);co=-(ra@m@rb+da@J@db)/2;lab=O.R.classify(float(co))
  if lab in('obtuse','angle_gap'):continue
  det=float((ra@m@ra)*(rb@m@rb)-(ra@m@rb)**2);pairs.append(dict(pair=[a['id'],b['id']],ambient_label=lab,intrinsic_finite=det>1e-8,projected_determinant=det))
 paths=[f'scripts/{TAG}.py',INPUT,f'scripts/{O.TAG}.py','results/fresh_h20_replaced_old_edge_profiles_20261017.json']
 return dict(status='BOUNDED_NUMERICAL_OPEN_PATCH_FINITE_EVENT_SCOUT',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},retained_ids=ret,focus_edges=FOCUS,censuses=censuses,statistics=dict(stats),candidates=candidates,compatible_pairs=pairs,at_most_two_bad_gates=near,event_trials=trials,scope='Fixed spherical endpoint profiles plus isolated singleton labels3..12, both conic branches, actual-old-facet anchor labels2..6. Every candidate is tested against all171 retained native gates and all17 finite target vertices. Numerical compatibility only; no all-label/continuous exclusion, exact promotion, or finite-volume claim. Newly closed old edges are recorded before any global replacement-boundary claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],'pairs',len(d['compatible_pairs']),flush=True)
