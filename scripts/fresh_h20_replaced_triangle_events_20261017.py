#!/usr/bin/env python3
"""Finite old-facet events after replacing only158 by the opposite-cap branch.

Numerical discovery only. Fixed endpoint profiles are inherited complete;
isolated I2 labels are explicitly bounded. Every native ambient gate is tested.
"""
import argparse,hashlib,json,math
from collections import Counter
import numpy as np
from scipy.linalg import null_space
import fresh_h20_opposite_prism_cap_20261017 as P
import fresh_h20_e8_partner_endpoints_20261017_bridges as B
import fresh_h20_prism_replacement_20261017 as R
ROOT=P.ROOT;TAG='fresh_h20_replaced_triangle_events_20261017'

def run():
 D=P.P.load(P.END);G,Y,CT,J,q,V,C=P.geometry();m=G['m'];T=G['T']
 repl=P.P.load('results/fresh_h20_prism_replacement_20261017.json')
 r0=np.array(repl['residual4'][1]);Y=np.r_[Y,[r0]];CT=np.r_[CT,np.zeros((1,17))];q=np.r_[q,2]
 retained=[j for j in range(169)if j!=158]
 pos=q-np.einsum('ij,jk,ik->i',CT,J,CT)
 oldfacets=[j for j in retained if pos[j]>1e-7 and j not in T]
 vnew=[]
 for cen in D['endpoint_censuses']:
  U=Y[cen['edge']+[168]];v=null_space(U@m)[:,0]
  assert v@m@v<0;v*=1 if v[0]>0 else -1;v/=math.sqrt(-(v@m@v));vnew.append(v)
 stats=Counter();out=[];near=[];trials=[]
 for ci,cen in enumerate(D['endpoint_censuses']):
  for row in cen['profiles']:
   for br in cen['families'][row['id']]['branches']:
    a=np.array(br['positive_k_coefficient']);b=np.array(br['inverse_k_coefficient']);p=np.array(br['particular']);dt=np.array(br['T_pairings']);offset=CT@J@dt
    aa=a@m@Y.T;bb=b@m@Y.T;cc=p@m@Y.T+offset
    for j in oldfacets:
     if j in cen['edge']:continue
     isolated=max(abs(dt))<1e-9 and max(abs(CT[j]))<1e-8
     labels=range(2,65)if isolated else range(2,6)
     for lab in labels:
      stats['equations']+=1;target=-math.sqrt(2*q[j])*(0 if lab==2 else math.cos(math.pi/lab))
      for k in B.roots(aa[j],cc[j]-target,bb[j]):
       stats['positive_roots']+=1;r=p+k*a+b/k
       if r@m@vnew[ci]>=-1e-9:stats['old_corner_wrong_side']+=1;continue
       ps=r@m@Y.T+offset;co=-ps/np.sqrt(2*q);native=[R.classify(z)for z in co]
       # The prescribed extra ambient contact must also be a genuine
       # finite intersection of target-face facets.
       h=2-dt@J@dt;proj=ps[j]-offset[j]
       if h*pos[j]-proj*proj<=1e-8:stats['extra_not_intrinsically_finite']+=1;continue
       bad=[u for u in retained if native[u]in('obtuse','angle_gap')]
       stats['intrinsic_events']+=1
       rec=dict(edge=ci,profile=row['id'],contacts=row['contacts'],branch=br['branch'],anchor=j,label=lab,k=k,residual4=r.tolist(),T_pairings=dt.tolist(),bad_wall_ids=bad,bad_gate_details=[dict(wall=u,kind=native[u],cosine=float(co[u]))for u in bad])
       if len(bad)<=3:
        rec['all_retained_labels']=[native[u]for u in retained];near.append(rec)
       if not bad:
        if any(max(abs(np.array(z['residual4'])-r))<1e-6 and max(abs(np.array(z['T_pairings'])-dt))<1e-6 for z in out):continue
        rec['face']=R.poly(np.r_[Y,[r]],retained+[169]);out.append(rec);stats['accepted']+=1
       trials.append([ci,row['id'],br['branch'],j,lab,float(k),bad])
  print('REPLACED_EDGE',ci,dict(stats),'near',len(near),'accepted',len(out),flush=True)
 paths=[f'scripts/{TAG}.py',P.END,'results/fresh_h20_prism_replacement_20261017.json',f'scripts/{P.TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{R.TAG}.py']
 return dict(status='BOUNDED_REPLACED_TRIANGLE_FINITE_EVENT_DISCOVERY',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},base=dict(old_wall_count=168,removed=[158],replacement_recipe_branch=1,replacement_internal_id=168),actual_T_facets=oldfacets,new_base_finite_corners=[v.tolist()for v in vnew],statistics=dict(stats),candidates=out,at_most_three_bad_gates=near,all_intrinsic_event_trials=trials,scope='All159 fixed endpoint profiles, both branches, with one additional intrinsic finite old-facet event. Nonisolated labels2..5; isolated labels2..64. The singleton151 endpoint family, unbounded isolated labels, and intermediate facets not terminating these three edges remain outside scope. All retained native gates are freely classified as right/finite/parallel/dotted; no dotted-other assumption. Entire calculation is numerical discovery, not exact promotion or a global obstruction.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],flush=True)
