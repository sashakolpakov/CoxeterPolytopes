#!/usr/bin/env python3
"""One complete fixed-profile generation on the three new nonright-T cap edges.

Numerical, conditional on the explicit one-wall replacement branch and seed83.
No unbounded label or recursive termination claim.
"""
import argparse,hashlib,itertools,json,math
from collections import Counter
import numpy as np
import fresh_h20_replaced_old_edge_profiles_20261017 as O
P,R,B,ROOT=O.P,O.R,O.B,O.ROOT
TAG='fresh_h20_h4_edge_continuation_20261017';INPUT='results/fresh_h20_replaced_old_edge_profiles_20261017.json';SEED=83
def run():
 G,Y,CT,J,q,V,C=P.geometry();m=G['m'];T=G['T'];rep=P.P.load('results/fresh_h20_prism_replacement_20261017.json');D=P.P.load(INPUT);seed=D['candidates'][SEED]
 assert seed['contacts']==[[61,3],[142,5]]and seed['anchor']==152 and seed['label']==5
 Y=np.r_[Y,[rep['residual4'][1],seed['residual4']]];CT=np.r_[CT,np.zeros((1,17)),[seed['T_pairings']]];q=np.r_[q,2,2];ret=[i for i in range(170)if i!=158]
 full=Y@m@Y.T+CT@J@CT.T;pos=q-np.einsum('ij,jk,ik->i',CT,J,CT);facets=[i for i in ret if i not in T and pos[i]>1e-7]
 face=R.poly(Y,ret);V=np.array([[1.]+v['point']for v in face['vertices']if v['radius_squared']<1-1e-6]);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
 censuses=[]
 for j in facets:
  if j==169:continue
  W=T+[j,169]
  if min(np.linalg.eigvalsh(full[np.ix_(W,W)]))<1e-8:continue
  c=O.census(W,full[np.ix_(W,W)]);c.update(edge=[j,169],finite_endpoint_indices=[k for k,v in enumerate(V)if max(abs(Y[[j,169]]@m@v))<1e-7]);censuses.append(c)
  print('NEW_EDGE',j,'profiles',len(c['profiles']),'components',c['components'],'finiteends',c['finite_endpoint_indices'],flush=True)
 stats=Counter();out=[];near=[];trials=[]
 for ci,c in enumerate(censuses):
  for row in c['profiles']:
   for sg,p,a,b,dt,offset in O.conic(c,row,G,Y,CT,J,q):
    aa=a@m@Y.T;bb=b@m@Y.T;cc=p@m@Y.T+offset;h=2-dt@J@dt
    for j in facets:
     if j in c['edge']:continue
     isolated=max(abs(dt))<1e-9 and max(abs(CT[j]))<1e-8;labs=range(2,65)if isolated else range(2,6)
     for lab in labs:
      stats['equations']+=1;target=-math.sqrt(2*q[j])*(0 if lab==2 else math.cos(math.pi/lab))
      for k in B.roots(aa[j],cc[j]-target,bb[j]):
       stats['real_branches']+=1;r=p+k*a+b/k
       if max(r@m@V.T)>=-1e-8:stats['old_vertex_not_strict']+=1;continue
       ps=r@m@Y.T+offset;co=-ps/np.sqrt(2*q);native=[R.classify(x)for x in co]
       if h*pos[j]-(ps[j]-offset[j])**2<=1e-8:stats['extra_not_finite']+=1;continue
       bad=[u for u in ret if native[u]in('obtuse','angle_gap')];stats['intrinsic_events']+=1;trials.append([ci,row['id'],sg,j,lab,k,bad])
       rec=dict(census=ci,edge=c['edge'],profile=row['id'],contacts=row['contacts'],branch=sg,anchor=j,label=lab,k=k,residual4=r.tolist(),T_pairings=dt.tolist(),bad_wall_ids=bad,bad_gate_details=[dict(wall=u,kind=native[u],cosine=float(co[u]))for u in bad])
       if len(bad)<=3:near.append(rec)
       if not bad:
        if any(max(abs(np.array(z['residual4'])-r))<1e-6 and max(abs(np.array(z['T_pairings'])-dt))<1e-6 for z in out):continue
        rec['all_retained_labels']=[native[u]for u in ret];rec['face']=R.poly(np.r_[Y,[r]],ret+[170]);out.append(rec);stats['accepted']+=1
  print('NEW_EDGE_SOLVED',ci,dict(stats),'accepted',len(out),'near',len(near),flush=True)
 paths=[f'scripts/{TAG}.py',INPUT,'results/fresh_h20_prism_replacement_20261017.json',f'scripts/{O.TAG}.py',f'scripts/{P.TAG}.py',f'scripts/{R.TAG}.py',f'scripts/{B.TAG}.py']
 return dict(status='BOUNDED_NONRIGHT_T_H4_EDGE_CONTINUATION',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},base=dict(removed_old_ids=[158],replacement_internal_id=168,seed_internal_id=169,seed_source=INPUT,seed_index=SEED,seed_contacts=seed['contacts'],seed_anchor=[seed['anchor'],seed['label']],wall_count=169),base_face=face,endpoint_censuses=censuses,statistics=dict(stats),candidates=out,at_most_three_bad_gates=near,all_intrinsic_event_trials=trials,scope='All fixed finite deletion profiles on every actual finite edge of the chosen nonright-T seed facet. Old labels are recognized numerically; the deletion matching and symbolic Schur residuals are exact conditional on those labels. Additional finite old-facet events use nonisolated labels2..5 and isolated2..64. Explicit singleton endpoint tails, parabolic endpoints, and caps with no further old finite contact remain open. This is one edge-driven discovery generation, not an exact ambient construction or complete continuous obstruction.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],flush=True)
