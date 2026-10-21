#!/usr/bin/env python3
"""Numerical endpoint families on the nine open edges of the five-cap seed."""
import argparse,hashlib,itertools,json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
from fresh_h20_two_anchor_families_20261016 import charts,allowed_intervals,constant_gate,point,samples
from fresh_h20_one_anchor_dotted_families_20261016 import cells
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158,gate_mask,face_data,poly
from fresh_h20_release144_cap_20261016_residual_scout import close
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_seed_edge_families_20261016';BASE='results/fresh_h20_open_edge_closure_20261016.json'
def run():
 D=json.loads((ROOT/BASE).read_text());G=geometry158();M=G['metric'];H=np.array(D['seed']['inherited_H']);seed=np.array(D['seed']['roots'])
 X=np.r_[G['X'],seed];q=np.r_[G['norms'],[2.]*5];Z=np.array([v['z']for v in D['seed']['face']['vertices']if v['radius2']>1+1e-6])
 families=[];pool=[];seen=set();stats=Counter();labels=[2,3,4,5,6,8,10]
 for edge in D['seed']['remaining_open_edges']:
  frame=[151,*edge['edge']];U=X[frame];B=U@M@U.T;assert min(np.linalg.eigvalsh(B))>1e-7
  inv=np.linalg.inv(B);K=null_space(U@M)
  for labs in itertools.product(labels,repeat=2):
   cs=np.array([0.,*[0 if m==2 else np.cos(np.pi/m)for m in labs]]);target=-np.sqrt(2*q[frame])*cs
   p=target@inv@U;w=2-p@M@p;stats['patterns']+=1
   if w<=1e-7:stats['not_finite_endpoint']+=1;continue
   for f in charts(p,K,M):
    num=np.array(f['numerator']);den=np.array(f['denominator']);di=np.argmax(abs(den));co=-(num@M@X.T).T/np.sqrt(2*q)[:,None]
    P=[];who=[];fixed=[];bad=False
    for j,cc in enumerate(co):
     v=cc[di]/den[di]
     if j in frame or max(abs(cc-v*den))<1e-7:
      kind=constant_gate(v)
      if kind is None:bad=True;break
      fixed.append(dict(wall=j,cosine=float(v),kind=kind))
     else:P.append(cc-den);who.append(j)
    if bad:stats['bad_fixed_contact']+=1;continue
    P.extend([-num@M@H,num[:,0]]);P=np.array(P);who+=['center','positive_time'];ranges=allowed_intervals(P,f['positive'])
    if not ranges:stats['empty_intervals']+=1;continue
    f.update(id=len(families),frame=frame,labels=[2,*labs],intervals=ranges,finite_endpoint_margin=float(w),
       fixed_contacts=fixed,inequality_coefficients=P.tolist(),inequality_rows=who)
    cc,pp=cells(f,Z);f['coverage_cells']=cc;f['point_pairing_coefficients']=pp.tolist();families.append(f)
    for cellid,c in enumerate(cc):
     for t in samples(c['interval']):
      r=point(f,t)
      if abs(r@M@r-2)>2e-4 or r@M@H>=0 or r[0]<=0 or not gate_mask((-(r@M@X.T)/np.sqrt(2*q))[None,:])[0]:stats['numerical_sample_reject']+=1;continue
      key=tuple(np.round(r,7))
      if key in seen:continue
      seen.add(key);ar=-r[1:]/r[0];cuts=np.where(Z@ar>1+1e-7)[0].tolist()
      if not cuts:continue
      pool.append(dict(id=len(pool),family=f['id'],cell=cellid,parameter=t,root=r.tolist(),constraint=ar.tolist(),cuts=cuts,finite_corner=frame))
  print('EDGE',edge['edge'],'families',len(families),'samples',len(pool),flush=True)
 src=[f'scripts/{TAG}.py',BASE,'scripts/fresh_h20_two_anchor_families_20261016.py','scripts/fresh_h20_one_anchor_dotted_families_20261016.py','scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py']
 return dict(status='PASS_BOUNDED_SEED_OPEN_EDGE_ENDPOINT_FAMILIES',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},label_domain=labels,stats=dict(stats),families=families,candidates=pool,exterior_points=Z.tolist(),
  scope='Numerical conic intervals for the nine open seed edges, right151 and finite labels2,3,4,5,6,8,10 to both edge sides. Positive endpoint Schur margin required. Every other variable contact is dotted; isolated third finite contacts are outside this family domain. No exact or continuous unrestricted exclusion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert close(d,json.loads(p.read_text()))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']),len(d['candidates']),flush=True)
