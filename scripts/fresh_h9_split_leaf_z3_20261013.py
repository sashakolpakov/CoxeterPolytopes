#!/usr/bin/env python3
"""Lazy numerical compatibility+escape search; emits explicit SAT conflicts.
An UNSAT result needs exact replays of every premise before geometric use.
"""
import itertools,json,os,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection,QhullError
import z3
from scripts.fresh_h9_split_leaf_caps_20261013 import seed,coeff,LABELS
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT/'data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_scout.json';OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_z3_20261013_exact_source_scout.json'

def run():
 Path('/tmp/fresh_h9_split_leaf_z3_20261013.pid').write_text(str(os.getpid()))
 old=json.loads(SOURCE.read_text());a=old['a'];b=old['b'];G=seed(a,b)
 import mpmath as mp
 mp.mp.dps=70;phi=(1+mp.sqrt(5))/2
 def kval(v):return mp.mpf(v['a'])+mp.mpf(v['b'])*phi
 exact=json.loads((ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json').read_text())
 R=np.array([[float(kval(v[0])+kval(v[1])*mp.sqrt(max(mp.mpf(0),kval(row['radicand'])))) for v in row['normal']] for row in exact['caps']]);n=len(R)
 eig,Q=np.linalg.eigh(G);idx=list(np.where(eig< -1e-7)[0])+list(np.where(eig>1e-7)[0]);N=Q[:,idx]*np.sqrt(abs(eig[idx]));N[:,0]*=np.sign(sum(N[:,0]));caps=R@N
 H=R@G@R.T;C=H<=-2+1e-8
 for v in LABELS:C|=abs(H+v)<1e-8
 np.fill_diagonal(C,False)
 solver=z3.Solver();solver.set(timeout=20000);xs=[z3.Bool('r'+str(i)) for i in range(n)];witnesses=[];seen=set();history=[];conflicts=[];positive=None;unsat_core=[];rng=np.random.default_rng(101307)
 def add_witness(v,origin):
  v=np.asarray(v)*(1-1e-6)
  assert v@v>1+1e-7 and max(N[:,1:]@v-N[:,0])<1e-7
  cut=tuple(np.flatnonzero(caps[:,1:]@v-caps[:,0]>1e-7).tolist())
  if cut in seen:return False
  seen.add(cut);tag='w'+str(len(witnesses));solver.assert_and_track(z3.Or([xs[i] for i in cut]),tag);witnesses.append({'klein_point':v.tolist(),'square':float(v@v),'cutters':list(cut),'origin':origin});return True
 def save():
  OUT.write_text(json.dumps({'scope':'Numerical only: exact normals are separately available, but these approximate nonedges and escape signs are not certified. Explicit Boolean UNSAT core may be independently promoted only after exact premise checks.','source':str(SOURCE.relative_to(ROOT)),'normal_source':'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json','history':history,'witnesses':witnesses,'incompatibilities':conflicts,'unsat_core':unsat_core,'positive':positive},indent=2)+'\n')
 for i,rec in enumerate(old['records']):
  if 'worst_klein_vertex' in rec['geometry']:add_witness(rec['geometry']['worst_klein_vertex'],'previous-greedy-'+str(i))
 print({'caps':n,'initial_escape_clauses':len(witnesses)},flush=True)
 prior=ROOT/'data/enumerations/fresh_h9_split_leaf_z3_20261013_scout.json'
 if prior.exists():
  for w in json.loads(prior.read_text())['witnesses']:add_witness(w['klein_point'],'prior-z3')
 for iteration in range(10000):
  start=time.time();status=solver.check();entry={'iteration':iteration,'status':str(status),'seconds':time.time()-start,'escapes':len(witnesses),'conflicts':len(conflicts)}
  if status==z3.unsat:
   unsat_core=[str(x) for x in solver.unsat_core()];entry['core_size']=len(unsat_core);history.append(entry);save();print(entry,flush=True);break
  if status!=z3.sat:
   entry['reason']=solver.reason_unknown();history.append(entry);save();print(entry,flush=True);break
  model=solver.model();chosen=[i for i,x in enumerate(xs) if z3.is_true(model.eval(x,model_completion=True))];bad=[(i,j) for p,i in enumerate(chosen) for j in chosen[p+1:] if not C[i,j]]
  if bad:
   for i,j in bad:
    tag='e'+str(len(conflicts));conflicts.append([i,j]);solver.assert_and_track(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])),tag)
   entry['new_conflicts']=len(bad);entry['selected_count']=len(chosen);history.append(entry)
   if iteration%100==0:print(entry,flush=True);save()
   continue
  rem=np.all(C[chosen],axis=0) if chosen else np.ones(n,dtype=bool);selected=chosen[:]
  while np.any(rem):
   avail=np.flatnonzero(rem);j=int(avail[np.argmax(C[np.ix_(avail,avail)].sum(axis=1))]);selected.append(j);rem&=C[j]
  W=np.vstack([N,caps[selected]]);new=0;maxq=0
  for direction in rng.normal(size=(60,9)):
   ans=linprog(-direction,A_ub=W[:,1:],b_ub=W[:,0],bounds=[(None,None)]*9,method='highs');assert ans.success;q=ans.x@ans.x;maxq=max(maxq,q)
   if q>1+1e-5:new+=add_witness(ans.x,'z3-'+str(iteration))
   if new>=12:break
  entry.update({'facets':16+len(selected),'selected':selected,'new_escape_clauses':new,'sampled_max_square':float(maxq)})
  if not new:
   try:V=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9)).intersections
   except QhullError:V=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9),qhull_options='Qs').intersections
   q=np.sum(V*V,axis=1);outside=np.flatnonzero(q>1+1e-6);entry.update({'vertices':len(V),'outside':len(outside),'max_square':float(max(q))})
   if not len(outside):positive=entry
   else:
    for j in outside[np.argsort(q[outside])[::-1]][:60]:new+=add_witness(V[j],'z3-full-'+str(iteration))
    assert new>0
  history.append(entry);save();print({k:v for k,v in entry.items() if k!='selected'},flush=True)
  if positive:break
 save();print({'done':True,'last':history[-1]['status'],'core':len(unsat_core),'positive':positive is not None},flush=True)
if __name__=='__main__':run()
