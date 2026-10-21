#!/usr/bin/env python3
"""Counterexample-guided compatible completion of the 3040-cap numeric bank.
All geometry/pair filtering is numerical; no exact theorem is claimed here.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.spatial import HalfspaceIntersection,QhullError
from ortools.sat.python import cp_model
from scripts.fresh_h9_split_leaf_caps_20261013 import seed,coeff,LABELS
ROOT=Path(__file__).resolve().parents[1];SOURCE=ROOT/'data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_scout.json';OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_scout.json'

def run(lazy=False):
 old=json.loads(SOURCE.read_text());a=old['a'];b=old['b'];G=seed(a,b);R=np.array([coeff(a,b,np.array(x),np.array(y)) for x,y in old['all_profiles']]);n=len(R)
 eig,Q=np.linalg.eigh(G);idx=list(np.where(eig< -1e-7)[0])+list(np.where(eig>1e-7)[0]);N=Q[:,idx]*np.sqrt(abs(eig[idx]));N[:,0]*=np.sign(sum(N[:,0]));caps=R@N
 H=R@G@R.T;C=H<=-2+1e-8
 for v in LABELS:C|=abs(H+v)<1e-8
 np.fill_diagonal(C,False)
 model=cp_model.CpModel();xs=[model.new_bool_var('r'+str(i)) for i in range(n)];neg=[~x for x in xs]
 if not lazy:
  for i in range(n):
   for j in np.flatnonzero(~C[i,i+1:])+i+1:model.add_bool_or([neg[i],neg[int(j)]])
 else:
  model.minimize(sum(xs))
 print({'caps':n,'incompatibility_clauses':len(model.proto.constraints)},flush=True)
 witnesses=[];seen=set();history=[];positive=None
 def add_witness(v,origin):
  v=np.asarray(v)*(1-1e-6)
  assert v@v>1+1e-7 and max(N[:,1:]@v-N[:,0])<1e-7
  cut=tuple(np.flatnonzero(caps[:,1:]@v-caps[:,0]>1e-7).tolist())
  if cut in seen:return False
  seen.add(cut);model.add_bool_or([xs[i] for i in cut]);witnesses.append({'klein_point':v.tolist(),'square':float(v@v),'cutters':list(cut),'origin':origin});return True
 for i,rec in enumerate(old['records']):
  if 'worst_klein_vertex' in rec['geometry']:add_witness(rec['geometry']['worst_klein_vertex'],'previous-greedy-'+str(i))
 print({'initial_escape_clauses':len(witnesses)},flush=True)
 solver=cp_model.CpSolver();solver.parameters.num_search_workers=1;solver.parameters.max_time_in_seconds=20
 for iteration in range(80):
  start=time.time();status=solver.solve(model);name=solver.status_name(status);entry={'iteration':iteration,'solver_status':name,'seconds':time.time()-start,'escape_clauses':len(witnesses)}
  print({'solver_returned':name,'iteration':iteration,'seconds':entry['seconds']},flush=True)
  if status in (cp_model.FEASIBLE,cp_model.OPTIMAL):
   chosen=np.array([i for i,x in enumerate(xs) if solver.value(x)],dtype=int)
   bad=[(int(i),int(j)) for p,i in enumerate(chosen) for j in chosen[p+1:] if not C[i,j]]
   if bad:
    assert lazy
    for i,j in bad:model.add_bool_or([neg[i],neg[j]])
    entry['new_incompatibilities']=len(bad);entry['chosen']=chosen.tolist();history.append(entry);print(entry,flush=True)
    OUT.write_text(json.dumps({'scope':'Numerical lazy-compatible SAT completion; no exact exclusion claim.','history':history,'witnesses':witnesses,'positive':None},indent=2)+'\n')
    continue
   # Add compatible unused roots: this can only improve volume while preserving
   # every SAT constraint. Prefer roots cutting previous sampled escapes.
   rem=np.all(C[chosen],axis=0) if len(chosen) else np.ones(n,dtype=bool)
   selected=chosen.tolist()
   while np.any(rem):
    avail=np.flatnonzero(rem);scores=C[np.ix_(avail,avail)].sum(axis=1);j=int(avail[np.argmax(scores)]);selected.append(j);rem&=C[j]
   print({'geometry_start':iteration,'facets':16+len(selected)},flush=True)
   W=np.vstack([N,caps[selected]])
   try:V=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9)).intersections
   except QhullError:V=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9),qhull_options='Qs').intersections
   q=np.sum(V*V,axis=1);outside=np.flatnonzero(q>1+1e-6);entry.update({'facets':16+len(selected),'selected':selected,'vertices':len(V),'outside':len(outside),'max_square':float(max(q))})
   if not len(outside):positive=entry;history.append(entry);print(entry,flush=True);break
   # Diverse most distant vertices; stop after 32 distinct new cutter clauses.
   added=0
   for j in outside[np.argsort(q[outside])[::-1]]:
    added+=add_witness(V[j],f'iteration{iteration}-vertex{int(j)}')
    if added>=32:break
   entry['new_escape_clauses']=added
  history.append(entry);print({k:v for k,v in entry.items() if k!='selected'},flush=True)
  OUT.write_text(json.dumps({'scope':'Numerical counterexample-guided SAT completion of the frozen3040cap bank. Approximate angles and numeric strict escapes; SAT infeasibility is not an exact geometric theorem without separate replay.','source':str(SOURCE.relative_to(ROOT)),'history':history,'witnesses':witnesses,'positive':positive},indent=2)+'\n')
  if status not in (cp_model.FEASIBLE,cp_model.OPTIMAL):break
 OUT.write_text(json.dumps({'scope':'Numerical counterexample-guided SAT completion of the frozen3040cap bank. Approximate angles and numeric strict escapes; SAT infeasibility is not an exact geometric theorem without separate replay.','source':str(SOURCE.relative_to(ROOT)),'history':history,'witnesses':witnesses,'positive':positive},indent=2)+'\n')
 print({'done':True,'last':history[-1]['solver_status'],'positive':positive is not None},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--lazy',action='store_true');args=p.parse_args()
 if args.lazy:OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_lazy_scout.json'
 run(args.lazy)
