#!/usr/bin/env python3
"""Completion search using certified necessary Boolean constraints.

Every graph nonedge and every escape clause has an exact arithmetic proof.
SAT means only that these necessary conditions remain consistent. UNSAT would
exclude finite-volume completions from the declared fixed mirror bank.
"""
import argparse,hashlib,json,os,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection,QhullError
import z3
from scripts.fresh_h9_split_leaf_escape_constraints_20261013 import setup,certify,DEN
from scripts.fresh_h9_split_leaf_caps_20261013 import seed
ROOT=Path(__file__).resolve().parents[1]
NORMALS=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
ESCAPES=ROOT/'data/enumerations/fresh_h9_split_leaf_escape_constraints_20261013.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_certified_sat_20261013.json'

def indices(mask):
 while mask:
  l=mask&-mask;yield l.bit_length()-1;mask-=l

def run(graph_path,steps):
 Path('/tmp/fresh_h9_split_leaf_certified_sat_20261013.pid').write_text(str(os.getpid()))
 graph=json.loads(graph_path.read_text());saved=json.loads(ESCAPES.read_text());raw=NORMALS.read_bytes();digest=hashlib.sha256(raw).hexdigest();assert graph['source_sha256']==saved['normal_sha256']==digest;normal=json.loads(raw);n=len(normal['caps']);adj=[int(x,16) for x in graph['adjacency_hex']];assert len(adj)==n
 A,B,pats,dints=setup();constraints=saved['constraints'][:]
 for c in constraints:assert certify(c['integer_coefficients'],A,B,pats,dints)==c
 print({'certified_initial_escapes':len(constraints),'graph_possible_edges':sum(m.bit_count() for m in adj)//2},flush=True)
 scratch=Path('/tmp/fresh_h9_split_leaf_certified_sat_20261013.smt2');edge_clauses=0
 with scratch.open('w') as f:
  for i in range(n):f.write(f'(declare-const r{i} Bool)\n')
  allbits=(1<<n)-1
  for i in range(n):
   nonedges=(allbits^adj[i])&~((1<<(i+1))-1)
   for j in indices(nonedges):f.write(f'(assert (or (not r{i}) (not r{j})))\n');edge_clauses+=1
  for j,c in enumerate(constraints):
   f.write(f'(declare-const w{j} Bool)\n(assert (=> w{j} (or '+ ' '.join('r'+str(i) for i in indices(int(c['possible_cutters_hex'],16)))+')))\n')
 solver=z3.Solver();solver.set(timeout=45000);start=time.time();solver.from_file(str(scratch));print({'loaded_certified_clauses':edge_clauses+len(constraints),'seconds':time.time()-start},flush=True)
 xs=[z3.Bool('r'+str(i)) for i in range(n)];assumptions=[z3.Bool('w'+str(i)) for i in range(len(constraints))];seen={r['possible_cutters_hex'] for r in constraints}
 import mpmath as mp
 mp.mp.dps=70;phi=(1+mp.sqrt(5))/2
 def kval(v):return mp.mpf(v['a'])+mp.mpf(v['b'])*phi
 R=np.array([[float(kval(v[0])+kval(v[1])*mp.sqrt(max(mp.mpf(0),kval(row['radicand'])))) for v in row['normal']] for row in normal['caps']]);p=(1+5**.5)/2;a=p+.5;b=-a/(12*a+1);G=seed(a,b);ev,Q=np.linalg.eigh(G);ix=list(np.where(ev< -1e-7)[0])+list(np.where(ev>1e-7)[0]);N=Q[:,ix]*np.sqrt(abs(ev[ix]));N[:,0]*=np.sign(sum(N[:,0]));caps=R@N;NI=np.linalg.pinv(N.T)
 history=[];core=[];candidate=None;rng=np.random.default_rng(6101307)
 def add_point(v):
  v=np.asarray(v)*(1-1e-6);w=[int(round(float(x)*DEN)) for x in NI@np.r_[1,v]]
  c=certify(w,A,B,pats,dints)
  if c['possible_cutters_hex'] in seen:return False
  seen.add(c['possible_cutters_hex']);j=len(constraints);constraints.append(c);tag=z3.Bool('w'+str(j));assumptions.append(tag);solver.add(z3.Implies(tag,z3.Or([xs[i] for i in indices(int(c['possible_cutters_hex'],16))])));return True
 def save(status):
  OUT.write_text(json.dumps({'status':status,'scope':'Exact necessary conditions for finite-volume completion of the frozen3034distinctmirrorbank. SAT is not a finite-volume certificate; an UNSAT result requires Boolean replay as well as the separate exact premises.','normal_source':str(NORMALS.relative_to(ROOT)),'normal_sha256':digest,'graph_source':str(graph_path.relative_to(ROOT)),'graph_sha256':hashlib.sha256(graph_path.read_bytes()).hexdigest(),'initial_escape_source':str(ESCAPES.relative_to(ROOT)),'initial_escape_sha256':hashlib.sha256(ESCAPES.read_bytes()).hexdigest(),'pair_conflict_count':edge_clauses,'history':history,'constraints':constraints,'unsat_core':core,'numerical_candidate':candidate},indent=2)+'\n')
 for iteration in range(steps):
  start=time.time();status=solver.check(*assumptions);entry={'iteration':iteration,'solver_status':str(status),'seconds':time.time()-start,'escape_constraints':len(constraints)};print(entry,flush=True)
  if status==z3.unsat:
   core=[int(str(t)[1:]) for t in solver.unsat_core()];entry['core_size']=len(core);history.append(entry);save('UNSAT_AWAITING_INDEPENDENT_BOOLEAN_REPLAY');break
  if status!=z3.sat:
   entry['reason']=solver.reason_unknown();history.append(entry);save('OPEN_SOLVER_UNKNOWN');break
  model=solver.model();selected=[i for i,x in enumerate(xs) if z3.is_true(model.eval(x,model_completion=True))];remaining=(1<<n)-1
  for i in selected:remaining&=adj[i]
  while remaining:
   options=list(indices(remaining));j=max(options,key=lambda i:(adj[i]&remaining).bit_count());selected.append(j);remaining&=adj[j]
  W=np.vstack([N,caps[selected]]);added=0;maxq=0
  for direction in rng.normal(size=(80,9)):
   ans=linprog(-direction,A_ub=W[:,1:],b_ub=W[:,0],bounds=[(None,None)]*9,method='highs');assert ans.success;q=float(ans.x@ans.x);maxq=max(maxq,q)
   if q>1+1e-5:added+=add_point(ans.x)
   if added>=24:break
  entry.update({'facets_if_actual_compatible':16+len(selected),'selected_entries':selected,'new_escape_constraints':added,'sampled_max_square':maxq})
  if not added:
   try:V=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9)).intersections
   except QhullError:V=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9),qhull_options='Qs').intersections
   q=np.sum(V*V,axis=1);outside=np.flatnonzero(q>1+1e-6);entry.update({'vertices':len(V),'outside':len(outside),'max_square':float(max(q))})
   if not len(outside):candidate=entry
   else:
    for j in outside[np.argsort(q[outside])[::-1]][:80]:added+=add_point(V[j])
    assert added
  history.append(entry);save('NUMERICAL_CANDIDATE_REQUIRES_EXACT_GEOMETRY' if candidate else 'OPEN_CERTIFIED_NECESSARY_CONDITIONS')
  print({k:v for k,v in entry.items() if k!='selected_entries'},flush=True)
  if candidate:break
 print({'finished':True,'constraints':len(constraints),'core':len(core),'candidate':candidate is not None},flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--graph',default='data/enumerations/fresh_h9_split_leaf_3034_safe_graph_20261013.json');p.add_argument('--steps',type=int,default=100);args=p.parse_args();run(ROOT/args.graph,args.steps)
