#!/usr/bin/env python3
"""Numerical compact-completion search beyond the closed split-seed bank.

Exact normal archives are inputs; cap-pair classification and geometry here
are numerical discovery only. No numerical UNSAT is an exact exclusion.
"""
import argparse,hashlib,json,os,time
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection,QhullError
import z3
from scripts.fresh_h9_split_leaf_caps_20261013 import seed,PHI
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
NEW=ROOT/'data/enumerations/fresh_h9_split_endpoint_families_20261013_exact.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_endpoint_selection_20261013_scout.json'

def bits(m):
 while m:
  l=m&-m;yield l.bit_length()-1;m-=l

def run(steps):
 Path('/tmp/fresh_h9_split_endpoint_selection_20261013.pid').write_text(str(os.getpid()))
 mp.mp.dps=60;p=(1+mp.sqrt(5))/2;k=lambda x:mp.mpf(x['a'])+mp.mpf(x['b'])*p
 rows=[];refs=[];sources=[]
 for path in [OLD,NEW]:
  data=json.loads(path.read_text());sources.append({'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
  for row in data['caps']:
   sd=mp.sqrt(max(mp.mpf(0),k(row['radicand'])));cov=np.array([float(k(x[0])+k(x[1])*sd)for x in row['covector']])
   if min(abs(cov+2))<1e-7:continue
   if any(sum(abs(cov[s:s+8]+1)<1e-7)>=4 for s in [0,8]):continue
   r=[float(k(x[0])+k(x[1])*sd)for x in row['normal']];rows.append(r);refs.append([path.name,row['index']])
 R=np.array(rows);a=PHI+.5;b=-a/(12*a+1);G=seed(a,b);n=len(R);pairs=R@G@R.T
 alphabet=[2,3,4,5,6,8,10,12,15,16,17,20,24,30]
 allowed=pairs<-2-1e-7
 for m in alphabet:allowed|=abs(pairs+2*np.cos(np.pi/m))<1e-7
 np.fill_diagonal(allowed,False);adj=[sum(1<<int(j)for j in np.flatnonzero(row))for row in allowed]
 del pairs
 scratch=Path('/tmp/fresh_h9_split_endpoint_selection_20261013.smt2');conflicts=0
 with scratch.open('w')as f:
  for i in range(n):f.write(f'(declare-const r{i} Bool)\n')
  for i in range(n):
   for j in np.flatnonzero(~allowed[i,i+1:])+i+1:f.write(f'(assert (or (not r{i}) (not r{j})))\n');conflicts+=1
 s=z3.Solver();s.set(timeout=45000);s.from_file(str(scratch));xs=[z3.Bool('r'+str(i))for i in range(n)]
 ev,Q=np.linalg.eigh(G);ii=list(np.where(ev< -1e-7)[0])+list(np.where(ev>1e-7)[0]);N=Q[:,ii]*np.sqrt(abs(ev[ii]));N[:,0]*=np.sign(sum(N[:,0]));W=R@N;eps=np.tile(np.r_[np.ones(4),-np.ones(4)],2);d=45+74*PHI;clauses=[]
 for side in [0,1]:
  P=np.ones(16);P[side*8:side*8+8]+=4*a;T=eps.copy();T[side*8:side*8+8]+=4*b*eps[side*8:side*8+8]
  for sg in [-1,1]:
   endpoint=P+sg*np.sqrt(d)*T;ids=np.flatnonzero(R@G@endpoint>1e-7).tolist();s.add(z3.Or([xs[i]for i in ids]));clauses.append({'kind':'seed_edge_endpoint','side':side,'sign':sg,'cutters':ids})
 v=np.array([8,8,8,1,-3,-3,-3,-3,-3,-3,-3,5,5,5,5,5]);ids=np.flatnonzero(R@G@v>1e-7).tolist();s.add(z3.Or([xs[i]for i in ids]));clauses.append({'kind':'old_universal_escape','cutters':ids})
 history=[];seen=set();rng=np.random.default_rng(2026101317);candidate=None
 def save(status):OUT.write_text(json.dumps({'status':status,'scope':'Numerical discovery using exactly reconstructed normals; no exact graph, complete all-affine test, or finite-volume/compactness certificate is asserted.','sources':sources,'references':refs,'pair_conflicts':conflicts,'necessary_clauses':clauses,'history':history,'candidate':candidate},indent=2)+'\n')
 def point_clause(v):
  v=np.asarray(v)*(1-1e-7);mask=W[:,1:]@v-W[:,0]>1e-8;ids=np.flatnonzero(mask).tolist();key=tuple(ids)
  if key in seen:return False
  seen.add(key);s.add(z3.Or([xs[i]for i in ids]));clauses.append({'kind':'numerical_spacelike_point','klein_point':v.tolist(),'square':float(v@v),'cutters':ids});return True
 print({'entries':n,'pair_conflicts':conflicts,'endpoint_clause_lengths':[len(c['cutters'])for c in clauses]},flush=True)
 for it in range(steps):
  start=time.time();status=s.check();row={'iteration':it,'status':str(status),'seconds':time.time()-start,'clauses':len(clauses)}
  if status!=z3.sat:
   if status==z3.unknown:row['reason']=s.reason_unknown()
   history.append(row);save('NUMERICAL_UNSAT_REQUIRES_EXACT_REPLAY'if status==z3.unsat else'NUMERICAL_SOLVER_UNKNOWN');print(row,flush=True);break
  model=s.model();chosen=[i for i in range(n)if z3.is_true(model.eval(xs[i],model_completion=True))];rem=(1<<n)-1
  for i in chosen:rem&=adj[i]
  while rem:
   opts=list(bits(rem));j=max(opts,key=lambda i:(adj[i]&rem).bit_count());chosen.append(j);rem&=adj[j]
  H=np.vstack([N,W[chosen]]);added=0;maxq=0
  for direction in rng.normal(size=(96,9)):
   ans=linprog(-direction,A_ub=H[:,1:],b_ub=H[:,0],bounds=[(None,None)]*9,method='highs');assert ans.success;q=float(ans.x@ans.x);maxq=max(maxq,q)
   if q>1+1e-6:added+=point_clause(ans.x)
   if added>=24:break
  if not added:
   try:verts=HalfspaceIntersection(np.c_[H[:,1:],-H[:,0]],np.zeros(9)).intersections
   except QhullError:verts=HalfspaceIntersection(np.c_[H[:,1:],-H[:,0]],np.zeros(9),qhull_options='Qs').intersections
   squares=np.sum(verts*verts,axis=1);maxq=float(max(squares));outside=np.flatnonzero(squares>1+1e-6);ideal=np.flatnonzero(abs(squares-1)<=1e-6)
   row.update({'vertices':len(verts),'outside':len(outside),'near_ideal':len(ideal)})
   if len(outside):
    for j in outside[np.argsort(squares[outside])[::-1]][:48]:added+=point_clause(verts[j])
   elif len(ideal):
    # A numerical affine no-good is discovery only; exact replay required.
    for j in ideal[:24]:
     active=[chosen[i]for i in range(len(chosen))if abs(W[chosen[i],1:]@verts[j]-W[chosen[i],0])<1e-6];assert active;s.add(z3.Or([~xs[i]for i in active]));clauses.append({'kind':'numerical_ideal_vertex_no_good','active_caps':active,'klein_point':verts[j].tolist()});added+=1
   else:candidate=chosen
  row.update({'selected':chosen,'facets':16+len(chosen),'new_constraints':added,'sampled_max_square':maxq});history.append(row);save('NUMERICAL_COMPACT_CANDIDATE_REQUIRES_EXACT_CERTIFICATE'if candidate else'OPEN_NUMERICAL_COMPACT_COMPLETION');print({k:v for k,v in row.items()if k!='selected'},flush=True)
  if candidate:break
 print({'finished':True,'candidate':candidate,'constraints':len(clauses)},flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--steps',type=int,default=40);a=p.parse_args();run(a.steps)
