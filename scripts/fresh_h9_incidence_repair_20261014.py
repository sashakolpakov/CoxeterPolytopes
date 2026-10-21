#!/usr/bin/env python3
"""Change old right incidences, guided by the complete spherical face complex.

The SAT models are proposals unless separately certified; no compactness
claim follows from a satisfying diagram alone. Old parallel distances and
new dotted distances are free. No previous producer arithmetic is imported.
"""
from pathlib import Path
from itertools import combinations
from fractions import Fraction as F
import argparse,hashlib,json,time
import numpy as np
import z3
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_h9_incidence_surgery_20261013.json'
PATTERNS=ROOT/'results/fresh_h9_incidence_surgery_20261013_dotted.json'
OUT=ROOT/'results/fresh_h9_incidence_repair_20261014.json'
TREE_CACHE={}

def pd(A):
 A=[r[:]for r in A]
 for k in range(len(A)):
  if A[k][k]<=0:return False
  for i in range(k+1,len(A)):
   for j in range(i,len(A)):
    A[j][i]=A[i][j]=A[i][j]-A[i][k]*A[k][j]/A[k][k]
 return True

def source():
 d=json.loads(SOURCE.read_text());r=next(r for r in d['cases']if r['id']=='A16_14_1_C01')
 G=[[F(x)for x in row]for row in r['source_positive_gram']]
 pats=next(r for r in json.loads(PATTERNS.read_text())['cases']if r['id']=='A16_14_1_C01')['fixed_principal_witnesses']
 return r,G,pats

def family(row,G,pr):
 S=[tuple(p)for p in row['finite_edge_cover']['minimum_supports'][pr['support_index']]]
 H=[r[:]for r in G];fixed={tuple(p)for p in row['parallel_pairs']}
 for k,(i,j)in enumerate(S):
  if pr['new_dotted_mask']>>k&1:fixed.add((i,j))
  else:H[i][j]=H[j][i]=F(0)
 # Every originally right pair is genuinely new incidence freedom.
 choices=[(i,j)for i,j in combinations(range(15),2)if not G[i][j]]
 return H,fixed,choices

def diagram_complex(H,fixed):
 # Exact incremental PD census: every subset of a PD matrix is PD.
 levels={0:[()]}
 for k in range(1,11):
  out=[]
  for B in levels[k-1]:
   for j in range(B[-1]+1 if B else 0,15):
    if any((i,j)in fixed for i in B):continue
    C=B+(j,)
    if pd([[H[a][b]for b in C]for a in C]):out.append(C)
  levels[k]=out
 return levels

def finite_labels(H):
 return [[1 if i==j else 2 if not H[i][j] else {F(1,4):3,F(1,2):4,F(1):0}[H[i][j]*H[i][j]/(H[i][i]*H[j][j])]for j in range(len(H))]for i in range(len(H))]

def spherical(labels,B):
 """Exact positive definiteness for Coxeter labels2,3,4.

 A cycle has a nonpositive all-ones quadratic form on its cycle vertices.
 A forest is diagonally congruent to an integer Gram: crossing label4
 switches squared scale1/2, whereas label3 preserves it.
 """
 B=set(B);adj={i:[j for j in B if i!=j and labels[i][j]>2]for i in B};seen=set()
 for root in B:
  if root in seen:continue
  todo=[root];C={root};parent={root:None};parity={root:0};cyclic=False
  while todo:
   i=todo.pop()
   for j in adj[i]:
    if j==parent[i]:continue
    if j in C:cyclic=True;break
    C.add(j);parent[j]=i;parity[j]=parity[i]^(labels[i][j]==4);todo.append(j)
   if cyclic:break
  if cyclic:return False
  seen|=C
  def code(i,p):return tuple(sorted((labels[i][j],code(j,i))for j in adj[i]if j!=p))
  key=min(code(i,None)for i in C)
  if key not in TREE_CACHE:
   Ns=sorted(C);A=[]
   for i in Ns:
    row=[]
    for j in Ns:
     if i==j:v=2*(2**parity[i])
     elif labels[i][j]==2:v=0
     elif labels[i][j]==3:assert parity[i]==parity[j];v=-(2**parity[i])
     else:assert labels[i][j]==4 and parity[i]!=parity[j];v=-2
     row.append(F(v))
    A.append(row)
   TREE_CACHE[key]=pd(A)
  if not TREE_CACHE[key]:return False
 return True

def diagram_complex_labels(labels,fixed):
 levels={0:[()]}
 for k in range(1,11):
  out=[]
  for B in levels[k-1]:
   for j in range(B[-1]+1 if B else 0,15):
    if any((i,j)in fixed for i in B):continue
    C=B+(j,)
    if spherical(labels,C):out.append(C)
  levels[k]=out
 return levels

def search(pr,ordinal,timeout_ms=5000,model_limit=4,additions=(),release_old_disjointness=False):
 row,G,_=source();H,fixed,choices=family(row,G,pr);labels=finite_labels(H)
 if release_old_disjointness:
  for i,j in fixed:labels[i][j]=labels[j][i]=2
  fixed=set()
 for i,j,m in additions:assert labels[i][j]==2;labels[i][j]=labels[j][i]=m
 choices=[(i,j)for i,j in combinations(range(15),2)if labels[i][j]==2]
 levels=diagram_complex_labels(labels,fixed)
 X={p:z3.Bool(f'x_{i}_{j}')for p in choices for i,j in[p]};s=z3.Tactic('sat').solver();s.set(timeout=timeout_ms)
 alive={}
 def face(B):
  if B not in alive:alive[B]=z3.And([z3.Not(X[p])for p in combinations(B,2)if p in X])
  return alive[B]
 # A compact H9 acute chamber has no elliptic10 and every elliptic face
 # extends to an elliptic9 vertex. Every elliptic8 ridge has two endpoints.
 for B in levels[10]:s.add(z3.Not(face(B)))
 contain={B:[]for k in range(1,9)for B in levels[k]}
 for V in levels[9]:
  for k in range(1,9):
   for B in combinations(V,k):contain[B].append(V)
 pending=[]
 for B,Vs in contain.items():
  if len(B)==8:s.add(z3.Implies(face(B),z3.PbEq([(face(V),1)for V in Vs],2)if Vs else z3.BoolVal(False)))
  elif len(B)==1:s.add(z3.Implies(face(B),z3.Or([face(V)for V in Vs])))
  else:pending.append((B,Vs))
 models=[];rounds=[];refinements=[]
 for attempt in range(100):
  res=s.check();rounds.append(str(res))
  if res!=z3.sat:break
  M=s.model();D=[p for p in choices if z3.is_true(M.eval(X[p],model_completion=True))]
  Vs=[list(V)for V in levels[9]if all(p not in D for p in combinations(V,2))]
  bad=[(B,Ws)for B,Ws in pending if not any(p in D for p in combinations(B,2))and not any(all(p not in D for p in combinations(W,2))for W in Ws)]
  if bad:
   for B,Ws in bad:s.add(z3.Implies(face(B),z3.Or([face(V)for V in Ws])))
   refinements.append({'attempt':attempt,'unextendible_faces':[list(B)for B,_ in bad]});continue
  models.append({'extra_dotted_pairs':[list(p)for p in D],'vertices':Vs})
  if len(models)>=model_limit:break
  s.add(z3.Or([z3.Not(X[p])if p in D else X[p]for p in choices]))
 out={'pattern_index':ordinal,'support_index':pr['support_index'],'new_dotted_mask':pr['new_dotted_mask'],'released_old_disjointness':release_old_disjointness,'added_finite_contacts':[list(e)for e in additions],'fixed_free_pairs':[list(p)for p in sorted(fixed)],'new_right_pair_variables':len(choices),'spherical_counts':{str(k):len(v)for k,v in levels.items()},'solver_rounds':rounds,'refinements':refinements,'models':models}
 print(ordinal,'support',pr['support_index'],'fixed',len(fixed),'PD9/10',len(levels[9]),len(levels[10]),rounds,'models',[(len(m['extra_dotted_pairs']),len(m['vertices']))for m in models],flush=True)
 return out

def main():
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=72);p.add_argument('--start',type=int,default=0);p.add_argument('--timeout-ms',type=int,default=5000);p.add_argument('--models',type=int,default=4);a=p.parse_args();row,G,ps=source();rows=[]
 for i in range(a.start,min(a.start+a.limit,len(ps))):
  rows.append(search(ps[i],i,a.timeout_ms,a.models))
  data={'status':'EXACT_SPHERICAL_COMPLEX_WITH_BOUNDED_BOOLEAN_SOLVER','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[SOURCE,PATTERNS,Path(__file__).resolve()]},'case':row['id'],'cases':rows,'scope':'72 previous minimal right/dotted patterns; arbitrary additional originally right pairs may become dotted. Necessary compact spherical-face extension/ridge conditions only. Solver outcomes need independent replay or proof before a family theorem; SAT is not a metric realization.'}
  OUT.write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':main()
