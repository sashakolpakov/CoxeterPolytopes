#!/usr/bin/env python3
"""Separate Python-integer/Fraction audit of incidence normals and their orbit.

Does not import the construction, NumPy, or its classifier. Full retained-bank
gates remain the separately replayed producer and the reviewed orbit argument.
"""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
CORE='results/fresh_h22_incidence_core_20261013_full.json'
BINARY='results/fresh_h22_incidence_binary_20261013.json'
SYMS='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
OUT='results/fresh_h22_incidence_core_20261013_review.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(M,v):return [dot(r,v)for r in M]
def kind(p,a,b):
 assert p<=0
 if p==0:return 'right'
 t=F(p*p,a*b)
 if t==F(1,4):return 'pi/3'
 if t==F(1,2):return 'pi/4'
 if t==F(3,4):return 'pi/6'
 assert t>=1
 return 'parallel'if t==1 else'dotted'
def run():
 raw=read(RAW)['coxeter_input'];Q=[[-int(x)for x in row]for row in raw['basis_gram']];W=[[int(x)for x in row]for row in raw['wall_coordinates_as_columns_transposed']]
 core=read(CORE);binary=read(BINARY);C=binary['accepted_old_projections'];wc=[mv(Q,w)for w in W];norms=[dot(w,q)for w,q in zip(W,wc,strict=True)]
 parts=core['incidence_parts'];assert sorted(parts[0]+parts[1])==[i for i,n in enumerate(norms)if n==2]
 assert all(dot(W[i],wc[j])==0 for part in parts for i,j in combinations(part,2))
 A=[[-dot(W[i],wc[j])for j in parts[1]]for i in parts[0]]
 assert len(A)==21 and all(len(r)==21 and set(r)<={0,1}and sum(r)==5 for r in A)
 assert all(dot(A[i],A[j])==4*(i==j)+1 for i in range(21)for j in range(21))
 a,b=[[sum(W[i][k]for i in part)for k in range(22)]for part in parts]
 V=[[F(2*x+5*y,21)for x,y in zip(a,b,strict=True)],[F(5*x+2*y,21)for x,y in zip(a,b,strict=True)]]
 assert all(x.denominator==1 for v in V for x in v);V=[[int(x)for x in v]for v in V]
 assert V==core['old_projections'];vc=[mv(Q,v)for v in V];h=[x+y for x,y in zip(*V,strict=True)];hc=mv(Q,h)
 assert h==core['center']and dot(h,hc)==-14 and [dot(v,q)for v,q in zip(V,vc,strict=True)]==[-2,-2]and dot(V[0],vc[1])==-5
 assert max(dot(w,hc)for w in W)<0 and all(dot(v,hc)==-7 for v in V)
 corecounts=Counter()
 R=[(v,z)for v in V for z in(2,-2)]
 for v,z in R:
  for w,q,n in zip(W,wc,norms,strict=True):corecounts[kind(dot(v,q),2,n)]+=1
 for (v,z),(w,t)in combinations(R,2):corecounts[kind(dot(v,mv(Q,w))+z*t,2,2)]+=1
 assert dict(corecounts)==core['new_pair_counts']
 assert len(C)==len({tuple(r)for r in C})==280
 cc=[mv(Q,c)for c in C]
 assert all(dot(c,q)==-2 and dot(c,hc)==-8 for c,q in zip(C,cc,strict=True))
 assert all(dot(c,q)==-4 for c in C for q in vc)
 oldcounts=Counter();newcounts=Counter()
 for c in C:
  for q,n in zip(wc,norms,strict=True):oldcounts[kind(dot(c,q),2,n)]+=2
 for i,j in combinations(range(280),2):newcounts[kind(dot(C[i],cc[j])+4,2,2)]+=2
 for c in C:
  for q in cc:newcounts[kind(dot(c,q)-4,2,2)]+=1
 assert dict(oldcounts)==binary['new_old210_pair_counts'];assert dict(newcounts)==binary['new_new_pair_counts']
 assert sum(newcounts.values())==560*559//2
 lookup={tuple(r):i for i,r in enumerate(C)};perms=[]
 for m in read(SYMS)['ambient_symmetry_generators']:
  g=m['matrix'];assert mv(g,h)==h
  cols=list(map(list,zip(*g,strict=True)));qc=[mv(Q,c)for c in cols]
  assert [[dot(c,q)for q in qc]for c in cols]==Q
  assert [mv(g,w)for w in W]==[W[i]for i in m['wall_permutation']]
  perm=[lookup[tuple(mv(g,c))]for c in C];assert sorted(perm)==list(range(280));perms.append(perm)
 assert perms==binary['symmetry_orbit']['generator_permutations']
 seen={0};todo=[0]
 while todo:
  i=todo.pop()
  for p in perms:
   if p[i]not in seen:seen.add(p[i]);todo.append(p[i])
 assert len(seen)==280
 esc=read('results/fresh_h22_new_cusp_completion_20261013_lower_escape.json');X=esc['numerator'];den=esc['denominator'];assert den>0
 xp=mv(Q,X[:22]);cuts=[{'index':i,'pairing':str(F(dot(c,xp)+2*abs(X[22]),den))}for i,c in enumerate(C)if dot(c,xp)+2*abs(X[22])>0]
 assert cuts==binary['prior_escape_cuts']and len(cuts)==8
 return {'status':'PASS_SEPARATE_INTEGER_INCIDENCE_AND_560_ORBIT_REVIEW','sources':{p:sha(p)for p in(RAW,CORE,BINARY,SYMS)},'core_new_pairs':dict(corecounts),'binary_new_old210_pairs':dict(oldcounts),'binary_new_new_pairs':dict(newcounts),'orbit_size':len(seen),'old_escape_cutters':cuts,'scope':'Separate exact integer/Fraction verification of the four core normals,280 projections,560-wall pair gates against each other and old210/core4, strict center, exact transitive orbit, and eight old-escape cutters. Full old-bank representative gates and bank invariance are source-reviewed and separately replayed, not reimplemented here. No whole-volume claim.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
