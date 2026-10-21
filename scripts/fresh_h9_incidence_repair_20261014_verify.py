#!/usr/bin/env python3
"""Stdlib exact flag obstruction and complete right-contact ridge peeling.

No floating-point or Boolean solver predicates enter these two conclusions.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from collections import Counter
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_h9_incidence_surgery_20261013.json'
PATTERNS=ROOT/'results/fresh_h9_incidence_surgery_20261013_dotted.json'
OUT=ROOT/'results/fresh_h9_incidence_repair_20261014_verify.json'
CACHE={}

def pivots(A):
 A=[list(map(F,r))for r in A];p=[]
 for k in range(len(A)):
  q=A[k][k]
  if q<=0:return False,p+[q]
  p.append(q)
  for i in range(k+1,len(A)):
   for j in range(i,len(A)):A[j][i]=A[i][j]=A[i][j]-A[i][k]*A[k][j]/q
 return True,p

def spherical(L,B):
 adj={i:[j for j in B if i!=j and L[i][j]>2]for i in B};seen=set()
 for r in B:
  if r in seen:continue
  C={r};queue=[r];parent={r:None};scale={r:1}
  while queue:
   i=queue.pop()
   for j in adj[i]:
    if j==parent[i]:continue
    if j in C:return False  # cycle: all-ones form on cycle vertices <=0
    C.add(j);parent[j]=i;scale[j]=scale[i]if L[i][j]==3 else 3-scale[i];queue.append(j)
  seen|=C
  def code(i,p):return tuple(sorted((L[i][j],code(j,i))for j in adj[i]if j!=p))
  key=min(code(i,None)for i in C)
  if key not in CACHE:
   A=[]
   for i in sorted(C):
    row=[]
    for j in sorted(C):
     if i==j:v=2*scale[i]
     elif L[i][j]==2:v=0
     elif L[i][j]==3:assert scale[i]==scale[j];v=-scale[i]
     else:assert L[i][j]==4 and scale[i]!=scale[j];v=-2
     row.append(v)
    A.append(row)
   CACHE[key]=pivots(A)[0]
  if not CACHE[key]:return False
 return True

def right_covers(row,G):
 P=list(map(tuple,row['parallel_pairs']));C=[set(a['nodes'])for a in row['affine_circuits']if not any(i in a['nodes']and j in a['nodes']for i,j in P)]
 Z=[p for p in combinations(range(15),2)if not G[p[0]][p[1]]and any(set(p)<=c for c in C)]
 req=[{p for p in Z if set(p)<=c}for c in C];forced=set(next(iter(z))for z in req if len(z)==1);rest=[p for p in Z if p not in forced];sizes=[]
 for k in range(len(forced),8):
  covers=[]
  for T in combinations(rest,k-len(forced)):
   S=forced|set(T)
   if all(S&r for r in req):covers.append(sorted(S))
  sizes.append({'size':k,'covers':len(covers)})
  if covers:return covers,sorted(forced),sizes
 raise AssertionError('no7-cover')

def peel(V):
 V=set(V);layers=[]
 while V:
  counts=Counter(B for W in V for B in combinations(W,8))
  bad=sorted(W for W in V if any(counts[B]<2 for B in combinations(W,8)))
  assert bad,'unexcluded ridge core'
  proof=[]
  for W in bad:
   B=next(B for B in combinations(W,8)if counts[B]<2);assert counts[B]==1
   proof.append({'vertex':list(W),'unique_ridge':list(B)})
  layers.append(proof);V-=set(bad)
 return layers

def verify():
 d=json.loads(SOURCE.read_text());row=next(r for r in d['cases']if r['id']=='A16_14_1_C01');G=[[F(x)for x in r]for r in row['source_positive_gram']];P=set(map(tuple,row['parallel_pairs']));pd=json.loads(PATTERNS.read_text());assert pd['source_sha256']==hashlib.sha256(SOURCE.read_bytes()).hexdigest()
 ps=next(r for r in pd['cases']if r['id']==row['id'])['fixed_principal_witnesses'];flags=[]
 for pr in ps:
  H=[r[:]for r in G];S=row['finite_edge_cover']['minimum_supports'][pr['support_index']]
  for i,j in list(P)+S:H[i][j]=H[j][i]=F(0)
  ok,p=pivots(H);assert ok and len(p)==15
  flags.append({'support_index':pr['support_index'],'new_dotted_mask':pr['new_dotted_mask'],'positive_pivots':[str(x)for x in p]})
 assert len(flags)==72
 S,forced,minimum=right_covers(row,G);assert len(S)==11 and len(S[0])==7
 L=[[1 if i==j else {F(0):2,F(1,4):3,F(1,2):4,F(1):0}[G[i][j]**2/(G[i][i]*G[j][j])]for j in range(15)]for i in range(15)]
 Bs=[B for B in combinations(range(15),9)if not any(i in B and j in B for i,j in P)];assert len(Bs)==631
 rows=[]
 for si,edges in enumerate(S):
  for labels in product([3,4],repeat=7):
   H=[r[:]for r in L]
   for(i,j),m in zip(edges,labels):H[i][j]=H[j][i]=m
   V=[B for B in Bs if spherical(H,B)];layers=peel(V)
   rows.append({'support_index':si,'labels':list(labels),'spherical9':len(V),'vertices':[list(B)for B in V],'peeling':layers})
  print('right cover',si,'128 profiles exact, maxvertices',max(r['spherical9']for r in rows if r['support_index']==si),flush=True)
 assert len(rows)==1408
 return {'status':'EXACT_FLAG_FAMILY_AND_RIGHT_CONTACT_FAMILY_COMPACTNESS_OBSTRUCTIONS','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[SOURCE,PATTERNS,Path(__file__).resolve()]},'case':row['id'],'flag_patterns':flags,'flag_proof':'The finite-only15-node matrix is positive definite. All nondotted cliques are therefore spherical; arbitrary additional right-to-dotted mutations keep the nerve flag. A flag8-sphere requires at least18 vertices: the9 ridges of one facet require distinct outside vertices, since a repeated vertex would form a forbidden10-clique. Fifteen facets cannot suffice.','right_only_covers':[[list(p)for p in s]for s in S],'forced_right_pairs':[list(p)for p in forced],'minimum_cover_census':minimum,'right_contact_profiles':rows,'tree_types_checked':len(CACHE),'right_contact_proof':'Each of1408 profiles retains original finite labels, changes seven original right contacts to labels3/4, and allows arbitrary dotted old parallel entries. Optional additional right-to-dotted changes only delete spherical9 blocks. Exact unique-ridge peeling removes every potential vertex, so no closed compact8-sphere nerve can survive.','scope':'No claim about other finite additions, changed old parallel-to-finite incidences, larger changed supports, other labels, or compactH9 in general.'}

def main():
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=verify()
 if a.verify:assert d==json.loads(OUT.read_text());print('PASS exact saved replay')
 else:OUT.write_text(json.dumps(d,indent=2)+'\n');print('WROTE',OUT)
if __name__=='__main__':main()
