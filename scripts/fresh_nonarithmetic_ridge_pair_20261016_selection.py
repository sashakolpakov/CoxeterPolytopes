#!/usr/bin/env python3
"""Numerical candidate selection with an exact tentative-diagram affine filter."""
import argparse,itertools,json,math
from collections import Counter
from hashlib import sha256
from functools import lru_cache
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_ridge_pair_20261016_polar as P
import fresh_h9_persistent_cusp_affine_20261016 as A
ROOT=P.ROOT;TAG='fresh_nonarithmetic_ridge_pair_20261016_selection';SC='results/fresh_nonarithmetic_ridge_pair_20261016_cutter_scout.json';INV=A.INPUT
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
@lru_cache(None)
def affine_check(gramlabels):
 matrix=tuple(tuple(A.scale(A.ONE,2)if i==j else A.VALUES[m]for j,m in enumerate(row))for i,row in enumerate(gramlabels))
 return A.classify(matrix)
def run():
 Q,N,H,E=P.geometry();s=P.P.B.read(SC);cl=P.P.B.read(P.CL);cone=P.P.B.read(P.CO);old=P.P.B.read(INV)
 L=[r+[0]*3 for r in old['pair_matrix']]+[[0]*52 for _ in range(3)]
 for c in cl['caps']:
  j=c['index'];L[j][j]=2
  for r in c['relations']:L[j][r['wall']]=L[r['wall']][j]=r.get('label',0)
 records=[];normals=[];basis=[0,1,2,3,4,5,6,7,8,11];M=mp.matrix([list(Q*N[i])for i in basis]);Mi=M**-1
 for index,r in enumerate(s['candidates']):
  vals=r['normalized_wall_pairings'];labels=[]
  for g in vals:
   if abs(g)<1e-6:labels.append(2)
   elif g<-2-1e-6:labels.append(0)
   else:
    m=math.pi/math.acos(-g/2);assert abs(m-round(m))<1e-5;labels.append(round(m))
  T=[row+[labels[i]]for i,row in enumerate(L)]+[labels+[2]];ids=[i for i,m in enumerate(labels)if m in A.VALUES]
  queue=[(52,)];seen=set(queue);aff=[];counts=Counter()
  while queue:
   ss=queue.pop();mat=tuple(tuple(T[i][j]for j in ss)for i in ss);kind,rank,nullity=affine_check(mat);counts[kind]+=1
   if kind=='positive_semidefinite':aff.append(list(ss));continue
   if kind!='positive_definite':continue
   for new in ids:
    if new in ss:continue
    ms=[T[new][j]for j in ss]
    if any(m not in A.VALUES for m in ms)or all(m==2 for m in ms):continue
    nxt=tuple(sorted(ss+(new,)))
    if nxt not in seen:seen.add(nxt);queue.append(nxt)
  blockers=[]
  for ss in aff:
   orth=[i for i in range(53)if i not in ss and all(T[i][j]==2 for j in ss)]
   for i,j in itertools.combinations(orth,2):
    if T[i][j]not in[0,2,3,4,6]:blockers.append(dict(affine=ss,forbidden_pair=[i,j],label=T[i][j]))
  point=Mi*mp.matrix([mp.mpf(vals[i])for i in basis]);normals.append(point)
  records.append(dict(index=index,recipe={k:v for k,v in r.items()if k not in['normalized_wall_pairings','bad_finite','near_parallel']},finite_contacts=[[i,m]for i,m in enumerate(labels)if m],affine_circuits=aff,permanent_cusp_blockers=blockers,affine_subset_counts=dict(counts)))
 print('affine filtered',dict(Counter('blocked'if r['permanent_cusp_blockers']else'affine_nonblocked'if r['affine_circuits']else'no_affine'for r in records)),flush=True)
 cov=[Q*n for n in N];hcov=Q*H;rays=[]
 for ids in cone['complete_census']['hyperideal_rays']:
  v=mp.matrix(P.P.null([list(cov[j])for j in ids]));v/=-(v.T*hcov)[0];rays.append([float(x)for x in v])
 V=np.array(rays);Qf=np.array(Q.tolist(),dtype=float);RV=np.array([[float(x)for x in r]for r in normals]);coverage=RV@Qf@V.T
 for r,g in zip(records,coverage):
  r['numerically_cut_ray_indices']=np.flatnonzero(g>1e-7).tolist();r['numerically_cut_ray_count']=len(r['numerically_cut_ray_indices'])
 keep=[r['index']for r in records if not r['affine_circuits']];pair=RV@Qf@RV.T;edges=[]
 for i,j in itertools.combinations(keep,2):
  g=pair[i,j]
  if g<-2-1e-6:edges.append(dict(indices=[i,j],kind='dotted',pairing=float(format(g,'.12g'))))
  elif abs(g)<1e-6:edges.append(dict(indices=[i,j],kind='right',pairing=float(format(g,'.12g'))))
  elif -2+1e-6<g<-1e-6:
   m=math.pi/math.acos(-g/2)
   if abs(m-round(m))<1e-5:edges.append(dict(indices=[i,j],kind='finite',label=round(m),pairing=float(format(g,'.12g'))))
 order=sorted(keep,key=lambda i:(-records[i]['numerically_cut_ray_count'],-records[i]['recipe']['escape_pairing']))
 print('top',[(i,records[i]['numerically_cut_ray_count'],records[i]['recipe']['escape_pairing'])for i in order[:20]],'pair edges',len(edges),flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',f'scripts/{A.TAG}.py',SC,INV,P.CO,P.CL,P.P.BATCH]
 return dict(status='NUMERICAL_CUTTER_SELECTION_WITH_EXACT_TENTATIVE_AFFINE_DIAGRAM_CHECK',sources={p:sha(p)for p in paths},records=records,no_affine_candidate_indices=keep,ranking=order,numerical_compatibility=edges,scope='The new normals and angle assignments remain numerical proposals. Conditional on each tentative literal contact graph, the connected affine-subdiagram scan uses exact Fraction arithmetic in Q(sqrt2,sqrt3). Coverage and mutual compatibility use binary64; all promoted caps, angle relations and cusp conclusions require exact reconstruction.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved candidate selection')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['no_affine_candidate_indices']),len(d['numerical_compatibility']))
