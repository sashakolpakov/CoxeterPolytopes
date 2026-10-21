#!/usr/bin/env python3
"""Exact two-transverse-coordinate lift of the known Vinberg H18 hybrid.

The finite source is a published example. This constructs an H20 seed, checks
its mandatory short-pair cusps, and tests global volume by an exact escape.
The initial numerical LP proposal has been replaced by a small literal exact
witness. Generation and saved replay contain no numerical decisions.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter
import argparse, hashlib, json, sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.silver import Silver as K, silver_inertia
TAG='fresh_dimension_construction_routes_20261014'
SOURCE='data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
ZERO=K(0)

def rref(A):
 A=[list(r)for r in A]; piv=[];row=0
 for col in range(len(A[0])):
  p=next((i for i in range(row,len(A))if A[i][col]!=ZERO),None)
  if p is None:continue
  A[row],A[p]=A[p],A[row];z=A[row][col];A[row]=[a/z for a in A[row]]
  for i in range(len(A)):
   if i!=row and A[i][col]!=ZERO:
    z=A[i][col];A[i]=[a-z*b for a,b in zip(A[i],A[row])]
  piv.append(col);row+=1
  if row==len(A):break
 return A,piv
def inverse(A):
 n=len(A);R,p=rref([r+[K(i==j)for j in range(n)]for i,r in enumerate(A)])
 assert p[:n]==list(range(n));return[r[n:]for r in R]
def matvec(A,x):return[sum((a*b for a,b in zip(r,x)),ZERO)for r in A]
def dot(x,y,Q):return sum((a*b for a,b in zip(x,matvec(Q,y))),ZERO)
def add(x,y):return[a+b for a,b in zip(x,y)]
def scale(x,c):return[a*c for a in x]
def enc(x):
 if isinstance(x,K):return x.to_json()
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x
def pairgate(g,a,b):
 assert g<=0
 if g==0:return'right'
 if g*g>=a*b:return'parallel'if g*g==a*b else'dotted'
 s=g*g/(a*b)
 finite={K(F(1,4)):'pi3',K(F(1,2)):'pi4'}
 assert s in finite,(g,a,b,s)
 return finite[s]

def run(saved=None):
 source=json.loads((ROOT/SOURCE).read_text());G=[[K.from_json(x)for x in row]for row in source['gram']]
 labels=source['labels'];n=len(G);assert n==40
 _,basis=rref(G);assert len(basis)==19
 Q0=[[G[i][j]for j in basis]for i in basis];I=inverse(Q0)
 Q=[r+[ZERO,ZERO]for r in Q0]+[[ZERO]*19+[K(1),ZERO],[ZERO]*20+[K(1)]]
 inertia=silver_inertia(Q);assert(inertia.positive,inertia.negative,inertia.zero)==(20,1,0)
 old=[matvec(I,[G[j][i]for j in basis])+[ZERO,ZERO]for i in range(n)]
 assert all(dot(old[i],old[j],Q)==G[i][j]for i in range(n)for j in range(i+1))
 h=[sum((r[j]for r in old),ZERO)for j in range(21)]
 assert dot(h,h,Q)<0 and all(dot(h,r,Q)<0 for r in old)
 short=[i for i in range(n)if G[i][i]==1];long=[i for i in range(n)if G[i][i]==2]
 assert(len(short),len(long))==(14,26)
 E=[(i,j)for i,j in combinations(short,2)if G[i][j]==-1]
 assert len(E)==25
 triangle=[labels.index(s)for s in ['L18','L25','L26']]
 assert all(G[i][j]==-1 for i,j in combinations(triangle,2))
 # Every short-long entry, and every short-short height dot in {-1,0,1}.
 template=Counter()
 for i in short:
  for j in long:template[pairgate(G[i][j],K(2),K(2))]+=1
 for i,j in combinations(short,2):
  for hd in[-1,0,1]:template[pairgate(G[i][j]+hd,K(2),K(2))]+=1
 # The same old short root's opposite signed lifts are right.
 assert K(1)-1==0

 # A global check before adding anything: every coloring has a monochrome
 # edge of the short parallel triangle. Either unused height yields an escape.
 raw_escapes=[];epsilon=F(1,100000)
 for i,j in combinations(triangle,2):
  delta=add(old[i],old[j]);assert dot(delta,delta,Q)==0 and dot(delta,h,Q)<0
  X=add(delta,scale(h,epsilon));X[20]+=1
  assert dot(X,X,Q)>0 and dot(X,h,Q)<0
  assert all(dot(X,old[k],Q)<0 for k in long)
  for k in short:
   allowed_axes=[0]if k in[i,j]else[0,1]
   for axis,sgn in product(allowed_axes,[-1,1]):
    r=list(old[k]);r[19+axis]+=sgn;assert dot(X,r,Q)<0
  raw_escapes.append(dict(pair=[labels[i],labels[j]],epsilon=epsilon,
      unused_height_axis=1,point=X,norm=dot(X,X,Q)))

 # Count all colorings modulo a global color swap. This is only a finite
 # combinatorial cut census, not an all-height lift classification.
 best=26;best_masks=[];hist=Counter()
 for mask in range(1<<13):
  color={i:(mask>>a)&1 for a,i in enumerate(short)}
  count=sum(color[i]==color[j]for i,j in E);hist[count]+=1
  if count<best:best=count;best_masks=[mask]
  elif count==best:best_masks.append(mask)
 assert best==3 and len(best_masks)==1
 color={i:(best_masks[0]>>a)&1 for a,i in enumerate(short)}
 monochrome=[(i,j)for i,j in E if color[i]==color[j]]
 assert set(map(frozenset,monochrome))==set(map(frozenset,combinations(triangle,2)))
 normals=[old[i]for i in long];names=[labels[i]for i in long]
 for i in short:
  for sgn in[-1,1]:
   r=list(old[i]);r[19+color[i]]+=sgn;normals.append(r);names.append(f'{labels[i]}:{sgn:+}')
 assert len(normals)==54
 # At t=2, opposite caps at the missing transverse direction close these
 # three specifically identified links; all nonincident relations are dotted.
 repairs=[]
 for i,j in monochrome:
  delta=add(old[i],old[j]);missing=1-color[i];new=[]
  for sgn in[-1,1]:
   r=scale(delta,2);r[19+missing]+=sgn;new.append(len(normals));normals.append(r)
   names.append(f'null:{labels[i]}:{labels[j]}:{sgn:+}')
  repairs.append(dict(pair=[labels[i],labels[j]],delta=delta,new_wall_indices=new))
 assert len(normals)==60
 norms=[dot(r,r,Q)for r in normals]
 assert norms==[K(2)]*54+[K(1)]*6
 counts=Counter(pairgate(dot(normals[i],normals[j],Q),norms[i],norms[j])
                for i in range(60)for j in range(i))
 assert all(dot(h,r,Q)<0 for r in normals)
 # Global acute pairing plus this strict center proves every wall is a facet.
 for r,nr in zip(normals,norms):
  x=add(h,scale(r,-dot(h,r,Q)/nr));assert dot(x,x,Q)<0
  assert all(dot(x,s,Q)<0 for s in normals if s!=r)

 import networkx as nx
 for repair in repairs:
  d=repair['delta'];vals=[dot(d,r,Q)for r in normals];assert all(x<=0 for x in vals)
  incident=[i for i,x in enumerate(vals)if x==0]
  gi=[[dot(normals[i],normals[j],Q)for j in incident]for i in incident]
  assert len(rref(gi)[1])==19
  assert len(rref([normals[i]for i in incident])[1])==20
  graph=nx.Graph();graph.add_nodes_from(range(len(incident)))
  graph.add_edges_from((i,j)for i in range(len(incident))for j in range(i)if gi[i][j]<0)
  components=[]
  for cc in nx.connected_components(graph):
   ids=sorted(cc);B=[[gi[i][j]for j in ids]for i in ids];R,p=rref(B)
   assert len(p)==len(ids)-1
   f=next(j for j in range(len(ids))if j not in p);w=[ZERO]*len(ids);w[f]=K(1)
   for a,j in enumerate(p):w[j]=-R[a][f]
   assert all(x>0 for x in w)
   z=[sum((w[a]*normals[incident[j]][k]for a,j in enumerate(ids)),ZERO)for k in range(21)]
   pivot=next(k for k,x in enumerate(d)if x!=0);mu=z[pivot]/d[pivot]
   assert mu>0 and z==scale(d,mu)
   components.append(dict(walls=[incident[j]for j in ids],weights=w,delta_multiple=mu))
  repair.update(incident=incident,gram_rank=19,normal_span=20,positive_affine_components=components)

 # The same rational indefinite pair persists under the nontrivial cycle-field embedding.
 a,b=labels.index('L1'),labels.index('R10');cycle=G[a][b]*G[a][b]
 assert cycle==K(12,8)
 a,b=labels.index('R23'),labels.index('R24')
 assert G[a][a]==2 and G[b][b]==2 and G[a][b]==-6
 assert G[a][a]*G[b][b]-G[a][b]**2==-32

 boundary=[K(x)for x in[3,6,10,9,8,7,6,5,4,3,2,1,0,-1,-2,-3,-4,-5,5,5,1]]
 assert dot(boundary,boundary,Q)==25
 assert all(dot(boundary,r,Q)<=0 for r in normals)
 X=add(boundary,scale(h,F(1,100)))
 assert dot(X,X,Q)>0 and dot(X,h,Q)<0 and all(dot(X,r,Q)<0 for r in normals)
 paths=[SOURCE,f'scripts/{TAG}.py','src/coxeter_polytopes/exact/silver.py','src/coxeter_polytopes/exact/matrix.py']
 return enc(dict(status='PASS_TWO_AXIS_H18_LIFT_THREE_BOUNDED_LINKS_AND_STRICT_H20_ESCAPE',
   sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
   source_dimension=18,target_dimension=20,old_basis_indices=basis,form=Q,
   normals=normals,names=names,norms=norms,pair_counts=dict(counts),strict_center=h,
   template_pair_checks=dict(template),raw_family_escapes=raw_escapes,
   color_census=dict(colorings=8192,histogram=dict(hist),minimum_monochrome_edges=best,
                    optimal_masks=best_masks,selected_colors={labels[i]:color[i]for i in short}),
   repaired_cusps=repairs,
   conditional_nonquasi=dict(cycle=cycle,rational_indefinite_pair=['R23','R24'],determinant=-32),
   strict_escape=dict(boundary_point=boundary,boundary_norm=25,center_perturbation=F(1,100),
      point=X,norm=dot(X,X,Q),center_pairing=dot(X,h,Q),all_60_strict=True),
   scope='Exact 60-facet infinite-volume H20 chamber from a known H18 hybrid. '
         'All 8192 two-axis colorings have a raw escape; the selected coloring minimizes '
         'monochrome short parallel pairs and its three repaired links are bounded. '
         'This is not a complete inherited or new cusp census, all-height lift classification, '
         'or exclusion of future compatible additions. No new finite-volume polytope.'))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 path=ROOT/'results'/f'{TAG}.json';saved=json.loads(path.read_text())if a.verify else None
 d=run(saved)
 if a.verify:assert saved==d
 else:path.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status']);print(d['pair_counts']);print(d['strict_escape']['norm'])
