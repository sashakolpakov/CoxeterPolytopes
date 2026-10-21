#!/usr/bin/env python3
"""Exact geometric comparison of affine ADE product alcoves of rank nineteen."""
import argparse,hashlib,json
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_intrinsic_cusp_benchmark_20261014'
def inv(A):
 n=len(A);M=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  k=next(i for i in range(j,n)if M[i][j]);M[j],M[k]=M[k],M[j]
  q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(n):
   if i!=j:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def pd(A):
 A=[[F(x)for x in r]for r in A];p=[]
 while A:
  q=A[0][0];assert q>0;p.append(q)
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,len(A))]for i in range(1,len(A))]
 return p
def cartan(kind,n):
 if kind=='A':edges=[(i,i+1)for i in range(n-1)]
 else:
  edges=[];j=1
  for length in ([1,1,n-3]if kind=='D'else[1,2,n-4]):
   old=0
   for _ in range(length):edges.append((old,j));old=j;j+=1
  assert j==n
 C=[[2*int(i==j)for j in range(n)]for i in range(n)]
 for i,j in edges:C[i][j]=C[j][i]=-1
 return C
def factor(kind,n):
 C=cartan(kind,n);piv=pd(C)
 roots={tuple(int(i==j)for i in range(n))for j in range(n)};todo=list(roots)
 while todo:
  a=todo.pop()
  for i in range(n):
   b=list(a);b[i]-=sum(x*y for x,y in zip(C[i],a));b=tuple(b)
   if b not in roots:roots.add(b);todo.append(b)
 expected=n*(n+1)if kind=='A'else 2*n*(n-1)if kind=='D'else{6:72,7:126,8:240}[n]
 assert len(roots)==expected
 assert all(sum(a[i]*C[i][j]*a[j]for i in range(n)for j in range(n))==2 for a in roots)
 height=max(map(sum,roots));highest=[r for r in roots if sum(r)==height];assert len(highest)==1
 theta=highest[0];assert all(x>0 for x in theta)
 b=[-sum(C[i][j]*theta[j]for j in range(n))for i in range(n)]
 assert all(x<=0 for x in b)
 affine=[[2]+b]+[[b[i]]+C[i]for i in range(n)];marks=(1,)+theta;h=sum(marks)
 assert all(sum(x*y for x,y in zip(row,marks))==0 for row in affine)
 I=[[2*x for x in r]for r in inv(C)]
 rhs=[[F(1)for _ in range(n)]]+[[F(1)-(F(h,theta[j])if i==j else 0)for i in range(n)]for j in range(n)]
 dot=lambda a,b:sum(a[i]*I[i][j]*b[j]for i in range(n)for j in range(n))
 # Facet0 is minus the highest-root combination of the finite normals.
 for omit,a in enumerate(rhs):
  assert -sum(theta[i]*a[i]for i in range(n))==(1-F(h)if omit==0 else 1)
  for i in range(n):assert a[i]==(1-F(h,theta[i])if omit==i+1 else 1)
 ns=[dot(a,a)for a in rhs]
 distances=[dot([x-y for x,y in zip(a,b)],[x-y for x,y in zip(a,b)])for i,a in enumerate(rhs)for b in rhs[:i]]
 assert min(distances)>0
 return dict(name=f'{kind}{n}',rank=n,finite_Cartan=C,finite_LDL_pivots=piv,
             finite_root_count=len(roots),affine_Cartan=affine,positive_marks=marks,
             unit_root_basis_inverse=I,vertex_pairings_with_finite_unit_roots=rhs,
             vertex_squares=ns,minimum_vertex_square=min(ns),maximum_vertex_square=max(ns),
             minimum_separation_square=min(distances),minimum_omitted_facet_distance=F(h,max(marks)))
def summary(chosen):
 mn=sum(r['minimum_vertex_square']for r in chosen);mx=sum(r['maximum_vertex_square']for r in chosen)
 sep=min(r['minimum_separation_square']for r in chosen)
 return dict(factors=[r['name']for r in chosen],factor_count=len(chosen),seed_facet_count=19+len(chosen),
             product_vertex_count=prod(r['rank']+1 for r in chosen),minimum_vertex_square=mn,
             maximum_vertex_square=mx,minimum_separation_square=sep,
             necessary_equal_cap_radius_squared_bound=sep/2,
             parallel_nearest_neighbor_radius_squared=sep/4,
             worst_vertex_square_over_acute_radius_bound=2*mx/sep)
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def run():
 rows=[factor(k,n)for k,ns in [('A',range(1,20)),('D',range(4,20)),('E',[6,7,8])]for n in ns]
 at={r['name']:r for r in rows};parts=[]
 def rec(start,left,chosen):
  if not left:parts.append(summary(chosen));return
  if len(chosen)==4:return
  for i in range(start,len(rows)):
   r=rows[i]
   if r['rank']<=left:rec(i,left-r['rank'],chosen+[r])
 rec(0,19,[])
 p3=[p for p in parts if p['factor_count']<=3]
 assert len(rows)==38 and len(p3)==229 and len(parts)==607
 best3=min(p3,key=lambda p:p['maximum_vertex_square']);best4=min(parts,key=lambda p:p['maximum_vertex_square'])
 assert best3['factors']==['A6','A6','A7']and best3['maximum_vertex_square']==196
 assert best4['factors']==['A4','A5','A5','A5']and best4['maximum_vertex_square']==125
 assert sum(p['maximum_vertex_square']==196 for p in p3)==1
 assert sum(p['maximum_vertex_square']==125 for p in parts)==1
 # A simple prospective parallel layer for the selected balanced seed.
 # These are geometric inequalities, not a full chamber/volume certificate.
 layer_R2=best3['parallel_nearest_neighbor_radius_squared'];assert layer_R2==21
 assert all(at[name]['minimum_omitted_facet_distance']**2>layer_R2 for name in best3['factors'])
 central_gap=best3['minimum_vertex_square']-4-layer_R2
 assert central_gap>0 and central_gap**2>16*layer_R2
 chosen=[['A19'],['A9','A10'],['A6','A6','A7'],['A4','A5','A5','A5'],
         ['A11','A2','E6'],['E8','E8','A3'],['E8','E7','D4'],['E7','E6','E6'],['A1']*19]
 p=str(Path(__file__).resolve().relative_to(ROOT))
 return enc(dict(status='PASS_EXACT_ADE_ALCOVE_GEOMETRY_AND_BOUNDED_PARTITION_BENCHMARK',
                 sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()},factors=rows,
                 rank19_partitions_at_most4_factors=parts,counts=dict(factor_types=38,partitions_at_most3=229,partitions_at_most4=607),
                 best_maximum_vertex_square_at_most3=best3,best_maximum_vertex_square_at_most4=best4,
                 chosen_parallel_layer_necessary_geometry=dict(radius_squared=layer_R2,
                    all_omitted_seed_relations_dotted=True,central_relation_dotted=True,
                    pair_separation_squared_lower_bound=best3['minimum_separation_square']),
                 selected_comparisons=[summary([at[x]for x in names])for names in chosen],
                 scope='Exact unit-inradius Euclidean geometry, not arithmetic root enumeration or a hyperbolic completion. Radius bound is necessary acuteness only; it need not yield Coxeter cap angles or coverage. Exhaustive partition scope is at most4 factors; A1^19 is one extra comparison.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['counts'])
