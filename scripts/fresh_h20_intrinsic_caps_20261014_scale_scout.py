#!/usr/bin/env python3
"""Exact small factor-distance benchmark; no unrestricted cap classification."""
import argparse, hashlib, json
from fractions import Fraction as F
from itertools import product
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_intrinsic_caps_20261014_scale_scout'
SRC='results/fresh_h20_intrinsic_caps_20261014_endpoint.json'
def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[0],x[1]-y[1]
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def sc(x,a):return x[0]*a,x[1]*a
def sign(x):
 a,b=x[0]+x[1]/2,x[1]/2
 if not b:return (a>0)-(a<0)
 if not a:return (b>0)-(b<0)
 if (a>0)==(b>0):return (a>0)-(a<0)
 d=a*a-5*b*b;assert d
 return ((a>0)-(a<0))if d>0 else((b>0)-(b<0))
def dot(x,y):return sum(a*b for a,b in zip(x,y,strict=True))
def inv(A):
 n=len(A);M=[list(r)+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  k=next(i for i in range(j,n)if M[i][j]);M[j],M[k]=M[k],M[j]
  q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(n):
   if i!=j:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(tuple,list)):return [enc(v)for v in x]
 if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
 return x
def run():
 d=json.loads((ROOT/SRC).read_text());ids=d['seed_order'];at={w:i for i,w in enumerate(ids)}
 G=[[F(x,2)for x in r]for r in d['normalized_seed_Gram']]
 marks=dict(zip(ids,d['positive_affine_marks']))
 factors=[[4,12,14,21,25,28,30,31,39,48,49,63],[23,50,200],[22,24,26,33,34,40,45]]
 records=[]
 for f in factors:
  I=inv([[G[at[i]][at[j]]for j in f[1:]]for i in f[1:]])
  total=sum(marks[i]for i in f)
  vs=[[F(1)-(F(total,marks[o])if w==o else 0)for w in f[1:]]for o in f]
  ip=lambda x,y:dot(x,[dot(r,y)for r in I])
  ns=[ip(x,x)for x in vs]
  ds=[ip([a-b for a,b in zip(x,y)],[a-b for a,b in zip(x,y)])for i,x in enumerate(vs)for y in vs[:i]]
  target=[x/4 for x in vs[0]]
  qdist={}
  for name,scale in [('vertices',F(1)),('half_vertices',F(1,2))]:
   qdist[name]=min(ip([a-scale*b for a,b in zip(target,y)],[a-scale*b for a,b in zip(target,y)])for y in vs)
  records.append(dict(walls=f,vertex_squares=ns,minimum_separation_square=min(ds),quarter_target_square=ip(target,target),quarter_target_distances=qdist))
 assert [r['minimum_separation_square']for r in records]==[264,12,48]
 assert [min(r['vertex_squares'])for r in records]==[286,4,12]
 assert all(r['quarter_target_distances']['half_vertices']==r['quarter_target_square']for r in records)
 assert all(r['quarter_target_distances']['vertices']==9*r['quarter_target_square']for r in records)
 sq={3:(F(1),F(0)),4:(F(2),F(0)),5:(F(1),F(1)),6:(F(3),F(0))}
 rows=[];best=None
 for labels in product(sq,repeat=3):
  radii=[sq[m]for m in labels]
  summ=lambda coeff:sumq([sc(x,c)for x,c in zip(radii,coeff)])
  n=summ([286,4,12]);rho2=sc(radii[1],6);parallel2=sc(radii[1],3);q0=sc(n,F(1,16))
  assert sign(sub(sc(radii[0],264),sc(radii[1],12)))>0
  assert sign(sub(sc(radii[2],48),sc(radii[1],12)))>0
  powers=dict(central=sub(q0,(F(4),F(0))),half_vertex=sub(q0,(F(1),F(0))),vertex=sub(sc(q0,9),rho2))
  assert all(sign(sub(x,(F(1,8),F(0))))>0 for x in powers.values())
  rows.append(dict(central_labels=labels,inradii_squared=radii,minimum_vertex_square=n,
                   necessary_equal_vertex_radius_squared_upper_bound=rho2,
                   explicit_parallel_vertex_layer_radius_squared=parallel2,quarter_ball_powers=powers))
  if best is None or sign(sub(mul(n,best[1]),mul(best[0],rho2)))<0:best=(n,rho2,labels)
 assert best[2]==(3,6,3)
 chosen=next(r for r in rows if r['central_labels']==best[2])
 assert chosen['quarter_ball_powers']==dict(central=(F(123,8),F(0)),half_vertex=(F(147,8),F(0)),vertex=(F(1251,8),F(0)))
 paths=[SRC,str(Path(__file__).resolve().relative_to(ROOT))]
 return enc(dict(status='PASS_EXACT_FACTOR_SCALE_BENCHMARK_64_CHOICES',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},factor_geometry=records,choices=rows,best_by_minimum_vertex_square_over_acute_radius_upper_bound=chosen,scope='Small geometric benchmark for central labels3,4,5,6 at radius2. All64 center/equal-radius-vertex/radius-one-half-vertex templates retain a quarter-vertex strict escape even at the necessary acute radius bound; that bound need not be Coxeter feasible. Parallel-layer radius is recorded separately and is not a maximum. No exclusion of other cap centers, varying-radius selections, or unrestricted completions.'))
def sumq(xs):
 out=(F(0),F(0))
 for x in xs:out=add(out,x)
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
