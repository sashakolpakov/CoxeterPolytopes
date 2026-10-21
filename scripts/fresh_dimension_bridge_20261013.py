#!/usr/bin/env python3
"""Construct the coupled Vinberg-H19 lift and simultaneous inherited-cusp repairs.

This producer builds exact 70-wall data over Q(sqrt(6)).  The numerical LP is
used only to propose an escape, which is accepted using Fraction arithmetic.
"""
from __future__ import annotations
import argparse, hashlib, json
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/vinberg_x3_h19.json'
OUT='data/enumerations/fresh_dimension_bridge_20261013.json'
Z=(F(0),F(0)); O=(F(1),F(0))
def k(a=0,b=0):return F(a),F(b)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def mul(x,y):return x[0]*y[0]+6*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def sc(x,a):return x[0]*a,x[1]*a
def sign(x):
 a,b=x;sa=(a>0)-(a<0);sb=(b>0)-(b<0)
 if not sb:return sa
 if not sa or sa==sb:return sb
 v=a*a-6*b*b;sv=(v>0)-(v<0)
 return 0 if not sv else sa if sv>0 else sb
def total(xs):
 z=Z
 for x in xs:z=add(z,x)
 return z
def dot(x,y):return total(mul(a,b)for a,b in zip(x,y,strict=True))
def mv(M,x):return tuple(dot(row,x)for row in M)
def pair(x,y,Q):return dot(x,mv(Q,y))
def enc(x):return list(map(str,x))
def dec(x):return tuple(map(F,x))
def encv(x):return list(map(enc,x))
def parsev(x):return tuple(map(dec,x))
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def value(x):return float(x[0])+6**.5*float(x[1])
def classify(p,qa,qb):
 assert sign(p)<=0
 if not sign(p):return '2'
 n=mul(p,p);d=mul(qa,qb)
 if sign(sub(n,d))==0:return 'parallel'
 if sign(sub(n,d))>0:return 'dotted'
 for m,f in ((3,F(1,4)),(4,F(1,2)),(6,F(3,4))):
  if n==sc(d,f):return str(m)
 raise AssertionError(('forbidden',enc(p),enc(qa),enc(qb)))

def build():
 import sympy as sp
 import numpy as np
 from scipy.optimize import linprog
 raw=json.loads((ROOT/SOURCE).read_text());G=-sp.Matrix(raw['wall_gram'])
 B=G[:20,:20];I=B.inv();W=(I*G[:20,:]).T
 assert W*B*W.T==G
 adj={i:[j for j in range(24,36)if j!=i and G[i,j]==-6]for i in range(24,36)}
 eps={24:1};todo=deque([24])
 while todo:
  i=todo.popleft()
  for j in adj[i]:
   if j not in eps:eps[j]=-eps[i];todo.append(j)
   assert eps[j]==-eps[i]
 assert len(eps)==12
 def rational(x):return F(int(x.p),int(x.q))
 Q=tuple(tuple(k(rational(x))for x in row)for row in sp.diag(B,3).tolist())
 old=[tuple(k(rational(x))for x in W[i,:])+(Z,)for i in range(36)]
 h=tuple(total(row[j]for row in old)for j in range(21))
 R=list(old[:24]);names=['short'+str(i)for i in range(24)];split={}
 for i in range(24,36):
  split[i]=[]
  for tag,height in [('A',2*eps[i]),('B',-eps[i])]:
   split[i].append(len(R));R.append(old[i][:-1]+(k(height),));names.append('long'+str(i)+tag)
 cusp_rows=[];T=k(2,F(1,2));missing=[]
 for ci,cusp in enumerate(raw['vertices_and_cusps']['cusps']):
  coeff=[F(0)]*36;delta=None
  for component in cusp['component_facets_one_based']:
   J=[j-1 for j in component];ker=G.extract(J,J).nullspace()[0]
   denominator=sp.ilcm(*[x.q for x in ker]);ker*=denominator
   if ker[0]<0:ker=-ker
   ker/=sp.igcd(*[abs(x)for x in ker]);assert all(x>0 for x in ker)
   v=sum((ker[i]*W[j,:]for i,j in enumerate(J)),sp.zeros(1,20))
   if delta is None:delta=v
   first=next(i for i,x in enumerate(v)if x);ratio=delta[first]/v[first];assert ratio>0 and ratio*v==delta
   for i,j in enumerate(J):coeff[j]+=rational(ratio*ker[i])
  d=tuple(k(rational(x))for x in delta)+(Z,)
  assert pair(d,d,Q)==Z and sign(pair(d,h,Q))<0
  incident=[j-1 for j in cusp['facets_one_based']]
  assert [i for i in range(36)if pair(d,old[i],Q)==Z]==incident
  weights={i:k(coeff[i])for i in incident if i<24}
  for i in incident:
   if i>=24:
    a,b=split[i];weights[a]=k(coeff[i]/3);weights[b]=k(2*coeff[i]/3)
  if max(incident)<24:
   missing.append(ci)
   for s in(1,-1):
    idx=len(R);R.append(tuple(mul(T,x)for x in d[:-1])+(k(s),));names.append('cusp'+str(ci)+('plus'if s>0 else'minus'));weights[idx]=O
  cusp_rows.append({'source_cusp_index':ci,'source_type':cusp['component_types'],'ray':encv(d),
                    'old_incident_walls':incident,'positive_null_weights':[[i,enc(v)]for i,v in sorted(weights.items())]})
 assert len(R)==70 and len(missing)==11
 norms=[pair(r,r,Q)for r in R]
 checks=[]
 for i,j in combinations(range(70),2):checks.append([i,j,classify(pair(R[i],R[j],Q),norms[i],norms[j])])
 assert sign(pair(h,h,Q))<0 and all(sign(pair(r,h,Q))<0 for r in R)
 for row in cusp_rows:
  d=parsev(row['ray']);inc=[i for i,r in enumerate(R)if pair(d,r,Q)==Z]
  assert all(sign(pair(d,r,Q))<=0 for r in R)
  assert inc==[i for i,x in row['positive_null_weights']]
  v=tuple(total(mul(dec(w),R[i][j])for i,w in row['positive_null_weights'])for j in range(21))
  assert all(sign(dec(w))>0 for i,w in row['positive_null_weights'])
  assert pair(v,v,Q)==Z and sign(pair(v,h,Q))<0
  # Numeric-free rank: every original cusp has rank19 including its null ray.
  # A changed long pair (or the new +/- cusp pair) supplies the new coordinate.
  original=[old[i]for i in row['old_incident_walls']]
  assert sp.Matrix([[x[0]for x in r]for r in original]).rank()==19
  assert any(R[i][-1]!=Z for i in inc)
  row['incident_walls']=inc;row['incident_normal_span_rank']=20;row['euclidean_rank']=19
  row['positive_null_vector']=encv(v)
 qn=np.array([[value(x)for x in row]for row in Q]);rn=np.array([[value(x)for x in row]for row in R]);hn=np.array([value(x)for x in h]);hn/=np.sqrt(-hn@qn@hn)
 A=rn@qn;rng=np.random.default_rng(19020);escape=None
 for objective in range(60):
  lp=linprog(rng.normal(size=21),A_ub=A,b_ub=np.full(70,-1e-7),A_eq=(hn@qn)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if lp.status or lp.x@qn@lp.x<=1e-6:continue
  v=tuple(add(k(F(round(float(x)*10**8),10**8)),sc(y,F(1,10**6)))for x,y in zip(lp.x,h))
  q=pair(v,v,Q);hp=pair(v,h,Q);ps=[pair(v,r,Q)for r in R]
  if sign(q)>0 and sign(hp)<0 and all(sign(z)<0 for z in ps):
   escape={'vector':encv(v),'norm':enc(q),'center_pairing':enc(hp),'all_wall_pairings':list(map(enc,ps)),'proposal_objective':objective};break
 assert escape
 data={'status':'EXACT_COXETER_H20_70_WALL_CONE_WITH_ALL_107_INHERITED_CUSPS_COMPLETE_BUT_STRICT_ESCAPE',
       'source':SOURCE,'source_sha256':sha(SOURCE),'dimension':20,'field':'Q(sqrt(6)), sqrt(6)>0',
       'basis':'Original first20 norm-two walls, followed by z with z^2=3',
       'form':[[enc(x)for x in row]for row in Q],'old_normals':list(map(encv,old)),
       'height_signs':[[i,eps[i]]for i in sorted(eps)],'normal_names':names,'normals':list(map(encv,R)),
       'norms':list(map(enc,norms)),'strict_center':encv(h),'center_norm':enc(pair(h,h,Q)),
       'pair_checks':checks,'pair_type_counts':dict(sorted(Counter(row[2]for row in checks).items())),
       'cusp_parameter':enc(T),'cusp_rows':cusp_rows,'repaired_cusps':missing,
       'strict_escape':escape,
       'scope':'70 actual acute Coxeter facets with strict interior; 107 inherited cusps have full rank19 Euclidean links. The strict escape proves infinite volume. This is not all cusps or a complete edge census; further compatible mirrors remain unrestricted.'}
 (ROOT/OUT).write_text(json.dumps(data,indent=2)+'\n');return data
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.parse_args();d=build();print(json.dumps({k:d[k]for k in ['status','dimension','pair_type_counts','center_norm']}))
