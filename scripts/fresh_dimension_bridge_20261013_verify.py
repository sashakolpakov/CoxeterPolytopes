#!/usr/bin/env python3
"""Producer-free exact reconstruction of the H19-to-H20 coupled lift.

Uses a separate quadratic-field class and rational row reduction. No producer,
SymPy, NumPy, LP solver, graph package, or numerical input is imported.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from collections import Counter,deque
from itertools import combinations
import argparse,hashlib,json,time,copy
ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_20261013.json'
OUT='data/enumerations/fresh_dimension_bridge_20261013_verify.json'
@dataclass(frozen=True)
class E:
 a:F=F(0)
 b:F=F(0)
 def __add__(x,y):return E(x.a+y.a,x.b+y.b)
 def __neg__(x):return E(-x.a,-x.b)
 def __sub__(x,y):return x+-y
 def __mul__(x,y):return E(x.a*y.a+6*x.b*y.b,x.a*y.b+x.b*y.a)
 def __truediv__(x,y):
  n=y.a*y.a-6*y.b*y.b;assert n
  return x*E(y.a/n,-y.b/n)
 def sign(x):
  if not x.b:return (x.a>0)-(x.a<0)
  if not x.a:return (x.b>0)-(x.b<0)
  if (x.a>0)==(x.b>0):return 1 if x.a>0 else-1
  if x.a*x.a==6*x.b*x.b:return 0
  larger=x.a if x.a*x.a>6*x.b*x.b else x.b
  return 1 if larger>0 else-1
 def enc(x):return [str(x.a),str(x.b)]
Z=E();O=E(F(1))
def parse(v):return E(*map(F,v))
def vec(v):return tuple(map(parse,v))
def scalar(n):return E(F(n))
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),Z)
def mv(M,x):return tuple(dot(row,x)for row in M)
def rank(M):
 a=[list(map(F,row))for row in M];n=len(a);width=len(a[0]);r=0
 for j in range(width):
  pivot=next((i for i in range(r,n)if a[i][j]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r];p=a[r][j];a[r]=[x/p for x in a[r]]
  for i in range(r+1,n):
   p=a[i][j]
   if p:a[i]=[x-p*y for x,y in zip(a[i],a[r])]
  r+=1
  if r==n:break
 return r

def inverse(M):
 n=len(M);a=[list(map(F,row))+[F(i==j)for j in range(n)]for i,row in enumerate(M)]
 for j in range(n):
  p=next(i for i in range(j,n)if a[i][j]);a[j],a[p]=a[p],a[j]
  p=a[j][j];a[j]=[x/p for x in a[j]]
  for i in range(n):
   if i==j:continue
   p=a[i][j]
   if p:a[i]=[x-p*y for x,y in zip(a[i],a[j])]
 return [row[n:]for row in a]

def inertia(M):
 a=[list(map(F,row))for row in M];pos=neg=zero=0
 while a:
  p=next((i for i in range(len(a))if a[i][i]),None)
  if p is None:
   off=next(((i,j)for i,j in combinations(range(len(a)),2)if a[i][j]),None)
   if off is None:zero+=len(a);break
   i,j=off
   # Replace e_i by e_i+e_j to obtain a nonzero diagonal without changing rank.
   old=[row[:]for row in a]
   for k in range(len(a)):a[i][k]=old[i][k]+old[j][k];a[k][i]=old[k][i]+old[k][j]
   a[i][i]=old[i][i]+2*old[i][j]+old[j][j];p=i
  a[0],a[p]=a[p],a[0]
  for row in a:row[0],row[p]=row[p],row[0]
  v=a[0][0];pos+=v>0;neg+=v<0
  a=[[a[i][j]-a[i][0]*a[0][j]/v for j in range(1,len(a))]for i in range(1,len(a))]
 return pos,neg,zero

def label(p,q,r):
 assert p.sign()<=0
 if p==Z:return'2'
 diff=p*p-q*r
 if diff==Z:return'parallel'
 if diff.sign()>0:return'dotted'
 for m,c in [(3,F(1,4)),(4,F(1,2)),(6,F(3,4))]:
  if p*p==scalar(c)*q*r:return str(m)
 raise AssertionError('forbidden pairing')

def check(d):
 rawpath=ROOT/d['source'];assert hashlib.sha256(rawpath.read_bytes()).hexdigest()==d['source_sha256']
 raw=json.loads(rawpath.read_text());G=[[-F(x)for x in row]for row in raw['wall_gram']]
 B=[row[:20]for row in G[:20]];I=inverse(B);assert inertia(B)==(19,1,0)
 q=[[scalar(B[i][j])if i<20 and j<20 else scalar(3)if i==j==20 else Z for j in range(21)]for i in range(21)]
 assert [[x.enc()for x in row]for row in q]==d['form']
 old=[tuple(scalar(sum((I[i][j]*G[j][n]for j in range(20)),F(0)))for i in range(20))+(Z,)for n in range(36)]
 assert old==list(map(vec,d['old_normals']))
 oldcov=[mv(q,r)for r in old]
 assert all(dot(old[i],oldcov[j])==scalar(G[i][j])for i in range(36)for j in range(36))
 eps={24:1};todo=deque([24])
 while todo:
  i=todo.popleft()
  for j in range(24,36):
   if G[i][j]!=-6:continue
   if j not in eps:eps[j]=-eps[i];todo.append(j)
   assert eps[j]==-eps[i]
 assert [[i,eps[i]]for i in sorted(eps)]==d['height_signs']
 roots=list(old[:24]);names=['short'+str(i)for i in range(24)]
 for i in range(24,36):
  for suffix,z in [('A',2*eps[i]),('B',-eps[i])]:roots.append(old[i][:-1]+(scalar(z),));names.append('long'+str(i)+suffix)
 h=tuple(sum((r[j]for r in old),Z)for j in range(21));assert vec(d['strict_center'])==h
 hc=mv(q,h);assert dot(h,hc)==E(F(-1518))==parse(d['center_norm'])
 rows=d['cusp_rows'];assert len(rows)==len(raw['vertices_and_cusps']['cusps'])==107
 t=E(F(2),F(1,2));assert parse(d['cusp_parameter'])==t;missing=[]
 for ci,(row,original)in enumerate(zip(rows,raw['vertices_and_cusps']['cusps'],strict=True)):
  assert ci==row['source_cusp_index']and row['source_type']==original['component_types']
  ds=vec(row['ray']);assert all(x.b==0 for x in ds)and ds[-1]==Z
  dp=mv(q,ds);assert dot(ds,dp)==Z and dot(h,dp).sign()<0
  inc=[j for j in range(36)if dot(old[j],dp)==Z]
  assert inc==row['old_incident_walls']==[j-1 for j in original['facets_one_based']]
  assert all(dot(r,dp).sign()<=0 for r in old)
  assert rank([[x.a for x in old[i]]for i in inc])==19
  if max(inc)<24:
   missing.append(ci)
   for s in(1,-1):roots.append(tuple(t*x for x in ds[:-1])+(scalar(s),));names.append('cusp'+str(ci)+('plus'if s>0 else'minus'))
 assert len(roots)==70 and roots==list(map(vec,d['normals']))and names==d['normal_names']and missing==d['repaired_cusps']
 cov=[mv(q,r)for r in roots];norms=[dot(r,p)for r,p in zip(roots,cov)]
 assert norms==list(map(parse,d['norms']))
 assert Counter(norms)==Counter({scalar(2):24,scalar(18):12,scalar(9):12,scalar(3):22})
 assert all(dot(r,hc).sign()<0 for r in roots)
 actual=[[i,j,label(dot(roots[i],cov[j]),norms[i],norms[j])]for i,j in combinations(range(70),2)]
 assert actual==d['pair_checks']and dict(sorted(Counter(row[2]for row in actual).items()))==d['pair_type_counts']
 cycle=scalar(4)*dot(roots[6],cov[48])*dot(roots[6],cov[48])/(norms[6]*norms[48])
 assert cycle==E(F(11,3),F(4,3))and (2*cycle.a).denominator!=1
 for row in rows:
  ds=vec(row['ray']);dp=mv(q,ds);ps=[dot(r,dp)for r in roots];inc=[i for i,p in enumerate(ps)if p==Z]
  assert all(p.sign()<=0 for p in ps)
  assert inc==row['incident_walls']==[i for i,v in row['positive_null_weights']]
  assert all(parse(w).sign()>0 for i,w in row['positive_null_weights'])
  v=tuple(sum((parse(w)*roots[i][j]for i,w in row['positive_null_weights']),Z)for j in range(21))
  assert v==vec(row['positive_null_vector'])
  k=next(i for i,x in enumerate(ds)if x!=Z);mu=v[k]/ds[k]
  assert mu.sign()>0 and tuple(mu*x for x in ds)==v
  # The old rank19 span is contained in the new incident span: split pairs
  # recover every old incident long root and supply z. For short-only cusps,
  # the opposite added pair supplies z instead. Thus actual span rank is20.
  for oldi in row['old_incident_walls']:
   if oldi>=24:
    a=24+2*(oldi-24);b=a+1
    assert a in inc and b in inc
    assert tuple((roots[a][j]+scalar(2)*roots[b][j])/scalar(3)for j in range(21))==old[oldi]
  assert any(roots[i][-1]!=Z for i in inc)
  assert row['incident_normal_span_rank']==20 and row['euclidean_rank']==19
 y=vec(d['strict_escape']['vector']);yc=mv(q,y);n=dot(y,yc);ch=dot(y,hc)
 assert n.sign()>0 and ch.sign()<0
 assert n==parse(d['strict_escape']['norm'])and ch==parse(d['strict_escape']['center_pairing'])
 ps=[dot(r,yc)for r in roots];assert all(p.sign()<0 for p in ps)
 assert ps==list(map(parse,d['strict_escape']['all_wall_pairings']))
 return {'status':'PASS_PRODUCER_FREE_H20_COUPLED_LIFT_AND_107_CUSP_REPAIRS','source':SOURCE,
         'source_sha256':hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),'walls':70,'dimension':20,
         'form_inertia':[20,1,0],'all_pairs':len(actual),'inherited_complete_cusps':107,
         'strict_center_norm':dot(h,hc).enc(),'strict_escape_norm':n.enc(),
         'nonintegral_length_two_cycle':{'walls':[6,48],'value':cycle.enc(),'field_trace':str(2*cycle.a),'conditional_consequence':'Any finite-volume completion retaining these facets is nonarithmetic.'},
         'actual_facets_reason':'Every pair is nonpositive and the center is strict timelike, so its projection to each wall is timelike and strict on all other walls.',
         'cusp_reason':'At each archived ray, every incident normal occurs in a strictly positive combination equal to a positive multiple of the ray; their actual span has rank20. The Euclidean link therefore has full rank19 and is bounded.',
         'scope':'All70 normal identities and Coxeter pairs, strict center, all107 inherited bounded cusp links, and strict spacelike escape replayed. This is not an exhaustive cusp/edge census or a no-completion theorem.'}

def main():
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();start=time.monotonic()
 data=json.loads((ROOT/SOURCE).read_text());out=check(data)
 bad=copy.deepcopy(data);bad['normals'][-1][-1]=['2','0']
 try:check(bad)
 except AssertionError:pass
 else:raise AssertionError('forged root accepted')
 bad=copy.deepcopy(data);bad['strict_escape']['vector']=data['strict_center']
 try:check(bad)
 except AssertionError:pass
 else:raise AssertionError('timelike escape accepted')
 out['forged_normal_and_timelike_escape_rejected']=True
 if a.verify:assert out==json.loads((ROOT/OUT).read_text())
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out),flush=True);print('seconds',round(time.monotonic()-start,3),flush=True)
if __name__=='__main__':main()
