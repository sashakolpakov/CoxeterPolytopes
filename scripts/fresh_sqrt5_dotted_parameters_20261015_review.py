#!/usr/bin/env python3
"""Independent complete dotted parameter box, using tuple-Fraction Q(sqrt5).
No producer imports: the previously independent radical-gauge helper is reused.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_quad_filters_20261015_review as R
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_sqrt5_dotted_parameters_20261015';TAG=BASE+'_review'
Z,ONE,TWO=R.Z,R.ONE,R.TWO
PHI=(F(1,2),F(1,2));FOUR=(F(4),F(0))
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dec(x):return tuple(map(F,x))
def q(x):return F(x),F(0)
def dot(x,y):return sumfield(R.mul(a,b)for a,b in zip(x,y))
def sumfield(values):
 s=Z
 for x in values:s=R.add(s,x)
 return s
def edges(row):return dict(zip(combinations(range(row['nodes']),2),row['labels']))
def finite_row(row):return dict(nodes=row['nodes'],white=row['white'],labels=[2 if x==0 else x for x in row['labels']])
def component(row,start):
 g=edges(row);seen={start};todo=[start]
 for i in todo:
  for j in range(row['nodes']):
   if j!=i and g[tuple(sorted((i,j)))] not in(0,2)and j not in seen:seen.add(j);todo.append(j)
 return seen
def matrix(row,ratio):
 r=finite_row(row);colors=R.coloring(r)
 if colors is None:return None,'REJECT_FINITE_CYCLE_FIELD'
 M=R.matrix(r,colors);z=R.mul(ratio,R.mul(M[0][0],M[2][2]));I=component(row,0)
 if 2 in I:
  try:x=R.square_root(z)
  except AssertionError:return None,'REJECT_DOTTED_CYCLE_FIELD'
  M[0][2]=M[2][0]=R.neg(x)
 else:
  if R.sign(R.conj(z))<=0:return None,'REJECT_DOTTED_CONJUGATE_SQUARE'
  for i in I:
   for j in I:M[i][j]=R.mul(M[i][j],z)
  M[0][2]=M[2][0]=R.neg(z)
 assert R.div(R.mul(FOUR,R.mul(M[0][2],M[0][2])),R.mul(M[0][0],M[2][2]))==R.mul(FOUR,ratio)
 return M,None

def canonical_cycle(c):
 c=list(c);i=c.index(min(c));a=c[i:]+c[:i];b=[a[0]]+list(reversed(a[1:]));return tuple(min(a,b))
def integral(x):return (2*x[1]).denominator==1 and (x[0]-x[1]).denominator==1

def cycles(M):
 n=len(M);adj=[set(j for j in range(n)if j!=i and M[i][j]!=Z)for i in range(n)]
 out={tuple((i,j))for i,j in combinations(range(n),2)if j in adj[i]}
 for s in range(n):
  def visit(path):
   for j in adj[path[-1]]:
    if j==s and len(path)>=3 and path[1]<path[-1]:out.add(tuple(path))
    elif j>s and j not in path:visit(path+[j])
  visit([s])
 result={}
 for c in sorted(out):
  v=ONE
  for i,j in zip(c,c[1:]+c[:1]):v=R.mul(v,R.div(R.mul(TWO,M[i][j]),M[i][i]))
  result[canonical_cycle(c)]=(v,integral(v))
 return result

def manifest(paths):
 seen=set()
 def check(p):
  if p in seen:return
  seen.add(p)
  if p.endswith('.json'):
   for name,h in read(p).get('sources',{}).items():assert sha(name)==h,(p,name);check(name)
 for p in paths:check(p)
 return {p:sha(p)for p in paths}

def run():
 primary=read('results/'+BASE+'.json');old=read('results/fresh_sqrt5_quad_contacts_20261015.json');unions=read('results/fresh_sqrt5_quad_exhaustion_20261015.json')['labelled_union_supports']
 source=[r for r in old['records']if r['status']=='UNPAID_DOTTED_CONTACT_DOMAIN'];assert len(source)==len(primary['records'])==24
 counts=Counter();records=[];cycle_count=0
 for origin,rec in zip(source,primary['records']):
  row=unions[origin['union_index']];assert row==rec['row']and origin['class_index']==rec['source_class_index']and origin['union_index']==rec['union_index']
  n=row['nodes'];assert [p for p,m in edges(row).items()if m==0]==[(0,2)]
  r=finite_row(row);M=R.matrix(r,R.coloring(r));I=[i for i in range(n)if i not in(0,2)];B=[[M[i][j]for j in I]for i in I];J=R.inverse(B)
  assert R.inertia(B)==[n-3,1,0]
  u=[M[0][i]for i in I];v=[M[2][i]for i in I];Ju=[dot(z,u)for z in J];Jv=[dot(z,v)for z in J]
  a=R.sub(M[0][0],dot(u,Ju));b=R.sub(M[2][2],dot(v,Jv));c=dot(u,Jv);assert R.sign(a)>=0 and R.sign(b)>=0
  bound=R.div(R.mul(q(8),R.add(R.mul(a,b),R.mul(c,c))),R.mul(M[0][0],M[2][2]));assert bound==dec(rec['upper_bound_expression'])
  U=rec['integer_upper_bound'];assert R.sign(R.sub(q(U),bound))>=0 and R.sign(R.sub(q(U-1),bound))<0
  # Scan a deliberately larger box and both signs of B independently.
  candidates={}
  for A in range(-U,U+1):
   for Bb in range(-U,U+1):
    t=R.add(q(A),R.mul(q(Bb),PHI));tc=R.conj(t)
    if R.sign(R.sub(t,FOUR))>0 and R.sign(R.sub(q(U),t))>=0 and R.sign(tc)>=0 and R.sign(R.sub(FOUR,tc))>=0:candidates[A,Bb]=t
  saved={(z['A'],z['B']):z for z in rec['trials']};assert len(saved)==len(rec['trials'])and set(saved)==set(candidates)
  accepted=[];local=Counter();review_trials=[]
  for key,t in sorted(candidates.items()):
   s=saved[key];assert t==dec(s['twice_gram_cycle']);G,error=matrix(row,R.div(t,FOUR));detail={}
   if error:status=error
   else:
    p=R.inertia(G);co=R.inertia([[R.conj(x)for x in z]for z in G]);assert p==s['physical']and co==s['conjugate'];detail=dict(physical=p,conjugate=co)
    if p[1]>1:status='REJECT_TWO_NEGATIVE'
    elif co[1]:status='REJECT_CONJUGATE'
    else:
     cc=cycles(G);cycle_count+=len(cc);bad=next((c for c,(v,ok)in cc.items()if not ok),None)
     if bad:status='REJECT_CYCLIC_INTEGRALITY'
     else:
      assert p[1]==1 and p[2]in(0,1)and co==[n-p[2],0,p[2]]
      status='ACCEPT_SINGULAR'if p[2]else'ACCEPT_NONSINGULAR'
      archived={canonical_cycle(z['cycle']):(dec(z['value']),z['integral'])for z in s['cycles']};assert cc==archived
      accepted.append(key);detail['cycle_count']=len(cc)
   assert status==s['status'],(origin['class_index'],key,status,s['status']);counts[status]+=1;local[status]+=1
   review_trials.append(dict(A=key[0],B=key[1],status=status,**detail))
  archived_accepted=[(rec['trials'][z['trial_index']]['A'],rec['trials'][z['trial_index']]['B'])for z in rec['accepted']];assert accepted==archived_accepted
  records.append(dict(source_class_index=origin['class_index'],union_index=origin['union_index'],upper_bound=R.enc(bound),integer_upper_bound=U,counts=dict(local),trials=review_trials))
 assert dict(counts)==primary['counts']and sum(counts.values())==532 and len([z for r in records for z in r['trials']if z['status'].startswith('ACCEPT')])==32
 paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'results/{BASE}.md','scripts/fresh_sqrt5_quad_filters_20261015_review.py','results/fresh_sqrt5_quad_filters_20261015_review.json','results/fresh_sqrt5_quad_contacts_20261015.json','results/fresh_sqrt5_quad_contacts_20261015_review.json','results/fresh_sqrt5_quad_exhaustion_20261015.json']
 return dict(status='PASS_INDEPENDENT_COMPLETE_DOTTED_PARAMETER_BOX',sources=manifest(paths),cores=24,integral_box_candidates=532,counts=dict(counts),accepted_cycle_products=cycle_count,records=records,method='Separate tuple-Fraction Qsqrt5 arithmetic; independent minimal radical gauge, full principal inverse/inertia, symmetric integer box and elementary simple-cycle enumeration. No producer imports.',scope='Complete necessary integral parameter box on the24 frozen unpaid minimal dotted cores, yielding32 principal seeds. The Lorentz Schur bound holds in every larger extension. Arbitrary attachment coverage, geometric realization, payment and dimension improvements are separate.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();path=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==read(str(path.relative_to(ROOT)))
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'])
