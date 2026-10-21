#!/usr/bin/env python3
"""Shared single-own-root high-label family, exact affine Schur gate."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product,combinations
from sympy import symbols,expand
from fresh_h9_star_zero_20261014_own_three import K,z,o,h,co,sign,enc,det3
from independent_h9_14_record30_a1_labels import finite_code
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_shared_high_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';GATES='results/fresh_h9_star_zero_20261014_gates.json'
S=[3,4,6,7,9,10,12,13];P=[2,5,8,11];C=[14,15,16,17]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 paths=[GEOM,GATES,'scripts/fresh_h9_star_zero_20261014_own_three.py','scripts/independent_h9_14_record30_a1_labels.py'];seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p not in seen:
   seen.add(p)
   for q,hh in d.get('sources',{}).items():
    assert sha(q)==hh,q
    if q.endswith('.json'):read(q)
  return d
 geo,gates=[read(p)for p in paths[:2]];g=geo['representatives'][0]['geometry'];V=[v['incident']for v in g['vertices']];sets=list(map(set,V));missing={frozenset(m)for m in g['minimal_nonfaces']}
 def face(s):
  support=sorted(s);vi=next(i for i,v in enumerate(sets)if set(support)<=v)
  return dict(support=support,actual_vertex=vi)
 own={tuple(sorted((leaf,a))):3 for i in range(3)for leaf in(P[i],C[i])for a in S[2*i:2*i+2]}
 own|={(11,12):2,(11,13):3,(12,17):2,(13,17):3}
 assert all(frozenset(S[2*i:2*i+2]+[leaf])in missing for i in range(4)for leaf in(P[i],C[i]))
 assert F(1,2)+F(1,3)+F(1,6)==1
 swap=list(range(18));swap[12],swap[13]=13,12
 assert {tuple(sorted(swap[a]for a in v))for v in V}==set(map(tuple,V))
 assert {frozenset(swap[a]for a in v)for v in missing}==missing
 high=[]
 for i in range(3):
  lab=own|{tuple(S[2*j:2*j+2]):(6 if i==j else 7 if j==3 else 4)for j in range(4)}
  hits=[f['id']for f in gates['families']if f['origin']=='two_dense'and not f['fixed']and all(lab.get((a,b),0)>=m for a,b,m in f['bounds'])]
  assert hits;high.append(dict(high_block=i,reviewed_gate=hits[0]))
 cross=[]
 for i,j in combinations(range(4),2):
  for a in S[2*i:2*i+2]:
   for b in S[2*j:2*j+2]:cross.append(dict(pair=sorted((a,b)),**face(S[2*i:2*i+2]+S[2*j:2*j+2]),relaxed_path=[4,3,4],determinant='0'))
 assert 1-F(1,2)-F(1,4)-F(1,2)+F(1,4)==0
 budgets=[]
 for leaf,drop in[(11,12),(17,12),(14,3)]:budgets.append(dict(leaf=leaf,**face((set(S)-{drop})|{leaf})))
 assert F(1,4)+2*2*F(1,4)>1 and F(1,4)+2*F(1,2)>1
 high_third=face([12,13,14]);mutual=face([11,13,17]);assert 1-F(3,4)-F(1,4)==0
 def failure(lab,required=None):
  for vi,v in enumerate(V):
   if required is not None and not set(required)<=sets[vi]:continue
   code=tuple(lab.get(tuple(sorted(p)),2)for p in combinations(v,2))
   if finite_code(9,code):continue
   B=[[o if a==b else-co[lab.get(tuple(sorted((a,b))),2)]for b in v]for a in v];piv=[]
   for k in range(9):
    pivot=B[k][k];piv.append(enc(pivot))
    if sign(pivot)<=0:return dict(actual_vertex=vi,prefix_order=k+1,pivots=piv)
    for a in range(k+1,9):
     for b in range(k+1,9):B[a][b]-=B[a][k]*B[k][b]/pivot
   raise AssertionError(('unverified classifier witness',vi))
  return None
 base=own|{tuple(S[2*i:2*i+2]):4 for i in range(4)}
 upper=[]
 for pair in((11,14),(14,17)):
  f=failure(base|{pair:4},pair);assert f
  upper.append(dict(pair=pair,test_label=4,relaxed_internal_labels=[4,4,4,4],**f))
 # A rank-one D contribution makes the residual determinant affine.
 a,b,d,e,f,q,D=symbols('a b d e f q D')
 symbolic=det3(a-D/4,b-D/4,e,d-D/4,f,q)
 intercept=det3(a,b,e,d,f,q);slope=((e-f)**2-(a+d-2*b)*q)/4
 assert expand(symbolic-intercept-D*slope)==0
 ar=[]
 for extra in [None]+list(range(6)):
  v=[z]*6
  if extra is not None:v[extra]=-h
  ar.append((extra,v))
 tr=[]
 for extra in[None,2,3,4,5]:
  v=[-h,-h]+[z]*4
  if extra is not None:v[extra]=-h
  tr.append((extra,v))
 rows=[];pos=0;bad=0;minA=None;minBound=None
 for blocks in product((4,5),repeat=3):
  inv=[o/(o-co[m]**2)for m in blocks]
  def pdot(x,y):return sum((inv[i]*(x[2*i]*y[2*i]+x[2*i+1]*y[2*i+1]+co[blocks[i]]*(x[2*i]*y[2*i+1]+x[2*i+1]*y[2*i]))for i in range(3)),z)
  for (ia,x),(ib,y),(ic,t)in product(ar,ar,tr):
   aa=o-pdot(x,x);bb=o-pdot(y,y);ab=-pdot(x,y);tt=o-pdot(t,t)
   for ma,mb in product((2,3),repeat=2):
    at=-co[ma]-pdot(x,t);bt=-co[mb]-pdot(y,t)
    B=det3(aa,ab,at,bb,bt,tt);A=((at-bt)**2-(aa+bb-2*ab)*tt)/4;bound=B+4*A
    assert sign(A)>0
    assert det3(aa-o,ab-o,at,bb-o,bt,tt)==bound
    if minA is None or sign(A-minA)<0:minA=A
    row=dict(block_labels=blocks,foreign_contacts=[ia,ib,ic],third_labels=[ma,mb],slope=enc(A),intercept=enc(B),boundary_D4=enc(bound))
    if sign(bound)>0:
     row['status']='STRICT_POSITIVE_FOR_ALL_D_GREATER_THAN_FOUR';pos+=1
     if minBound is None or sign(bound-minBound)<0:minBound=bound
    else:
     lab=base|{tuple(S[2*i:2*i+2]):blocks[i]for i in range(3)}|{(11,14):ma,(14,17):mb}
     for leaf,extra in zip((11,17,14),(ia,ib,ic)):
      if extra is not None:lab[tuple(sorted((leaf,S[extra])))]=3
     f=failure(lab);assert f,(blocks,ia,ib,ic,ma,mb)
     row['status']='NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR';row['vertex_failure']=f;bad+=1
    rows.append(row)
  print('blocks',blocks,'positive',pos,'vertex_rejected',bad,flush=True)
 assert len(rows)==7840 and pos==7392 and bad==448
 paths+=[f'scripts/{TAG}.py']
 return dict(status='PASS_SHARED_SINGLE_OWN_ROOT_ALL_HIGH_LABELS_EXCLUDED',sources={p:sha(p)for p in paths},case=0,fixed_own_contacts=sorted([list(p)+[m]for p,m in own.items()]),internal_lower_bounds=[S[2*i:2*i+2]+[4 if i<3 else 7]for i in range(4)],literal_twin_swap=swap,first_three_high_reductions=high,cross_right_guards=cross,projection_budgets=budgets,high_block_third_right_guard=high_third,mutual_right_guard=mutual,third_upper_guards=upper,affine_identity='det([[a-D/4,b-D/4,e],[b-D/4,d-D/4,f],[e,f,q]])=det([[a,b,e],[b,d,f],[e,f,q]])+D*((e-f)^2-(a+d-2b)*q)/4',tail_argument='m>=7 implies cos(pi/m)>cos(pi/6)=sqrt3/2, hence D=1/sin^2(pi/m)>4; no upper label bound or coefficient field assumption on D.',cases=rows,positive_cases=pos,actual_vertex_rejections=bad,minimum_slope=enc(minA),minimum_positive_boundary=enc(minBound),field='Coefficient field Q(sqrt2,sqrt5), descending ANP at sqrt2+sqrt5; D arbitrary real>4',scope='Star0 first3blocks all-own3 and arbitrary internal labels. Last two leaves have the same single own pi3 contact and a right contact to the other last root. Its missing triples force every finite internal label>=7. Both root orientations checked; omitted contacts/distances unrestricted. Not a whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['positive_cases'],d['actual_vertex_rejections'],flush=True)
