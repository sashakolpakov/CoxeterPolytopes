#!/usr/bin/env python3
"""Opposite single-own-root pi3 family, exact continuous quadratic gate."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product,combinations
from sympy import symbols,expand,chebyshevt
from fresh_h9_star_zero_20261014_own_three import K,z,o,h,co,sign,enc,det3
from independent_h9_14_record30_a1_labels import finite_code
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_opposite_high_20261014'
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
 own|={(11,12):2,(11,13):3,(12,17):3,(13,17):2}
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
 for leaf,drop in[(11,12),(17,13),(14,3)]:budgets.append(dict(leaf=leaf,**face((set(S)-{drop})|{leaf})))
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
 for pair,test in[((11,14),4),((14,17),5)]:
  f=failure(base|{pair:test},pair);assert f
  upper.append(dict(pair=pair,test_label=test,relaxed_internal_labels=[4,4,4,4],**f))
 # Clearing 1-c² makes the residual determinant quadratic in c.
 a,b,d,e,f,q,D,c=symbols('a b d e f q D c')
 symbolic=det3(a-D/4,b-c*D/4,e,d-D/4,f,q)
 B0=det3(a,b,e,d,f,q);A0=-(a+d)*q/4+(e*e+f*f)/4+q/16;C0=(b*q-e*f)/2
 assert expand((symbolic-B0-D*(A0+c*C0)).subs(D,1/(1-c*c))).cancel()==0
 x=symbols('x');p7=8*x**3-4*x**2-4*x+1
 assert expand(chebyshevt(7,x)+1-(x+1)*p7**2)==0
 assert p7.subs(x,F(9,10))==-F(1,125)
 assert (24*x*x-8*x-4).subs(x,F(3,4))==F(7,2)
 assert (48*x-8).subs(x,F(3,4))==28
 assert F(3,4)>F(36,49) and F(3,4)>F(9,16)
 def bern(B,A,C,lo):
  l=K.convert(lo);b0=B*(o-l*l)+A+l*C
  return [b0,b0+(o-l)*(-2*B*l+C)/2,A+C]
 def bern_ok(bs):return all(sign(v)>=0 for v in bs)and any(v!=z for v in bs)
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
 rows=[];pos=0;bad=0;refined=0;zero_coefficients=0
 for blocks in product((4,5),repeat=3):
  inv=[o/(o-co[m]**2)for m in blocks]
  def pdot(x,y):return sum((inv[i]*(x[2*i]*y[2*i]+x[2*i+1]*y[2*i+1]+co[blocks[i]]*(x[2*i]*y[2*i+1]+x[2*i+1]*y[2*i]))for i in range(3)),z)
  for (ia,x),(ib,y),(ic,t)in product(ar,ar,tr):
   aa=o-pdot(x,x);bb=o-pdot(y,y);tt=o-pdot(t,t)
   for mutual_label,ma,mb in product((2,3,4,5),(2,3),(2,3,4)):
    ab=-co[mutual_label]-pdot(x,y);at=-co[ma]-pdot(x,t);bt=-co[mb]-pdot(y,t)
    B=det3(aa,ab,at,bb,bt,tt);A=-(aa+bb)*tt/4+(at*at+bt*bt)/4+tt/16;Cq=(ab*tt-at*bt)/2
    bs=bern(B,A,Cq,F(6,7))
    row=dict(block_labels=blocks,foreign_contacts=[ia,ib,ic],mutual_label=mutual_label,third_labels=[ma,mb],B=enc(B),A=enc(A),C=enc(Cq),bernstein_6_7=[enc(v)for v in bs])
    if bern_ok(bs):
     row['status']='STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE';pos+=1;zero_coefficients+=sum(v==z for v in bs)
    else:
     lab=base|{tuple(S[2*i:2*i+2]):blocks[i]for i in range(3)}|{(11,17):mutual_label,(11,14):ma,(14,17):mb}
     for leaf,extra in zip((11,17,14),(ia,ib,ic)):
      if extra is not None:lab[tuple(sorted((leaf,S[extra])))]=3
     ff=failure(lab)
     if ff is not None:
      row['status']='NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR';row['vertex_failure']=ff;bad+=1
     else:
      br=bern(B,A,Cq,F(9,10));assert all(sign(v)>0 for v in br)
      row['status']='STRICT_BERNSTEIN_9_10_TO_ONE';row['bernstein_9_10']=[enc(v)for v in br];refined+=1
    rows.append(row)
  print('blocks',blocks,'positive',pos,'vertex_rejected',bad,'refined',refined,flush=True)
 assert len(rows)==47040 and pos==43412 and bad==3612 and refined==16
 paths+=[f'scripts/{TAG}.py']
 return dict(status='PASS_OPPOSITE_SINGLE_PI3_ALL_HIGH_LABELS_EXCLUDED',sources={p:sha(p)for p in paths},case=0,fixed_own_contacts=sorted([list(p)+[m]for p,m in own.items()]),internal_lower_bounds=[S[2*i:2*i+2]+[4 if i<3 else 7]for i in range(4)],literal_twin_swap=swap,first_three_high_reductions=high,cross_right_guards=cross,projection_budgets=budgets,high_block_third_right_guard=high_third,mutual_label_six_guard=mutual,third_upper_guards=upper,quadratic_identity='(1-c^2)detR=B*(1-c^2)+A+c*C, R=[[a-D/4,b-cD/4,e],[b-cD/4,d-D/4,f],[e,f,q]], D=1/(1-c^2); A=-(a+d)q/4+(e^2+f^2)/4+q/16; C=(bq-ef)/2',tail_bounds=dict(coarse='cos(pi/m)>sqrt3/2>6/7 for m>=7',refined='cos(pi/7)>9/10 from T7(x)+1=(x+1)*(8x^3-4x^2-4x+1)^2 and strict cubic monotonicity on[3/4,1]',cubic_at_9_10='-1/125',derivative_at_3_4='7/2',second_derivative_at_3_4='28'),cases=rows,coarse_bernstein_cases=pos,actual_vertex_rejections=bad,refined_bernstein_cases=refined,zero_bernstein_coefficients=zero_coefficients,field='Coefficient field Q(sqrt2,sqrt5), descending ANP at sqrt2+sqrt5; c arbitrary real in the proved intervals',scope='Star0 first3blocks all-own3 and arbitrary internal labels. Last two leaves have opposite single own pi3 contacts. Its missing triples force every finite internal label>=7. Both root orientations checked; omitted contacts/distances unrestricted. Not a whole-incidence exclusion.')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['coarse_bernstein_cases'],d['actual_vertex_rejections'],d['refined_bernstein_cases'],flush=True)
