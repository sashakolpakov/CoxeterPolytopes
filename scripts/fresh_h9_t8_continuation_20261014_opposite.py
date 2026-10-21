#!/usr/bin/env python3
"""Exact continuous rank gate for opposite-root T8 anchors.

Own nonright contacts may have arbitrary finite labels >=3.  The finite
400-case table treats their minimal value; explicit derivative inequalities
cover every larger value.  No dotted distance or tail label is substituted.
"""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product,combinations
from sympy import QQ,sqrt,symbols,expand,Matrix,diff
import fresh_h9_star_cube_metric_20261014 as G
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_continuation_20261014_opposite'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
STATE='results/fresh_h9_t8_continuation_20261014_case_08.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 geo,state=[json.loads((ROOT/p).read_text())for p in[GEOM,STATE]]
 for obj in[geo,state]:
  for p,h in obj['sources'].items():assert sha(p)==h,p
 assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
 V=[v['incident']for v in geo['representatives'][8]['geometry']['vertices']]
 lab={(a,b):m for a,b,m in state['solution']['pair_labels']}
 expected={(0,1):2,(0,2):3,(1,2):7,(1,14):2,(1,15):3,(1,16):3,(2,14):3,(2,15):2,(2,16):2,
  (3,4):3,(3,5):3,(4,5):4,(4,14):3,(5,14):3,(6,7):3,(6,8):3,(7,8):4,(7,15):3,(8,15):3,
  (11,12):2,(11,13):3,(12,13):7,(12,17):3,(13,17):2}
 assert lab==expected
 S=[1,2,4,5,7,8,12,13];blocks=[S[2*i:2*i+2]for i in range(4)]
 def actual(s):return next(i for i,v in enumerate(V)if set(s)<=set(v))
 cross=[]
 for i,j in combinations(range(4),2):
  for a,b in product(blocks[i],blocks[j]):
   if i==0:
    support=[next(v for v in blocks[0]if v!=a),a,b]
    cross.append(dict(pair=[a,b],support=support,actual_vertex=actual(support),
     witness='3-path: determinant 3/4-t<0 for t=cos(pi/m)^2>3/4'))
   else:
    support=[next(v for v in blocks[i]if v!=a),a,b,next(v for v in blocks[j]if v!=b)]
    M=[[G.O if x==y else G.Z for y in range(4)]for x in range(4)]
    for x,c0 in enumerate([G.A,G.H,G.A]):M[x][x+1]=M[x+1][x]=G.neg(c0)
    w=[G.O,(F(0),F(1)),(F(0),F(1)),G.O]
    assert all(G.total(G.mul(M[x][y],w[y])for y in range(4))==G.Z for x in range(4))
    cross.append(dict(pair=[a,b],support=support,actual_vertex=actual(support),
     relaxed_gram=M,positive_null_vector=w,relaxation='Both internal labels lowered to4'))
 assert len(cross)==24
 anchor_budgets=[]
 for leaf,kept,removed in[(11,13,12),(17,12,13)]:
  support=sorted((set(S)-{removed})|{leaf})
  anchor_budgets.append(dict(leaf=leaf,retained_own_root=kept,support=support,actual_vertex=actual(support),
   at_most_one_foreign='owncost>=1/4; each B2 nonright contact costs>=1/2; I2(m) contact costs>1',
   foreign_label_bound='a label>=4 costs>=1, plus positive owncost',
   own_label_with_foreign='ownlabel>=4 costs>=1/2, plus foreigncost>=1/2; thus ownlabel=3'))
 third_budgets=[]
 for kept in[4,5]:
  support=sorted((set(S)-{1}-({4,5}-{kept}))|{14})
  third_budgets.append(dict(leaf=14,retained_own_roots=[2,kept],support=support,actual_vertex=actual(support)))
 assert F(1,4)+2*F(1,2)>1 and F(1,2)+F(1,2)==1
 assert F(1,4)+F(1,4)+F(1,2)==1
 # The anchor pair is not forced right.  Its first forbidden label is6.
 path=[11,12,17];assert lab[11,12]==2 and lab[12,17]==3
 actual_path=actual(path);assert 1-F(3,4)-F(1,4)==0
 # Literal upper witnesses use minimal own contacts3; all unassigned
 # contacts are right. Increasing any specified contact preserves failure.
 co2={2:G.Z,3:G.H,4:G.A};upper=[]
 for anchor in[11,17]:
  pair=tuple(sorted((anchor,14)))
  for vi,v in enumerate(V):
   if not set(pair)<=set(v):continue
   assert not {1,2}<=set(v) and not {12,13}<=set(v)
   labels=lab|{pair:4}
   M=[[G.O if a==b else G.neg(co2[labels.get(tuple(sorted((a,b))),2)])for b in v]for a in v]
   pp=G.pivots(M)
   if G.sign(pp[-1])<=0:
    upper.append(dict(pair=pair,test_label=4,actual_vertex=vi,vertex=v,gram=M,
     leading_order=len(pp),last_pivot=pp[-1]));break
  else:raise AssertionError(('no fresh upper witness',pair))
 # cos^2(pi/7)>81/100: U6 and an increasing cubic above3/4.
 u0=[1];u1=[0,2]
 def ipadd(a,b):return[(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))]
 for _ in range(2,7):u0,u1=u1,ipadd([0]+[2*v for v in u1],[-v for v in u0])
 assert u1==[-1,0,24,0,-80,0,64]
 p=lambda t:64*t**3-80*t*t+24*t-1
 assert p(F(81,100))==F(-559,15625)
 assert 192*F(3,4)**2-160*F(3,4)+24==12 and 384*F(3,4)-160==128
 K=QQ.algebraic_field(sqrt(2),sqrt(5));z,o=K.zero,K.one
 k=lambda q:K.convert(QQ(q.numerator,q.denominator))if isinstance(q,F)else K.convert(q)
 h=k(F(1,2));r2=K.from_sympy(sqrt(2));r5=K.from_sympy(sqrt(5))
 def enc(x):return[str(v)for v in x.to_list()]
 def sign2(a,b):
  sg=lambda v:int(v>0)-int(v<0)
  if not b:return sg(a)
  if not a:return sg(b)
  if sg(a)==sg(b):return sg(a)
  return sg(a)*sg(a*a-2*b*b)
 def sign(x):
  pp=[F(str(v))for v in x.to_list()];pp=[F(0)]*(4-len(pp))+pp
  a3,a2,a1,a0=pp;a=a0+7*a2;b=a1+17*a3;c=a1+11*a3;d=2*a2
  su,sv=sign2(a,b),sign2(c,d)
  if not sv:return su
  if not su:return sv
  if su==sv:return su
  return su*sign2(a*a+2*b*b-5*c*c-10*d*d,2*a*b-10*c*d)
 assert K.mod.to_list()==[1,0,-14,0,9] and sign(r2)>0 and sign(r5)>0
 def pa(a,b):return[(a[i]if i<len(a)else z)+(b[i]if i<len(b)else z)for i in range(max(len(a),len(b)))]
 def ps(a,b):return pa(a,[-v for v in b])
 def pm(a,b):
  out=[z]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):out[i+j]+=x*y
  return out
 def sc(a,b):return[v*b for v in a]
 def trim(a):
  a=a[:]
  while len(a)>1 and a[-1]==z:a.pop()
  return a
 def eq(a,b):assert trim(a)==trim(b)
 def bern(a):
  l=k(F(9,10));v=a[0]+a[1]*l+a[2]*l*l
  return[v,v+k(F(1,20))*(a[1]+2*a[2]*l),sum(a,z)]
 Q=r2/2+k(F(25,19));q=-Q;assert sign(Q-2)>0
 cos={2:z,3:h,4:r2/2,5:(o+r5)/4};B=[o,z,-o];rows=[]
 patterns=[None,4,5,7,8]
 for ia,ib,mab,ma,mb in product(patterns,patterns,[2,3,4,5],[2,3],[2,3]):
  u=z if ia is None else h;v=z if ib is None else h
  pairproj=z if ia is None or ib is None or ((ia in[4,5])!=(ib in[4,5]))else h if ia==ib else r2/4
  d=cos[mab]+pairproj
  ac=-cos[ma]-(h+r2/4 if ia in[4,5]else z)
  bc=-cos[mb]-(h+r2/4 if ib in[4,5]else z)
  C=(o-u)*(o-v)-d*d
  f=q*C-(o-u)*bc*bc-(o-v)*ac*ac-2*d*ac*bc
  N=[f+(-q*(k(F(7,4))-u-v)+bc*bc+ac*ac)/4,(-2*q*d-2*ac*bc)/4,-f]
  A=[C-(k(F(7,4))-u-v)/4,-d/2,-C]
  qa=ps(sc(B,o-u),[k(F(1,4))]);qb=ps(sc(B,o-v),[k(F(1,4))]);ab=ps(sc(B,-d),[z,k(F(1,4))])
  acp,bcp,qp=sc(B,ac),sc(B,bc),sc(B,q)
  anchor=ps(pm(qa,qb),pm(ab,ab));eq(anchor,pm(A,B))
  det=pa(ps(pm(qa,ps(pm(qb,qp),pm(bcp,bcp))),pm(ab,ps(pm(ab,qp),pm(acp,bcp)))),pm(acp,ps(pm(ab,bcp),pm(qb,acp))))
  eq(det,pm(N,pm(B,B)))
  nb,abern=bern(N),bern(A)
  assert all(sign(x)>0 for x in nb) and all(sign(x)<0 for x in abern)
  rows.append(dict(foreign_contacts=[ia,ib],anchor_mutual_label=mab,third_contact_labels=[ma,mb],
   anchor_numerator=list(map(enc,A)),residual_numerator=list(map(enc,N)),
   anchor_bernstein=list(map(enc,abern)),residual_bernstein=list(map(enc,nb))))
 assert len(rows)==400
 # Symbolic identities for the entire own-contact extension, separate from
 # the fixed-contact polynomial calculation above.
 aa,bb,c,u,v,d,ac,bc,Qs=symbols('alpha beta c u v d ac bc Q',real=True)
 Bs=1-c*c
 An=((1-u)*(1-v)-d*d)*Bs-aa*aa*(1-v)-bb*bb*(1-u)-2*d*c*aa*bb+aa*aa*bb*bb
 N=-Qs*An-Bs*((1-u)*bc*bc+(1-v)*ac*ac+2*d*ac*bc)+aa*aa*bc*bc+bb*bb*ac*ac-2*c*aa*bb*ac*bc
 M=Matrix([[1-u-aa*aa/Bs,-d-c*aa*bb/Bs,ac],[-d-c*aa*bb/Bs,1-v-bb*bb/Bs,bc],[ac,bc,-Qs]])
 from sympy import cancel
 assert cancel(M.det()*Bs-N)==0
 assert expand(diff(N,aa)/2-(aa*(Qs*(1-v-bb*bb)+bc*bc)+c*bb*(Qs*d-ac*bc)))==0
 assert expand(diff(An,aa)+2*aa*(1-v-bb*bb)+2*d*c*bb)==0
 special=Qs*(c*c-(1-aa*aa)*(1-bb*bb))+(aa*bc-bb*ac)**2+2*aa*bb*(1-c)*ac*bc-Bs*(ac*ac+bc*bc)
 assert expand(N.subs({u:0,v:0,d:0})-special)==0
 assert F(2,8)-F(1,32)==F(7,32)>0
 assert 2*(F(81,100)-F(9,16))-F(19,100)/2==F(2,5)>0
 # In the one-free/one-foreign case z=|bc| and epsilon=|ac|<=1/2:
 # N_alpha/2 >= Q/8+z^2/2-z/4 = Q/8-1/32+(z-1/4)^2/2 >0.
 # If both free and delta>=1/2, Q*delta-ac*bc>1-1/4>0.
 paths=[GEOM,STATE,f'scripts/{TAG}.py','scripts/fresh_h9_star_cube_metric_20261014.py']
 return G.enc(dict(status='PASS_T8_OPPOSITE_ROOT_ALL_OWN_LABELS_FAMILY_EXCLUDED',sources={p:sha(p)for p in paths},
  field='Q(sqrt2,sqrt5), primitive element sqrt2+sqrt5; ANP descending coefficients',
  source_mandatory_labels=state['solution']['pair_labels'],free_mandatory_pairs=[[1,2],[12,13],[11,13],[12,17]],
  parameter_domains='m,n finite integers>=7; alpha=cos(pi/a),beta=cos(pi/b), a,b finite integers>=3; t=cos(pi/m)^2>81/100, c=cos(pi/n)>9/10',
  positive_support=S,positive_support_actual_vertex=actual(S),positive_support_determinant='(1-t)*(1-c^2)/4>0',
  cross_right_guards=cross,anchor_budgets=anchor_budgets,third_budgets=third_budgets,
  anchor_mutual_upper_guard=dict(support=path,actual_vertex=actual_path,test_label=6,determinant='1-3/4-1/4=0',allowed_labels=[2,3,4,5]),
  third_contact_upper_guards=upper,chebyshev_U6=u1,threshold_cubic_value=F(-559,15625),
  Q_at_t_lower=enc(Q),Q_minus_two=enc(Q-2),polynomial_interval=[F(9,10),F(1)],cases=rows,
  all_own_extension=dict(symbolic_determinant_identity=True,symbolic_derivative_identities=True,
   anchor_derivative='A_alpha=-2alpha(1-v-beta^2)-2delta*c*beta<0 for every free alpha',
   residual_derivative='N_alpha/2=alpha[Q(1-v-beta^2)+bc^2]+c*beta[Q*delta-ac*bc]',
   one_free_one_foreign_lower=F(7,32),both_free_positive_mutual_bracket_lower=F(3,4),
   both_free_right_mutual_lower=F(2,5),foreign_forces_own_label=3,
   argument='Only foreign-free anchors can have alpha>1/2. Reduce any such variable by positive derivatives except when both anchors have no foreign contact and mutual label2; the explicit 2/5 lower bound handles that case. Anchor minors decrease under every free increase. The residual determinant increases with t because derivative=-A/[4B(1-t)^2]>0, where A is the anchor numerator and B=1-c^2.'),
  scope='Compact H9 Coxeter realizations of the saved opposite-root mandatory family, allowing both I2 labels to be arbitrary finite >=7 and both own nonright contacts arbitrary finite >=3. Other saved mandatory labels remain fixed; omitted meeting contacts and dotted distances remain unrestricted subject to necessary actual-vertex gates. The whole incidence is not excluded. The PD8 Schur complement would have rank at most2, but its selected 3x3 determinant is strictly positive.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['cases']),flush=True)
