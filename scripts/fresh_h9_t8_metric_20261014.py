#!/usr/bin/env python3
"""Continuous rank obstruction for T-case8's saved m>=7 assignment."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product,combinations
import fresh_h9_star_cube_metric_20261014 as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_metric_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';STATE='results/fresh_h9_tree_cube_wedges_20261014_case_08.json'
Z,O,H,A=K.Z,K.O,K.H,K.A
COS={2:Z,3:H,4:A}
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 geo,state=[json.loads((ROOT/p).read_text())for p in[GEOM,STATE]]
 for d in[geo,state]:
  for p,h in d['sources'].items():assert sha(p)==h,p
 assert state['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR'
 g=geo['representatives'][8]['geometry'];V=[v['incident']for v in g['vertices']]
 labels={(a,b):m for a,b,m in state['solution']['pair_labels']}
 expected={(0,1):2,(0,2):3,(1,2):7,(1,14):2,(1,15):3,(1,16):3,(2,14):3,(2,15):2,(2,16):2,
  (3,4):3,(3,5):3,(4,5):4,(4,14):3,(5,14):3,(6,7):3,(6,8):3,(7,8):4,(7,15):3,(8,15):3,
  (11,12):3,(11,13):3,(12,13):4,(12,17):3,(13,17):3}
 assert labels==expected
 S=[1,2,4,5,7,8,12,13];blocks=[S[2*i:2*i+2]for i in range(4)]
 def actual(s):return next(i for i,v in enumerate(V)if set(s)<=set(v))
 svertex=actual(S)
 cross=[]
 for i,j in combinations(range(4),2):
  for a,b in product(blocks[i],blocks[j]):
   if i==0:
    support=[next(v for v in blocks[0]if v!=a),a,b]
    cross.append(dict(pair=[a,b],support=support,actual_vertex=actual(support),
     witness='three-path determinant 3/4-t<0, t=cos(pi/m)^2>3/4'))
   else:
    support=[next(v for v in blocks[i]if v!=a),a,b,next(v for v in blocks[j]if v!=b)]
    # Set the explicit 4-3-4 chain, independent of all leaf contacts.
    M=[[O if x==y else Z for y in range(4)]for x in range(4)]
    for x,c in enumerate([A,H,A]):M[x][x+1]=M[x+1][x]=K.neg(c)
    w=[O,(F(0),F(1)),(F(0),F(1)),O]
    assert all(K.total(K.mul(M[x][y],w[y])for y in range(4))==Z for x in range(4))
    cross.append(dict(pair=[a,b],support=support,actual_vertex=actual(support),gram=M,null_vector=w))
 assert len(cross)==24
 own=[11,12,17];assert all(labels.get(tuple(sorted(p)),2)==3 for p in[(11,12),(12,17)])
 ownM=[[O if i==j else K.neg(H)for j in range(3)]for i in range(3)]
 assert all(K.total(row)==Z for row in ownM)
 own_guard=dict(pair=[11,17],support=own,actual_vertex=actual(own),test_label=3,gram=ownM,null_vector=[1,1,1])
 budgets=[]
 for leaf in[11,17]:
  for kept in[12,13]:
   support=sorted((set(S)-({12,13}-{kept}))|{leaf})
   budgets.append(dict(leaf=leaf,retained_own_root=kept,support=support,actual_vertex=actual(support)))
 third_budgets=[]
 for kept in[4,5]:
  support=sorted((set(S)-{1}-({4,5}-{kept}))|{14})
  third_budgets.append(dict(leaf=14,retained_own_roots=[2,kept],support=support,actual_vertex=actual(support)))
 assert F(1,4)+2*F(1,4)*2>1 and F(1,4)+2*F(1,2)>1
 assert F(1,4)+F(1,4)+2*F(1,4)==1
 # t>3/4 implies the I2 inverse diagonal 1/(1-t)>4. A minimal
 # nonright anchor/I2 contact costs >1, already impossible.
 def gram(v,lab):
  assert not {1,2}<=set(v),'parameter edge must be absent from this finite witness'
  return[[O if a==b else K.neg(COS[lab.get(tuple(sorted((a,b))),2)])for b in v]for a in v]
 upper=[]
 for anchor in[11,17]:
  pair=(anchor,14)if anchor<14 else(14,anchor)
  for vi,v in enumerate(V):
   if not set(pair)<=set(v):continue
   pp=K.pivots(gram(v,labels|{pair:4}))
   if K.sign(pp[-1])<=0:
    upper.append(dict(pair=pair,test_label=4,actual_vertex=vi,leading_order=len(pp),last_pivot=pp[-1]));break
  else:raise AssertionError(('no upper bound',pair))
 def dot(x,y):
  # The two anchors have zero I2 coordinates; only three B2 blocks
  # contribute to any cross product used below.
  return K.total(K.add(K.scale(K.add(K.mul(x[k],y[k]),K.mul(x[k+1],y[k+1])),2),
   K.mul((F(0),F(1)),K.add(K.mul(x[k],y[k+1]),K.mul(x[k+1],y[k]))))for k in[2,4,6])
 def pattern(extra):
  p=[Z]*8;p[6]=p[7]=K.neg(H)
  if extra is not None:p[extra]=K.neg(H)
  return p
 third=[Z]*8;third[1]=third[2]=third[3]=K.neg(H)
 rows=[];anchors=[];boundary_values=[]
 for ia,ib in product([None,2,3,4,5],repeat=2):
  x,y=pattern(ia),pattern(ib)
  qa=K.sub(O,dot(x,x));qb=K.sub(O,dot(y,y));ab=K.neg(dot(x,y))
  ad=K.sub(K.mul(qa,qb),K.mul(ab,ab));assert K.sign(ad)<0
  anchors.append(dict(foreign_contacts=[ia,ib],residual_gram=[[qa,ab],[ab,qb]],determinant=ad))
  for ma,mb in product([2,3],repeat=2):
   ac=K.sub(K.neg(COS[ma]),dot(x,third));bc=K.sub(K.neg(COS[mb]),dot(y,third))
   D0=K.determinant3([[qa,ab,ac],[ab,qb,bc],[ac,bc,K.neg(A)]])
   boundary=K.sub(D0,ad);assert K.sign(boundary)>0
   slope=K.scale(K.neg(ad),F(1,4));assert K.sign(slope)>0
   # D(t)=D0+slope/(1-t). It is strictly increasing for 3/4<t<1.
   assert K.add(D0,K.scale(slope,4))==boundary
   boundary_values.append(boundary)
   rows.append(dict(foreign_contacts=[ia,ib],contact_labels=[ma,mb],constant=D0,
    positive_inverse_denominator_coefficient=slope,boundary_at_three_quarters=boundary))
 minimum=(F(3,4),F(5,8));assert minimum in boundary_values
 assert all(K.sign(K.sub(v,minimum))>=0 for v in boundary_values)
 assert len(rows)==100 and len(anchors)==25
 assert minimum==(F(3,4),F(5,8))
 paths=[GEOM,STATE,f'scripts/{TAG}.py','scripts/fresh_h9_star_cube_metric_20261014.py']
 return K.enc(dict(status='PASS_T8_SAVED_ASSIGNMENT_ALL_M_GE_7_RANK_ELEVEN',sources={p:sha(p)for p in paths},
  positive_support=S,positive_support_actual_vertex=svertex,positive_determinant='(1-t)/8>0',
  parameter='t=cos(pi/m)^2 lies strictly between3/4 and1 for every finite integer m>=7',
  forced_cross_right_witnesses=cross,anchor_pair_right_witness=own_guard,anchor_budget_witnesses=budgets,
  third_leaf_budget_witnesses=third_budgets,anchor_contact_upper_witnesses=upper,
  third_projection='pairings at S[1],S[2],S[3] are -1/2; all other S pairings zero',
  third_residual_norm='-sqrt2/2 -1/(4(1-t))',anchors=anchors,residual_cases=rows,minimum_boundary=minimum,
  scope='The literal saved T-case8 mandatory assignment, with category7 meaning every finite integer m>=7. The PD8 block leaves residual dimension at most2. All25 anchor patterns and100 third-contact cases have a strictly positive residual determinant throughout3/4<t<1, hence rank>=11. Literal actual-vertex and Schur-budget witnesses cover all omitted finite contacts used; arbitrary dotted distances remain unconstrained. Other assignments and the whole T incidence remain open.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['residual_cases']),flush=True)
