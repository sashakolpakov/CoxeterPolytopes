#!/usr/bin/env python3
"""Continuous rank obstruction for T-case8's saved m>=7 assignment."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from itertools import product,combinations
import fresh_h9_star_cube_metric_20261014 as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_metric_20261014_two_high'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json';STATE='results/fresh_h9_t8_metric_20261014_family_search_case_08.json';TWO='results/fresh_h9_t8_metric_20261014_two_labels.json'
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
  (11,12):2,(11,13):3,(12,13):7,(12,17):2,(13,17):3}
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
    cross.append(dict(pair=[a,b],support=support,actual_vertex=actual(support),gram=M,null_vector=w,relaxation='Internal label at least4 lowered to4; positive null vector remains nonpositive when restored'))
 assert len(cross)==24
 own=[11,13,17];assert all(labels.get(tuple(sorted(p)),2)==3 for p in[(11,13),(13,17)])
 ownM=[[O if i==j else K.neg(H)for j in range(3)]for i in range(3)]
 assert all(K.total(row)==Z for row in ownM)
 own_guard=dict(pair=[11,17],support=own,actual_vertex=actual(own),test_label=3,gram=ownM,null_vector=[1,1,1])
 budgets=[]
 for leaf in[11,17]:
  for kept in[13]:
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
  assert not {1,2}<=set(v) and not {12,13}<=set(v),'parameter edges must be absent from this finite witness'
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
 two=json.loads((ROOT/TWO).read_text())
 for p,h in two['sources'].items():assert sha(p)==h,p
 assert two['status']=='PASS_T8_TWO_LABEL_FAMILY_M_GE_7_N_GE_4_EXCLUDED'
 parse=lambda x:tuple(F(v)for v in x)
 k0=K.add(O,A);kmin=(F(5,4),F(0));rows=[];anchors=[]
 for rec in two['anchors']:
  slope=parse(rec['negative_k_slope']);old=parse(rec['anchor_determinant_at_k0'])
  value=K.sub(old,K.mul(K.sub(k0,kmin),slope));assert K.sign(value)<0 and K.sign(slope)<0
  anchors.append(dict(foreign_contacts=rec['foreign_contacts'],determinant_at_five_quarters=value,negative_k_slope=slope))
 vals=[]
 for rec in two['cases']:
  slope=parse(rec['positive_k_slope_at_boundary']);old=parse(rec['positive_minimum'])
  value=K.sub(old,K.mul(K.sub(k0,kmin),slope));assert K.sign(value)>0 and K.sign(slope)>0
  vals.append(value)
  rows.append(dict(foreign_contacts=rec['foreign_contacts'],contact_labels=rec['contact_labels'],
   negative_anchor_slope=parse(rec['negative_anchor_slope']),nonnegative_square=parse(rec['nonnegative_square']),
   positive_k_slope_at_boundary=slope,positive_minimum=value))
 minimum=(F(1,8),F(1,8));assert minimum in vals and all(K.sign(K.sub(v,minimum))>=0 for v in vals)
 # U6(cos(pi/7))=0, by sin(7theta)/sin(theta). Generate U6 exactly.
 def padd(a,b):return[(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))]
 u0=[1];u1=[0,2]
 for _ in range(2,7):u0,u1=u1,padd([0]+[2*v for v in u1],[-v for v in u0])
 assert u1==[-1,0,24,0,-80,0,64]
 p=lambda t:64*t**3-80*t*t+24*t-1
 assert p(F(4,5))==F(-29,125)
 assert 192*F(3,4)**2-160*F(3,4)+24==12
 assert 384*F(3,4)-160==128
 # Thus p is strictly increasing above3/4; its root cos²(pi/7)>3/4
 # lies above4/5. For n>=7, k=1/(4(1-cos²(pi/n)))>5/4.
 paths=[GEOM,STATE,TWO,f'scripts/{TAG}.py','scripts/fresh_h9_star_cube_metric_20261014.py']
 return K.enc(dict(status='PASS_T8_TWO_HIGH_LABELS_M_N_GE_7_EXCLUDED',sources={p:sha(p)for p in paths},
  positive_support=S,positive_support_actual_vertex=svertex,positive_determinant='(1-t)*(1-u)/4>0',
  parameters='t=cos(pi/m)^2>3/4; u=cos(pi/n)^2>4/5; k=1/(4(1-u))>5/4; m,n finite integers>=7',
  forced_cross_right_witnesses=cross,anchor_pair_right_witness=own_guard,anchor_budget_witnesses=budgets,
  third_leaf_budget_witnesses=third_budgets,anchor_contact_upper_witnesses=upper,
  chebyshev_second_kind_degree6=u1,cosine_square_cubic=[-1,24,-80,64],cubic_at_four_fifths=F(-29,125),
  derivative_at_three_quarters=12,second_derivative_at_three_quarters=128,
  minimum_boundary=minimum,anchors=anchors,cases=rows,
  scope='Exactly the new saved mandatory assignment, with both category7 labels covering all finite integers>=7. Anchors11/17 have own contacts(2,3), so their common projection cost is k=1/(4(1-u)), not the earlier two-nonright cost. Changed literal actual-vertex guards are recomputed. All100 residual determinants remain positive above the new real boundary t=3/4,k=5/4. Other mandatory assignments and whole T incidence remain open.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['cases']),flush=True)
