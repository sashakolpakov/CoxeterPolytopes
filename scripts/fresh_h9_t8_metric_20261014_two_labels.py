#!/usr/bin/env python3
"""Two-label continuous extension: m>=7 and the last block label n>=4."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from itertools import product
import fresh_h9_star_cube_metric_20261014 as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_t8_metric_20261014'
BASE=f'results/{TAG}.json'
Z,O,H,A=K.Z,K.O,K.H,K.A

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def parse(v):return tuple(F(t)for t in v)
def run():
 base=json.loads((ROOT/BASE).read_text())
 for p,h in base['sources'].items():assert sha(p)==h,p
 assert base['status']=='PASS_T8_SAVED_ASSIGNMENT_ALL_M_GE_7_RANK_ELEVEN'
 k0=K.add(O,A)
 def fdot(a,b):
  if a is None or b is None:return Z
  if a==b:return H
  if a//2==b//2:return K.scale(A,F(1,2))
  return Z
 def third_dot(a):return K.add(H,K.scale(A,F(1,2)))if a in[2,3]else Z
 rows=[];anchors={}
 for row in base['residual_cases']:
  ia,ib=row['foreign_contacts'];ma,mb=row['contact_labels']
  u=Z if ia is None else H;v=Z if ib is None else H;p=fdot(ia,ib)
  a1=K.sub(K.add(u,v),K.add(K.scale(O,2),K.scale(p,2)));assert K.sign(a1)<0
  ac=K.sub(Z if ma==2 else K.neg(H),third_dot(ia))
  bc=K.sub(Z if mb==2 else K.neg(H),third_dot(ib))
  square=K.mul(K.sub(ac,bc),K.sub(ac,bc));assert K.sign(square)>=0
  def D(k,q):
   qa=K.sub(K.sub(O,k),u);qb=K.sub(K.sub(O,k),v);ab=K.neg(K.add(k,p))
   return K.determinant3([[qa,ab,ac],[ab,qb,bc],[ac,bc,q]])
  # Algebraic coefficient checks: D is affine in k, and its coefficient
  # equals a1*q+(ac-bc)^2 for every q.
  for q in[Z,O]:
   d0,d1,d2=[D(k,q)for k in[Z,O,K.scale(O,2)]]
   assert K.sub(K.sub(d2,K.scale(d1,2)),K.neg(d0))==Z
   assert K.sub(d1,d0)==K.add(K.mul(a1,q),square)
  qbnd=K.neg(K.add(A,O))
  boundary=D(k0,qbnd);assert boundary==parse(row['boundary_at_three_quarters']) and K.sign(boundary)>0
  slope_boundary=K.add(K.mul(a1,qbnd),square);assert K.sign(slope_boundary)>0
  qa=K.sub(K.sub(O,k0),u);qb=K.sub(K.sub(O,k0),v);ab=K.neg(K.add(k0,p))
  ad=K.sub(K.mul(qa,qb),K.mul(ab,ab));assert K.sign(ad)<0
  anchors[ia,ib]=dict(foreign_contacts=[ia,ib],anchor_determinant_at_k0=ad,negative_k_slope=a1)
  rows.append(dict(foreign_contacts=[ia,ib],contact_labels=[ma,mb],negative_anchor_slope=a1,
   nonnegative_square=square,positive_k_slope_at_boundary=slope_boundary,positive_minimum=boundary))
 assert len(rows)==100 and len(anchors)==25
 paths=[BASE,f'scripts/{TAG}.py',f'scripts/{TAG}_two_labels.py','scripts/fresh_h9_star_cube_metric_20261014.py']
 return K.enc(dict(status='PASS_T8_TWO_LABEL_FAMILY_M_GE_7_N_GE_4_EXCLUDED',sources={p:sha(p)for p in paths},
  parameter_pairs={'m':[1,2],'n':[12,13]},k0=k0,
  parameters='t=cos(pi/m)^2 in(3/4,1); k=1/(2(1-cos(pi/n)))>=1+sqrt2/2',
  positive_support_determinant='(1-t)*(1-s^2)/4>0, s=cos(pi/n)<1',
  anchor_own_projection_norm='((1/4)+(1/4)+s/2)/(1-s^2)=1/(2(1-s))=k',
  inherited_guards='Only the internal12--13 label increases from4 to n>=4. All nonpositive actual-vertex Gram witnesses persist by Z-matrix monotonicity. The own-root deletion budgets remove this parameter edge; third14 foreign-contact costs use the inverse diagonal>=2 of I2(n). Thus anchors still have five patterns and third14 still has its fixed projection.',
  determinant_identity='D_k=a1*q14+(ac-bc)^2>0; a1<0 and q14=-sqrt2/2-1/(4(1-t))<0. Anchor determinant decreases in k, remainsnegative, so D_t>0 too.',
  anchors=list(anchors.values()),cases=rows,
  scope='All other mandatory labels retain their source T8 values; omitted contacts remain free subject to the proved necessary bounds. The two internal finite labels vary over m>=7,n>=4. Every one of the100 residual cases is positive on this whole real parameter region. This excludes this two-label metric family only, not the incidence or other mandatory assignments.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_two_labels.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['cases']),flush=True)
