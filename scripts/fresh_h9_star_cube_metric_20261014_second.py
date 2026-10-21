#!/usr/bin/env python3
"""Exact metric obstruction for the second saved star mandatory state."""
import argparse,json
from itertools import product
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from fresh_h9_star_cube_metric_20261014 import Z,O,H,A,add,neg,sub,mul,scale,sign,enc,pdot,determinant3
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_cube_metric_20261014'
FIRST=f'results/{TAG}.json';OLD='results/fresh_h9_tree_cube_wedges_20261014_case_00.json';STATE=f'results/{TAG}_search_case_00.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 first,old,state=[json.loads((ROOT/p).read_text())for p in[FIRST,OLD,STATE]]
 for d in[first,old,state]:
  for p,h in d['sources'].items():assert sha(p)==h,p
 assert first['status']=='PASS_FIRST_STAR_MANDATORY_ASSIGNMENT_RANK_AT_LEAST_ELEVEN'
 before={(i,j):m for i,j,m in old['solution']['pair_labels']};labels={(i,j):m for i,j,m in state['solution']['pair_labels']}
 assert set(before)==set(labels)and[(p,before[p],labels[p])for p in labels if before[p]!=labels[p]]==[((13,17),3,4)]
 assert all(labels[p]>=m for p,m in before.items())
 # All30 first-state contact gates persist under this monotone increase.
 S=first['fixed_positive_support'];COS={2:Z,3:H,4:A}
 def pattern(v,group,extra):
  x=[Z]*8
  for j in range(2):x[2*group+j]=neg(COS[labels[tuple(sorted((v,S[2*group+j])))]] )
  if extra is not None:x[extra]=neg(H)
  return tuple(x)
 aa=[pattern(11,3,e)for e in[None,0,1,2,3,4,5]]
 bb=[pattern(17,3,None)]
 cc=[pattern(14,0,e)for e in[None,2,3,4,5,6,7]]
 assert len(aa)==7 and len(bb)==1 and len(cc)==7
 assert sub(O,pdot(bb[0],bb[0]))==(F(-3,2),F(0))
 # The saved C3+7 support retains own root13, now a label4 contact.
 budget=next(w for w in first['foreign_contact_budget_witnesses']if w['leaf']==17)
 assert 13 in budget['support']and 12 not in budget['support']
 assert F(1,2)+2*F(1,4)==1
 # Strict spherical positivity prohibits even one foreign label3.
 rows=[];anchors=[];counts={-1:0,1:0}
 for ia,x in enumerate(aa):
  y=bb[0];qx=sub(O,pdot(x,x));qy=sub(O,pdot(y,y));xy=neg(pdot(x,y));ad=sub(mul(qx,qy),mul(xy,xy));assert sign(ad)<0
  anchors.append(dict(pattern=ia,determinant=ad))
  for ic,z in enumerate(cc):
   qz=sub(O,pdot(z,z))
   for ma,mb in product([2,3],repeat=2):
    xz=sub(neg(COS[ma]),pdot(x,z));yz=sub(neg(COS[mb]),pdot(y,z))
    det=determinant3([[qx,xy,xz],[xy,qy,yz],[xz,yz,qz]]);assert sign(det)!=0;counts[sign(det)]+=1
    rows.append(dict(anchor_pattern=ia,third_pattern=ic,contact_labels=[ma,mb],residual_determinant=det))
 assert len(rows)==196
 paths=[FIRST,OLD,STATE,f'scripts/{TAG}.py',f'scripts/{TAG}_second.py']
 return enc(dict(status='PASS_SECOND_STAR_MANDATORY_ASSIGNMENT_RANK_AT_LEAST_ELEVEN',sources={p:sha(p)for p in paths},
  changed_contact=[13,17,3,4],inherited_monotone_contact_witnesses=30,
  stronger_C3_foreign_contact_budget=budget,projection_patterns=dict(P3=aa,C3=bb,C0=cc),
  anchors=anchors,rank_obstructions=rows,determinant_sign_counts=counts,
  scope='Exactly the second saved assignment. Every first-state forced-right and finite-contact upper bound remains valid by nonnegative-matrix monotonicity. The increased own contact forces zero foreign C3 contacts by a strict spherical Schur budget. All7 independent anchor patterns and196 third-vector cases give rank11; no central or dotted-distance choices remain relevant to this assignment.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_second.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['rank_obstructions']),d['determinant_sign_counts'])
