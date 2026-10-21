#!/usr/bin/env python3
"""Exact deficient cusp created by the radius-four replacement."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from fresh_h20_intrinsic_caps_20261014_golden_seed import ROOT,ZERO,ONE,q,add,sub,mul,scale,sign,enc,sha,pd
TAG='fresh_h20_intrinsic_caps_20261014_parallel_cusp_gap'
BASE='results/fresh_h20_intrinsic_caps_20261014_'
def rank(A):
 A=[r[:]for r in A];row=0
 for j in range(len(A[0])):
  ix=next((i for i in range(row,len(A))if A[i][j]),None)
  if ix is None:continue
  A[row],A[ix]=A[ix],A[row];pivot=A[row][j];A[row]=[v/pivot for v in A[row]]
  for i in range(row+1,len(A)):
   a=A[i][j]
   if a:A[i]=[v-a*w for v,w in zip(A[i],A[row])]
  row+=1
  if row==len(A):break
 return row
def run():
 paths=[BASE+s+'.json'for s in ['golden_seed','golden_extension','quarter_face_exact','radius4_replacement']]
 seed,ext,quarter,rep=[json.loads((ROOT/p).read_text())for p in paths]
 for data in(seed,ext,quarter,rep):
  for p,h in data['sources'].items():assert sha(p)==h,p
 dec=lambda x:tuple(map(F,x));B=[[F(x)for x in row]for row in seed['spatial_Gram']];pd(B)
 walls=[[dec(x)for x in r]for r in seed['walls']+ext['new_walls']+quarter['new_walls']]
 assert len(walls)==587 and rep['replacement_index']==23;walls[23]=list(map(dec,rep['new_normal']))
 def ld(a,b):
  out=ZERO
  for i,row in enumerate(B):
   for j,x in enumerate(row):
    if x:out=add(out,scale(mul(a[i],b[j]),x))
  return sub(sub(out,mul(a[19],b[20])),mul(a[20],b[19]))
 assert seed['source_wall_order'][15]==22
 assert ld(walls[15],walls[15])==ONE and ld(walls[23],walls[23])==q(16)and ld(walls[15],walls[23])==q(-4)
 delta=[add(a,scale(b,4))for a,b in zip(walls[23],walls[15])]
 center=list(map(dec,seed['common_strict_center']['point']))
 assert ld(center,center)==q(-9)and all(sign(ld(center,w))<0 for w in walls)
 assert ld(delta,delta)==ZERO and ld(delta,center)==q(F(-709,6),-2)
 pairs=[ld(delta,w)for w in walls];assert all(sign(v)<=0 for v in pairs)
 incident=[i for i,p in enumerate(pairs)if p==ZERO]
 assert incident==[2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,19,23]
 gram=[[ld(walls[i],walls[j])for j in incident]for i in incident]
 assert all(p[1]==0 for row in gram for p in row)
 rational=[[p[0]for p in row]for row in gram];assert rank(rational)==16
 assert len(pairs)-len(incident)==570
 paths+=['scripts/fresh_h20_intrinsic_caps_20261014_golden_seed.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return enc({'status':'PASS_EXACT_RADIUS4_PARALLEL_CUSP_HAS_THREE_DIMENSIONAL_LINK_LINEALITY','sources':{p:sha(p)for p in paths},'field':seed['field'],'Lorentz_signature':[20,1],'ray':delta,'positive_null_relation':{'replacement_index':23,'coefficient':1,'old_seed_index':15,'old_source_wall':22,'old_seed_coefficient':4},'ray_norm':ZERO,'ray_center_pairing':ld(delta,center),'all587_pairings':pairs,'incident_indices':incident,'strict_nonincident_count':570,'incident_Gram':rational,'incident_Gram_rank':16,'horospherical_dimension':19,'lineality_dimension':3,'strict_interior_point':center,'proof':'The ambient Lorentz signature is(20,1); delta is nonzero future null and feasible. Its perpendicular quotient is Euclidean19. All incident normals project there and their rational Gram has rank16, so the orthogonal complement has dimension3 and is contained in both signs of the horospherical recession cone. The common strict timelike point gives link interior. All570 other inequalities are strict at delta and disappear sufficiently deep in the cusp. Hence the actual cusp link is unbounded.','scope':'This is one explicit cusp of the frozen radius4 replacement, not a complete cusp census. Rank16 is the Euclidean quotient rank, not the ambient rank. No finite-volume claim.'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
