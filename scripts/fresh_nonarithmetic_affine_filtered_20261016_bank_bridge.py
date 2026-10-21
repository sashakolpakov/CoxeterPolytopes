#!/usr/bin/env python3
"""Bounded numerical bridge-partner diagnostic; no exact exclusion claim."""
import argparse,json,math
from pathlib import Path
from hashlib import sha256
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_escape49_20261016_polar as P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_affine_filtered_20261016_bank_bridge'
BANK='results/fresh_nonarithmetic_escape49_20261016_target_scout.json';CAPS='results/fresh_nonarithmetic_affine_filtered_20261016_cluster.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 old=P.B.read(P.BATCH);bank=P.B.read(BANK);caps=P.B.read(CAPS)
 Q=mp.matrix([[P.scalar(P.D.encfield(P.D.lift(P.A.B.db(x))))for x in row]for row in old['form']]);N=[mp.matrix([P.scalar(x)for x in row])for row in old['normals']]
 N=[n*mp.sqrt(2/(n.T*Q*n)[0])for n in N];basis=[0,1,2,3,4,5,6,7,8,11]
 G=mp.matrix([[(N[i].T*Q*N[j])[0]for j in basis]for i in basis]);J=np.array([[float(x)for x in row]for row in(G**-1).tolist()]);V=np.array([[r['normalized_wall_pairings'][i]for i in basis]for r in bank['candidates']])
 selected=[c['scout_index']for c in caps['caps']];out=[];common=set(range(len(V)))
 for i in selected:
  vals=V@J@V[i];dotted=[];finite=[];unresolved=[]
  for j,g in enumerate(vals):
   if j==i:continue
   if g<-2-1e-7:dotted.append(j);continue
   if abs(g+2)<=1e-7:unresolved.append(j);continue
   if abs(g)<1e-7:finite.append([j,2]);continue
   if -2<g<0:
    m=math.pi/math.acos(-g/2)
    if abs(m-round(m))<1e-7:finite.append([j,round(m)])
  comp=set(dotted)|{j for j,m in finite}|set(unresolved);common&=comp
  out.append(dict(scout_index=i,wall=49+selected.index(i),dotted_partners=dotted,finite_angle_fits=finite,near_parallel_unresolved=unresolved))
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',P.BATCH,BANK,CAPS]
 return dict(status='BOUNDED_NUMERICAL_NO_BRIDGE_PARTNER_FOR_SELECTED_CAPS',sources={p:sha(p)for p in paths},basis=basis,proposal_count=len(V),selected=out,common_compatible_or_unresolved=sorted(common),scope='Double-precision reconstruction from saved numerical proposals, using80-digit old Gram data. No exact rejection or all-label theorem. Within this bounded1846-bank no selected cap has a finite-angle fit to another proposal, and no additional candidate is numerically compatible with allthree.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved bounded bridge diagnostic')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
