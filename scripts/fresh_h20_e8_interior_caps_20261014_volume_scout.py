#!/usr/bin/env python3
"""Bounded numerical global witness search for the chosen coupled pair."""
import json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_interior_caps_20261014_volume_scout'
BASE='results/fresh_h20_e8_cusp_completion_20261014.json';CAPS='results/fresh_h20_e8_interior_caps_20261014_scout.json'
def run():
 d,c=[json.loads((ROOT/p).read_text())for p in [BASE,CAPS]];Q=np.array(d['form'],float);rad=np.sqrt([1,2,79,158])
 def nums(row):return np.array([[float(F(v))for v in x]for x in row])@rad
 N=np.vstack([np.array([nums(r)for r in d['normals']]),np.array([c['candidates'][i]['normal']for i in[12,15]])]);H=nums(d['strict_center']);h=H/np.sqrt(-H@Q@H)
 B=null_space((h@Q)[None,:]);ev,U=np.linalg.eigh(B.T@Q@B);assert ev.min()>0;P=B@U@np.diag(1/np.sqrt(ev))
 assert np.max(abs(P.T@Q@P-np.eye(20)))<1e-7 and np.max(abs(h@Q@P))<1e-7
 A=N@Q@P;b=-N@Q@h;assert b.min()>0
 rng=np.random.default_rng(1979);directions=np.vstack([np.eye(20),-np.eye(20),rng.normal(size=(32,20))]);escapes=[];seen=set()
 for direction in directions:
  for step in range(4):
   lp=linprog(-direction,A_ub=A,b_ub=b,bounds=[(None,None)]*20,method='highs')
   assert lp.success,lp.message
   z=lp.x;direction=z;key=tuple(np.flatnonzero(abs(A@z-b)<1e-6))
   if key not in seen:
    seen.add(key);zz=.999*z;x=h+P@zz;norm=x@Q@x
    if norm>1e-7:escapes.append(dict(point=x.tolist(),norm=float(norm),active=list(map(int,key)),max_pair=float(np.max(N@Q@x))))
 escapes.sort(key=lambda x:-x['norm'])
 return dict(status='NUMERICAL_81_WALL_GLOBAL_ESCAPE_SCOUT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in [BASE,CAPS,f'scripts/{TAG}.py']},selected_candidates=[12,15],LP_calls=288,distinct_proposals=len(seen),escapes=escapes,scope='Bounded numerical proposal search only; no finite-volume or whole-ray census. Exact checker must reconstruct and verify the chosen strict witness.')
if __name__=='__main__':
 d=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(d,indent=2)+'\n');print(d['status'],len(d['escapes']),d['escapes'][0]['norm']if d['escapes']else None)
