#!/usr/bin/env python3
"""Bounded numerical volume scout for the changed three-anchor family."""
import json
import numpy as np
from scipy.optimize import linprog
from fresh_h20_four_cusp_completion_20261013 import state,OUT as SOURCE,ROOT,read,vec,B
from fresh_h20_coupled_completion_20261013 import num
BANK='data/enumerations/fresh_h20_four_cusp_completion_20261013_gap_scout.json'
OUT='data/enumerations/fresh_h20_four_cusp_completion_20261013_gap_volume_scout.json'
Q,h,hc,roots=state();roots.extend((B,vec(c['normal']))for c in read(SOURCE)['new_caps'])
Q=np.array([[float(x.a)for x in row]for row in Q]);h=np.array([float(x.a)for x in h]);R=np.array([[num(f,x)for x in r]for f,r in roots]);caps=np.array([c['normal']for c in read(BANK)['new_candidates']]);hn=h/np.sqrt(-h@Q@h)
def escape(U,objectives=120):
 ns=np.einsum('ij,jk,ik->i',U,Q,U);A=U@Q/np.sqrt(ns[:,None]);rng=np.random.default_rng(202610139)
 for k in range(objectives):
  ans=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if ans.success and ans.x@Q@ans.x>1e-7:return {'vector':ans.x.tolist(),'norm':float(ans.x@Q@ans.x),'active_indices':np.flatnonzero(abs(A@ans.x)<1e-7).tolist(),'attempt':k}
 return None
out={'status':'NUMERICAL_ONLY','source':SOURCE,'bank':BANK,'bank_count':len(caps),'tests':[]}
for ids in([0],list(range(len(caps)))):
 e=escape(np.vstack([R,caps[ids]]));out['tests'].append({'cap_ids':ids,'escape':e});print('selection',len(ids),'escape',e,flush=True)
(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
