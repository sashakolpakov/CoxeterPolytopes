#!/usr/bin/env python3
"""Numerical whole-bank and coherent-selection volume scouts."""
import json,hashlib,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup

ROOT=Path(__file__).resolve().parents[1]
BANK='data/enumerations/fresh_h20_targeted_escape_20261013_vertex_bank.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_bank_geometry.json'

def run():
 bank=json.loads((ROOT/BANK).read_text());Q,R,h,_,_,_,_=extended_setup()
 C=np.array([c['normal']for c in bank['caps']]);n=len(C);cov=R@Q;ns=np.einsum('ij,jk,ik->i',R,Q,R)
 center=h/np.sqrt(-h@Q@h);equ=(center@Q)[None,:]
 cutoff=np.cos(np.pi/np.arange(2,121));dots=C@Q@C.T/2
 ratios=-dots;ix=np.clip(np.searchsorted(cutoff,ratios),1,len(cutoff)-1)
 errors=np.minimum(abs(ratios-cutoff[ix]),abs(ratios-cutoff[ix-1]));compat=(ratios>=1-1e-7)|(errors<1e-7)
 np.fill_diagonal(compat,False);assert np.all(compat==compat.T)
 rng=np.random.default_rng(204687);objectives=rng.normal(size=(128,21));checks=[]
 def volume(ids,limit):
  walls=np.vstack([R,C[ids]]);norms=np.r_[ns,np.ones(len(ids))*2];A=walls@Q/np.sqrt(norms[:,None])
  best=None
  for j,obj in enumerate(objectives[:limit]):
   lp=linprog(obj,A_ub=A,b_ub=np.zeros(len(A)),A_eq=equ,b_eq=[-1],bounds=[(None,None)]*21,method='highs')
   if not lp.success:return {'status':'LP_UNRESOLVED','lp_status':int(lp.status),'objective':j}
   q=lp.x@Q@lp.x
   if best is None or q>best['norm']:best={'norm':float(q),'point':lp.x.tolist(),'objective':j,'active_indices':np.flatnonzero(abs(A@lp.x)<1e-7).tolist()}
   if q>1e-7:
    best['status']='NUMERICAL_SPACELIKE_ESCAPE';return best
  best['status']='NO_ESCAPE_IN_BOUNDED_OBJECTIVES';return best
 full=volume(list(range(n)),128);print({'full_bank':full['status'],'norm':full.get('norm')},flush=True)
 seeds=list(bank['escape_cutters'])+[None]*18
 for iteration,seed in enumerate(seeds):
  allowed=np.ones(n,dtype=bool);selected=[]
  if seed is not None:selected=[seed];allowed&=compat[seed];allowed[seed]=False
  while np.any(allowed):
   ids=np.flatnonzero(allowed)
   degrees=np.sum(compat[np.ix_(ids,ids)],axis=1)
   scores=degrees+(rng.random(len(ids))*2 if iteration>=6 else 0)
   j=int(ids[np.argmax(scores)]);selected.append(j);allowed&=compat[j];allowed[j]=False
  assert all(compat[i,j]for k,i in enumerate(selected)for j in selected[:k])
  geom=volume(selected,64);checks.append({'seed':seed,'cap_indices':selected,'new_caps':len(selected),'walls':len(R)+len(selected),'geometry':geom})
  print({'selection':iteration,'caps':len(selected),'geometry':geom['status'],'norm':geom.get('norm')},flush=True)
 out={'status':'NUMERICAL_H20_SHORT_VERTEX_BANK_GEOMETRY','source':BANK,
  'source_sha256':hashlib.sha256((ROOT/BANK).read_bytes()).hexdigest(),
  'bank_entries':n,'approximate_compatible_pairs':int(np.sum(compat)//2),
  'whole_bank_geometry':full,'coherent_selections':checks,
  'scope':'Floating fullbank relaxation and24greedy maximal compatible subsets usinglabels2..120 or dotted. LP objectives are bounded. Spacelike proposals need exact certification; no finite-volume conclusion from sampling.'}
 (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':run()
