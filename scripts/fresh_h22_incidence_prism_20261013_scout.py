#!/usr/bin/env python3
"""Numerical finite/remote-vertex and polar-truncation scout for new812 walls."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fresh_h22_borcherds_extension_20261013 import prepare
P=Path(__file__).resolve().parents[1]
PRE='fresh_h22_incidence_prism_20261013'
def run():
 Q,W,G,h,X,A,b,Qn,_=prepare()
 data=json.loads((P/'data/enumerations/fresh_h22_incidence_general_binary_20261013_bank.json').read_text())['candidates']
 chosen=[c for c in data if c['twice_height_square']==4];C=np.array([c['twice_old_projection']for c in chosen],dtype=float)/2
 normals=np.vstack([np.column_stack([W,np.zeros(len(W))]),np.column_stack([C,np.full(len(C),np.sqrt(2))]),np.column_stack([C,np.full(len(C),-np.sqrt(2))])])
 metric=np.zeros((23,23));metric[:22,:22]=Q;metric[-1,-1]=1
 transform=np.zeros((23,22));transform[:22,:21]=X;transform[-1,-1]=1
 center=np.r_[h/np.sqrt(14),0];limits=-(normals@metric@center);H=(normals@metric@transform)/limits[:,None]
 norms=np.einsum('ij,ij->i',normals,normals@metric);cosines=np.cos(np.pi/np.arange(2,121))
 rng=np.random.default_rng(22602);runs=[];seen=set();polars=[]
 for trial in range(300):
  direction=rng.normal(size=22);lp=linprog(direction,A_ub=H,b_ub=np.ones(len(H)),bounds=[(None,None)]*22,method='highs');assert lp.success
  key=tuple(np.flatnonzero(abs(H@lp.x-1)<1e-7));
  if key in seen:continue
  seen.add(key);q=float(lp.x@lp.x-1);x=center+transform@lp.x
  record={'trial':trial,'norm':q,'point':x.tolist(),'active':key}
  record['active']=list(map(int,key))
  if q>1e-7:
   c=-(normals@metric@x)/np.sqrt(q*norms);finite=(c<1-1e-7);invalid=(c<-1e-7)|(finite&(np.min(abs(c[:,None]-cosines),axis=1)>1e-6));record['polar_invalid_contacts']=list(map(int,np.flatnonzero(invalid)));record['polar_invalid_count']=int(invalid.sum())
   if not np.any(invalid):polars.append(record)
  runs.append(record)
  if trial%25==0:print('trial',trial,'vertices',len(runs),'spacelike',sum(r['norm']>1e-7 for r in runs),'max',max(r['norm']for r in runs),'polars',len(polars),flush=True)
 out={'status':'NUMERICAL_NEW812_VERTEX_POLAR_SCOUT','facets':812,'runs':runs,'polar_candidates':polars,'scope':'LP proposals only; no exact or exhaustive vertex census.'}
 (P/('data/enumerations/'+PRE+'_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE vertices',len(runs),'finite',sum(r['norm']<-1e-7 for r in runs),'spacelike',sum(r['norm']>1e-7 for r in runs),'null',sum(abs(r['norm'])<1e-7 for r in runs),'polars',len(polars),'badpolarcounts',sorted(set(r.get('polar_invalid_count',-1)for r in runs)),flush=True)
if __name__=='__main__':run()
