import sys,json,numpy as np
sys.path[:0]=['scripts','src']
from fresh_h20_gap_next_20261013 import setup,read
from scipy.optimize import linprog
Q,R,h,w=setup();d1=read('data/enumerations/fresh_h20_gap_next_20261013_coupled.json');d=read('data/enumerations/fresh_h20_gap_next_20261013_asymmetric.json');U=np.vstack([R,np.array(d1['new_normal'])]);hn=h/np.sqrt(-h@Q@h);out=[]
for i,c in enumerate(d['caps']):
 V=np.vstack([U,np.array(c['normal'])]);ns=np.einsum('ij,jk,ik->i',V,Q,V);A=V@Q/np.sqrt(ns[:,None]);rng=np.random.default_rng(20272003);e=None
 for k in range(100):
  v=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
  if v.success and v.x@Q@v.x>1e-7:
   e={'attempt':k,'vector':v.x.tolist(),'norm':float(v.x@Q@v.x),'active_indices':np.flatnonzero(abs(A@v.x)<1e-7).tolist()};break
 out.append({'candidate':i,'labels':c['labels'],'escape':e});print(i,c['labels'],None if e is None else(k,e['norm']),flush=True)
json.dump({'status':'BOUNDED_22_ALTERNATIVE_98_WALL_VOLUME_SCOUT','source':'data/enumerations/fresh_h20_gap_next_20261013_asymmetric.json','attempt_limit_per_candidate':100,'candidates':out},open('data/enumerations/fresh_h20_gap_next_20261013_asymmetric_volume.json','w'),indent=2)
