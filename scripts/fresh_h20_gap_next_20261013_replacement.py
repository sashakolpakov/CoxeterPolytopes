import sys,json,numpy as np
sys.path[:0]=['scripts','src']
from fresh_h20_gap_next_20261013 import setup,read,TARGET
from scipy.optimize import linprog
Q,R,h,w=setup();ids=read(TARGET)['certificate']['numeric_active_indices'];_,s,vh=np.linalg.svd(R[ids]@Q,full_matrices=True);p=vh[-1];p*=np.sign(-h@Q@p);p*=np.sqrt(2/(p@Q@p))
keep=[i for i in range(96)if i not in(3,42)];U=np.vstack([R[keep],p]);ns=np.einsum('ij,jk,ik->i',U,Q,U);A=U@Q/np.sqrt(ns[:,None]);hn=h/np.sqrt(-h@Q@h);rng=np.random.default_rng(20272003)
for k in range(200):
 v=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
 if v.success and v.x@Q@v.x>1e-7:
  print(k,v.x@Q@v.x,'olddeleted',R[[3,42]]@Q@v.x,'active',np.flatnonzero(abs(A@v.x)<1e-7).tolist(),flush=True)
  json.dump({'status':'NUMERICAL_95_WALL_POLAR_REPLACEMENT_ESCAPE','removed':[3,42],'new_normal':p.tolist(),'target':w.tolist(),'escape':v.x.tolist(),'norm':float(v.x@Q@v.x),'active_indices':np.flatnonzero(abs(A@v.x)<1e-7).tolist()},open('data/enumerations/fresh_h20_gap_next_20261013_replacement.json','w'),indent=2);break
else:print('NOESCAPE IN200')
