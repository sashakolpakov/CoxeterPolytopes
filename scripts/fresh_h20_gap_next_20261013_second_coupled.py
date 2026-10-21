import sys,json,numpy as np
sys.path[:0]=['scripts','src']
from fresh_h20_gap_next_20261013 import setup,read
from fresh_h20_edge_family_completion_20261013_bank import gate
from scipy.optimize import linprog
Q,R,h,w=setup();d=read('data/enumerations/fresh_h20_gap_next_20261013_coupled.json');U=np.vstack([R,np.array(d['new_normal'])]);ids=d['escape']['active_indices'];_,s,vh=np.linalg.svd(U[ids]@Q,full_matrices=True);p=vh[-1];p*=np.sign(-h@Q@p);p*=np.sqrt(2/(p@Q@p));r=5*p+np.sqrt(2)*(R[5]+R[14]+R[36]/3)
V=np.vstack([U,r]);ns=np.einsum('ij,jk,ik->i',V,Q,V);ratios=-(V@Q@V.T)/np.sqrt(ns[:,None]*ns[None,:]);good=gate(ratios);np.fill_diagonal(good,True)
print('norm',ns[-1],'bad',np.argwhere(~good).tolist(),'cut',r@Q@np.array(d['escape']['vector']),'finite',[(i,x)for i,x in enumerate(ratios[-1,:-1])if abs(x)>1e-7 and x<1+1e-7],flush=True)
A=V@Q/np.sqrt(ns[:,None]);hn=h/np.sqrt(-h@Q@h);rng=np.random.default_rng(20272003);trials=[];escape=None
for k in range(250):
 v=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
 if not v.success:continue
 q=float(v.x@Q@v.x);trials.append({'attempt':k,'norm':q})
 if q>1e-7:
  print('escape',k,q,'active',np.flatnonzero(abs(A@v.x)<1e-7).tolist(),flush=True)
  escape={'vector':v.x.tolist(),'norm':q,'active_indices':np.flatnonzero(abs(A@v.x)<1e-7).tolist()};break
out={'status':'NUMERICAL_COHERENT_98_WALL_SECOND_COUPLED_POLAR_REPAIR','source':'data/enumerations/fresh_h20_gap_next_20261013_coupled.json','polar':p.tolist(),'new_normal':r.tolist(),'norm':float(ns[-1]),'target':d['escape']['vector'],'target_pairing':float(r@Q@np.array(d['escape']['vector'])),'escape':escape,'trials':trials}
json.dump(out,open('data/enumerations/fresh_h20_gap_next_20261013_second_coupled.json','w'),indent=2)
