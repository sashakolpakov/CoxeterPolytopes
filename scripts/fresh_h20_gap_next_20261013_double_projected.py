import sys,json,numpy as np
sys.path[:0]=['scripts','src']
from fresh_h20_gap_next_20261013 import setup,read,TARGET
from fresh_h20_edge_family_completion_20261013_bank import gate
from scipy.optimize import linprog
Q,R,h,w=setup();ids=read(TARGET)['certificate']['numeric_active_indices'];_,s,vh=np.linalg.svd(R[ids]@Q,full_matrices=True);p=vh[-1];p*=np.sign(-h@Q@p);p*=np.sqrt(2/(p@Q@p))
a=p-(R[3]@Q@p)/2*R[3];b=p-(R[42]@Q@p)/18*R[42];U=np.vstack([R,a,b]);ns=np.einsum('ij,jk,ik->i',U,Q,U);ratios=-(U@Q@U.T)/np.sqrt(ns[:,None]*ns[None,:]);good=gate(ratios);np.fill_diagonal(good,True)
print('norms',ns[-2:],'bad',np.argwhere(~good).tolist(),'cross',a@Q@b,'cuts',U[-2:]@Q@w,flush=True)
delta=a+b;pp=U@Q@delta;inc=np.flatnonzero(abs(pp)<1e-7);print('cusp',delta@Q@delta,'inc',inc.tolist(),'rank',np.linalg.matrix_rank(U[inc]@Q@U[inc].T,tol=1e-7),flush=True)
A=U@Q/np.sqrt(ns[:,None]);hn=h/np.sqrt(-h@Q@h);rng=np.random.default_rng(20272003);trials=[];escape=None
for k in range(200):
 v=linprog(rng.normal(size=21),A_ub=A,b_ub=np.zeros(len(A)),A_eq=(hn@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
 if not v.success:continue
 q=float(v.x@Q@v.x);trials.append({'attempt':k,'norm':q})
 if q>1e-7:
  print('escape',k,q,'active',np.flatnonzero(abs(A@v.x)<1e-7).tolist(),flush=True)
  escape={'vector':v.x.tolist(),'norm':q,'active_indices':np.flatnonzero(abs(A@v.x)<1e-7).tolist()};break
out={'status':'NUMERICAL_COHERENT_98_WALL_DOUBLE_PROJECTION_REPAIR','source':TARGET,'polar':p.tolist(),'new_normals':[a.tolist(),b.tolist()],'new_norms':ns[-2:].tolist(),'target':w.tolist(),'escape':escape,'trials':trials,'new_cusp':{'ray':delta.tolist(),'incident_indices':inc.tolist()}}
json.dump(out,open('data/enumerations/fresh_h20_gap_next_20261013_double_projected.json','w'),indent=2)
