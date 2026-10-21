#!/usr/bin/env python3
"""Choose mutually Coxeter-compatible interior caps at quarter-vertex centers."""
import json,math,time
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import lil_matrix
from fresh_h20_intrinsic_caps_20261014 import ROOT,TAG,data,sha
s=(1+math.sqrt(5))/2;D=data(s);V=np.array([v['point']for v in D['vertices']]);C=V/4;R=1.5
DD=np.maximum(np.sum(C*C,axis=1)[:,None]+np.sum(C*C,axis=1)[None,:]-2*C@C.T,0)
co=DD/(2*R*R)-1;compatible=(co>=1-1e-8)|(abs(co)<1e-8)|(abs(co-.5)<1e-8)|(abs(co-s/2)<1e-8)
np.fill_diagonal(compatible,True)
cover=DD<=R*R+1e-8
conflicts=np.array(np.where(np.triu(~compatible,1))).T
A=lil_matrix((len(conflicts)+252,252));lo=np.full(A.shape[0],-np.inf);hi=np.ones(A.shape[0])
for k,(i,j)in enumerate(conflicts):A[k,i]=A[k,j]=1
for k in range(252):A[len(conflicts)+k,np.flatnonzero(cover[k])]=1;lo[len(conflicts)+k]=1;hi[len(conflicts)+k]=np.inf
print('conflicts',len(conflicts),'coverage counts',np.unique(cover.sum(axis=1),return_counts=True),flush=True)
t=time.time();res=milp(np.ones(252),integrality=np.ones(252),bounds=Bounds(0,1),constraints=LinearConstraint(A.tocsr(),lo,hi),options={'time_limit':60})
selected=[]
if res.x is not None:
 selected=np.flatnonzero(res.x>.5).tolist();assert compatible[np.ix_(selected,selected)].all()and cover[:,selected].any(axis=1).all()
out={'status':'NUMERICAL_QUARTER_CAP_COMPATIBILITY_COVER_SEARCH','sources':{p:sha(p)for p in[f'scripts/{TAG}.py',f'scripts/{TAG}_quarter_scout.py']},'caps':252,'radius':'3/2','center_scale':'1/4','conflict_pairs':len(conflicts),'solver_status':int(res.status),'solver_message':res.message,'selected':selected,'elapsed_seconds':time.time()-t,'scope':'Finite numerical binary compatibility/coverage model; covers only252 quarter-vertex probes. A feasible solution would need exact pair checks and global volume certification. Infeasibility does not exclude other caps.'}
(ROOT/f'results/{TAG}_quarter_scout.json').write_text(json.dumps(out,indent=2)+'\n');print(out,flush=True)
