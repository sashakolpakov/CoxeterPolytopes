#!/usr/bin/env python3
"""Global H5 halfspace relaxation of saved candidate banks; numerical only."""
import json,sys,hashlib,time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fresh_h9_h5_closure_next_20261013_scout as S
P=S.P;TAG='fresh_h9_global_h5_20261013'
BANKS=[
 ('fresh_h9_higher_face_20261013_scout','proposals'),
 ('fresh_h9_higher_face_20261013_branch_scout','bank'),
 ('fresh_h9_higher_face_20261013_profile_scout','candidates'),
 ('fresh_h9_higher_face_20261013_finite_truncation_scout','candidates'),
 ('fresh_h9_higher_face_20261013_truncation_scout','candidates'),
 ('fresh_h9_h5_closure_20261013_conic_scout','candidates'),
 ('fresh_h9_h5_closure_20261013_conic_scout_stage1','candidates'),
 ('fresh_h9_h5_closure_20261013_conic_scout_stage2','candidates'),
 ('fresh_h9_h5_closure_20261013_endpoint_scout','accepted'),
 ('fresh_h9_h5_closure_next_20261013_scout','candidates'),
 ('fresh_h9_h5_closure_next_20261013_stage1_scout','candidates'),
 ('fresh_h9_h5_closure_next_20261013_stage2_scout','candidates'),
]
def load():
 R,G,c,names,doc,p,C,J,RF,phys,target=S.load()
 path='data/enumerations/fresh_h9_h5_closure_next_20261013_verify.json';d=S.S.B.V.read(path)
 for row in d['new_caps']:
  R=np.vstack([R,np.array([[(int(a)+int(b))/(2*10**80)for a,b in row['normal']]])]);names.append(row['specification']['name'])
 assert len(R)==125
 roots=[];rows=[];sources={path:S.S.B.V.sha(path)};counts={}
 for tag,key in BANKS:
  path='data/enumerations/'+tag+'.json';data=json.loads((P/path).read_text());sources[path]=S.S.B.V.sha(path);counts[tag]=len(data[key])
  for i,row in enumerate(data[key]):
   assert len(row['normal'])==10;roots.append(row['normal']);rows.append({'source':tag,'row':i,'array':key})
 return R,G,c,names,d,p,C,np.array(roots),rows,sources,counts
def run():
 R,G,c,names,doc,p,C,NR,rows,sources,counts=load();W=np.vstack([R,NR]);A=W@G@C;b=-W@G@p;keep=np.linalg.norm(A,axis=1)>1e-7
 assert min(b[keep])>0;aa=A[keep]/b[keep,None];rhs=np.ones(len(aa));old=S.S.geometry(R,G,p,C);old.sort(key=lambda v:v['norm'],reverse=True)
 checks=[];best=None;start=time.time()
 # First use directions of every original projective vertex, then fixed random
 # objectives. One strict spacelike feasible point suffices for a negative bank.
 dirs=[C.T@G@np.array(v['point'])for v in old]
 dirs+=list(np.random.default_rng(913).normal(size=(100,5)))
 for i,direction in enumerate(dirs):
  ans=linprog(-direction,A_ub=aa,b_ub=rhs,bounds=[(None,None)]*5,method='highs')
  if not ans.success:checks.append({'objective':i,'status':int(ans.status)});continue
  z=ans.x;q=float(z@z-1);checks.append({'objective':i,'norm':q})
  if best is None or q>best['norm']:
   x=p+C@z;best={'norm':q,'point':x.tolist(),'slice_coordinates':z.tolist(),'maximum_normalized_residual':float(max(aa@z-1)),'active_entries':np.flatnonzero(abs(W@G@x)<2e-6).tolist()}
   print('best',i,q,'seconds',time.time()-start,flush=True)
  if q>1e-4:
   # Move a controlled amount toward the common H5 center to make every
   # nonfixed inequality strict, preserving spacelikeness.
   epsilon=min(.001,(1-1/np.sqrt(q+1))/4);zz=(1-epsilon)*z;xx=p+C@zz
   best['strict_point']=xx.tolist();best['strict_norm']=float(zz@zz-1);best['strict_maximum_normalized_residual']=float(max(aa@zz-1));break
 out={'status':'NUMERICAL_ALL_SAVED_BANK_HALFSPACE_RELAXATION','sources':sources,'bank_counts':counts,'literal_candidates':len(NR),'retained_walls':125,'best':best,'lp_checks':checks,'scope':'Numerical common-halfspace test only, ignoring all candidate conflicts. Candidate formulas are not certified here. A strict common spacelike point, if exactified, rules out every subset of this finite literal bank.'}
 (P/'data/enumerations'/f'{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n');print(out['literal_candidates'],out['best']['norm'],flush=True)
if __name__=='__main__':run()
