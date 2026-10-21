#!/usr/bin/env python3
"""Numerical whole-chamber search after the 560 incidence-profile walls.

Derived from the previous constraint-generation scout; exact escape checks
are separate. The redundant old central inequalities are retained in the LP."""
import argparse,json,numpy as np
from pathlib import Path
from scipy.optimize import linprog
from fresh_h22_borcherds_extension_20261013 import prepare
from fresh_h22_new_cusp_completion_20261013_verify import setup,families,DEN,SHIFT,PRE
P=Path(__file__).resolve().parents[1]
def run(tight=False,lower=True):
 Q,W,G,h,X,A,b,Qn,_=prepare();_,_,_,_,_,_,D,V,B,AA,ids,_=setup()
 z=np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz');K=z['old_null_projections'];four=z['four_heights']
 U=np.load(P/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections']
 R=np.load(P/'data/enumerations/fresh_h22_gap_face_orbit_20261013.npz')['six_old_projections']
 blocks=[]
 def add(name,old,z):
  lim=-old@(Q@h)/np.sqrt(14);assert min(lim)>0
  blocks.append((name,(old@Q@X)/lim[:,None],np.broadcast_to(z,len(old))/lim))
 add('U',U/2,np.sqrt(2.5));add('D',K,four/4);add('R',R/6,np.sqrt(11))
 if tight or lower:
  from fresh_h22_new_cusp_completion_20261013_tight import shifted,DEN as DT
  sizes=np.load(P/('data/enumerations/'+PRE+'_highest_roots.npz'))['component_sizes']
  for name,old,t,ix,index in shifted(D,V,B,AA,ids,sizes,Q,lower):add(name,old/DT,t*np.sqrt(11)/10)
 else:
  for name,base,ix in families(D,V,B,AA,ids):add(name,base/DEN+SHIFT*D[ix]/6,2*SHIFT*np.sqrt(11))
 base=np.vstack([np.column_stack([A/b[:,None],np.zeros(210)]),[[0]*21+[4/np.sqrt(14)],[0]*21+[-4/np.sqrt(14)]]])
 from fresh_h22_incidence_core_20261013 import read
 core=read('results/fresh_h22_incidence_core_20261013.json')
 binary=read('results/fresh_h22_incidence_binary_20261013.json')
 add('incidence_core',np.array(core['old_projections'],dtype=float),2.)
 add('incidence_binary',np.array(binary['accepted_old_projections'],dtype=float),2.)
 # The old m=3 center inequalities are redundant positive sums of the new
 # incidence normals and may remain in the numerical LP.
 print('prepared',sum(len(p)for _,p,_ in blocks),'positive caps',flush=True)
 active=set();rng=np.random.default_rng(221350);attempts=[]
 for trial in range(30):
  direction=rng.normal(size=22)
  for step in range(150):
   chosen=sorted(active);extra=np.array([np.r_[blocks[g][1][i],s*blocks[g][2][i]]for g,i,s in chosen]).reshape((-1,22))
   lp=linprog(direction,A_ub=np.vstack([base,extra]),b_ub=np.ones(len(base)+len(extra)),bounds=[(None,None)]*22,method='highs')
   assert lp.success,lp.message
   worst=-np.inf;new=[]
   for g,(name,proj,z)in enumerate(blocks):
    vals=proj@lp.x[:21]+z*abs(lp.x[21])-1;worst=max(worst,float(max(vals)))
    wrong=np.flatnonzero(vals>2e-8)
    if len(wrong):
     wrong=wrong[np.argsort(vals[wrong])[-100:]];new.extend((g,int(i),1 if lp.x[21]>=0 else -1)for i in wrong)
   if worst<2e-8:break
   before=len(active);active.update(new);assert len(active)>before
   if step%10==0:print('trial',trial,'step',step,'active',len(active),'max',worst,flush=True)
  else:raise RuntimeError('constraint generation exceeded stated bound')
  square=float(lp.x@lp.x-1);point=np.r_[h/np.sqrt(14)+X@lp.x[:21],lp.x[21]]
  incident=[]
  for g,(name,proj,z)in enumerate(blocks):
   vals=proj@lp.x[:21]+z*abs(lp.x[21])-1
   incident.extend([name,int(i),1 if lp.x[21]>=0 else -1]for i in np.flatnonzero(abs(vals)<1e-7))
  row={'trial':trial,'point':point.tolist(),'norm':square,'old_active':np.flatnonzero(abs(base[:210]@lp.x-1)<1e-7).tolist(),'cap_active':incident,'active_constraints_used':len(active),'iterations':step+1,'max_numerical_violation':worst};attempts.append(row)
  print('trial',trial,'norm',square,'incidence',row['old_active'],incident,flush=True)
  out={'status':'NUMERICAL_FULL_COMPLETED_CUSP_ORBIT_VOLUME_SCREEN','runs':attempts,'scope':'Floating-point proposal; exact rational strict-escape certificate required.'}
  (P/('data/enumerations/fresh_h22_incidence_binary_20261013_volume_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
  if square>1e-6:break
if __name__=='__main__': run(lower=True)
