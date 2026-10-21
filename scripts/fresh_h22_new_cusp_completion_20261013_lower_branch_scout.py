import sys,json,numpy as np
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fresh_h22_new_cusp_completion_20261013_verify import setup,DEN,PRE,P,read
p,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup();c=read('results/'+PRE+'_verify.json');bs=np.array(c['base_old_numerators'])/DEN;qn=np.array([float(F(x))for x in c['new_cap_norms']]);dp=W@(Q@D[0])/6;vp=(W@Q@bs.T).T;wn=np.diag(W@Q@W.T)
U=np.load(P/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections']/2
zz=np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz');K=zz['old_null_projections'];kt=zz['four_heights']/4
R=np.load(P/'data/enumerations/fresh_h22_gap_face_orbit_20261013.npz')['six_old_projections']/6
cos=np.cos(np.pi/np.arange(2,101));rows=[]
for j,ms in [(0,[3]),(4,[2]),(5,[2]),(6,[4,6,7]),(10,[3,4,5,6])]:
 vals=np.where(dp<0,(vp[j]+np.sqrt(qn[j]*wn))/np.maximum(-dp,1e-10),-np.inf);i=int(np.argmax(vals));a=vp[j,i];b=-dp[i];r=np.sqrt(qn[j]*wn[i])
 for m in ms:
  t=(a+r*np.cos(np.pi/m))/b;old=bs[j]+t*D[0]/6;height=2*t*np.sqrt(11);q=qn[j];groups=[('old',W,0,wn),('U',U,np.sqrt(2.5),1),('D',K,kt,kt*kt),('R',R,np.sqrt(11),2/3),('center',h[None,:],4,2)];bad=[];counts={}
  for name,rs,z,n in groups:
   vals=-(rs@(Q@old)+z*height)/np.sqrt(q*n);finite=np.flatnonzero(vals<1-1e-7)
   err=np.min(abs(vals[finite,None]-cos),axis=1);wrong=finite[err>1e-6];counts[name]=len(wrong)
   bad.extend((name,int(i),float(vals[i]))for i in wrong[:3])
  rec={'seed':j,'anchor':i,'label':m,'shift':t,'bad_counts':counts,'first_bad':bad[:4]};rows.append(rec);print(rec,flush=True)
(P/('data/enumerations/'+PRE+'_lower_branch_scout.json')).write_text(json.dumps(rows,indent=2)+'\n')
