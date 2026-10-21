#!/usr/bin/env python3
"""Exact rational strict escape for the complete H22 new-cusp repair."""
import argparse,json
from fractions import Fraction as F
import numpy as np
from fresh_h22_new_cusp_completion_20261013_verify import setup,families,DEN,SHIFT,PRE,P,read,sha
OUT='results/'+PRE+'_escape.json';SOURCE='results/'+PRE+'_verify.json'

def check(nums,den):
 prior,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup();n=np.array(nums[:22],dtype=object);z=abs(int(nums[22]));cov=Q.astype(object)@n
 square=int(n@cov)+z*z;hp=int(h.astype(object)@cov)
 assert square>0 and hp<0 and hp+4*z<0
 assert max(W.astype(object)@cov)<0
 bound=22*max(int(abs(x).max())for x in(D,V,B,A,W))*max(abs(int(x))for x in cov)
 assert bound<2**63
 cov=np.array(cov,dtype=np.int64)
 U=np.load(P/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections'];v=int(max(U@cov));assert v<0 and v*v>10*z*z
 zz=np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz');K=zz['old_null_projections'];four=zz['four_heights'];kp=K@cov
 assert np.all(4*kp+four*z<0)
 R=np.load(P/'data/enumerations/fresh_h22_gap_face_orbit_20261013.npz')['six_old_projections'];rp=int(max(R@cov));assert rp<0 and rp*rp>396*z*z
 oldmax={'old210':int(max(W@cov)),'twice_timelike':v,'four_null_signed':int(max(4*kp+four*z)),'six_k6_projection':rp,'center_at_m3':hp+4*z}
 maxima={}
 for name,base,ix in families(D,V,B,A,ids):
  roots=base+SHIFT*(DEN//6)*D[ix]
  assert 22*int(abs(roots).max())*int(abs(cov).max())<2**63
  maximum=int(max(roots@cov));assert maximum<0 and maximum**2>44*SHIFT**2*DEN**2*z*z
  maxima[name]=maximum
 return {'status':'PASS_EXACT_STRICT_ESCAPE_AFTER_COMPLETE_NEW_CUSP_ORBIT_REPAIR','source':SOURCE,'source_sha256':sha(SOURCE),'numerator':list(map(int,nums)),'denominator':den,'square':str(F(square,den*den)),'center_pairing':str(F(hp,den)),'old_family_integer_maxima':oldmax,'new_family_old_pair_integer_maxima':maxima,'scope':'All5,397,548 wall inequalities are strict, uniformly for every old center label m>=3 at new label k=6. Positive square proves infinite volume. All square-root comparisons are exact integer squares; both height signs are covered by absolute added coordinate.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 if a.verify:
  s=read(OUT);out=check(s['numerator'],s['denominator']);assert out==s
 else:
  s=read('data/enumerations/'+PRE+'_volume_scout.json')['runs'][0]
  _,Q,W,h,*_=setup();x=np.array(s['point']);x*=.999;x[:22]+=.001*h/np.sqrt(14)
  den=10**6;num=np.rint(x*den).astype(np.int64).tolist();out=check(num,den)
  (P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['square'])
