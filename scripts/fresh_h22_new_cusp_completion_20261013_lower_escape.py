#!/usr/bin/env python3
"""Exact target cut and new strict escape for the finite-label cusp placements."""
import argparse,json
import numpy as np
from fresh_h22_new_cusp_completion_20261013_verify import setup,PRE,P,read,sha
from fresh_h22_new_cusp_completion_20261013_escape import check
from fresh_h22_new_cusp_completion_20261013_tight import shifted,DEN
OUT='results/'+PRE+'_lower_escape.json';SOURCE='results/'+PRE+'_lower.json';OLD='results/'+PRE+'_escape.json'
def run(nums,den):
 base=check(nums,den)
 prior,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup();cov=Q.astype(object)@np.array(nums[:22],dtype=object);z=abs(nums[22]);cov=np.array(cov,dtype=np.int64)
 old=read(OLD);ocov=Q.astype(object)@np.array(old['numerator'][:22],dtype=object);oz=old['numerator'][22];assert oz>0;ocov=np.array(ocov,dtype=np.int64)
 maxima={};cut=None
 for name,roots,t,ix,index in shifted(D,V,B,A,ids,sizes,Q,True):
  assert 22*int(abs(roots).max())*max(int(abs(cov).max()),int(abs(ocov).max()))<2**63
  pair=roots@cov;op=roots@ocov;maxima[name]={}
  for value in np.unique(t):
   v=int(value);subset=np.flatnonzero(t==v);m=int(max(pair[subset]));assert m<0 and 100*m*m>11*v*v*DEN*DEN*z*z;maxima[name][str(v)]=m
   k=int(subset[np.argmax(op[subset])]);a=10*int(op[k]);b=v*oz*DEN
   if cut is None and(a>=0 or 11*b*b>a*a):
    cut={'family':name,'entry':k,'cusp_index':int(ix[k]),'old_normal_numerator':roots[k].tolist(),'old_normal_denominator':DEN,'twenty_shift':v,'pairing_numerator_rational':str(a),'pairing_numerator_sqrt11':str(b),'pairing_denominator':10*DEN*old['denominator'],'positive_sign_reason':'B>0 and either A>=0 or11B^2>A^2'}
 assert cut is not None
 return {'status':'PASS_EXACT_LOWER_REPAIR_TARGET_CUT_AND_NEW_STRICT_ESCAPE','source':SOURCE,'source_sha256':sha(SOURCE),'previous_escape_source':OLD,'previous_escape_sha256':sha(OLD),'old_target_cut':cut,'numerator':list(map(int,nums)),'denominator':den,'square':base['square'],'center_pairing':base['center_pairing'],'old_family_integer_maxima':base['old_family_integer_maxima'],'new_integer_maxima_by_twenty_shift':maxima,'scope':'The lower-placement chamber cuts the exact preceding uniform/tight remote escape. Another rational point is strict on all5,397,548 facets and spacelike, proving infinite volume for every old center label m>=3. The uniform check used here is a stronger explicit replay dependency for the inherited old constraints; no numerical inequality is accepted.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 if a.verify:s=read(OUT);d=run(s['numerator'],s['denominator']);assert d==s
 else:
  row=read('data/enumerations/'+PRE+'_lower_volume_scout.json')['runs'][0];_,Q,W,h,*_=setup();x=np.array(row['point']);x*=.999;x[:22]+=.001*h/np.sqrt(14);den=10**6;nums=np.rint(x*den).astype(np.int64).tolist();d=run(nums,den);(P/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['square'])
