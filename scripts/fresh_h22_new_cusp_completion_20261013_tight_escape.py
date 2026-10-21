#!/usr/bin/env python3
"""Recheck the exact remote escape against the closer cusp-repair placement."""
import argparse,json
import numpy as np
from fresh_h22_new_cusp_completion_20261013_verify import setup,PRE,P,read,sha
from fresh_h22_new_cusp_completion_20261013_escape import check
from fresh_h22_new_cusp_completion_20261013_tight import shifted,DEN
OUT='results/'+PRE+'_tight_escape.json';SOURCE='results/'+PRE+'_tight.json'
def run():
 s=read('results/'+PRE+'_escape.json');base=check(s['numerator'],s['denominator'])
 assert base==s
 prior,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup();cov=Q.astype(object)@np.array(s['numerator'][:22],dtype=object);z=abs(s['numerator'][22]);cov=np.array(cov,dtype=np.int64);maximums={}
 for name,roots,t,ix,index in shifted(D,V,B,A,ids,sizes):
  assert 22*int(abs(roots).max())*int(abs(cov).max())<2**63
  pair=roots@cov;maximums[name]={}
  for value in np.unique(t):
   v=int(value);m=int(max(pair[t==v]));assert m<0 and 100*m*m>11*v*v*DEN*DEN*z*z
   maximums[name][str(v)]=m
 return {'status':'PASS_EXACT_STRICT_ESCAPE_FOR_TIGHT_CUSP_REPAIR','source':SOURCE,'source_sha256':sha(SOURCE),'uniform_escape_source':'results/'+PRE+'_escape.json','uniform_escape_sha256':sha('results/'+PRE+'_escape.json'),'numerator':s['numerator'],'denominator':s['denominator'],'square':s['square'],'new_integer_maxima_by_twenty_shift':maximums,'scope':'All5,397,548 tight-placement inequalities are strict; the same rational escape works for all old center labels m>=3. The uniform and tight chambers are separately certified constructions, not their union.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
 if a.verify:assert d==read(OUT)
 else:(P/OUT).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['square'])
