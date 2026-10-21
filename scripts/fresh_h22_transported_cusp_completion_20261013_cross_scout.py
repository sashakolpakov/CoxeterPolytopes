#!/usr/bin/env python3
"""Bounded numerical common-scale screen for complete equivariant cusp repairs."""
import json,numpy as np
from fractions import Fraction as F
from fresh_h22_new_cusp_completion_20261013_verify import P,read,RAW
PRE='fresh_h22_transported_cusp_completion_20261013'
def run():
    raw=read(RAW)['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64)
    z=np.load(P/('data/enumerations/'+PRE+'_frames.npz'));B=z['old_base_numerators'];H=z['height_numerators'];D=z['six_delta_old'];den=int(z['denominator'])
    meta=read('data/enumerations/'+PRE+'_frames.json');seed=meta['seed'];N=len(D)
    placement=read('data/enumerations/'+PRE+'_placements.json');ts=np.floor(placement['thresholds']).astype(np.int64)+1
    spec=[(i,1)for i in range(5)]+[(i,s)for i in range(5,9)for s in(1,-1)]
    q0=[F(x)for x in meta['base_norms']];q=np.array([float(q0[j])for j,s in spec])
    bases=np.array([s*B[seed,j]/den for j,s in spec]);heights=np.array([s*H[j]/den for j,s in spec]);delta=D[seed]/6
    results=[]
    for scale in(1,1.5,2,3,4,6,8):
        t=ts*scale;C=bases+t[:,None]*delta;Z=heights+2*t;cv=Q@C.T
        worst=-np.inf;where=None;bad=0
        for k,(j,sgn)in enumerate(spec):
            roots=sgn*B[:,j]/den+t[k]*D/6;zz=sgn*H[j]/den+2*t[k]
            pair=roots@cv+11*zz*Z;measure=pair+np.sqrt(q[k]*q);measure[seed]=-np.inf
            maximum=float(measure.max());bad+=int(np.sum(measure>1e-6))
            if maximum>worst:
                worst=maximum;ii=np.unravel_index(np.argmax(measure),measure.shape);where=[k,int(ii[0]),int(ii[1]),1]
            opposite=pair-22*zz*Z;maximum=float((opposite+np.sqrt(q[k]*q)).max());bad+=int(np.sum(opposite+np.sqrt(q[k]*q)>1e-6))
            if maximum>worst:worst=maximum;where=[k,-1,-1,-1]
        rec={'scale':scale,'integer_starting_shifts':ts.tolist(),'shifts':t.tolist(),'nonlocal_dotted_excess':worst,'failed_dotted_pairs':bad,'worst_location':where};results.append(rec);print(rec,flush=True)
        (P/('data/enumerations/'+PRE+'_cross_scout.json')).write_text(json.dumps({'status':'NUMERICAL_COMMON_SCALE_CROSS_SCREEN','runs':results},indent=2)+'\n')
        if worst<-1e-6:break
if __name__=='__main__':run()
