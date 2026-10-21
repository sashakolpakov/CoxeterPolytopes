#!/usr/bin/env python3
"""One16-short residual5-space; three active long anchors plus a fourth."""
import itertools,json
from collections import Counter
import numpy as np
from fresh_h20_edge_family_completion_20261013_bank import ROOT,PREFIX,gate,read,sha
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup
TARGET='data/enumerations/'+PREFIX+'_bank7.json'
OUT='data/enumerations/'+PREFIX+'_four_angles.json'
LABELS=[0,2,3,4,6]

def run():
    target=read(TARGET)['whole_bank_geometry'];Q,R,h,*_=extended_setup();cov=R@Q;ns=np.einsum('ij,jk,ik->i',R,Q,R);w=np.array(target['point'])
    S=[i for i in target['active_indices']if i<24];fixed=[i for i in target['active_indices']if 24<=i<48]
    assert len(S)==16 and fixed==[24,30,38]and np.linalg.eigvalsh(R[S]@Q@R[S].T)[0]>1e-7
    _,_,vh=np.linalg.svd(cov[S],full_matrices=True);N=vh[16:].T;H=N.T@Q@N;T=cov@N;wc=w@Q@N;hc=h@Q@N
    profiles=list(itertools.product(LABELS,repeat=4));cos=np.array([[1 if m==0 else np.cos(np.pi/m)for m in p]for p in profiles]);seen=set();caps=[];counts=Counter()
    for fourth in range(87):
        if fourth in S or fourth in fixed:continue
        ids=fixed+[fourth];A=T[ids];_,sing,vh=np.linalg.svd(A,full_matrices=True)
        if sing[-1]<1e-8:continue
        d=vh[-1];aa=d@H@d;pinv=A.T@np.linalg.inv(A@A.T);rhs=-cos*np.sqrt(2*ns[ids]);initial=rhs@pinv.T
        bb=2*initial@H@d;cc=np.einsum('ij,jk,ik->i',initial,H,initial)-2;disc=bb*bb-4*aa*cc;counts['profiles']+=len(profiles)
        for branch in((0,)if abs(aa)<1e-10 else(-1,1)):
            rows=np.flatnonzero(abs(bb)>1e-10)if branch==0 else np.flatnonzero(disc>=-1e-8)
            t=-cc[rows]/bb[rows]if branch==0 else(-bb[rows]+branch*np.sqrt(np.maximum(0,disc[rows])))/(2*aa)
            u=initial[rows]+t[:,None]*d;good=(u@wc>1e-7)&(u@hc<-1e-7);rows=rows[good];u=u[good]
            if not len(rows):continue
            counts['oriented_cutters']+=len(rows);ratios=-(u@T.T)/np.sqrt(2*ns)
            for k in np.flatnonzero(np.all(gate(ratios),axis=1)):
                r=N@u[k]
                if abs(r@Q@r-2)>1e-5:continue
                key=tuple(np.round(r/np.linalg.norm(r),9))
                if key in seen:continue
                seen.add(key);caps.append({'short_orthogonal_support':S,'anchors':ids,'labels':list(profiles[rows[k]]),'normal':r.tolist(),'target_pairing':float(r@Q@w),'ratios':ratios[k].tolist()})
    out={'status':'NUMERICAL_TARGETED_16_SHORT_FOUR_ANGLE_SCOUT','target_source':TARGET,'target_sha256':sha(TARGET),'short_orthogonal_support':S,'fixed_long_anchors':fixed,
         'labels':LABELS,'parallel_label':0,'counts':dict(counts),'caps':caps,
         'scope':'One residual5-space of16 active oldshortwalls; fixed three active lifted-long anchors plus every remaining retained87-wall fourth anchor. Two norm-quadratic branches/linear case; bounded labels; numeric all87 pair gate and target cut. No completeness or exact discovery claim.'}
    (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print({'counts':dict(counts),'caps':len(caps)});return out

if __name__=='__main__':run()
