#!/usr/bin/env python3
"""Two prescribed angles in a residual (2,1) space for the H22 scout."""
import itertools,json
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h22_borcherds_extension_20261013 import prepare
ROOT=Path(__file__).resolve().parents[1]
TARGET='data/enumerations/fresh_h22_borcherds_extension_20261013_scout.json'
OUT='data/enumerations/fresh_h22_two_anchor_20261013.json'
LABELS=[0,2,3,4,5,6,8,10,12]

def run():
    oldQ,W,G,h,X,A,b,Q,R=prepare();g=json.loads((ROOT/TARGET).read_text())
    profiles=list(itertools.product(LABELS,repeat=2));cos=np.array([[1 if m==0 else np.cos(np.pi/m) for m in row] for row in profiles]);finite=np.cos(np.pi/np.arange(2,121))
    counts=Counter();runs=[];caps=[]
    for record in g['runs']:
        m=record['center_short_label'];T=record['center_height'];roots=np.vstack([R,np.r_[h,T],np.r_[h,-T]])
        center=np.r_[h,0];cov=roots@Q;ns=np.einsum('ij,jk,ik->i',roots,Q,roots)
        target=record['escape'];w=np.array(target['point']);active=target['active_old']
        supports=list(itertools.combinations(active,20)) if len(active)>=20 else [tuple(sorted(active+list(extra))) for extra in itertools.combinations([i for i in range(42) if i not in active],20-len(active))]
        # Bound the broad lower-incidence branch explicitly.
        supports=supports[:40];start=len(caps);accepted=[]
        for S in supports:
            M=roots[list(S)]@Q@roots[list(S)].T
            if np.linalg.eigvalsh(M)[0]<1e-7:continue
            _,_,vh=np.linalg.svd(cov[list(S)],full_matrices=True);N=vh[20:].T;H=N.T@Q@N;P=cov@N;wc=w@Q@N;hc=center@Q@N
            accepted.append(list(S));others=[i for i in range(212) if i not in S]
            seen=set()
            for a,bb in itertools.combinations(others,2):
                B=P[[a,bb]];d=np.cross(B[0],B[1]);length=np.linalg.norm(d)
                if length<1e-9:continue
                d/=length;aa=d@H@d;pinv=B.T@np.linalg.inv(B@B.T)
                rhs=-cos*np.sqrt(2*ns[[a,bb]]);initial=rhs@pinv.T
                linear=2*initial@H@d;cc=np.einsum('ij,jk,ik->i',initial,H,initial)-2;disc=linear*linear-4*aa*cc;counts['quadratics']+=len(profiles)
                for branch in ((0,) if abs(aa)<1e-10 else(-1,1)):
                    rows=np.flatnonzero(abs(linear)>1e-10) if branch==0 else np.flatnonzero(disc>=-1e-8)
                    t=-cc[rows]/linear[rows] if branch==0 else(-linear[rows]+branch*np.sqrt(np.maximum(0,disc[rows])))/(2*aa)
                    u=initial[rows]+t[:,None]*d;good=(u@hc<-1e-7)&(u@wc>1e-7);rows=rows[good];u=u[good]
                    if not len(rows):continue
                    counts['oriented_cutters']+=len(rows);ratios=-(u@P.T)/np.sqrt(2*ns)
                    ix=np.clip(np.searchsorted(finite,ratios),1,len(finite)-1);err=np.minimum(abs(ratios-finite[ix]),abs(ratios-finite[ix-1]));err=np.where(ratios>=1-1e-7,0,err)
                    for j in np.flatnonzero(np.max(err,axis=1)<1e-7):
                        r=N@u[j];key=tuple(np.round(r/np.linalg.norm(r),9))
                        if key in seen:continue
                        seen.add(key);assert abs(r@Q@r-2)<1e-5
                        caps.append({'center_label':m,'orthogonal_support':list(S),'anchors':[a,bb],'labels':list(profiles[rows[j]]),'normal':r.tolist(),'target_pairing':float(r@Q@w)})
            print('m',m,'support',len(accepted),'caps',len(caps)-start,'quadratics',counts['quadratics'],flush=True)
        runs.append({'center_label':m,'supports':accepted,'new_caps':len(caps)-start})
        (ROOT/OUT).write_text(json.dumps({'status':'NUMERICAL_BOUNDED_H22_TWO_ANCHOR_SCOUT','target':TARGET,'labels':LABELS,'counts':dict(counts),'runs':runs,'caps':caps,'scope':'At most40 lexicographic20-wall supports per saved target; prescribed two finite/parallel angles, norm2, all212 pair gates. Max numerical label120. Not a complete angle or direction search.'},indent=2)+'\n')
    print(dict(counts),len(caps),flush=True)

if __name__=='__main__':run()
