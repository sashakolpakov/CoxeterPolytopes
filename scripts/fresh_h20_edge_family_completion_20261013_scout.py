#!/usr/bin/env python3
"""Targeted old-edge directions, with two angle constraints and norm2."""
import argparse,hashlib,itertools,json
from pathlib import Path
from collections import Counter
import numpy as np
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup

ROOT=Path(__file__).resolve().parents[1]
PREFIX='fresh_h20_edge_family_completion_20261013'
LABELS=[0,2,3,4,5,6,8,10,12]

def run(target_path,out_path):
    data=json.loads((ROOT/target_path).read_text());target=data['whole_bank_geometry']
    Q,R,h,*_=extended_setup();cov=R@Q;ns=np.einsum('ij,jk,ik->i',R,Q,R);w=np.array(target['point'])
    active=[i for i in target['active_indices']if i<24]
    if len(active)>=18:supports=list(itertools.combinations(active,18))
    else:supports=[tuple(sorted(active+list(extra)))for extra in itertools.combinations([i for i in range(24)if i not in active],18-len(active))]
    assert len(supports)<=200,'Target needs a wider support design, not an unbounded expansion.'
    profiles=list(itertools.product(LABELS,repeat=2));cos=np.array([[1 if m==0 else np.cos(np.pi/m)for m in row]for row in profiles])
    finite=np.cos(np.pi/np.arange(2,121));counts=Counter();caps=[];seen=set();good_supports=[]
    for S in supports:
        S=list(S);M=R[S]@Q@R[S].T
        if np.linalg.eigvalsh(M)[0]<1e-7:continue
        good_supports.append(S);_,_,vh=np.linalg.svd(cov[S],full_matrices=True);N=vh[18:].T;H=N.T@Q@N
        T=cov@N;wc=w@Q@N;hc=h@Q@N;other=[i for i in range(87)if i not in S];start=len(caps)
        for a,b in itertools.combinations(other,2):
            A=T[[a,b]];d=np.cross(A[0],A[1]);length=np.linalg.norm(d)
            if length<1e-9:continue
            d/=length;aa=d@H@d;pinv=A.T@np.linalg.inv(A@A.T)
            rhs=-cos*np.sqrt(2*ns[[a,b]]);initial=rhs@pinv.T
            bb=2*initial@H@d;cc=np.einsum('ij,jk,ik->i',initial,H,initial)-2
            counts['profiles']+=len(profiles);disc=bb*bb-4*aa*cc
            for branch in((0,)if abs(aa)<1e-10 else(-1,1)):
                rows=np.flatnonzero(abs(bb)>1e-10)if branch==0 else np.flatnonzero(disc>=-1e-8)
                t=-cc[rows]/bb[rows]if branch==0 else(-bb[rows]+branch*np.sqrt(np.maximum(0,disc[rows])))/(2*aa)
                u=initial[rows]+t[:,None]*d;good=(u@hc<-1e-7)&(u@wc>1e-7);rows=rows[good];u=u[good]
                if not len(rows):continue
                counts['oriented_target_cutters']+=len(rows)
                ratios=-(u@T.T)/np.sqrt(2*ns);ix=np.clip(np.searchsorted(finite,ratios),1,len(finite)-1)
                err=np.minimum(abs(ratios-finite[ix]),abs(ratios-finite[ix-1]));err=np.where(ratios>=1-1e-7,0,err)
                for j in np.flatnonzero(np.max(err,axis=1)<1e-7):
                    r=N@u[j];assert abs(r@Q@r-2)<1e-5;key=tuple(np.round(r/np.linalg.norm(r),9))
                    if key in seen:continue
                    seen.add(key);caps.append({'short_orthogonal_support':S,'anchors':[a,b],'labels':list(profiles[rows[j]]),'normal':r.tolist(),
                                             'target_pairing':float(r@Q@w),'ratios':ratios[j].tolist()})
        print({'support':S,'new_caps':len(caps)-start,'caps':len(caps),'profiles':counts['profiles']},flush=True)
    out={'status':'NUMERICAL_TARGETED_18_SHORT_TWO_ANGLE_SCOUT','target_source':target_path,
         'target_sha256':hashlib.sha256((ROOT/target_path).read_bytes()).hexdigest(),'labels':LABELS,'parallel_label':0,
         'active_short_walls':active,'candidate_supports':len(supports),'supports':good_supports,'counts':dict(counts),'caps':caps,
         'scope':'All18-short supports obtained by dropping from or extending the current active short set, within explicit200-support cap. Two anchor equations against retained87 walls, both norm-quadratic branches/linear cases. Angle checks use labels2..120 or parallel/dotted numerically. Every saved proposal cuts the one target. No exact geometry/completeness claim.'}
    (ROOT/out_path).write_text(json.dumps(out,indent=2)+'\n');return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--target',default='data/enumerations/'+PREFIX+'_bank.json');p.add_argument('--output',default='data/enumerations/'+PREFIX+'_scout0.json');a=p.parse_args();run(a.target,a.output)
