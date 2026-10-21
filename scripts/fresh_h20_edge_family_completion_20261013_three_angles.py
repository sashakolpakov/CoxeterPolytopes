#!/usr/bin/env python3
"""A target-driven17-short residual4-space with three angle constraints."""
import argparse,hashlib,itertools,json
from collections import Counter
import numpy as np
from fresh_h20_edge_family_completion_20261013_bank import ROOT,PREFIX,gate,read
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup
LABELS=[0,2,3,4,6]

def run(target_path,out_path):
    data=read(target_path);target=data['whole_bank_geometry'];Q,R,h,*_=extended_setup();cov=R@Q;ns=np.einsum('ij,jk,ik->i',R,Q,R);w=np.array(target['point'])
    active=[i for i in target['active_indices']if i<24]
    if len(active)>=17:supports=list(itertools.combinations(active,17))
    else:supports=[tuple(sorted(active+list(e)))for e in itertools.combinations([i for i in range(24)if i not in active],17-len(active))]
    assert len(supports)<=40
    profiles=list(itertools.product(LABELS,repeat=3));cos=np.array([[1 if m==0 else np.cos(np.pi/m)for m in p]for p in profiles])
    caps=[];seen=set();counts=Counter();passed=[]
    for S in supports:
        S=list(S)
        if np.linalg.eigvalsh(R[S]@Q@R[S].T)[0]<1e-7:continue
        passed.append(S);_,_,vh=np.linalg.svd(cov[S],full_matrices=True);N=vh[17:].T;H=N.T@Q@N;T=cov@N;wc=w@Q@N;hc=h@Q@N
        other=[i for i in range(87)if i not in S];start=len(caps)
        for ai,a in enumerate(other):
            for b,c in itertools.combinations(other[ai+1:],2):
                A=T[[a,b,c]];_,sing,vh=np.linalg.svd(A,full_matrices=True)
                if sing[-1]<1e-8:continue
                d=vh[-1];aa=d@H@d;pinv=A.T@np.linalg.inv(A@A.T);rhs=-cos*np.sqrt(2*ns[[a,b,c]])
                initial=rhs@pinv.T;bb=2*initial@H@d;cc=np.einsum('ij,jk,ik->i',initial,H,initial)-2
                counts['profiles']+=len(profiles);disc=bb*bb-4*aa*cc
                for branch in((0,)if abs(aa)<1e-10 else(-1,1)):
                    rows=np.flatnonzero(abs(bb)>1e-10)if branch==0 else np.flatnonzero(disc>=-1e-8)
                    t=-cc[rows]/bb[rows]if branch==0 else(-bb[rows]+branch*np.sqrt(np.maximum(0,disc[rows])))/(2*aa)
                    u=initial[rows]+t[:,None]*d;good=(u@wc>1e-7)&(u@hc<-1e-7);rows=rows[good];u=u[good]
                    if not len(rows):continue
                    counts['oriented_cutters']+=len(rows);ratios=-(u@T.T)/np.sqrt(2*ns)
                    for k in np.flatnonzero(np.all(gate(ratios),axis=1)):
                        r=N@u[k]
                        if abs(r@Q@r-2)>1e-5:counts['ill_conditioned_rejected']+=1;continue
                        key=tuple(np.round(r/np.linalg.norm(r),9))
                        if key in seen:continue
                        seen.add(key);caps.append({'short_orthogonal_support':S,'anchors':[a,b,c],'labels':list(profiles[rows[k]]),'normal':r.tolist(),'target_pairing':float(r@Q@w),'ratios':ratios[k].tolist()})
            if ai%20==0:print({'support':S,'anchor_prefix':ai,'caps':len(caps),'profiles':counts['profiles']},flush=True)
        print({'completed_support':S,'new_caps':len(caps)-start,'caps':len(caps)},flush=True)
    out={'status':'NUMERICAL_TARGETED_17_SHORT_THREE_ANGLE_SCOUT','target_source':target_path,'target_sha256':hashlib.sha256((ROOT/target_path).read_bytes()).hexdigest(),
         'labels':LABELS,'parallel_label':0,'active_short_walls':active,'supports':passed,'counts':dict(counts),'caps':caps,
         'scope':'All17-short supports obtained by extending/dropping current active shorts within explicit40-support cap. Three retained87-wall anchor equations; norm2 quadratic both branches or linear case; numerical labels2..120 or dotted gate. New direction dimension; no exact positive volume claim or exhaustive continuous-angle theorem.'}
    (ROOT/out_path).write_text(json.dumps(out,indent=2)+'\n');return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--target',default='data/enumerations/'+PREFIX+'_bank4.json');p.add_argument('--output',default='data/enumerations/'+PREFIX+'_three_angles.json');a=p.parse_args();run(a.target,a.output)
