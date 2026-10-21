#!/usr/bin/env python3
"""Changed local H20 search: one finite anchor and a free conic parameter.

Numerical proposal only. For each fixed profile the remaining nonconstant
pairings are required strictly dotted throughout an open parameter interval.
"""
import argparse,itertools,json,math
from collections import Counter
from pathlib import Path
import numpy as np
from fresh_h20_edge_family_completion_20261013_four_volume import setup as old_setup,read,sha
from fresh_h20_targeted_escape_20261013 import ev
from fresh_h20_edge_family_completion_20261013_bank import gate
ROOT=Path(__file__).resolve().parents[1]
BASE='data/enumerations/fresh_h20_edge_family_completion_20261013_cusp_repairs.json'
TARGET='data/enumerations/fresh_h20_edge_family_completion_20261013_four_volume.json'
OUT='data/enumerations/fresh_h20_gap_next_20261013.json'
LABELS=[2,3,4,6]

def setup():
    Q,R,h,*_=old_setup();data=read(BASE)
    R=np.vstack([R,np.array([[ev(x)for x in c['normal']]for c in data['new_caps']])])
    w=np.array([ev(x)for x in data['strict_escape']['vector']])
    assert len(R)==96
    return Q,R,h,w

def less_zero(a,b,c):
    """Open positive-x intervals where a*x*x+b*x+c<0, floating scout only."""
    scale=max(abs(a),abs(b),abs(c),1.);tol=1e-11*scale
    if abs(a)<tol:
        if abs(b)<tol:return[(0.,math.inf)]if c<-tol else[]
        z=-c/b
        if b>0:return[(0.,z)]if z>0 else[]
        return[(max(0.,z),math.inf)]
    disc=b*b-4*a*c
    if disc<=tol*scale:
        return[(0.,math.inf)]if a<0 else[]
    d=math.sqrt(disc);z1,z2=sorted(((-b-d)/(2*a),(-b+d)/(2*a)))
    if a>0:return[(max(0.,z1),z2)]if z2>max(0.,z1)else[]
    out=[]
    if z1>0:out.append((0.,z1))
    out.append((max(0.,z2),math.inf));return out

def intersect(I,J):
    out=[]
    for a,b in I:
        for c,d in J:
            lo,hi=max(a,c),min(b,d)
            if hi>lo+1e-9*max(1.,lo):out.append((lo,hi))
    return out

def samples(lo,hi):
    if math.isinf(hi):return[max(1.,lo*2+1),max(2.,lo*4+2)]
    if lo==0:return[hi/2,hi/4]
    return[(lo+hi)/2,math.sqrt(lo*hi)]

def run():
    Q,R,h,w=setup();cov=R@Q;ns=np.einsum('ij,jk,ik->i',R,Q,R);old=read(TARGET)
    active=[i for i in old['certificate']['numeric_active_indices']if i<24]
    assert len(active)==19
    caps=[];seen=set();count=Counter();families=[]
    for omitted in active:
        S=[i for i in active if i!=omitted]
        M=R[S]@Q@R[S].T
        if np.linalg.eigvalsh(M)[0]<=1e-8:continue
        for a in range(82):
            if a in S:continue
            J=S+[a];G=R[J]@Q@R[J].T
            if np.linalg.eigvalsh(G)[0]<=1e-8:continue
            _,_,vh=np.linalg.svd(cov[J],full_matrices=True);N=vh[19:].T
            vals,V=np.linalg.eigh(N.T@Q@N)
            assert vals[0]<-1e-8 and vals[1]>1e-8
            u=N@V[:,1]/math.sqrt(vals[1]);v=N@V[:,0]/math.sqrt(-vals[0])
            pu,pv=cov@u,cov@v
            fixed=(abs(pu)+abs(pv)<1e-7)
            for m in LABELS:
                rhs=np.zeros(19);rhs[-1]=-math.sqrt(2*ns[a])*math.cos(math.pi/m)
                r0=np.linalg.solve(G,rhs)@R[J];p0=cov@r0;kappa=2-r0@Q@r0;count['profiles']+=1
                if np.max(abs(p0[J]-rhs))>1e-6:count['ill_conditioned']+=1;continue
                if not np.all(gate(-p0[fixed]/np.sqrt(2*ns[fixed]))):count['constant_forbidden']+=1;continue
                if abs(kappa)<1e-8:count['null_residual_skipped']+=1;continue
                for sign in(-1,1):
                    if kappa>0:B=sign*math.sqrt(kappa)*u;C=math.sqrt(kappa)*v
                    else:B=sign*math.sqrt(-kappa)*v;C=math.sqrt(-kappa)*u
                    pb,pc=cov@B,cov@C
                    # First orient/cut, then force every variable relation dotted.
                    conditions=[(-((B+C)@Q@w),-2*(r0@Q@w),-((B-C)@Q@w)),
                                ((B+C)@Q@h,2*(r0@Q@h),(B-C)@Q@h)]
                    conditions += [(pb[i]+pc[i],2*(p0[i]+math.sqrt(2*ns[i])),pb[i]-pc[i])for i in range(96)if not fixed[i]]
                    intervals=[(0.,math.inf)]
                    for aa,bb,cc in conditions:
                        intervals=intersect(intervals,less_zero(aa,bb,cc))
                        if not intervals:break
                    if not intervals:continue
                    count['open_interval_families']+=1
                    families.append({'short_support':S,'anchor':a,'label':m,'branch':sign,'kappa':float(kappa),'intervals':[[lo,None if math.isinf(hi)else hi]for lo,hi in intervals]})
                    for lo,hi in intervals:
                        for x in samples(lo,hi):
                            r=r0+B*(x+1/x)/2+C*(x-1/x)/2
                            ratios=-(cov@r)/np.sqrt(2*ns)
                            if abs(r@Q@r-2)>1e-5 or r@Q@w<=1e-8 or r@Q@h>=-1e-8 or not np.all(gate(ratios)):continue
                            key=tuple(np.round(r/np.linalg.norm(r),9))
                            if key in seen:continue
                            seen.add(key);caps.append({'short_support':S,'anchor':a,'label':m,'branch':sign,'parameter':x,'normal':r.tolist(),'norm':float(r@Q@r),'target_pairing':float(r@Q@w),'center_pairing':float(r@Q@h),'constant_indices':np.flatnonzero(fixed).tolist(),'ratios':ratios.tolist()})
        print({'omitted':omitted,'profiles':count['profiles'],'families':len(families),'caps':len(caps)},flush=True)
    out={'status':'NUMERICAL_FREE_DOTTED_ANCHOR_CONIC_SCOUT','source':BASE,'source_sha256':sha(BASE),'target':TARGET,'target_sha256':sha(TARGET),'active_short_walls':active,'labels':LABELS,'counts':dict(count),'families':families,'caps':caps,
         'scope':'All18-short supports dropping one current active short wall, old82 single anchors whose19-wall span isPD, labels2/3/4/6, two nonnull norm-conic branches. Remaining variable pairings strictly dotted; numerical polynomial interval intersections cover open intervals in each tested profile. Null residuals omitted. No exact or unrestricted nonexistence assertion.'}
    (ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n');print(out['counts']);return out
if __name__=='__main__':run()
