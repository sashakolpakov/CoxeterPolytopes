#!/usr/bin/env python3
"""Bounded numerical proposal bank for joint Makarov H5+H4 completion.

Every recipe is five exact old-wall finite pairings, a transverse finite
profile, and the norm-two equation. Floating point predicates scout only.
No exact compatibility, exhaustion or finite-volume claim is made here.
"""
from pathlib import Path
from fractions import Fraction
from itertools import combinations,product
import argparse,hashlib,json,time
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/fresh_h9_makarov_face_extension_20261013.json'
OUT=ROOT/'results/fresh_h9_makarov_joint_20261013.json'

def read():
    data=json.loads(SOURCE.read_text())
    def f(a):return float(Fraction(a[0]))+np.sqrt(5)*float(Fraction(a[1]))
    Q=np.array([[f(a)for a in row]for row in data['form']])
    R=np.array([[f(a)for a in row]for row in data['seed_normals']])
    center=np.array([f(a)for a in data['strict_center']])
    return data,Q,R,center

def build_bank(max_nonzero=2,mixed=False):
    start=time.time();data,Q,R,center=read();G=Q[:6,:6];H=Q[6:,6:]
    W=R[:7,:6];J=np.linalg.inv(H);WN=np.einsum('ij,jk,ik->i',W,G,W)
    labels=[3,4,5,6,10];value={m:2*np.cos(np.pi/m)for m in labels}
    accepted_labels=[2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]
    finite=np.array([-2*np.cos(np.pi/m)for m in accepted_labels])
    profiles=[([2]*5,np.zeros(5))]
    for n in range(1,max_nonzero+1):
        for ids in combinations(range(5),n):
            for ms in product(labels,repeat=n):
                a=np.zeros(5);labs=[2]*5
                for i,m in zip(ids,ms):a[i]=-value[m];labs[i]=m
                profiles.append((labs,a))
    transverse=[]
    for j in range(4):
        for m in labels:
            cov=np.zeros(4);cov[j]=-value[m]
            transverse.append(([j],[m],cov,cov@J))
    if mixed:
        for i,j in combinations(range(4),2):
            for a,b in product([3,4,5],repeat=2):
                cov=np.zeros(4);cov[i]=-value[a];cov[j]=-value[b]
                transverse.append(([i,j],[a,b],cov,cov@J))
    bank=[];keys={};counts={'rank5_supports':0,'quadratic_roots':0,'acute':0,'old_angle_pass':0}
    for ids in combinations(range(7),5):
        A=W[list(ids)]@G
        if np.linalg.matrix_rank(A,tol=1e-9)!=5:continue
        counts['rank5_supports']+=1
        P=A.T@np.linalg.inv(A@A.T);z=null_space(A)[:,0];za=float(z@G@z)
        if abs(za)<1e-9:continue
        rhs=np.array([p[1]*np.sqrt(WN[list(ids)]/2)for p in profiles])
        old0=rhs@P.T
        qb=2*old0@G@z;qc=np.einsum('ij,jk,ik->i',old0,G,old0)
        for ti,(tids,tlabels,tcov,y)in enumerate(transverse):
            residual=2-y@H@y;disc=qb*qb-4*za*(qc-residual)
            valid=np.flatnonzero(disc>1e-9)
            for sign in [-1,1]:
                roots=old0[valid]+((-qb[valid]+sign*np.sqrt(disc[valid]))/(2*za))[:,None]*z
                full=np.hstack([roots,np.repeat(y[None,:],len(valid),axis=0)])
                counts['quadratic_roots']+=len(full)
                ps=roots@G@W.T
                mask=(ps.max(axis=1)<1e-7)&(full@Q@center<0)
                full=full[mask];ps=ps[mask];vis=valid[mask]
                counts['acute']+=len(full)
                normalized=ps/np.sqrt(WN[None,:]/2)
                gap=np.min(np.abs(normalized[:,:,None]-finite[None,None,:]),axis=2)
                mask=((normalized<-2-1e-7)|(gap<2e-7)).all(axis=1)
                counts['old_angle_pass']+=int(mask.sum())
                for r,p,pi in zip(full[mask],ps[mask],vis[mask]):
                    assert abs(r@Q@r-2)<1e-6
                    key=tuple(np.rint(r*1e7).astype('int64'))
                    if key in keys:continue
                    keys[key]=len(bank)
                    bank.append({'id':len(bank),'old_support':list(ids),'old_labels':profiles[pi][0],
                                 'transverse_support':tids,'transverse_labels':tlabels,
                                 'quadratic_sign':sign,'normal':r.tolist(),
                                 'old_pairings':p.tolist()})
    out={'status':'NUMERICAL_RECIPE_BANK_ONLY','source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
         'max_nonzero_old_contacts':max_nonzero,'mixed_transverse':mixed,
         'proposal_finite_labels':labels,'numerical_gate_labels':accepted_labels,
         'counts':counts,'distinct_rounded_entries':len(bank),'caps':bank,
         'elapsed_seconds':time.time()-start}
    OUT.write_text(json.dumps(out,indent=2)+'\n')
    print('bank',len(bank),counts,'seconds',time.time()-start,flush=True)
    return out

def volume_scout(bank,objectives=200):
    data,Q,R,c=read();caps=np.array([row['normal']for row in bank['caps']])
    R=np.vstack([R,caps]);A=R@Q;eq=c@Q
    rng=np.random.default_rng(906131);x0=c/(-c@Q@c)
    best=(-1e100,None);rows=[];unbounded=0
    for i in range(objectives):
        direction=rng.normal(size=10) if i%4 else Q@(best[1]if best[1]is not None else rng.normal(size=10))
        res=linprog(-direction,A_ub=A,b_ub=np.zeros(len(R)),A_eq=[eq],b_eq=[-1],
                    bounds=[(None,None)]*10,method='highs')
        if not res.success:
            unbounded+=int(res.status==3);continue
        q=float(res.x@Q@res.x)
        if q>best[0]:best=(q,res.x)
        if q>1e-8:rows.append({'norm':q,'point':res.x.tolist(),'max_pair':float(max(A@res.x))})
        if i%20==0:print('volume',i,'best',best[0],flush=True)
    out={'status':'NUMERICAL_WHOLE_BANK_VOLUME_SCOUT','entries':len(caps),'objectives':objectives,
         'strict_center_max_pair':float(max(A@x0)),'unbounded_LP_count':unbounded,
         'best_norm':best[0],'best_point':None if best[1]is None else best[1].tolist(),
         'positive_proposals':rows[:20],
         'scope':'Positive numerical proposals require exactification. No positive proposal is not a finite-volume certificate.'}
    p=OUT.with_name(OUT.stem+'_volume.json');p.write_text(json.dumps(out,indent=2)+'\n')
    print('WHOLE BANK',best[0],flush=True);return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--max-nonzero',type=int,default=2);p.add_argument('--mixed',action='store_true')
    p.add_argument('--volume-only',action='store_true');p.add_argument('--objectives',type=int,default=200);a=p.parse_args()
    bank=json.loads(OUT.read_text())if a.volume_only else build_bank(a.max_nonzero,a.mixed)
    volume_scout(bank,a.objectives)
