#!/usr/bin/env python3
"""Numerical continuous edge profiles with all nonedge walls strictly dotted.

The scalar x=exp(theta)>1 gives exact-definable norm-two recipes; numerical
interval endpoints and acceptance are proposals only.
"""
from fractions import Fraction
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np
from fresh_h9_makarov_joint_completion_20261013_profiles import setup,ROOT,P,CORE,FRONT
OUT=ROOT/f'results/{P}_continuous.json'
def intervals(a,b,c):
    # a*x^2+b*x+c <0, returning its open intervals.
    if abs(a)<1e-10:
        if abs(b)<1e-10:return[(-np.inf,np.inf)]if c<0 else[]
        r=-c/b;return[(-np.inf,r)]if b>0 else[(r,np.inf)]
    d=b*b-4*a*c
    if d<=0:return[(-np.inf,np.inf)]if a<0 else[]
    roots=sorted([(-b-np.sqrt(d))/(2*a),(-b+np.sqrt(d))/(2*a)])
    return[(roots[0],roots[1])]if a>0 else[(-np.inf,roots[0]),(roots[1],np.inf)]
def intersect(a,b):
    return[(max(x,z),min(y,w))for x,y in a for z,w in b if min(y,w)-max(x,z)>1e-8]
def build():
    start=time.time();src,Q,W,c,E=setup();front=json.loads(FRONT.read_text())
    f=lambda x:float(Fraction(x[0]))+np.sqrt(5)*float(Fraction(x[1]))
    ds=np.array([f(x)for x in front['square_root_basis_over_Qsqrt5']])
    def evalue(x):
        return sum(f(a)*np.prod([np.sqrt(ds[j])for j in range(len(ds))if mask>>j&1])for mask,a in x)
    vertices=np.array([[evalue(x)for x in row['point']]for row in front['finite_vertices']])
    target0=np.array(json.loads((ROOT/f'results/{P}_all_relaxed_volume.json').read_text())['best_point'])
    target0=.999*target0+.001*c/(-c@Q@c);assert target0@Q@target0>0
    target2=np.array(src['escape']['point'],dtype=float);target2/=-(target2@Q@c)
    targets=np.array([target0,E,target2]);regions=[];caps=[];seen=set();stats={'profiles':0,'feasible_profiles':0,'regions':0}
    for eid,packet in front['open_edge_spherical_profiles'].items():
        edge=front['edges'][int(eid)];support=edge['support'];vertex=edge['finite_vertices'][0]
        v=vertices[vertex];pv=v/np.sqrt(-(v@Q@v));G=W[support]@Q@W[support].T;Gi=np.linalg.inv(G)
        released=edge['released_at_start'];w=W[released];b=w-(W[support]@Q@w)@Gi@W[support]
        b/=np.sqrt(b@Q@b);assert abs(b@Q@pv)<1e-6
        outside=sorted(set(range(27))-set(support));AA=pv@Q@W.T;BB=b@Q@W.T
        profiles=[r['contacts']for r in packet['finite_profiles']]
        profiles += [[[r['wall'],m]]for r in packet['dihedral_tail_profiles']for m in[3,4,5,6,8,10,12,16,24]]
        for pi,contacts in enumerate(profiles):
            stats['profiles']+=1;p=np.zeros(8);lookup=dict(contacts)
            for i,wid in enumerate(support):
                if wid in lookup:p[i]=-2*np.cos(np.pi/lookup[wid])
            u=p@Gi@W[support];res=2-p@Gi@p
            if res<=1e-8:continue
            q=np.sqrt(res);dd=u@Q@W.T
            allowed=[(1.,np.inf)]
            for j in outside:
                allowed=intersect(allowed,intervals(q*(AA[j]-BB[j])/2,dd[j]+2,-q*(AA[j]+BB[j])/2))
                if not allowed:break
            if not allowed:continue
            stats['feasible_profiles']+=1
            for lo,hi in allowed:
                stats['regions']+=1;proposals=[lo+max(1e-5,abs(lo)*1e-5)]
                if np.isfinite(hi):proposals +=[(lo+hi)/2,hi-max(1e-5,abs(hi)*1e-5)]
                else:proposals +=[max(lo+1,lo*1.2),max(lo+2,lo*2)]
                for target in targets:
                    A=pv@Q@target;B=b@Q@target
                    if abs(A-B)>1e-12 and -(A+B)/(A-B)>0:proposals.append(np.sqrt(-(A+B)/(A-B)))
                region={'edge':int(eid),'support':support,'vertex':vertex,'released':released,'contacts':contacts,'residual_norm':float(res),'interval':[float(lo),None if not np.isfinite(hi)else float(hi)]}
                regions.append(region)
                for x in proposals:
                    if not lo<x<hi or x>1e6:continue
                    xrat=Fraction(float(x)).limit_denominator(1000000);x=float(xrat)
                    if not lo<x<hi:continue
                    r=u+q*((x-1/x)/2*pv-(x+1/x)/2*b)
                    ps=W@Q@r
                    if max(ps[outside])>=-2-1e-8 or r@Q@c>=0:continue
                    assert abs(r@Q@r-2)<1e-3
                    key=tuple(np.rint(r*1e6).astype('int64'))
                    if key in seen:continue
                    seen.add(key);cuts=r@Q@targets.T
                    caps.append({'id':len(caps),'region':len(regions)-1,'boost_exp_rational':str(xrat),'normal':r.tolist(),'target_pairings':cuts.tolist(),'all27_pairings':ps.tolist()})
    out={'status':'NUMERICAL_ALL_DOTTED_CONTINUOUS_PROFILE_REGIONS','source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[CORE,FRONT]},
         'targets':targets.tolist(),'regions':regions,'caps':caps,'counts':stats,'target_cut_counts':(np.array([r['target_pairings']for r in caps])>1e-8).sum(axis=0).tolist()if caps else[0]*3,'seconds':time.time()-start,
         'scope':'Numerical quadratic interval intersections and rational boost proposals. Finite profiles are complete input; isolated tails are sampled. All nonedge relations are required strictly dotted; no exact interval or global-volume claim.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('CONTINUOUS',stats,'caps',len(caps),'cuts',out['target_cut_counts'],'seconds',out['seconds'],flush=True)
if __name__=='__main__':build()
