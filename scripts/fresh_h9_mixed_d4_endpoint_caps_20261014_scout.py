#!/usr/bin/env python3
"""Numerical general endpoint proposals for the mixed-D4 finite-corner queue.

Exact local profiles are an input, not a claim that this metric scout is exact.
The explicit dihedral tails are sampled. Continuous all-dotted intervals and
single-anchor finite-angle events are tested; other events remain outside scope.
"""
import argparse
import json
import math
from collections import Counter
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_mixed_d4_endpoint_caps_20261014_scout'
FILES=['results/fresh_h9_mixed_d4_coupling_20261014.json',
       'results/fresh_h9_mixed_d4_coupling_20261014_frontier.json',
       'results/fresh_h9_mixed_d4_endpoint_profiles_20261014.json']
ANCHORS=[2,3,4,5,6,7,8,9,10,11,12,16,20,24,32,64]
TAILS=[6,7,8,10,12,16,24,32,64]
mp.mp.dps=100
def base(x):return mp.mpf(x[0])+mp.mpf(x[1])*mp.sqrt(5)
RAD=(16+4*mp.sqrt(5))/11
def scalar(x):return base(x[0])+base(x[1])*mp.sqrt(RAD)
def vec(x):return np.array([float(scalar(t)) for t in x])
def matrix(x):return np.array([vec(r) for r in x])

def data():
    d,f,p=[json.loads((ROOT/x).read_text()) for x in FILES]
    Q=matrix(d['form']);N=matrix(d['final_normals'])
    values,U=np.linalg.eigh(Q);assert values[0]<0 and min(values[1:])>0
    T=np.diag(np.sqrt(abs(values)))@U.T
    h=vec(d['center'])@T.T
    if h[0]<0:T[0]*=-1;h=vec(d['center'])@T.T
    eta=np.diag([-1.]+[1.]*9);assert np.max(abs(T.T@eta@T-Q))<1e-12
    W=N@T.T;norms=np.einsum('ij,jk,ik->i',W,eta,W)
    assert min(norms)>0
    W*=np.sqrt(2/norms)[:,None]
    G=W@eta@W.T;assert max(abs(np.diag(G)-2))<1e-10
    return d,f,p,Q,N,T,eta,h,W,G,vec(d['strict_escape']['vector'])@T.T

def interval_sample(coefs):
    """All roots partition T>1; return one sample per all-negative interval."""
    cuts=[1.]
    for a,b,c in coefs:
        if abs(a)<1e-12:
            if abs(b)>1e-12:
                z=-c/b
                if z>1:cuts.append(z)
        else:
            disc=b*b-4*a*c
            if disc>=0:
                for z in ((-b-math.sqrt(disc))/(2*a),(-b+math.sqrt(disc))/(2*a)):
                    if z>1 and math.isfinite(z):cuts.append(z)
    cuts=sorted(set(cuts));samples=[]
    for lo,hi in zip(cuts,cuts[1:]+[math.inf]):
        if hi-lo<1e-8:continue
        t=(lo+hi)/2 if math.isfinite(hi) else max(lo+1,2*lo)
        if all(a*t*t+b*t+c< -1e-8 for a,b,c in coefs):samples.append((lo,hi,t))
    return samples

def run():
    d,f,p,Q,N,T,eta,h,W,G,escape=data();invT=np.linalg.inv(T)
    fvals=-2*np.cos(np.pi/np.arange(2,257));counts=Counter();found={};edge_stats=[]
    priority=[18,19,20,21,45,46,47,48]
    order=priority+[int(k) for k in p['profiles'] if int(k) not in priority]
    def record(root,recipe,block,start,rpar,rho):
        vals=root@eta@W.T;npair=float(root@eta@start)
        if npair>=-1e-8:return
        residual=root-rpar;end=start-npair/rho*residual
        eq=float(end@eta@end)
        if eq>=-1e-8 or end[0]<=0:return
        ev=end@eta@W.T;outside=[i for i in range(12) if i not in block]
        if max(ev[outside])>=-1e-8:return
        err=np.min(abs(vals[:,None]-fvals[None,:]),axis=1)
        okay=(vals< -2-1e-7)|(err<1e-8)
        if not np.all(okay):return
        assert abs(root@eta@root-2)<1e-6
        key=tuple(np.round(vals,8));complexity=sum(m not in (2,3,5) for _,m in recipe['contacts'])
        complexity+=5*any(i==11 for i,_ in recipe['contacts'])
        complexity+=3*(recipe.get('anchor')==11)+int(recipe.get('anchor_label',2) not in (2,3,5))
        candidate=dict(**recipe,complexity=complexity,normal_physical=root.tolist(),
            normal_original_coordinates=(invT@root).tolist(),normalized_wall_pairings=vals.tolist(),
            endpoint_original_coordinates=(invT@end/end[0]).tolist(),
            cuts_escape=bool(root@eta@escape>1e-8),escape_pairing=float(root@eta@escape),
            start_pairing=npair,endpoint_norm=eq,endpoint_other_pairings=ev[outside].tolist())
        if key not in found or complexity<found[key]['complexity']:found[key]=candidate
    for eid in order:
        row=p['profiles'][str(eid)];block=row['support'];queue=f['queue'][eid]
        start=vec(queue['start_point'])@T.T;other=vec(d['final_complete_ray_census']['rays'][queue['other_cone_ray']]['vector'])@T.T
        U=start/math.sqrt(-start@eta@start);V=other+(other@eta@U)*U;V/=math.sqrt(V@eta@V)
        assert np.max(abs(W[block]@eta@U))<1e-8 and np.max(abs(W[block]@eta@V))<1e-8
        profiles=[r['contacts'] for r in row['finite_profiles']]
        profiles += [[[r['wall'],m]] for r in row['dihedral_tail_profiles'] for m in TAILS]
        B=G[np.ix_(block,block)];Bi=np.linalg.inv(B)
        P=np.zeros((len(profiles),8))
        for k,contacts in enumerate(profiles):
            for i,m in contacts:P[k,block.index(i)]=-2*math.cos(math.pi/m)
        rho=2-np.einsum('ij,jk,ik->i',P,Bi,P);assert min(rho)>0
        pars=P@Bi@W[block];outside=[i for i in range(12) if i not in block]
        before=len(found);all_dotted=0;acute=0
        for k in range(len(profiles)):
            a=math.sqrt(rho[k])/2*((U+V)@eta@W[outside].T)
            b=pars[k]@eta@W[outside].T
            c=math.sqrt(rho[k])/2*((V-U)@eta@W[outside].T)
            if interval_sample(list(zip(a,b,c))):acute+=1
            intervals=interval_sample(list(zip(a,b+2,c)))
            for lo,hi,z in intervals:
                root=pars[k]+math.sqrt(rho[k])*((z-1/z)/2*U+(z+1/z)/2*V)
                record(root,dict(edge_id=eid,contacts=profiles[k],method='all_dotted_interval',
                    exponential_parameter=z,interval=[lo,hi if math.isfinite(hi) else None]),block,start,pars[k],rho[k])
                all_dotted+=1
        vals_anchor=-2*np.cos(np.pi/np.array(ANCHORS));nangles=len(ANCHORS)
        for anchor in outside:
            A=W[block+[anchor]]@eta;uu,ss,vv=np.linalg.svd(A,full_matrices=True)
            if min(ss)<1e-10:counts['rank_deficient_anchor']+=1;continue
            direction=vv[-1];Ai=np.linalg.pinv(A)
            target=np.column_stack([np.repeat(P,nangles,axis=0),np.tile(vals_anchor,len(P))])
            X=target@Ai.T;aa=float(direction@eta@direction);bb=2*(X@eta@direction)
            cc=np.einsum('ij,jk,ik->i',X,eta,X)-2;dd=bb*bb-4*aa*cc
            for sign in (-1,1):
                ids=np.flatnonzero(dd>1e-10)
                if abs(aa)<1e-12:
                    if sign==1:continue
                    ids=ids[abs(bb[ids])>1e-10];ts=-cc[ids]/bb[ids]
                else:ts=(-bb[ids]+sign*np.sqrt(dd[ids]))/(2*aa)
                roots=X[ids]+ts[:,None]*direction;counts['real_anchor_roots']+=len(roots)
                vals=roots@eta@W.T;err=np.min(abs(vals[:,:,None]-fvals[None,None,:]),axis=2)
                good=np.all((vals< -2-1e-7)|(err<1e-8),axis=1)&((roots@eta@start)<-1e-8)
                for j in np.flatnonzero(good):
                    ix=int(ids[j]);k=ix//nangles
                    record(roots[j],dict(edge_id=eid,contacts=profiles[k],method='finite_anchor',anchor=anchor,
                        anchor_label=ANCHORS[ix%nangles],quadratic_branch=sign),block,start,pars[k],rho[k])
        counts['sampled_profiles']+=len(profiles);counts['acute_parameter_profiles']+=acute
        counts['all_dotted_parameter_intervals']+=all_dotted
        edge_stats.append(dict(edge_id=eid,profiles=len(profiles),acute_parameter_profiles=acute,
                               all_dotted_intervals=all_dotted,new_distinct_caps=len(found)-before))
        print(edge_stats[-1],'total',len(found),'escape_cutters',sum(x['cuts_escape'] for x in found.values()),flush=True)
    candidates=sorted(found.values(),key=lambda x:(not x['cuts_escape'],x['complexity'],x['edge_id']))
    files=FILES+[f'scripts/{TAG}.py']
    return dict(status='NUMERICAL_MIXED_D4_GENERAL_ENDPOINT_CAP_PROPOSALS',
        sources={s:sha256((ROOT/s).read_bytes()).hexdigest() for s in files},
        counts=dict(counts),finite_anchor_labels=ANCHORS,sampled_dihedral_tail_labels=TAILS,
        nonsupport_finite_angle_numeric_gate='integer labels 2 through 256, tolerance 1e-8; strict dotted gap 1e-7',
        edges=edge_stats,candidate_count=len(candidates),escape_cutter_count=sum(x['cuts_escape'] for x in candidates),
        candidates=candidates,scope='Numerical proposals only. Complete exact local finite profiles are input, but unbounded dihedral tails are sampled. Every all-dotted open parameter interval is sampled numerically, as are the stated finite-anchor events. No exact metric exclusion, complete all-parameter/all-label search, or finite-volume claim.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.parse_args();d=run()
    (ROOT/f'results/{TAG}.json').write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['candidate_count'],d['escape_cutter_count'],flush=True)
