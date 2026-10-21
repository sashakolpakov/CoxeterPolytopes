#!/usr/bin/env python3
"""Separate binary-program check of the saved 1227-root numerical core.

No imports from the incidence producer or either clique-search implementation.
This independently checks the remaining core, not the preceding 50314-root
peeling calculation, and does not convert floating-point data to exact geometry.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import Bounds,LinearConstraint,milp
from scipy.sparse import coo_matrix
from fresh_h20_quad_adjacent_h3_closure_20261016 import geometry158

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_incidence_core_milp_review_20261016'

def run():
    names=['fresh_h20_two_anchor_families_20261016',
           'fresh_h20_two_anchor_families_20261016_pairs',
           'fresh_h20_nonright151_h3_closure_20261016',
           'fresh_h20_two_anchor_families_20261016_coordinated_angles',
           'fresh_h20_intrinsic_facet_incidence_20261016',
           'fresh_h20_continuous_sample_joint_cover_20261016',
           'fresh_h20_incidence_core_joint_cover_20261016']
    paths=[f'results/{x}.json'for x in names]
    d,p,e,c,screen,cover,answer=[json.loads((ROOT/x).read_text())for x in paths]
    raw=[];labels=[]
    for r in d['candidates']:
        raw.append(r['orthonormal_chart']);labels.append(d['families'][r['family']]['labels'][0])
    for r in p['strong_pairs']:
        raw.append(r['partner_orthonormal_chart']);labels.append(d['families'][r['family']]['labels'][0])
    for r in e['candidates']:raw.append(r['orthonormal_chart']);labels.append(r['labels'][0])
    for r in c['successes']:
        for f,x in zip(r['family_ids'],r['roots']):raw.append(x);labels.append(d['families'][f]['labels'][0])
    assert len(raw)==screen['candidate_count']==50314
    keep=screen['surviving_degree_core'];R=np.array(raw)[keep];labs=np.array(labels)[keep];n=len(keep)
    G=geometry158();M=G['metric'];X=G['X'].astype(np.longdouble);u=X[151];qu=u@M@u
    h=-2*np.cos(np.pi/labs);h[labs==2]=0;qr=2-h*h/2
    oldh=X@M@u;oldq=np.sum((X@M)*X,axis=1)-oldh*oldh/qu
    oldids=np.flatnonzero(oldq>1e-7);oldip=R.astype(np.longdouble)@M@X[oldids].T-h[:,None]*oldh[oldids]/qu
    olddegree=(qr[:,None]*oldq[oldids]-oldip**2>1e-6).sum(axis=1)
    cosine=-(R@M@R.T)/2
    ii,jj=np.triu_indices(n,1);v=cosine[ii,jj]
    angles=np.pi/np.arccos(np.clip(v,-1+1e-14,1-1e-14))
    valid=(abs(v)<2e-6)|(v>=1-2e-6)|((v>0)&(v<1-2e-6)&(angles>=2.9999)&(abs(angles-np.rint(angles))<2e-4))
    intrinsic=qr[ii]*qr[jj]-(2*v+h[ii]*h[jj]/2)**2
    edge=valid&(intrinsic>1e-6)
    ei,ej=ii[edge],jj[edge];ci,cj=ii[~valid],jj[~valid]
    Z=np.array([s['z']for s in cover['samples']]);ar=-R[:,1:]/R[:,0,None]
    cuts=Z@ar.T>1+1e-7
    rr=[];cc=[];vv=[];lower=[];upper=[]
    for row in cuts:
        i=len(lower);ids=np.flatnonzero(row);rr.extend([i]*len(ids));cc.extend(ids.tolist());vv.extend([1.]*len(ids));lower.append(1.);upper.append(np.inf)
    start=len(lower)
    rr.extend(np.repeat(np.arange(start,start+len(ci)),2).tolist());cc.extend(np.c_[ci,cj].ravel().tolist());vv.extend([1.]*(2*len(ci)))
    lower.extend([-np.inf]*len(ci));upper.extend([1.]*len(ci))
    adj=[[]for _ in range(n)]
    for i,j in zip(ei,ej):adj[i].append(int(j));adj[j].append(int(i))
    for v in range(n):
        if olddegree[v]>=3:continue
        i=len(lower);rr.append(i);cc.append(v);vv.append(float(3-olddegree[v]))
        rr.extend([i]*len(adj[v]));cc.extend(adj[v]);vv.extend([-1.]*len(adj[v]));lower.append(-np.inf);upper.append(0.)
    A=coo_matrix((vv,(rr,cc)),shape=(len(lower),n)).tocsc()
    result=milp(np.zeros(n),integrality=np.ones(n),bounds=Bounds(np.zeros(n),np.ones(n)),
                constraints=LinearConstraint(A,np.array(lower),np.array(upper)),
                options=dict(time_limit=300.,mip_rel_gap=0.))
    assert result.status==2,(result.status,result.message)
    assert answer['termination']=='finite_bank_compact_degree_cover_infeasible'
    assert answer['steps'][-1]['search']['complete']
    src=[f'scripts/{TAG}.py','scripts/fresh_h20_quad_adjacent_h3_closure_20261016.py',*paths]
    return dict(status='PASS_INDEPENDENT_NUMERICAL_CORE_MILP_INFEASIBILITY',
                sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},
                variables=n,ideal_sample_constraints=len(Z),incompatible_pairs=len(ci),
                finite_intrinsic_edges=len(ei),degree_constraints=int(sum(olddegree<3)),
                rows=len(lower),nonzeros=len(vv),solver_status=int(result.status),
                scope='Independent reconstruction and binary feasibility check of saved1227-root core only. Does not independently certify prior all-bank peeling or exact geometric angle/incidence signs. Solver infeasibility agrees with complete52-node finite-graph search.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['variables'],d['rows'],flush=True)
