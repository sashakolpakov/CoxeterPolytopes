#!/usr/bin/env python3
"""Bounded numerical quotient scout; no exact exclusion or accepted cap claim."""
import json, math, argparse
from pathlib import Path
from fractions import Fraction
from hashlib import sha256
import numpy as np
import cvxpy as cp

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_null_link_complete_scout_20261016'
PATHS=['results/fresh_h20_dotted_arc_cap_20261016_scout.json',
       'results/fresh_h20_dotted_arc_cap_20261016.json',
       'results/fresh_h20_complementary_cap_20261016.json',
       'results/fresh_h20_complementary_frontier_20261016.json',
       'results/fresh_h20_affine_a17_completion_20261016.json']
def mid(a):return float(sum(map(Fraction,a))/2)
def run():
    scout,arc,cap,front,completion=[json.loads((ROOT/p).read_text())for p in PATHS]
    b=scout['inputs'];P=np.array(b['P'],float);W=np.array(b['W'],float)
    J=np.array(b['J'],float);S=b['support'];G=np.linalg.inv(J)
    Q=np.zeros((21,21));Q[:3,:3]=np.diag([-1,1,1]);Q[3:,3:]=J
    N=np.c_[-P[:,0],P[:,1:],(G@W).T]
    for record,idx,key in [(arc,31,'dual_coefficient_bounds'),(cap,5,None)]:
        y=np.array([mid(x)for x in record['definition']['frame_bounds']])
        d=np.zeros(18);d[S.index(idx)]=mid(record['definition'][key])if key else 1.
        N=np.r_[N,[np.r_[y,-d]]]
    nu=np.sqrt(3)*N[134]+2*N[0]+N[1]
    E=np.array(list(map(Fraction,cap['strict_escape']['frame_coefficients']+cap['strict_escape']['dual_coefficients'])),float)
    hp=np.array([mid(x)for x in cap['common_center']['all135_pair_bounds']]);hs=hp[S]
    hf=np.linalg.lstsq(P,hp[:133]-W.T@hs,rcond=None)[0];H=np.r_[hf,hs]
    white=[3,4,6,7,8,9,10,12,13,14,15,17,31]
    U=N[white+[0,1]];JU=np.linalg.inv(U@Q@U.T)
    orth=lambda v:v-(v@Q@U.T)@JU@U
    hbar=orth(H);k=hbar/(-nu@Q@H);eta=k+(k@Q@k)*nu/2
    free=[]
    for e in np.eye(21):
        v=orth(e);v+=(v@Q@eta)*nu+(v@Q@nu)*eta
        for w in free:v-=(v@Q@w)*w
        q=v@Q@v
        if q>1e-8:free.append(v/np.sqrt(q))
        if len(free)==4:break
    free=np.array(free);assert np.max(abs(free@Q@free.T-np.eye(4)))<1e-8
    incident=set(front['exact_null_frontier']['incident_walls']);outside=[i for i in range(135)if i not in incident]
    profiles=json.loads((ROOT/'results/fresh_h20_null_link_profiles_20261016.json').read_text())
    Jwhite=np.linalg.inv(N[white]@Q@N[white].T)
    allprofiles=profiles['profiles']+profiles['continuous_tails']
    source='results/fresh_h20_null_link_profiles_20261016.json'
    PATHS.append(source)
    p=cp.Parameter(21);radius=cp.Parameter(nonneg=True);axis=cp.Parameter(21)
    d=cp.Variable();z=cp.Variable(4);height=cp.Variable();dlo=cp.Parameter();dhi=cp.Parameter()
    r=p+d*axis+z@free+height*nu
    lo=cp.Parameter(4);hi=cp.Parameter(4)
    scale=np.maximum(1,np.linalg.norm(N[outside]@Q,axis=1))
    constraints=[cp.norm(cp.hstack([z,d/math.sqrt(2)]))<=radius,d>=dlo,d<=dhi,
                 (N[outside]@Q@r)/scale<=0,H@Q@r<=0,
                 N[[2,5,61,87]]@Q@r>=lo,N[[2,5,61,87]]@Q@r<=hi]
    problem=cp.Problem(cp.Maximize(E@Q@r),constraints)
    records=[];positive=[]
    for idx,pr in enumerate(allprofiles):
        tail='wall'in pr
        if tail:
            p.value=np.zeros(21);axis.value=-N[pr['wall']]/2;radius.value=math.sqrt(2)
            dlo.value=math.sqrt(3);dhi.value=2.
        else:
            contacts=np.zeros(13)
            for i,m in zip(pr['walls'],pr['labels']):contacts[white.index(i)]=2*math.cos(math.pi/m)
            cost=contacts@Jwhite@contacts
            p.value=-(contacts@Jwhite)@N[white];axis.value=np.zeros(21)
            radius.value=math.sqrt(max(0,2-cost));dlo.value=dhi.value=0.
        for mask in range(16):
            active=[bool(mask&(1<<j))for j in range(4)]
            lo.value=np.array([-1e7 if yes else 0. for yes in active])
            hi.value=np.array([-1. if yes else 0. for yes in active])
            try:problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=300)
            except cp.error.SolverError:
                records.append(dict(profile=idx,mask=mask,status='solver_error'));continue
            row=dict(profile=idx,mask=mask,status=problem.status)
            if problem.status in ['optimal','optimal_inaccurate']:
                rr=np.array(r.value);pairs=N@Q@rr
                row.update(max_E1=float(E@Q@rr),norm=float(rr@Q@rr),center=float(H@Q@rr),
                           max_acute=float(max(pairs)),axis_contact=float(d.value),height=float(height.value))
                if row['max_E1']>1e-6:
                    row['normal']=rr.tolist();row['pairs']=pairs.tolist();positive.append(row)
            records.append(row)
        if idx%10==0:print('profile',idx,'solves',len(records),'positive',len(positive),flush=True)
    return dict(status='NUMERICAL_COMPLETE_INCIDENT_DOMAIN_MINIMUM_CONTACT_RELAXATION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in PATHS+[f'scripts/{TAG}.py']},
        profile_count=len(allprofiles),profile_definitions=allprofiles,gate_walls=[2,5,61,87],
        total_solves=len(records),positive_relaxations=len(positive),records=records,
        scope='All-label incident domain is exact in the separate bound certificate. This is only a double-precision convex relaxation of norm2 and all other Coxeter contacts; no exact exclusion or accepted wall is asserted.')
if __name__=='__main__':
    data=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['total_solves'],data['positive_relaxations'],flush=True)
