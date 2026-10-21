#!/usr/bin/env python3
"""Bounded numerical quotient scout; no exact exclusion or accepted cap claim."""
import json, math, argparse
from pathlib import Path
from fractions import Fraction
from hashlib import sha256
import numpy as np
import cvxpy as cp

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_null_link_cap_scout_20261016'
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
    # Fixed incident profile: zero or one contact d at isolated root17.
    d=cp.Parameter(nonneg=True);radius=cp.Parameter(nonneg=True)
    z=cp.Variable(4);height=cp.Variable()
    r=-d*N[17]/2+z@free+height*nu
    lo=cp.Parameter(2);hi=cp.Parameter(2)
    scale=np.maximum(1,np.linalg.norm(N[outside]@Q,axis=1))
    constraints=[cp.norm(z)<=radius,(N[outside]@Q@r)/scale<=0,
                 H@Q@r<=-1e-8,N[[5,87]]@Q@r>=lo,N[[5,87]]@Q@r<=hi]
    problem=cp.Problem(cp.Maximize(E@Q@r),constraints)
    labels=list(range(2,13))+['parallel','dotted']
    value=lambda m:0. if m==2 else 2. if m=='parallel' else 2*math.cos(math.pi/m)
    records=[];success=[]
    for own in [2]+list(range(3,13))+['parallel']:
        ownvalue=value(own);d.value=ownvalue;radius.value=math.sqrt(max(0,2-ownvalue**2/2))
        for a in labels:
            for b in labels:
                lo.value=np.array([-1e7 if m=='dotted'else -value(m)for m in [a,b]])
                hi.value=np.array([-2. if m=='dotted'else -value(m)for m in [a,b]])
                try:problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9,max_iter=200)
                except cp.error.SolverError:
                    records.append(dict(own17=own,contact5=a,contact87=b,status='solver_error'));continue
                rec=dict(own17=own,contact5=a,contact87=b,status=problem.status)
                if problem.status in ['optimal','optimal_inaccurate']:
                    rr=np.asarray(r.value);pairs=N@Q@rr
                    rec.update(max_E1=float(E@Q@rr),norm=float(rr@Q@rr),free_norm=float(np.linalg.norm(z.value)),
                               radius=float(radius.value),height=float(height.value),max_acute=float(max(pairs)),
                               center=float(H@Q@rr))
                    if rec['max_E1']>1e-6:
                        rec['normal']=rr.tolist();rec['pairs']=pairs.tolist();success.append(rec)
                records.append(rec)
        print('profile',own,'solves',len(records),'positive relaxations',len(success),flush=True)
    return dict(status='BOUNDED_NUMERICAL_QUOTIENT_CONTACT_RELAXATION',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in PATHS+[f'scripts/{TAG}.py']},
        definitions=dict(incident_white=white,own17_profiles=[2]+list(range(3,13))+['parallel'],
                         nonincident_5_87_choices=labels,all_other_incident_contacts='right',
                         norm_constraint='at most2; equality separately reported',
                         all_other_old_contacts='acute only, not Coxeter'),
        total_solves=len(records),positive_relaxations=len(success),records=records,
        scope='Bounded double-precision SOCP scout only. Infeasibility is not an exact all-label exclusion. Positive rows are not accepted Coxeter caps; other old contacts and exact norm/strict signs remain unverified.')
if __name__=='__main__':
    argparse.ArgumentParser().parse_args();data=run()
    (ROOT/f'results/{TAG}.json').write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['total_solves'],data['positive_relaxations'],flush=True)
