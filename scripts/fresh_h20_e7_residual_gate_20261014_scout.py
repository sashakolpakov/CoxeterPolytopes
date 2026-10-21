#!/usr/bin/env python3
"""Numerical SOC upper bounds for fixed E7 contacts, all intrinsic residuals.

The timelike norm constraint is relaxed inward, so a positive optimum is
only a proposal, never a verified Coxeter wall or a completeness certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import cvxpy as cp
import numpy as np
from decimal import Decimal, localcontext
import fresh_h20_e7_residual_gate_20261014_exact_data as X

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_residual_gate_20261014_scout'
BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
PARENT='results/fresh_h20_e7_intrinsic_caps_20261014.json'
FACE='results/fresh_h20_e7_face_20261014.json'


def run():
    exact=X.load();S=exact['S'];A=exact['A']
    I=np.array(exact['I'],dtype=float);J=np.array(exact['J'],dtype=float)
    coords=np.array([[X.value(x)for x in row]for row in exact['coords']])
    times=np.array(list(map(X.value,exact['times'])))
    wp=np.array([[X.value(x)for x in row]for row in exact['wp']])
    a=np.array(list(map(X.value,exact['a'])));B=X.value(exact['B'])
    weightE=np.array(list(map(X.value,exact['weightE'])))
    assert np.count_nonzero(coords[:,S])==0 and np.count_nonzero(times[S+A])==0
    L=np.linalg.cholesky(J).T
    d=cp.Variable(13,nonneg=True);t=cp.Variable(nonneg=True)
    contact=cp.Parameter(7,nonneg=True);mass=cp.Parameter(nonneg=True)
    acute=-coords.T@d+times*t-wp.T@contact<=0
    norm=cp.norm(cp.hstack([L@d,mass]),2)<=np.sqrt(2)*t
    objective=a@d-B*t+weightE@contact
    prob=cp.Problem(cp.Maximize(objective),[acute,norm])
    profiles=[];skipped=[]
    for size in [1,2,3]:
        for support in combinations(range(7),size):
            labels=[3,4,6]if size<=2 else[3]
            for ls in product(labels,repeat=size):
                # Every contact is 1,sqrt2 or sqrt3. Compute its mass in
                # Q(sqrt2,sqrt3) exactly before high-precision evaluation.
                roots={3:1,4:2,6:3};terms={1:F(-2),2:F(0),3:F(0),6:F(0)}
                for i,li in zip(support,ls):
                    for j,lj in zip(support,ls):
                        r=roots[li]*roots[lj];factor=1
                        if r in [4,9]:factor={4:2,9:3}[r];r=1
                        terms[r]+=factor*exact['I'][i][j]
                with localcontext()as ctx:
                    ctx.prec=90
                    mass2=sum(Decimal(v.numerator)/Decimal(v.denominator)*Decimal(r).sqrt()for r,v in terms.items())
                    m2=float(mass2)
                if mass2<0:
                    skipped.append(dict(support=list(support),labels=list(ls),reason='negative_mass_outside_causal_model',mass_square=str(mass2)));continue
                c=np.zeros(7);c[list(support)]=[np.sqrt(roots[x])for x in ls]
                contact.value=c;mass.value=np.sqrt(m2)
                prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9)
                row=dict(support=list(support),walls=[S[i]for i in support],labels=list(ls),
                    status=prob.status,upper_bound=float(prob.value)if np.isfinite(prob.value)else None,mass_square=float(m2),
                    exact_mass_square={str(r):str(v)for r,v in terms.items()})
                if d.value is not None:
                    row.update(face_contacts=d.value.tolist(),time=float(t.value),
                        actual_projected_norm=float(d.value@J@d.value-2*t.value*t.value),
                        primal_max_acute_violation=float(max(-coords.T@d.value+times*t.value-wp.T@c)),
                        primal_soc_violation=float(np.linalg.norm(np.hstack([L@d.value,np.sqrt(m2)]))-np.sqrt(2)*t.value),
                        dual_oldwalls=acute.dual_value.tolist(),dual_norm=float(norm.dual_value))
                profiles.append(row)
                if row['upper_bound']is not None and row['upper_bound']>1e-6:
                    print('NUMERIC_POSITIVE_RELAXATION',support,ls,prob.value,
                          row.get('actual_projected_norm'),flush=True)
    paths=[BASE,PARENT,FACE,f'scripts/{TAG}.py','scripts/fresh_h20_e7_residual_gate_20261014_exact_data.py',
           'scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
    result=dict(status='NUMERICAL_ALL_RESIDUAL_FIXED_CONTACT_SOC_SCOUT_ONLY',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        E7_support=S,face_walls=A,profiles=profiles,skipped=skipped,
        exact_zero_counts=dict(face_wall_pairings=sum(x==X.K.Z for row in exact['coords']for x in row),
                               time_wall_pairings=sum(x==X.K.Z for x in exact['times']),
                               weight_wall_pairings=sum(x==X.K.Z for row in exact['wp']for x in row)),
        coefficient_method='Exact physical eight-field pairings and rational inverses; 90-digit Decimal evaluation only after exact-zero tests.',
        scope='Finite fixed-contact numerical upper bounds; norm inequality relaxes exact norm. No exact exclusion, no Coxeter realization, no finite-volume claim.')
    (ROOT/f'results/{TAG}.json').write_text(json.dumps(result,indent=2)+'\n')
    print('profiles',len(profiles),'positive',sum(v['upper_bound']is not None and v['upper_bound']>1e-6 for v in profiles),flush=True)


if __name__=='__main__':run()
