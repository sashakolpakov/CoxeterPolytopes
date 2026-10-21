#!/usr/bin/env python3
"""Joint determinant-class and maximal-simplex partial payment experiment."""
import argparse,hashlib,json,time
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_discriminant_simplex_lp_20261017'
INPUT='fresh_sqrt5_rank13_simplex_aggregate_lp_20261017'
GATE='fresh_sqrt5_rank13_triangle_discriminant_gate_20261017'

def domain():
    old=B.read(INPUT);gate=B.read(GATE);records={r['face']:r for r in gate['records']}
    cases=[]
    order=[None]+[r['index'] for r in gate['classes']]
    for cl in order:
        included=[];witnesses=[]
        for f in old['faces']:
            eligible=[p for p in f['provenance'] if not records[p['original_face']]['fixed_discriminant'] or
                records[p['original_face']]['discriminant_class']==cl]
            if eligible:
                included.append(f['index']);witnesses.append(dict(face=f['index'],original_triangle_face=eligible[0]['original_face']))
        cases.append(dict(discriminant_class=cl,included_faces=included,witnesses=witnesses))
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'results/{INPUT}.json',f'results/{GATE}.json']
    return dict(sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        vertices=old['vertices'],variables=old['variables'],faces=old['faces'],cases=cases,
        classes=gate['classes'],counts=dict(vertices=len(old['vertices']),variables=len(old['variables']),classes=len(cases)))

def submodel(d,case):
    faces=[]
    for old in case['included_faces']:
        f=dict(d['faces'][old]);f.update(index=len(faces),original_aggregate_face=old);faces.append(f)
    return dict(vertices=d['vertices'],variables=d['variables'],faces=faces,
        counts=dict(vertices=len(d['vertices']),variables=len(d['variables']),distinct_face_constraints=len(faces)))

def optimize(d,limit):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    N=len(d['variables']);rows=[];cols=[];values=[];rhs=[]
    for f in d['faces']:
        r=len(rhs);rhs.append(-f['required'])
        for j,c in f['coefficients']:rows.append(r);cols.append(j);values.append(-c)
    for v in d['vertices']:
        r=len(rhs);rhs.append(0)
        for j,m in v['variables']:rows.append(r);cols.append(j);values.append(m)
        rows.append(r);cols.append(N);values.append(-1)
    matrix=coo_matrix((values,(rows,cols)),shape=(len(rhs),N+1)).tocsr();objective=np.zeros(N+1);objective[-1]=1
    result=linprog(objective,A_ub=matrix,b_ub=np.array(rhs,dtype=float),bounds=(0,None),method='highs',options={'time_limit':limit})
    solver=dict(status=int(result.status),message=str(result.message),success=bool(result.success),objective=None if result.fun is None else float(result.fun))
    certificates=[]
    if result.success:
        primal=B.recover_primal(d,result.x[:-1])
        if primal is not None:certificates.append(primal)
        certificates.append(D.recover_unminimized(d,result.ineqlin.marginals[:len(d['faces'])]))
    for c in certificates:B.verify_certificate(d,c)
    return solver,certificates

def run(verify):
    d=domain();B.dual_bound=D.sparse_dual_bound;out=ROOT/f'results/{TAG}.json'
    if verify:
        saved=json.loads(out.read_text())
        for k,v in d.items():assert saved[k]==v,k
        for result in saved['results']:
            model=submodel(d,d['cases'][result['case_index']])
            for c in result['certificates']:B.verify_certificate(model,c)
        return saved
    results=[];start=time.monotonic()
    for ci,case in enumerate(d['cases']):
        remaining=300-(time.monotonic()-start)
        if remaining<=0:break
        print('CLASS_SIMPLEX_START',ci,case['discriminant_class'],len(case['included_faces']),flush=True)
        solver,certificates=optimize(submodel(d,case),min(60,remaining))
        results.append(dict(case_index=ci,discriminant_class=case['discriminant_class'],solver=solver,certificates=certificates))
        print('CLASS_SIMPLEX_DONE',ci,solver,[(c['kind'],c.get('maximum',c.get('bound'))) for c in certificates],flush=True)
        if not solver['success']:break
    complete=len(results)==len(d['cases'])
    feasible=sum(any(c.get('reaches_target12',False) for c in r['certificates']) for r in results)
    excluded=sum(any(c.get('excludes_target12',False) for c in r['certificates']) for r in results)
    d.update(status='PASS_EXACT_COMPLETE_PARTIAL_CLASS_SIMPLEX_BATCH' if complete else 'UNRESOLVED_TIME_BOUNDED_PARTIAL_CLASS_SIMPLEX_BATCH',
        results=results,completed_cases=len(results),cases_with_budget_at_most12=feasible,cases_with_exact_lower_bound_above12=excluded,
        scope='Partial domain only: a row is retained when one saved literal simplex representative witnesses the chosen discriminant class or has an odd component. All605 represented whole-forest budgets remain. Other representatives, other simplex stars, quadrilaterals, absent forests, full-form invariants and simultaneous incidence remain outside. No class-specific feasibility theorem or dimension bound follows.')
    out.write_text(json.dumps(d,indent=2)+'\n');return d

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run(a.verify)
    print(d['status'],d['completed_cases'],d['cases_with_budget_at_most12'],d['cases_with_exact_lower_bound_above12'],flush=True)
