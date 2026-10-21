#!/usr/bin/env python3
"""Exact diagnostic for three source-family rescalings of the retained rule.

Only inherited zero-corner quadrilaterals and one marked forest are tested.
This is not a new weighting rule or a complete face/vertex proof.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from scipy.optimize import linprog
from sympy import Matrix,Rational
import fresh_sqrt5_dimension22_h3_cross_20261016_review as H

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_marked_old_bonus_scaling_20261016'
R=H.R


def run():
    rule=H.Rule(H.read(H.PRIMARY+'_cap2'))
    types={kind:{}for kind in H.KINDS}
    for name,g in R.MODELS.items():
        for u in g:
            for kind in H.KINDS[1:]:types[kind][R.tree_key(g,[u])]=name
        for u in g:
            for v in g:
                if u<v:types['connected'][R.tree_key(g,[u,v])]=name
    selected=[('hubcross','A5'),('a1cross','I2(5)'),('connected','A5')]
    old=H.read('fresh_sqrt5_parameter_weight_stress_20261015');queue=H.read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    rows=[('prior',i,r['row'])for i,r in enumerate(old['faces'])if r['row']['nodes']<=23]
    rows += [('finite_queue',r['index'],r['row'])for r in queue['records']if r['result']['singular']and r['row']['nodes']<=23]
    constraints={};cases=0
    for family,index,row in rows:
        corners=H.corners(row);assert all(c['baseline']==0 for c in corners)
        for ea,eb,ec in H.contexts(23-row['nodes']):
            p=[F(0)]*4
            for corner in corners:
                a,b,c=corner['a']+ea,corner['b']+eb,corner['c']+ec
                for kind,mult in [('connected',1),('a1cross',a),('hubcross',2*b),('h3cross',3*c)]:
                    key=corner[kind]
                    if key not in rule.orbits[kind]:continue
                    assert mult>0
                    v=rule.bonus(kind,key,a+1,b,c)/mult
                    group=(kind,types[kind][key])
                    p[1+selected.index(group)if group in selected else 0]+=v
            assert sum(p)>=1
            constraints.setdefault(tuple(p),dict(family=family,index=index,external_counts=[ea,eb,ec]))
            cases+=1
    # Literal coefficients of marked7A1+2I25+2A5, independently counted.
    objective=[]
    for kind,name in selected:
        g=R.MODELS[name]
        if kind=='connected':
            value=sum(rule.bonus(kind,R.tree_key(g,[u,v]),8,2,0)for u in g for v in g if u<v)
        else:value=sum(rule.bonus(kind,R.tree_key(g,[u]),8,2,0)for u in g)
        if kind=='a1cross':value*=F(6,7)
        objective.append(2*value)
    assert F(16)+sum(objective)==F(1516843,35000)
    profiles=list(constraints)
    result=linprog(list(map(float,objective)),A_ub=[[-float(v)for v in p[1:]]for p in profiles],
                   b_ub=[float(p[0]-1)for p in profiles],bounds=[(0,None)]*3,method='highs')
    assert result.success
    active=[];rhs=[]
    for p in profiles:
        if abs(float(p[0])+sum(float(v)*x for v,x in zip(p[1:],result.x))-1)<1e-7:
            active.append(list(p[1:]));rhs.append(1-p[0])
    for j,x in enumerate(result.x):
        if abs(x)<1e-7:
            p=[F(0)]*3;p[j]=1;active.append(p);rhs.append(F(0))
    rat=lambda v:Rational(v.numerator,v.denominator)if isinstance(v,F)else Rational(v)
    A=Matrix([[rat(v)for v in row]for row in active]);pivots=list(A.T.rref()[1]);assert len(pivots)==3
    sol=A[pivots,:].inv()*Matrix([rat(rhs[i])for i in pivots]);x=[F(int(v.p),int(v.q))for v in sol]
    assert min(x)>=0 and all(p[0]+sum(v*y for v,y in zip(p[1:],x))>=1 for p in profiles)
    support=[i for i,v in enumerate(result.ineqlin.marginals)if v<-1e-8]
    positive=[j for j,y in enumerate(x)if y]
    da=Matrix([[rat(profiles[i][j+1])for i in support]for j in positive]);db=Matrix([rat(objective[j])for j in positive])
    ds=da.gauss_jordan_solve(db);assert not ds[1]
    dual=[F(int(v.p),int(v.q))for v in ds[0]];assert min(dual)>=0
    assert all(sum(mu*profiles[i][j+1]for i,mu in zip(support,dual))<=objective[j]for j in range(3))
    lower=sum(mu*(1-profiles[i][0])for i,mu in zip(support,dual))
    value=sum(a*b for a,b in zip(objective,x));assert value==lower
    paths=[f'scripts/{TAG}.py',f'scripts/{H.TAG}.py',f'results/{H.TAG}.json',
        f'results/{H.PRIMARY}_cap2.json','results/fresh_sqrt5_parameter_weight_stress_20261015.json',
        'results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json']
    return dict(status='PASS_EXACT_THREE_SOURCE_FAMILY_RESCALING_DIAGNOSTIC',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},variable_order=selected,
        rows=len(rows),common_context_cases=cases,distinct_profiles=len(profiles),
        objective_coefficients=list(map(str,objective)),constant='16',exact_primal=list(map(str,x)),
        exact_minimum_forest_cost=str(16+value),target='35',target_met=16+value<=35,
        exact_dual=[dict(multiplier=str(mu),profile=list(map(str,profiles[i])),witness=constraints[profiles[i]])for i,mu in zip(support,dual)],
        all_quad_profiles=[dict(profile=list(map(str,p)),witness=constraints[p])for p in profiles],
        scope='Only three complete source-family multipliers of the frozen shifted rule are varied. Other old bonuses stay fixed. The exact LP minimizes one marked forest subject to inherited zero-corner quadrilaterals. It does not establish a complete compensated rule or any dimension bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    result=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['exact_minimum_forest_cost'],result['target_met'],flush=True)
