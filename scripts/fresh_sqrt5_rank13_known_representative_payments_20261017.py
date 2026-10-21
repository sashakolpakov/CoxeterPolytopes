#!/usr/bin/env python3
"""Bind every known literal simplex representative to the58 saved exact payment rules."""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_sqrt5_rank13_discriminant_simplex_lp_20261017 as L
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_known_representative_payments_20261017'
INPUT='fresh_sqrt5_rank13_triangle_all_representatives_20261017'

def run():
    reps=B.read(INPUT);lp=B.read(L.TAG);B.dual_bound=D.sparse_dual_bound
    assert reps['counts']['new_stars']==0 and reps['counts']['all_unique_stars']==1408
    assert lp['completed_cases']==lp['cases_with_budget_at_most12']==58
    finite={B.C.code(g):name for name,g in B.C.models() if len(g)<=13}
    variable_lookup={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index'] for v in lp['variables']}
    face_lookup={(f['required'],tuple(map(tuple,f['coefficients']))):f['index'] for f in lp['faces']}
    stars=[]
    for star in reps['stars']:
        g={i:{} for i in range(14)}
        for i,j,m in star['edges']:g[i][j]=g[j][i]=m
        S=[s['omitted'] for s in star['spherical_deletions']];coeff=Counter();corners=[]
        for removed in S:
            end={i:{j:m for j,m in row.items() if j!=removed} for i,row in g.items() if i!=removed}
            pieces=B.C.pieces(end);forest=tuple(sorted(finite[B.C.code(p)] for p in pieces));vv=[]
            for pair in combinations([i for i in S if i!=removed],2):
                code=B.dump(sorted(B.C.code(p,pair) for p in pieces));index=variable_lookup[forest,code]
                coeff[index]+=1;vv.append(dict(pair=list(pair),variable=index))
            corners.append(dict(removed=removed,forest=list(forest),pairs=vv))
        required=2*comb(len(S),3);key=(required,tuple(sorted(coeff.items())));face=face_lookup[key]
        assert sum(coeff.values())==3*comb(len(S),3)
        stars.append(dict(star=star['index'],aggregate_face=face,required=required,coefficients=[list(x) for x in sorted(coeff.items())],
            simplex_dimension=len(S)-1,corners=corners,fixed_discriminant=star['fixed_discriminant'],
            discriminant_class=star.get('discriminant_class')))
    assert len({s['aggregate_face'] for s in stars})==len(stars)==len(lp['faces'])
    cases=[]
    for result in lp['results']:
        ci=result['case_index'];case=lp['cases'][ci];model=L.submodel(lp,case)
        primal=next(c for c in result['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
        B.verify_certificate(model,primal);assert F(primal['maximum'])<12
        weights={i:F(a) for i,a in primal['weights']};payments=[]
        for star in stars:
            compatible=not star['fixed_discriminant'] or star['discriminant_class']==case['discriminant_class']
            assert compatible==(star['aggregate_face'] in case['included_faces'])
            if not compatible:continue
            payment=sum((m*weights.get(i,F(0)) for i,m in star['coefficients']),F(0))
            assert payment>=star['required']
            payments.append(dict(star=star['star'],aggregate_face=star['aggregate_face'],payment=str(payment),required=star['required']))
        cases.append(dict(case_index=ci,discriminant_class=case['discriminant_class'],maximum_vertex_budget=primal['maximum'],
            paid_stars=len(payments),payments=payments))
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'scripts/{L.TAG}.py',f'results/{INPUT}.json',f'results/{L.TAG}.json']
    return dict(status='PASS_EXACT_ALL_KNOWN_REPRESENTATIVE_CLASS_PAYMENTS',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(literal_representatives=len(reps['occurrences']),unique_simplex_stars=len(stars),class_rules=len(cases),
            exact_class_star_payments=sum(c['paid_stars'] for c in cases),whole_forests=len(lp['vertices'])),
        worst_vertex_budget=str(max(F(c['maximum_vertex_budget']) for c in cases)),stars=stars,cases=cases,
        scope='This closes all46032 saved literal provenance records for the3649 retained triangle rows, via exact weighted graph maps and complete simplex-corner sums. It does not expand compressed old-ADE occurrences, include other triangle rows or other simplex stars, cover quadrilaterals, or supply budgets/coverage for absent finite forests. No dimension theorem follows.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['worst_vertex_budget'],flush=True)
