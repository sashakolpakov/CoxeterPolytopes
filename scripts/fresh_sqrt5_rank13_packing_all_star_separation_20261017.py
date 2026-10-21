#!/usr/bin/env python3
"""Exact nonoptimizing separation of the packing rule on all12420 saved stars."""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as D
import fresh_simplex_port_incidence_20261017 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_rank13_packing_all_star_separation_20261017'
RULE='fresh_sqrt5_rank13_simplex_packing_repair_20261017'
STARS='fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'


def run():
    rule=B.read(RULE);full=B.read(STARS);packing=B.read(P.TAG)
    primal=next(c for c in rule['certificates'] if c['kind']=='EXACT_PARTIAL_MODEL_PRIMAL')
    B.dual_bound=D.sparse_dual_bound;B.verify_certificate(rule,primal)
    weights={i:F(a) for i,a in primal['weights']}
    models=dict(B.C.models());finite={B.C.code(g):name for name,g in models.items() if len(g)<=13}
    types=[(name,g) for name,g in models.items() if len(g)<=13]
    @lru_cache(None)
    def forests(rank,first=0):
        if not rank:return ((),)
        return tuple(sorted(tuple(sorted((name,)+tail)) for j,(name,g) in enumerate(types[first:],first)
            if len(g)<=rank for tail in forests(rank-len(g),j)))
    universe=set(forests(13));assert len(universe)==6789
    present={tuple(v['forest']) for v in rule['vertices']};assert present<=universe
    pairs={(tuple(v['forest']),B.dump(v['marked_forest_code'])):v['index']
           for v in rule['variables'] if v.get('kind')!='maximal_simplex_clique_charge'}
    active_columns={};column_models={}
    for c in rule['packing_columns']:
        gamma=weights.get(c['variable'],F(0))*c['denominator']
        if not gamma:continue
        col=packing['successor_columns'][c['source_column']];forest=tuple(col['forest'])
        assert forest in present
        active_columns.setdefault(forest,[]).append((col,gamma))
        graph=P.forest_graph(forest,models)
        column_models[col['index']]=B.C.pieces(graph)
    block_cache={}
    def block_coefficient(col,ports,blockcode):
        key=(col['index'],blockcode)
        if key not in block_cache:
            size=len(ports);U=set(col['subset']);k=col['minimum_block_size']
            if size<k or size>len(U):a=F(0)
            else:
                pieces=column_models[col['index']]
                orbit=[x for x in combinations(range(13),size) if P.code(pieces,x)==blockcode]
                assert orbit
                a=F(sum(set(x)<=U for x in orbit),len(orbit))
                assert a==F(sum(set(orbit[0])<=set(v) for v in col['subset_orbit']),len(col['subset_orbit']))
            block_cache[key]=a
        return block_cache[key]
    evaluations=[];missing=set();counts=Counter();byclass=defaultdict(list);witnesses=[]
    for star in full['stars']:
        corners=P.star_corners(star,finite);required=2*comb(len(corners),3)
        pair_coeff=Counter();packing_coeff=defaultdict(F);absent=set()
        for omitted,forest,ports,pieces,blockcode in corners:
            assert forest in universe
            if forest not in present:absent.add(forest)
            for pair in combinations(ports,2):
                key=(forest,B.dump(P.code(pieces,pair)))
                if forest in present:
                    assert key in pairs;pair_coeff[pairs[key]]+=1
            for col,gamma in active_columns.get(forest,[]):
                packing_coeff[col['index']]+=block_coefficient(col,ports,blockcode)
        pp=sum((m*weights.get(i,F(0)) for i,m in pair_coeff.items()),F(0))
        gp=F(0)
        for cols in active_columns.values():
            for col,gamma in cols:gp+=gamma*packing_coeff[col['index']]
        total=pp+gp;passed=total>=required;missing.update(absent)
        row=dict(star=star['index'],fixed_discriminant=star['fixed_discriminant'],
            discriminant_class=star.get('discriminant_class'),required=required,
            pair_coefficients=[list(x) for x in sorted(pair_coeff.items())],
            packing_coefficients=[[i,str(a)] for i,a in sorted(packing_coeff.items()) if a],
            pair_payment=str(pp),packing_payment=str(gp),payment=str(total),slack=str(total-required),
            ratio=str(total/required),passed=passed,zero_extended_forests=[list(x) for x in sorted(absent)])
        evaluations.append(row)
        kind='fixed' if star['fixed_discriminant'] else 'flexible'
        counts[kind+'_stars']+=1;counts[kind+('_passed' if passed else '_failed')]+=1
        if absent:counts['stars_using_zero_extended_forests']+=1
        if star['fixed_discriminant']:byclass[star['discriminant_class']].append(row)
        else:assert passed
        if not passed:
            g=P.C.graph(range(14),{(i,j):m for i,j,m in star['edges']})
            core_sizes=[]
            for piece in B.C.pieces(g):
                finite_piece=(sum(map(len,piece.values()))==2*len(piece)-2 and B.C.code(piece) in finite)
                if not finite_piece:core_sizes.append(len(piece))
            assert len(core_sizes)==1
            witnesses.append((core_sizes[0],len(star['edges']),star['index']))
        if (star['index']+1)%2000==0:print('ALL_PACKING_STAR_PAYMENTS',star['index']+1,dict(counts),flush=True)
    assert len(evaluations)==12420 and counts['flexible_stars']==10185
    class_rows=[]
    for ci,rr in sorted(byclass.items()):
        worst=min(rr,key=lambda r:(F(r['ratio']),r['star']))
        class_rows.append(dict(discriminant_class=ci,stars=len(rr),passed=sum(r['passed'] for r in rr),
            failed=sum(not r['passed'] for r in rr),worst_star=worst['star'],worst_ratio=worst['ratio']))
    failed=[r for r in evaluations if not r['passed']]
    minimal=None;worst=None
    if failed:
        _,_,si=min(witnesses);minimal=dict(evaluation=evaluations[si],literal_star=full['stars'][si],
            selection='Smallest connected Lorentz-core order, then fewest edges, then star index among unpaid saved stars.')
        w=min(failed,key=lambda r:(F(r['ratio']),r['star']));worst=dict(evaluation=w,literal_star=full['stars'][w['star']])
    paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',f'scripts/{D.TAG}.py',f'scripts/{P.TAG}.py',
        f'results/{P.TAG}.json',f'results/{RULE}.json',f'results/{STARS}.json']
    return dict(status='PASS_EXACT_ALL_SAVED_STAR_SEPARATION',sources={p:P.sha(p) for p in paths},
        dimension=13,target=12,counts=dict(counts),vertex_extension=dict(all_finite_forests=6789,
            explicitly_budgeted_forests=len(present),additional_zero_forests=len(universe-present),
            maximum=primal['maximum'],definition='Every pair payment is zero on absent forests. Packing charges are supported only on their22 original forests. Thus every absent finite forest has total vertex budget zero; the old exact maximum bounds the complete finite forest universe.'),
        missing_forests_in_saved_stars=[list(x) for x in sorted(missing)],
        fixed_classes=class_rows,evaluations=evaluations,minimal_witness=minimal,worst_witness=worst,
        scope='Nonoptimizing exact evaluation of one generic packing rule on all12420 saved stars. Fixed-class failures are listed separately and do not disprove the rule on a discriminant class excluding those stars. This is not complete triangle/star coverage; quadrilaterals and other connected cores remain outside. Zero extension completes vertex budgets only, not face coverage.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['vertex_extension'],flush=True)
