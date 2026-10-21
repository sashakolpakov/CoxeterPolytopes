#!/usr/bin/env python3
"""Exact current ambient61 finite/cusp frontier; no dimension60 claim.

The finite part restricts frozen residual certificates to rank61. The cusp
part reconstructs the complete rank60 low9-only domain and exact convolution.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import fresh_finite_rank15_rank16_rank62_budget_20261018 as B

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ambient61_payment_frontier_20261018'
COMPLETED = 'fresh_rank15_rank16_dimension61_20261018'
IDEAL = 'fresh_order14_ambient62_ideal_margin_closure_20261018'
STEMS = (COMPLETED, B.TAG, B.OLD, IDEAL)
N, D = 61, 60


def large_domain():
    out = []
    def visit(r, left, cs):
        if r == 18:
            c = dict(zip(range(9,18),cs))
            key = (2*c[10]+3*sum(c[s] for s in range(11,17))+2*c[17],
                   c[11]+c[12],c[13],c[11],c[14],c[9],c[10])
            out.append((tuple(cs),N-left,key));return
        for count in range(left//r+1):visit(r+1,left-r*count,cs+[count])
    visit(9,N,[])
    return sorted(out)


def finite(docs,curvature):
    previous,old = docs[B.TAG],docs[B.OLD]
    assert previous['finite_maximum']=='434' and previous['geometric_inputs_final']
    tables = {tuple(t['key']):t for t in old['certificate']['tables']}
    digest = sha256(json.dumps(old['certificate']['tables'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert digest==previous['inherited_table_sha256']
    dom=large_domain();assert len(dom)==1132
    inherited={tuple(r['counts']):r for r in previous['alternatives']
               if sum(r*c for r,c in zip(range(9,18),r['counts']))<=N}
    assert set(inherited)=={r[0] for r in dom}
    keys={r[2] for r in dom};assert len(keys)==854 and keys<=set(tables)
    types=[t for t in B.type_list() if t[1]<=N]
    for a in types:
        assert B.upper_base(B.representative(a))>=B.actual_base(a)
        for b in types:
            assert B.cross(B.representative(a),B.representative(b))>=B.cross(a,b)
    records=[];item_checks=0
    for counts,used,key in dom:
        large=[('B' if r==9 else 'A' if r==17 else 'D',r)
               for r,c in zip(range(9,18),counts) for _ in range(c)]
        assert key==tuple(inherited[counts]['key'])
        assert all(B.cross(a,b)==0 for i,a in enumerate(large) for b in large[i+1:])
        for r in tuple(range(1,9))+tuple(range(17,N+1)):
            t=B.representative(('A' if r<=8 else 'D',r))
            value=B.upper_base(t)+sum(B.cross(t,b) for b in large)
            broad,extra,a4,d11,d14,b9,d10=key
            encoded=B.before_base(t)+(r*broad if r<=5 else 0)
            encoded+=(2*extra+d14 if r==1 else
                      4*d11+2*d14+2*b9+4*d10 if r==2 else
                      2*d10 if r==3 else 2*a4 if r==4 else 8*d11 if r==6 else 0)
            assert value==encoded;item_checks+=1
        tab=tables[key];rem=N-used
        assert 0<=rem<len(tab['values'])==len(tab['last'])
        upper=sum(B.upper_base(t) for t in large)+tab['values'][rem]
        forest=large[:]
        while rem:
            r=tab['last'][rem]
            assert r<=rem and (1<=r<=8 or r>=17)
            forest.append(B.representative(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==N
        assert B.cost(forest,B.upper_base)==upper
        actual=B.cost(forest,B.actual_base);assert actual<=upper
        records.append(dict(counts=list(counts),key=list(key),forest=sorted(k+str(r) for k,r in forest),
                            doubled_upper=upper,doubled_actual=actual))
    maximum=max(r['doubled_upper'] for r in records)
    assert maximum==854
    witnesses=[['D17','D44'],['A16','D45']]
    for names in witnesses:
        ts=[(s[0],int(s[1:])) for s in names]
        assert sum(r for k,r in ts)==N and B.cost(ts,B.actual_base)==maximum
    singleton=B.actual_base(('D',61));assert singleton==852
    above=[r for r in records if F(r['doubled_upper'],2)>curvature]
    # Actual formal witnesses are necessary targets, not a sufficient set of constraints.
    necessary=[]
    for names in (['D61'],*witnesses,['D15','D46'],['A16','D15','D15','D15'],
                  ['A16','A16','D29'],['A17','D44'],['A1','D15','D15','D15','D15']):
        ts=[(s[0],int(s[1:])) for s in names]
        assert sum(r for k,r in ts)==N
        value=F(B.cost(ts,B.actual_base),2)
        assert value>curvature
        necessary.append(dict(forest=names,actual_cost=str(value),required_total_saving=str(value-curvature)))
    return dict(exact_maximum=str(F(maximum,2)),curvature=str(curvature),finite_slack=str(curvature-F(maximum,2)),
        actual_attaining_forests=witnesses,complete_large_vectors=len(dom),interaction_keys=len(keys),
        simultaneous_type_pairs=len(types)**2,unchanged_residual_items=item_checks,
        inherited_table_sha256=digest,residual_rank_restricted_to=N,new_Bellman_optimization=False,
        alternatives=records,representative_alternatives_above_curvature=sorted(above,key=lambda r:(-r['doubled_upper'],r['forest'])),
        actual_singleton_D61_cost=str(F(singleton,2)),necessary_actual_witness_targets=necessary,
        representative_D16_upper='113',actual_D16_cost='225/2',scope='Complete current-rule rank61 upper and actual attainment. Conservative representative alternatives need not attain their named actual costs. Necessary witness targets do not assert sufficient reductions or new geometric rules.')


def high_parts(total,least=9):
    if not total:yield ();return
    for r in range(least,total+1):
        for tail in high_parts(total-r,r):yield (r,)+tail


def family_cost(kind,m):
    assert m>=18
    if kind=='A':return 5*m
    offset,cap,first={'B':(-6,13,8),'C':(-15,14,6),'D':(3,12,10)}[kind]
    return 5*m+offset+sum(min(m-k,cap) for k in range(first,first+4))


def envelope(r):
    if r<=9:return F(r,2)
    small={10:F(54,11),11:F(63,12),12:F(73,13),13:F(83,14),
           14:F(94,15),15:F(106,16),16:F(116,17)}
    if r in small:return small[r]
    return max(F(family_cost(k,r+1),r+1) for k in 'ABCD')


def scalar(rs):
    h=[1]
    for r in rs:
        out=[0]*(len(h)+r)
        for i,x in enumerate(h):
            for j in range(r+1):out[i+j]+=x
        h=out
    V=prod(r+1 for r in rs)
    assert sum(rs)==D and len(h)==D+1 and h==h[::-1] and sum(h)==V
    G=sum(x*min(i,D-i) for i,x in enumerate(h))
    assert 2*G==D*V-sum(x*abs(2*i-D) for i,x in enumerate(h))
    cumulative=0;deviation=F()
    for i in range(D):
        cumulative+=h[i]
        deviation+=F(cumulative,V) if i<D//2 else F(V-cumulative,V)
    assert F(G,V)==F(D,2)-deviation
    curvature=F(5*G-D*V,2);cost=sum(map(envelope,rs),F())
    margin=curvature/V-cost
    return dict(ranks=list(rs),link_vertices=V,h_polynomial=h,G=G,curvature=str(curvature),
                normalized_cost=str(cost),normalized_margin=str(margin),absolute_margin=str(V*margin))


def summarize(rows):
    normal=min(F(r['normalized_margin']) for r in rows)
    absolute=min(F(r['absolute_margin']) for r in rows)
    return dict(count=len(rows),negative_count=sum(F(r['normalized_margin'])<0 for r in rows),
        zero_count=sum(F(r['normalized_margin'])==0 for r in rows),
        normalized_minimum=str(normal),normalized_minimizers=[r['ranks'] for r in rows if F(r['normalized_margin'])==normal],
        absolute_minimum=str(absolute),absolute_minimizers=[r['ranks'] for r in rows if F(r['absolute_margin'])==absolute],
        negative_rows=[dict(ranks=r['ranks'],normalized_margin=r['normalized_margin'],absolute_margin=r['absolute_margin'])
                       for r in sorted(rows,key=lambda x:(F(x['normalized_margin']),x['ranks'])) if F(r['normalized_margin'])<0])


def conditional(rows,drops):
    out=[]
    for row in rows:
        gain=sum(F(row['ranks'].count(r),r+1)*delta for r,delta in drops.items())
        margin=F(row['normalized_margin'])+gain
        out.append(dict(ranks=row['ranks'],normalized_margin=str(margin),absolute_margin=str(row['link_vertices']*margin)))
    return dict(conditional_total_cost_savings_by_affine_rank={str(r):str(v) for r,v in drops.items()},
                summary=summarize(out),geometry_proved=False,additional_finite_coupons_counted=False)


def threshold(rows,varied_rank,fixed):
    bound=F();critical=[];unaffected=[]
    for row in rows:
        margin=F(row['normalized_margin'])+sum(F(row['ranks'].count(r),r+1)*delta for r,delta in fixed.items())
        slope=F(row['ranks'].count(varied_rank),varied_rank+1)
        if not slope:
            if margin<=0:unaffected.append(dict(ranks=row['ranks'],margin=str(margin)))
            continue
        need=-margin/slope
        if need>bound:bound=need;critical=[row['ranks']]
        elif need==bound:critical.append(row['ranks'])
    return dict(varied_affine_rank=varied_rank,fixed_total_cost_savings={str(r):str(v) for r,v in fixed.items()},
                strict_threshold=str(bound),critical_profiles=critical,unaffected_nonpositive=unaffected,
                sufficient_for_strict_terminal_positivity=not unaffected,geometry_proved=False)


def run():
    docs={}
    for s in STEMS:
        d=json.loads((ROOT/f'results/{s}.json').read_text());assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
        docs[s]=d
    assert docs[COMPLETED]['dimension_bound']==61
    assert docs[COMPLETED]['actual_component_costs']==dict(A15='104',B15='104',D15='104',A16='113',B16='104',D16='225/2',A17='115',B17='105',D17='120')
    prior=docs[IDEAL]['rows'];assert len(prior)==2641
    assert all(sum(map(envelope,r['ranks']),F())==F(r['normalized_cost']) for r in prior)
    finite_G=sum(min(j,D-j) for j in range(D+1))
    curvature=F(5*finite_G-D*(D+1),2)
    extremum=sum(2*j-D for j in range(D+1) if 2*j>D)
    assert finite_G==900 and curvature==420 and extremum==930
    profiles=set(high_parts(D))
    for low in range(1,9):profiles.update((low,)+p for p in high_parts(D-low))
    assert len(profiles)==2365
    rows=[scalar(rs) for rs in sorted(profiles)]
    summary=summarize(rows)
    assert summary['negative_count']==21 and summary['zero_count']==0
    scenarios=[conditional(rows,d) for d in ({10:F(1,4),15:F(1)}, {10:F(1,2),15:F(1,2)},
                {10:F(1,2),15:F(1)}, {10:F(1)}, {15:F(1)}, {15:F(2)})]
    thresholds=[threshold(rows,10,{15:F(1)}),threshold(rows,15,{10:F(1,4)}),
                threshold(rows,10,{}),threshold(rows,15,{})]
    paths=[Path(__file__),Path(B.__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS
                                           for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_CURRENT_AMBIENT61_PAYMENT_FRONTIER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        inherited_completed_dimension_bound=61,ambient_dimension=N,link_dimension=D,
        finite=finite(docs,curvature),ideal=dict(complete_profiles=len(rows),rows=rows,summary=summary,
            small_rank_envelopes={str(r):str(envelope(r)) for r in range(1,17)},
            exact_family_tail=[dict(order=m,costs={k:family_cost(k,m) for k in 'ABCD'},envelope=str(envelope(m-1))) for m in range(18,63)],
            prior_ambient62_cost_comparisons=len(prior),conditional_scenarios=scenarios,conditional_thresholds=thresholds,
            low9_concentration_only=True,high_factor_concentration_used=False,
            negative_terminal_absolute_margin_is_universal_bound=False),
        finite_curvature=str(curvature),finite_G=finite_G,extremum_difference=extremum,
        strict_global_constant=str(5*extremum),new_dimension_bound=None,
        scope='Exact current rank61 finite maximum and rank60 low9-only ideal scalar frontier under completed61 rule. Frozen residual certificates are restricted to smaller residual ranks with unchanged items; no Bellman or global optimization. Conditional ideal savings are arithmetic targets only, with no geometry or extra-coupon claim. Negative terminal absolute margins are not extended to all cusps.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(p.read_text())
    else:
        assert not p.exists();p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'ideal',out['ideal']['complete_profiles'],
          'negative',out['ideal']['summary']['negative_count'],'constant',out['strict_global_constant'],flush=True)
