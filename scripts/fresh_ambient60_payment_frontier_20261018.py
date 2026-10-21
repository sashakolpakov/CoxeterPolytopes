#!/usr/bin/env python3
"""Exact ambient60 frontier under the completed dimension60 payment rule.

New finite tables use item-first generation; replay checks every last-item
upper inequality and an attaining predecessor. Ideal rows use fresh exact
convolutions on the complete low9-only terminal partition domain.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import fresh_rank61_complete_finite420_budget_20261018 as B

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_ambient60_payment_frontier_20261018'
COMPLETE='fresh_order11_order16_dimension60_20261018'
PREVIOUS='fresh_ambient61_payment_frontier_20261018'
STEMS=(COMPLETE,B.TAG,B.B.TAG,PREVIOUS)
N,D=60,59
ITEMS=tuple(range(1,9))+tuple(range(17,N+1))


def current_cross(a,b):
    # The final BC theorem adds the optional A1/B10 half already absorbed
    # by the frozen finite420 representative allowance.
    return B.actual_cross(a,b)+(1 if {a,b}=={('A',1),('B',10)} else 0)


def actual_cost(ts):
    return sum(B.base(t) for t in ts)+sum(current_cross(a,b) for i,a in enumerate(ts) for b in ts[i+1:])


def large_domain():
    rows=[]
    def visit(r,left,cs):
        if r==18:
            c=dict(zip(range(9,18),cs))
            key=(2*c[10]+3*sum(c[s] for s in range(11,17))+2*c[17],
                 c[11]+c[12],c[13],c[11],c[14],c[9],c[10])
            rows.append((tuple(cs),N-left,key));return
        for count in range(left//r+1):visit(r+1,left-r*count,cs+[count])
    visit(9,N,[])
    return sorted(rows)


def generate_table(key):
    values=[0]+[-10**9]*N
    for r in ITEMS:
        w=B.item(r,key)
        for s in range(r,N+1):values[s]=max(values[s],values[s-r]+w)
    last=[0]+[next(r for r in ITEMS if r<=s and values[s]==values[s-r]+B.item(r,key)) for s in range(1,N+1)]
    return dict(key=list(key),values=values,last=last)


def check_table(tab,key):
    assert tab['key']==list(key)
    v,last=tab['values'],tab['last']
    assert len(v)==len(last)==N+1 and v[0]==last[0]==0
    checks=0
    for s in range(1,N+1):
        for r in ITEMS:
            if r<=s:assert v[s]>=v[s-r]+B.item(r,key);checks+=1
        r=last[s]
        assert r in ITEMS and r<=s and v[s]==v[s-r]+B.item(r,key)
    return checks


def finite(docs,curvature,saved=None):
    old=docs[B.TAG]
    assert old['finite_maximum']=='420' and old['geometric_inputs_final']
    assert old['actual_component_costs']==B.ACTUAL
    assert old['optional_A1_B10_half_to_unit_absorption'] and not old['optional_A1_B10_rule_assumed']
    assert current_cross(('A',1),('B',10))==B.upper_cross(('A',1),('D',10))==2
    types=[t for t in B.B.types() if t[1]<=N]
    for a in types:
        assert B.base(B.B.representative(a))>=B.base(a)
        for b in types:
            assert current_cross(a,b)==current_cross(b,a)
            assert B.upper_cross(B.B.representative(a),B.B.representative(b))>=current_cross(a,b)
    residual_types=[B.B.representative(('A' if r<=8 else 'D',r)) for r in ITEMS]
    assert all(B.upper_cross(a,b)==0 for a in residual_types for b in residual_types)
    dom=large_domain()
    assert len(dom)==1044
    prior={tuple(r['counts']) for r in old['alternatives'] if sum(i*c for i,c in zip(range(9,18),r['counts']))<=N}
    assert {c for c,u,k in dom}==prior
    keys=sorted({k for c,u,k in dom})
    assert len(keys)==792
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {k:generate_table(k) for k in keys}
    assert set(tables)==set(keys) and (not saved or len(saved['tables'])==len(keys))
    inequalities=sum(check_table(tables[k],k) for k in keys)
    old_tables={tuple(t['key']):t for t in old['tables']}
    for k in keys:
        # A comparison with the frozen rank61 prefixes checks the generated
        # values, but generation above did not import those values.
        assert tables[k]['values']==old_tables[k]['values'][:N+1]
    records=[];item_checks=0
    for cs,used,key in dom:
        large=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in large)==used
        assert all(B.upper_cross(a,b)==0 for i,a in enumerate(large) for b in large[i+1:])
        for r in ITEMS:
            t=B.B.representative(('A' if r<=8 else 'D',r))
            assert B.item(r,key)==B.base(t)+sum(B.upper_cross(t,b) for b in large)
            item_checks+=1
        rem=N-used;tab=tables[key]
        upper=sum(B.base(t) for t in large)+tab['values'][rem]
        forest=large[:]
        while rem:
            r=tab['last'][rem];forest.append(B.B.representative(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==N
        assert B.cost(forest,True)==upper
        actual=actual_cost(forest);assert actual<=upper
        records.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
                            doubled_upper=upper,doubled_actual=actual))
    maximum=max(r['doubled_upper'] for r in records)
    witnesses=[['D60'],['D17','D43'],['D17','D17','D26']]
    assert maximum==826
    for fs in witnesses:
        ts=[(name[0],int(name[1:])) for name in fs]
        assert sum(r for k,r in ts)==N and actual_cost(ts)==maximum
    assert witnesses,'An actual formal forest must attain the representative upper.'
    above=sorted((r for r in records if F(r['doubled_upper'],2)>curvature),key=lambda r:(-r['doubled_upper'],r['forest']))
    actual_above=[dict(forest=r['forest'],cost=str(F(r['doubled_actual'],2)),deficit=str(F(r['doubled_actual'],2)-curvature))
                  for r in above if F(r['doubled_actual'],2)>curvature]
    # This is an explicitly hypothetical monotone deletion, not a local
    # theorem input. Its unchanged attaining forest makes a new DP unnecessary.
    deep7_witness=[('D',17),('D',17),('D',26)]
    deep7_removed=sum(max(0,r-{'A':29,'B':30,'D':31}[k]) for k,r in deep7_witness)
    assert deep7_removed==0 and actual_cost(deep7_witness)==maximum
    return dict(exact_maximum=str(F(maximum,2)),curvature=str(curvature),finite_slack=str(curvature-F(maximum,2)),
        current_component_costs=B.ACTUAL,actual_attaining_forests=witnesses,
        actual_singleton_D60_cost=str(F(B.base(('D',60)),2)),actual_four_D15_cost=str(F(actual_cost([('D',15)]*4),2)),
        complete_large_vectors=len(dom),interaction_keys=len(keys),simultaneous_type_pairs=len(types)**2,
        residual_item_checks=item_checks,Bellman_upper_inequalities=inequalities,attaining_predecessors=len(keys)*N,
        tables=[tables[k] for k in keys],alternatives=records,representative_alternatives_above_curvature=above,
        actual_witnesses_above_curvature=actual_above,final_A1_B10_unit_counted=True,
        conditional_deep7_depth11_only=dict(geometry_proved=False,all_other_costs_and_coupons_unchanged=True,
            assumed_deleted_pair_counts={k:'max(0,n-'+str(t)+')' for k,t in dict(A=29,B=30,D=31).items()},
            unchanged_actual_attainer=['D17','D17','D26'],attainer_deleted_pairs=deep7_removed,
            exact_maximum_by_monotonicity_and_attainment=str(F(maximum,2)),
            scope='Conditional monotone deletion only. No local geometry assumed proved or new forest optimization used.'),
        new_Bellman_generation=True,replay_optimizing_recurrence=False,
        rank61_table_prefix_comparisons=len(keys),scope='Complete current-rule formal rank60 forest upper and actual attainment; no physical realization asserted. Representative alternatives may be conservative. Witness deficits are necessary arithmetic targets, not sufficient new geometric rules.')


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
    small={10:F(53,11),11:F(63,12),12:F(73,13),13:F(83,14),
           14:F(94,15),15:F(105,16),16:F(116,17)}
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
    # For odd D the two middle positions straddle D/2; each CDF term
    # contributes the exact mean absolute deviation on its side.
    cumulative=0;deviation=F()
    for i in range(D):
        cumulative+=h[i]
        if 2*i<D:deviation+=F(cumulative,V)
        else:deviation+=F(V-cumulative,V)
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


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[COMPLETE]['dimension_bound']==60
    assert docs[COMPLETE]['actual_component_costs']==B.ACTUAL
    assert docs[COMPLETE]['ideal_order11_costs']==dict(A='99/2',B='53',C='52',D='53')
    assert docs[COMPLETE]['ideal_order16_costs']==dict(A='80',B='105',C='103',D='105')
    assert all(sum(map(envelope,r['ranks']),F())==F(r['normalized_cost']) for r in docs[COMPLETE]['rows'])
    oldtail=docs[PREVIOUS]['ideal']['exact_family_tail']
    assert all(r['costs']=={k:family_cost(k,r['order']) for k in 'ABCD'} and F(r['envelope'])==envelope(r['order']-1) for r in oldtail)
    finite_G=sum(min(j,D-j) for j in range(D+1))
    curvature=F(5*finite_G-D*(D+1),2)
    extremum=sum(2*j-D for j in range(D+1) if 2*j>D)
    assert (finite_G,curvature,extremum)==(870,F(405),900)
    profiles=set(high_parts(D))
    for low in range(1,9):profiles.update((low,)+p for p in high_parts(D-low))
    rows=[scalar(rs) for rs in sorted(profiles)]
    summary=summarize(rows)
    assert len(rows)==2113 and summary['negative_count']==77 and summary['zero_count']==0
    assert summary['normalized_minimum']=='-1823/11560' and summary['absolute_minimum']=='-192250'
    active=sorted({r for row in rows if F(row['normalized_margin'])<=0 for r in row['ranks'] if r>=10})
    scenarios=[conditional(rows,{r:delta for r in active}) for delta in (F(1,4),F(1,2),F(1),F(2))]
    coordinated={11:F(1,2),12:F(1,2),13:F(1,4),14:F(2),15:F(5,4),16:F(3,2)}
    scenarios.append(conditional(rows,coordinated))
    assert scenarios[-1]['summary']['negative_count']==scenarios[-1]['summary']['zero_count']==0
    assert scenarios[-1]['summary']['normalized_minimum']=='73/104000'
    thresholds=[threshold(rows,r,{}) for r in active]
    coordinated_thresholds=[threshold(rows,r,{k:v for k,v in coordinated.items() if k!=r}) for r in coordinated]
    assert all(r['sufficient_for_strict_terminal_positivity'] for r in coordinated_thresholds)
    isolated=[]
    for row in rows:
        high={r for r in row['ranks'] if r>=10}
        if F(row['normalized_margin'])<0 and len(high)==1:
            r=next(iter(high))
            isolated.append(dict(ranks=row['ranks'],affine_rank=r,affine_order=r+1,
                strict_required_total_cost_saving=str(-F(row['normalized_margin'])*F(r+1,row['ranks'].count(r)))))
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths}
    if saved:assert saved['sources']==sources
    return dict(status='PASS_CURRENT_AMBIENT60_PAYMENT_FRONTIER',sources=sources,
        inherited_completed_dimension_bound=60,ambient_dimension=N,link_dimension=D,
        finite=finite(docs,curvature,saved['finite'] if saved else None),
        ideal=dict(complete_profiles=len(rows),rows=rows,summary=summary,
            small_rank_envelopes={str(r):str(envelope(r)) for r in range(1,17)},
            exact_family_tail=[dict(order=m,costs={k:family_cost(k,m) for k in 'ABCD'},envelope=str(envelope(m-1))) for m in range(18,61)],
            completed60_cost_comparisons=len(docs[COMPLETE]['rows']),conditional_scenarios=scenarios,
            conditional_thresholds=thresholds,coordinated_strict_thresholds=coordinated_thresholds,
            isolated_order_necessary_targets=isolated,
            conditional_coordinated_new_order_totals={str(r+1):str((r+1)*envelope(r)-v) for r,v in coordinated.items()},
            low9_concentration_only=True,high_factor_concentration_used=False,
            negative_terminal_absolute_margin_is_universal_bound=False),
        finite_curvature=str(curvature),finite_G=finite_G,extremum_difference=extremum,
        strict_global_constant=str(5*extremum),new_dimension_bound=None,
        scope='Exact current rank60 finite maximum and rank59 low9-only ideal scalar frontier under completed60 rule. Fresh finite DP generation with separate all-last-item/attainment replay; fresh ideal convolution and two exact moment identities. Conditional savings are arithmetic targets only, with no local geometry or extra-coupon assertion. Negative terminal absolute margins are not extended to arbitrary cusps.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';saved=json.loads(p.read_text()) if args.verify else None
    out=run(saved)
    if args.verify:assert out==saved
    else:
        assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'ideal',out['ideal']['complete_profiles'],
          'negative',out['ideal']['summary']['negative_count'],'constant',out['strict_global_constant'],flush=True)
