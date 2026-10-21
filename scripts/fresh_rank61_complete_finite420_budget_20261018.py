#!/usr/bin/env python3
"""Complete finite rank61 certificate after the simultaneous local transfers.

New residual tables are generated item first. Fresh replay checks every
last-item upper inequality and an attaining predecessor, without repeating
the optimizing recurrence. Geometric inputs must be frozen before running.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_depth8_d17_rank61_forest_budget_20261018 as B

ROOT=B.ROOT;N=61;TAG='fresh_rank61_complete_finite420_budget_20261018'
A15='fresh_finite_a15_endpoint_coupon_zero_20261018'
BD15='fresh_finite_b15_d15_existing_payment_transfer_20261018'
A16='fresh_finite_a16_endpoint_sparse_half_transfer_20261018'
D16='fresh_finite_d16_fork_exchange_20261018'
STEMS=(B.TAG,A15,BD15,A16,D16)
D16_FINAL_MD='8ac18b4643f4d3f81dd937f0a3ee03a05859250c603112d8f2d40490cc47c0f9'
ACTUAL=dict(A15='103',B15='103',D15='103',A16='111',B16='104',D16='111',A17='115',B17='105',D17='119')
# These are increments beyond the complete frozen425 counting rule.
EXTRA_CROSS={
    (('A',2),('A',13)):4,
    (('A',1),('A',13)):2,
    (('A',2),('A',14)):4,
    (('A',4),('A',12)):4,
}


def base(t):
    name=t[0]+str(t[1])
    return int(2*F(ACTUAL[name])) if name in ACTUAL else B.base(t)


def actual_cross(a,b):
    return B.cross(a,b)+EXTRA_CROSS.get((a,b),0)+EXTRA_CROSS.get((b,a),0)


def upper_cross(a,b):
    # The A15 coupon upgrade is dominated by an additional two units on
    # representative A2/D13. This is an allowance, not an actual new coupon.
    return B.cross(a,b)+(4 if {a,b}=={('A',2),('D',13)} else 0)


def cost(ts,upper=False):
    cross=upper_cross if upper else actual_cross
    return sum(base(t) for t in ts)+sum(cross(a,b) for i,a in enumerate(ts) for b in ts[i+1:])


def item(r,key):
    # D17 is a residual item; its restored two halves change every table.
    # The A2/D13 allowance changes the rank2 interaction in that table.
    return B.item(r,key)+(4*key[2] if r==2 else 2 if r==17 else 0)


def generate_table(key):
    values=[0]+[-10**9]*N
    for r in B.ITEMS:
        w=item(r,key)
        for s in range(r,N+1):values[s]=max(values[s],values[s-r]+w)
    last=[0]+[next(r for r in B.ITEMS if r<=s and values[s]==values[s-r]+item(r,key)) for s in range(1,N+1)]
    return dict(key=list(key),values=values,last=last)


def check_table(tab,key):
    assert tab['key']==list(key)
    v,last=tab['values'],tab['last'];assert len(v)==len(last)==N+1 and v[0]==last[0]==0
    checks=0
    for s in range(1,N+1):
        for r in B.ITEMS:
            if r<=s:assert v[s]>=v[s-r]+item(r,key);checks+=1
        r=last[s];assert r in B.ITEMS and r<=s and v[s]==v[s-r]+item(r,key)
    return checks


def interfaces(docs):
    assert D16_FINAL_MD!='PENDING','Await final D16 geometric report binding.'
    assert sha256((ROOT/f'results/{D16}.md').read_bytes()).hexdigest()==D16_FINAL_MD
    assert docs[B.TAG]['finite_maximum']=='425' and docs[B.TAG]['current_component_costs']==dict(A15='104',B15='104',D15='104',A16='113',D16='225/2',A17='115',D17='118')
    a=docs[A15];assert a['new_cost']=='103'
    assert a['zero_pairs']==[[0,13],[1,14]] and a['literal_debited_corners']==0
    c=a['upgraded_coupon']
    assert (c['small'],c['small_marks'],c['large'],c['large_marks'])==('whole A2',[0,1],'whole A13',[0,12])
    assert (c['old_weight'],c['new_weight'],c['component_pair_cost_increase'])==('1/2','1','2')
    assert docs[BD15]['new_costs']==dict(B15='103',D15='103')
    assert docs[BD15]['new_cross_component_coupons']==[] and docs[BD15]['literal_changed_corners']=={}
    a=docs[A16];assert a['new_cost']=='111' and a['literal_debited_corners']==0
    assert a['new_cross_component_coupons']==[
        dict(name='A1_A13end',small='whole A1',small_nodes=[0],large='whole A13',large_nodes=[0,12],weight='1/2',component_pair_cost='1'),
        dict(name='A2_A14third',small='whole A2',small_nodes=[0,1],large='whole A14',large_nodes=[2,11],weight='1/2',component_pair_cost='2')]
    d=docs[D16]
    assert d['new_costs']==dict(D16='111',D17='119')
    assert d['D16_half_pairs']==[[0,14],[0,15],[2,12]]
    assert d['D17_restored_half_pairs']==[[1,15],[1,16]]
    assert d['retained_rank15_alternatives']==dict(D16_pairs=[[1,14],[1,15]],weight='1/2')
    assert d['new_cross_component_coupons']==[
        dict(name='A4end_A12third',small='whole A4',small_nodes=[0,3],large='whole A12',large_nodes=[2,9],weight='1/2',component_pair_cost='2')]
    assert sum(F(r['weight']) for r in d['D16_weighted_pairs'])==111
    assert sum(F(r['weight']) for r in d['D17_weighted_pairs'])==119
    assert d['compatible_feature_depth']==8 and d['literal_configurations']==216
    for s in (A15,BD15,A16,D16):
        assert docs[s]['triangle_minimum']==2 and docs[s]['quadrilateral_minimum']==1
    return dict(A15_coupon_increase='2',A16_new_coupon_costs=['1','2'],D16_new_coupon_cost='2',
                D17_restoration_increase='1',all_local_inputs_final=True)


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    local=interfaces(docs)
    ts=B.types();pairchecks=0
    for a in ts:
        assert base(B.representative(a))>=base(a)
        for b in ts:
            assert actual_cross(a,b)==actual_cross(b,a)
            assert upper_cross(B.representative(a),B.representative(b))>=actual_cross(a,b)
            pairchecks+=1
    # A possible future A1/B10 half-to-unit increase is not used here. The
    # saved inequality only records that its representative allowance fits.
    optional_A1_B10_absorption=(actual_cross(('A',1),('B',10))+1<=upper_cross(('A',1),('D',10)))
    assert optional_A1_B10_absorption
    dom=B.domain();assert len(dom)==1132
    assert {c for c,u,k in dom}=={tuple(r['counts']) for r in docs[B.TAG]['alternatives']}
    keys=sorted({k for c,u,k in dom});assert len(keys)==854
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {k:generate_table(k) for k in keys}
    assert set(tables)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    inequalities=sum(check_table(tables[k],k) for k in keys)
    records=[];itemchecks=0
    for cs,used,key in dom:
        large=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in large)==used
        assert all(upper_cross(a,b)==0 for i,a in enumerate(large) for b in large[i+1:])
        for r in B.ITEMS:
            t=B.representative(('A' if r<=8 else 'D',r))
            assert item(r,key)==base(t)+sum(upper_cross(t,b) for b in large);itemchecks+=1
        rem=N-used;tab=tables[key];total=sum(base(t) for t in large)+tab['values'][rem]
        forest=large[:]
        while rem:
            r=tab['last'][rem];forest.append(B.representative(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==N and cost(forest,True)==total
        actual=cost(forest);assert actual<=total
        records.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),doubled_upper=total,doubled_actual=actual))
    maximum=max(r['doubled_upper'] for r in records);assert maximum==840
    witnesses=[[('D',61)],[('D',17),('D',17),('D',27)],[('D',15)]*3+[('A',16)],[('D',15)]*3+[('D',16)]]
    for forest in witnesses:assert sum(r for k,r in forest)==N and cost(forest)==maximum
    G=sum(min(j,N-1-j) for j in range(N));curvature=F(5*G-(N-1)*N,2)
    assert (G,curvature)==(900,F(420))
    sources={f'scripts/{TAG}.py':sha256(Path(__file__).read_bytes()).hexdigest()}
    for s in STEMS:
        for folder,ext in (('scripts','py'),('results','json'),('results','md')):
            p=f'{folder}/{s}.{ext}';sources[p]=sha256((ROOT/p).read_bytes()).hexdigest()
    if saved:assert saved['sources']==sources
    return dict(status='PASS_COMPLETE_RANK61_FINITE420_BUDGET',sources=sources,
        ambient_dimension=N,finite_maximum=str(F(maximum,2)),finite_curvature=str(curvature),finite_slack='0',
        finite_link_G=G,finite_link_vertices=N,actual_component_costs=ACTUAL,local_interfaces=local,
        actual_added_cross_costs=[dict(small=a[0]+str(a[1]),large=b[0]+str(b[1]),increment=str(F(v,2))) for (a,b),v in EXTRA_CROSS.items()],
        conservative_extra_A2_D13_allowance='2',optional_A1_B10_half_to_unit_absorption=optional_A1_B10_absorption,
        optional_A1_B10_rule_assumed=False,simultaneous_type_pairs=pairchecks,complete_large_vectors=len(dom),interaction_keys=len(keys),
        tables=[tables[k] for k in keys],Bellman_upper_inequalities=inequalities,attaining_predecessors=len(keys)*N,
        residual_item_checks=itemchecks,alternatives=records,
        actual_attaining_forests=[sorted(k+str(r) for k,r in fs) for fs in witnesses],
        upper_attaining_alternatives=[r for r in records if r['doubled_upper']==maximum],
        new_residual_tables=True,inherited_old_tables_reused=False,geometric_inputs_final=True,new_dimension_bound=None,
        scope='Complete finite rank61 upper420 and actual formal-forest attainment under the simultaneous bound local transfers. All added typed roles are counted; A2/D13 extra2 is conservative only. New residual tables include D17=119 and the changed A2 interaction. No physical realization or ambient60 theorem is asserted.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    result=run(saved)
    if args.verify:assert result==saved
    else:
        assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'maximum',result['finite_maximum'],'vectors',result['complete_large_vectors'],'tables',result['interaction_keys'],flush=True)
