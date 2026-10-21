#!/usr/bin/env python3
"""Conditional complete finite count for depth6 finite-distance7 deletion.

The geometric local-payment hypothesis is not proved by this certificate.
All other weights are those of the frozen sparse-D16/ideal17 rebalancing.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_small_d16_coupon_and_ideal17_rebalance_20261018 as S

ROOT=S.ROOT;N=60;Z=S.Z
TAG='fresh_depth6_sparse_ambient60_forest_frontier_20261018'
FINAL_MD='c0c10822e322c31689b21a53b2bc64b03522dd7077fbcac1758dfdf19ddacb28'
OLD={'A':23,'B':24,'D':25};NEW={'A':19,'B':20,'D':21}


def deletions(t):
    k,n=t
    return max(0,n-NEW[k])-max(0,n-OLD[k]) if k in NEW else 0


def actual_base(t):return S.actual_base(t)-8*deletions(t)
def upper_base(t):return S.upper_base(t)-8*deletions(t)
def actual_cross(a,b):return S.actual_cross(a,b)
def upper_cross(a,b):return S.upper_cross(a,b)
def item(r,key):return S.item(r,key)-8*deletions(('A' if r<=8 else 'D',r))


def value(forest,upper=False):
    base,cross=(upper_base,upper_cross) if upper else (actual_base,actual_cross)
    return sum(base(t) for t in forest)+sum(cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])


def table(key):
    ws={r:item(r,key) for r in Z.ITEMS};v=[0]+[-10**9]*N
    for r,w in ws.items():
        for n in range(r,N+1):v[n]=max(v[n],v[n-r]+w)
    last=[0]+[next(r for r,w in ws.items() if r<=n and v[n]==v[n-r]+w) for n in range(1,N+1)]
    return dict(key=list(key),values=v,last=last)


def check_table(t,key):
    assert t['key']==list(key)
    v,last=t['values'],t['last'];assert len(v)==len(last)==61 and v[0]==last[0]==0
    ws={r:item(r,key) for r in Z.ITEMS};checks=0
    for n in range(1,61):
        for r,w in ws.items():
            if r<=n:assert v[n]>=v[n-r]+w;checks+=1
        r=last[n];assert r in ws and r<=n and v[n]==v[n-r]+ws[r]
    return checks


def finite(previous,saved=None):
    rep=Z.B.B.representative;types=[t for t in Z.B.B.types() if t[1]<=60]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    dom=[(cs,used,(*key,cs[7])) for cs,used,key in Z.large_domain()]
    keys=sorted({key for cs,used,key in dom});assert len(dom)==1044 and len(keys)==1042
    assert set(keys)=={tuple(t['key']) for t in previous['tables']}
    tabs={tuple(t['key']):t for t in saved['tables']} if saved else {k:table(k) for k in keys}
    assert set(tabs)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tabs[k],k) for k in keys)
    residual=[rep(('A' if r<=8 else 'D',r)) for r in Z.ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    before={tuple(r['counts']):r for r in previous['alternatives']};rows=[];items=0
    assert set(before)=={cs for cs,used,key in dom}
    for cs,used,key in dom:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert all(deletions(t)==0 for t in fixed) and sum(r for k,r in fixed)==used
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.ITEMS,residual):
            assert item(r,key)==upper_base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=N-used;tab=tabs[key];v=sum(upper_base(t) for t in fixed)+tab['values'][rem];forest=fixed[:]
        while rem:
            r=tab['last'][rem];forest.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==60 and value(forest,True)==v and value(forest)<=v
        assert v<=before[cs]['octupled_upper']
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            octupled_upper=v,octupled_actual=value(forest)))
    maximum=max(r['octupled_upper'] for r in rows);witness=[('A',15)]*4
    assert maximum==3232==value(witness)
    prior_six=[r for r in previous['alternatives'] if r['octupled_upper']>3240];assert len(prior_six)==6
    tail=[]
    for r in prior_six:
        forest=[(s[0],int(s[1:])) for s in r['forest']]
        assert S.value(forest)==3260 and sum(deletions(t) for t in forest)==4 and value(forest)==3228
        tail.append(dict(forest=r['forest'],old='815/2',new='807/2'))
    return dict(exact_maximum='404',curvature='405',finite_slack='1',complete_fixed_vectors=1044,
        interaction_keys=1042,simultaneous_type_pairs=len(types)**2,residual_item_checks=items,
        Bellman_upper_inequalities=checks,attaining_predecessors=60*len(keys),
        tables=[tabs[k] for k in keys],alternatives=rows,above_curvature=0,
        actual_attainer=['A15']*4,actual_attaining_cost='404',previous_six_obstructions=tail,
        actual_rank15_costs=dict(A15='101',B15='97',D15='201/2'),rank15_representative_upper='101',
        finite_cost_unit='1/8',all_cross_component_costs_unchanged=True,
        new_residual_tables=True,replay_optimizes=False,formal_forest_attainment_not_geometric_realization=True)


def ideal(previous):
    domain=sorted((low,)+p if low else p for low in range(9) for p in Z.high_parts(59-low))
    assert len(domain)==2113 and domain==sorted(tuple(r['ranks']) for r in previous['rows'])
    for r in previous['rows']:
        h=r['h_polynomial'];assert len(h)==60 and h==h[::-1] and all(c>0 for c in h)
        v=r['link_vertices'];assert sum(h)==v
        g=sum(min(i,59-i)*c for i,c in enumerate(h));assert g==r['G']
        curvature=F(5*g-59*v,2);assert curvature==F(r['curvature'])
        margin=curvature/v-F(r['normalized_cost']);assert margin==F(r['normalized_margin'])>0
        assert v*margin==F(r['absolute_margin'])
    summary=Z.summarize(previous['rows']);assert summary==previous['summary']
    assert summary['normalized_minimum']=='31/6000' and summary['absolute_minimum']=='54'
    return dict(complete_profiles=2113,rows=previous['rows'],summary=summary,ideal_weights_unchanged=True,
        inherited_order17_cost='229/2',convolutions_replayed=False,low9_concentration_only=True,
        positive_terminal_absolute_surplus_extends_to_all_actual_cusps=True)


def run(saved=None):
    prior=json.loads((ROOT/f'results/{S.TAG}.json').read_text())
    assert prior['status']=='PASS_SMALL_D16_COUPON_IDEAL17_REBALANCE'
    for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert sha256((ROOT/f'results/{S.TAG}.md').read_bytes()).hexdigest()==FINAL_MD
    assert prior['finite']['exact_maximum']=='815/2' and prior['finite']['D5_D12_actual_cost']=='45/4'
    paths=[Path(__file__),Path(S.__file__),Path(Z.__file__),ROOT/f'results/{S.TAG}.json',ROOT/f'results/{S.TAG}.md']
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths}
    if saved:assert sources==saved['sources']
    return dict(status='PASS_CONDITIONAL_DEPTH6_SPARSE_AMBIENT60_FOREST_FRONTIER',sources=sources,
        completed_unrestricted_dimension_bound=60,new_unrestricted_dimension_bound=None,
        ambient_dimension=60,link_dimension=59,local_payment_geometry_proved_here=False,
        conditional_local_rule='Delete finite distance7 pairs when both actual marks have feature depth at least6; preserve all other current sparse-D16/ideal17 weights.',
        old_all_rank_deleted_pair_formulas={k:f'max(0,n-{v})' for k,v in OLD.items()},
        conditional_all_rank_deleted_pair_formulas={k:f'max(0,n-{v})' for k,v in NEW.items()},
        finite=finite(prior['finite'],saved['finite'] if saved else None),ideal=ideal(prior['ideal']),
        strict_global_constant='4500',constant_counted_once=True,
        conditional_global_contradiction='0 >= 1*v_f + 54*v_i + 4500 > 0',
        geometric_hypothesis_and_full_payment_composition_required=True,
        scope='Complete changed finite forest certificate and unchanged ideal positive comparison, conditional on the separately proposed depth6 local-payment theorem. No new geometric proof or unconditional dimension bound is asserted here.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    if not args.verify:assert not path.exists(),'Frozen certificate exists; use --verify.'
    out=run(saved)
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'ideal',out['ideal']['summary'],flush=True)
