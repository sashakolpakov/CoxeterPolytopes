#!/usr/bin/env python3
"""Saved-witness necessary payment constraints, not a new forest optimum."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
from pathlib import Path
import fresh_d12_type_split_ambient58_forest_frontier_20261018 as C
ROOT=C.ROOT;TAG='fresh_ambient58_obstruction_payment_sensitivities_20261018';enc=C.enc
BASE_PINS={
 'scripts/'+C.TAG+'.py':'e99a08189a5c04a2bacac470121ffee458bcb5729984ac0ed5d58096810e0c04',
 'results/'+C.TAG+'.json':'9e02ab6b26f0cc2e779eb8bac820e8272794c51afe0d65c6da55e051e822a70d',
 'results/'+C.TAG+'.md':'4e6bf11e259af23023e43a10f7525ec160d00b0e3266d55ccd35295b7c30e7f8'}

def name(t):return t[0]+str(t[1])
def display(fs):return ' + '.join((str(n)+' 'if n>1 else'')+name(t)for t,n in sorted(Counter(fs).items()))
def features(fs):
    co=Counter(fs);out={};definitions={}
    for t,n in co.items():
        if t[1]>=9:
            key='base:'+name(t);out[key]=n
            definitions[key]=dict(kind='component_base',types=[t],current_component_cost=F(C.actual_base(t),32))
    for a,b in combinations_with_replacement(sorted(co),2):
        n=co[a]*co[b]if a!=b else co[a]*(co[a]-1)//2
        if n and C.actual_cross(a,b):
            key='coupon:'+name(a)+'|'+name(b);out[key]=n
            definitions[key]=dict(kind='ordinary_component_pair_total',types=[a,b],current_component_cost=F(C.actual_cross(a,b),32))
    # Grouped alternatives expose ties hidden by selecting one type label.
    # They are evaluated separately; overlapping columns are not a proposed
    # simultaneous payment rule.
    n=co['A',15]+co['D',15]
    if n:
        out['base_group:A15_D15']=n
        definitions['base_group:A15_D15']=dict(kind='common_component_base_decrease',types=[('A',15),('D',15)],current_component_cost=F(100))
    n=co['A',2]*(co['A',13]+co['D',13])
    if n:
        assert C.actual_cross(('A',2),('A',13))==C.actual_cross(('A',2),('D',13))==32*F(39,8)
        key='coupon_group:A2|A13_D13';out[key]=n
        definitions[key]=dict(kind='common_component_pair_total_decrease',types=[('A',2),('A',13),('D',13)],current_component_cost=F(39,8))
    n15=co['D',15];marked=0
    for rank in (3,4,5):
        small=[t for t in C.TYPES if t[1]==rank]
        assert all(C.actual_cross(t,('D',15))==32*F(3*rank,2)for t in small)
        n=sum(co[t]for t in small)*n15
        if n:
            key=f'coupon_group:rank{rank}|D15';out[key]=n
            definitions[key]=dict(kind='common_component_pair_total_decrease',types=small+[('D',15)],current_component_cost=F(3*rank,2))
            marked+=rank*n
    if marked:
        key='proportional_coupon:rank3to5|D15';out[key]=marked
        definitions[key]=dict(kind='per_small_mark_decrease_across_rank3_to_5',types=[t for t in C.TYPES if 3<=t[1]<=5]+[('D',15)],
            current_component_cost=F(3,2),suggested_probe=F(1,8))
    return out,definitions

def run():
    for p,h in BASE_PINS.items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    doc=json.loads((ROOT/f'results/{C.TAG}.json').read_text());classes=doc['actual_typed_refinement']['classes']
    assert doc['finite']['certified_upper']=='1527/4'
    allrows=[];maxrows=[];defs={};classrows=[]
    for ci,c in enumerate(classes):
        if F(c['actual_maximum'])<=377:continue
        ids=[];mx=[]
        for ri,r in enumerate(c['rows']):
            if r['thirtyseconds_actual_maximum']<=32*377:continue
            fs=tuple(map(tuple,r['attaining_forest']));cost=F(C.value(fs),32)
            assert sum(t[1]for t in fs)==58 and cost==F(r['thirtyseconds_actual_maximum'],32)
            # The current above-curvature attaining witnesses have no active
            # contextual charge. This is not an inventory of every forest.
            assert C.context(fs)==0
            coeff,dd=features(fs)
            for k,v in dd.items():assert k not in defs or defs[k]==v;defs[k]=v
            row=dict(class_index=ci,fixed_type_row=ri,forest=fs,display=display(fs),cost=cost,
                excess=cost-377,coefficients=coeff,contextual_charge='0')
            ids.append(len(allrows));allrows.append(row)
            if cost==F(c['actual_maximum']):mx.append(len(allrows)-1);maxrows.append(row)
        assert mx and len(mx)==len(c['maximizers'])
        classrows.append(dict(class_index=ci,actual_maximum=c['actual_maximum'],witness_indices=ids,maximizing_indices=mx,
            representative_constraint=allrows[mx[0]]))
    assert len(classrows)==23 and len(allrows)==114 and len(maxrows)==64
    ranked=[]
    for key,d in sorted(defs.items()):
        hit=[r for r in allrows if r['coefficients'].get(key,0)]
        probe=d.get('suggested_probe',min(F(1),d['current_component_cost']))
        after=[r['cost']-probe*r['coefficients'].get(key,0)for r in allrows]
        classes_all_max_hit=sum(all(allrows[i]['coefficients'].get(key,0)for i in c['maximizing_indices'])for c in classrows)
        survivors=sum(any(after[i]>377 for i in c['witness_indices'])for c in classrows)
        max_after=max(after);threshold=max(r['excess']/r['coefficients'][key]for r in hit)
        ranked.append(dict(feature=key,**d,affected_saved_witnesses=len(hit),classes_all_maximizers_affected=classes_all_max_hit,
            necessary_decrease_to_clear_all_affected_saved_witnesses=threshold,
            admissible_cost_probe=probe,remaining_saved_witness_maximum=max_after,
            classes_with_a_remaining_saved_counterexample=survivors,
            witness_at_remaining_maximum=allrows[after.index(max_after)]['display'],
            single_feature_can_clear_all_saved_witnesses=len(hit)==len(allrows)and threshold<=d['current_component_cost']))
    ranked.sort(key=lambda r:(r['classes_with_a_remaining_saved_counterexample'],r['remaining_saved_witness_maximum'],-r['classes_all_maximizers_affected'],r['feature']))
    assert not any(r['single_feature_can_clear_all_saved_witnesses']for r in ranked)
    important=[]
    for c in sorted(classrows,key=lambda c:-F(c['actual_maximum'])):
        r=c['representative_constraint']
        important.append(dict(display=r['display'],excess=r['excess'],coefficients=r['coefficients'],
            single_feature_necessary_thresholds={k:r['excess']/v for k,v in r['coefficients'].items()}))
    paths={f'scripts/{TAG}.py',*BASE_PINS}
    return enc(dict(status='PASS_SAVED_AMBIENT58_PAYMENT_SENSITIVITIES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
        baseline=C.TAG,completed_unrestricted_bound=58,new_dimension_bound=None,
        actual_obstruction_classes=23,saved_above_curvature_witnesses=114,class_maximizing_witnesses=64,
        witnesses=allrows,class_constraints=important,candidate_ranking=ranked,
        no_single_feature_clears_all_saved_witnesses=True,
        feature_units='One unit of same-component base cost or total ordinary component-pair coupon cost; not automatically one marked role.',
        grouped_features_are_separate_overlapping_parameterizations=True,
        necessary_constraints_only=True,conditional_geometric_payment_changes_unproved=True,
        all_ordinary_and_contextual_costs_in_baseline=True,
        no_new_DP_no_local_geometry_no_optimizer=True,
        no_sufficiency_or_reoptimized_maximum_claim=True,
        scope='Exact lower obstructions and sensitivities from frozen actual attaining witnesses. Clearing these witnesses is necessary, not sufficient; unrecorded forests may become maximizers after a change. A proposed coupon compensation requires its own complete count.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    p=ROOT/f'results/{TAG}.json';assert args.verify or not p.exists();out=run()
    if args.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['actual_obstruction_classes'],out['saved_above_curvature_witnesses'],len(out['candidate_ranking']),flush=True)
