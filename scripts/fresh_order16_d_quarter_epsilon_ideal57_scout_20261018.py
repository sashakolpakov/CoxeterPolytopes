#!/usr/bin/env python3
"""Separate order16 sparse-D epsilon1/4 extension scout."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations_with_replacement,product
import fcntl
import json
from pathlib import Path
import signal
import fresh_order10_d_quarter_ideal57_scout_20261018 as Q

ROOT=Q.ROOT;C=Q.C;V=Q.V;U=Q.U
TAG='fresh_order16_d_quarter_epsilon_ideal57_scout_20261018'
LOCAL='fresh_order16_cd_sparse_sixteenth_transfer_20261018'


def increment(fs):
    co=Counter(fs);return 12*co['A',1]*co['D',15]


def value(fs):return Q.value(fs)+increment(fs)


def refine(rows,fibers,target):
    """Split every still costly class, retaining actual D15 and positive A1."""
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    residual=[t for t in C.TYPES if C.representative(t)in reps and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in residual for b in residual)
    classes=[];checks=preds=options=contexts=pairs=0
    for old in rows:
        if old['upper']<=target:continue
        assert old['counts'][6]==1 and old['counts'][1]==1 and old['k5']==0
        fixed=C.fixed_components(old['counts'],old['key']);groups=sorted(Counter(fixed).items())
        fixedsets=sorted(tuple(sorted(t for group in choice for t in group))for choice in product(*(
            combinations_with_replacement(fibers[t],n)for t,n in groups)))
        fixedsets=[fs for fs in fixedsets if('D',15)in fs];rr=[];rem=old['remaining_rank']
        for fs in fixedsets:
            tabs=[];best=None
            for a1 in range(1,rem+1):
                for a2 in range((rem-a1)//2+1):
                    anchor=list(fs)+[('A',1)]*a1+[('A',2)]*a2;cap=rem-a1-2*a2
                    bc=value(anchor);assert bc==C.value(anchor)+increment(anchor)
                    ts=[t for t in residual if t[1]<=cap]
                    ws=[value(anchor+[t])-bc for t in ts]
                    # Once actual rank10, positive A1 and full A2 counts are
                    # fixed, every remaining old context is linear in A4.
                    for t,w in zip(ts,ws):
                        assert value(anchor+[t])==C.value(anchor+[t])+increment(anchor+[t])
                    contexts+=len(ts)
                    if a1==1:
                        for i,a in enumerate(ts):
                            for j in range(i,len(ts)):
                                b=ts[j]
                                if a[1]+b[1]<=cap:
                                    assert value(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
                    tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b
                    tab.update(a1=a1,a2=a2);tabs.append(tab)
                    z=tab['values'][cap]
                    if z is not None:
                        options+=1
                        if best is None or bc+z>best[0]:best=(bc+z,anchor,tab)
            assert best is not None
            z,att,tab=best;att=list(att);left=len(tab['values'])-1
            while left:
                t=tab['types'][tab['last'][left]];att.append(t);left-=t[1]
            assert sum(t[1]for t in att)==57 and value(att)==z
            rr.append(dict(fixed_actual=fs,tables=tabs,maximum=F(z,64),attaining_forest=sorted(att)))
        classes.append(dict(counts=old['counts'],k5=old['k5'],a12=old['a12'],remaining_rank=rem,rows=rr,
            maximum=max(r['maximum']for r in rr)))
    assert len(classes)==4
    return dict(classes=classes,actual_fixed_fibers=sum(len(c['rows'])for c in classes),
        Bellman_inequalities=checks,attaining_predecessors=preds,explicit_options=options,
        singleton_context_checks=contexts,feasible_pair_additivity_checks=pairs,
        actual_refined_maximum=max(c['maximum']for c in classes))


def local():
    old=V.load(LOCAL);d=next(r for r in old['rules']if r['kind']=='D');changed=set(map(tuple,d['changed_pairs']))
    assert changed=={(1,14),(1,15),(10,12),(10,13)}
    weights={tuple(r['pair']):F(r['weight'])for r in d['weighted_pairs']};assert sum(weights.values())==F(411,4)
    assert all(weights[p]==F(7,16)for p in changed)
    new={p:F(1,4)if p in changed else w for p,w in weights.items()};assert sum(new.values())==102
    assert all(new[tuple(p)]==weights[tuple(p)]for p in d['retained_basic_pairs'])
    assert all(new[tuple(p)]==1 for p in d['retained_triangle_pairs'])
    counts=Counter();payments=[]
    for c in d['cases']:
        for r in c['completions']:counts[r['proof']['outcome']]+=1
        for p in c['current_payments']:
            cs=[]
            for z in p['corners']:
                if z['kind']=='affineD15_order16':pay=F(1,4)
                else:pay=sum((F(t['payment'])for t in z['existing_coupon']),F())+F(len(z['new_coupon']),4)
                cs.append(dict(marks=z['marks'],kind=z['kind'],payment=pay))
            total=sum(z['payment']for z in cs);assert total>=1
            assert total==(1 if c['case']in(2368,2512)else F(7,4))
            payments.append(dict(case=c['case'],corners=cs,total=total))
    assert counts==dict(excluded_projection_above_norm=32,excluded_lorentz_index=16,sealed_completion=32)
    return dict(changed_pairs=sorted(changed),old_weight='7/16',new_weight='1/4',old_epsilon='1/16',new_epsilon='1/4',
        new_family_totals=dict(A='80',B='102',C='101',D='102'),new_rank15_envelope='51/8',
        finite_A1_D15_node11_increment='3/16',new_node11_coupon='1/4',new_actual_A1_D15_total='7/4',
        complete_selected_patterns=394,complete_affected_maps=80,outcomes=dict(counts),whole_map_payments=payments,
        target_basic_and_triangle_pairs_unchanged=True,arbitrary_support_seals_inherited_unchanged=True)


def ideal():
    old=V.load(Q.TAG)['ideal'];assert old['new_negative_count']==18;rows=[]
    for r in old['rows']:
        saving=F(3*r['ranks'].count(15),64);m=F(r['normalized_margin'])+saving
        rows.append(dict(ranks=r['ranks'],old_margin=r['normalized_margin'],normalized_saving=saving,normalized_margin=m))
    bad=[r for r in rows if r['normalized_margin']<0];zero=[r for r in rows if not r['normalized_margin']]
    return dict(complete_profiles=len(rows),rows=rows,negative_profiles=bad,zero_profiles=zero,
        old_negative_count=18,new_negative_count=len(bad),normalized_minimum=min(r['normalized_margin']for r in rows),
        normalized_saving_per_rank15_factor='3/64',no_negative_surplus_extension=True)


def finite():
    rep=C.representative;rows=[];checks=preds=ordinary_options=quarter_options=0
    fibers={}
    for t in C.TYPES:fibers.setdefault(rep(t),[]).append(t)
    f15=fibers['D',15];assert f15==[('A',15),('B',15),('D',15)]
    base15=max(C.actual_base(t)+(32 if t in(('A',15),('D',15))else 0)for t in f15)
    base_gap=C.C.upper_base(('D',15))-base15;assert base_gap>=0
    pair_gaps={}
    for partner,ts in fibers.items():
        typed=max(C.actual_cross(a,b)+Q.increment([a,b])+increment([a,b])for a in f15 for b in ts)
        old=C.C.upper_cross(('D',15),partner)+(12 if partner==('A',1)else 0)
        assert typed<=old;pair_gaps[partner]=old-typed
    extra={}
    for a in C.TYPES:
        for b in C.TYPES:
            need=max(0,C.actual_cross(a,b)+Q.increment([a,b])+increment([a,b])-C.C.upper_cross(rep(a),rep(b)))
            if need:
                pair=tuple(sorted((rep(a),rep(b))));extra[pair]=max(extra.get(pair,0),need)
    expected={tuple(sorted((t,('B',9)))):24 for t in(('A',5),('E',6),('A',7),('A',8))}
    expected[tuple(sorted((('A',1),('D',15))))]=12;assert extra==expected
    for cs,used,key in C.fixed_domain():
        if not cs[6]:continue
        for k5 in range((57-used)//5+1):
            fixed=C.fixed_components(cs,key)+[('A',5)]*k5;rem=57-used-5*k5
            pair_saving=sum(pair_gaps[b]if a==('D',15)else pair_gaps[a]if b==('D',15)else 0 for i,a in enumerate(fixed)for b in fixed[i+1:])
            fc=C.upper_value(fixed)+24*cs[0]*k5-base_gap*cs[6]-pair_saving
            ws={r:C.item(r,key,k5)+(12*cs[6]if r==1 else 24*cs[0]if r in(6,7,8)else 0)-cs[6]*pair_gaps[rep(('A'if r<=8 else'D',r))]for r in C.ITEMS};tabs=[]
            for no2 in(False,True):
                ranks=[r for r in C.ITEMS if r<=rem and(not no2 or r!=2)]
                ts=[rep(('A'if r<=8 else'D',r))for r in ranks];tab,a,b=C.additive(ts,[ws[r]for r in ranks],rem)
                checks+=a;preds+=b;tab['exclude2']=no2;tabs.append(tab)
            opts=[tabs[0]['values'][rem]]
            for a1 in(1,2):
                og=a1==2 and bool(cs[1]);ng=bool(cs[0]and k5)
                if rem<a1 or not(og or ng)or(a1==2 and not og):continue
                for a2 in range((rem-a1)//2+1):
                    v=tabs[1]['values'][rem-a1-2*a2]
                    if v is not None:opts.append(a1*ws[1]+a2*ws[2]+v+(36*a2*(a2-1)if og else 0)+(64*a2*k5 if ng else 0))
            ordinary_options+=len(opts);order17=16*k5*key[9]+(8*k5*rem if cs[2]else 0)
            ordinary=fc+max(v for v in opts if v is not None)+order17;quarter=None
            if cs[1]and rem>=7:
                for a1 in range(1,rem-5):
                    for a4 in range(1,(rem-a1-2)//4+1):
                        for a2 in range(1,(rem-a1-4*a4)//2+1):
                            v=tabs[1]['values'][rem-a1-4*a4-2*a2]
                            if v is None:continue
                            z=fc+a1*ws[1]+a2*ws[2]+a4*ws[4]+v+order17+(36*a2*(a2-1)if a1>=2 else 0)+(64*a2*k5 if cs[0]else 0)+64*a2*a4
                            quarter=z if quarter is None else max(quarter,z);quarter_options+=1
            rows.append(dict(counts=cs,k5=k5,a12=key[9],key=key,remaining_rank=rem,tables=tabs,
                ordinary_upper=F(ordinary,64),quarter_upper=None if quarter is None else F(quarter,64),
                upper=F(max(ordinary,quarter if quarter is not None else -1),64)))
    coarse_maximum=max(r['upper']for r in rows);target=F(1483,4);refined=refine(rows,fibers,target)
    maximum=max(max(r['upper']for r in rows if r['upper']<=target),refined['actual_refined_maximum'])
    guarded_D17_overlap=F(365)+F(33,16)
    return dict(rows=rows,classes=len(rows),Bellman_inequalities=checks,attaining_predecessors=preds,
        ordinary_context_options=ordinary_options,quarter_context_options=quarter_options,
        base_order17_and_quarterD10_D15_family_upper=maximum,coarse_upper_before_actual_refinement=coarse_maximum,
        actual_refinement=refined,guarded_D17_overlap_upper=guarded_D17_overlap,
        rank15_actual_fiber=f15,rank15_typed_base_with_context_padding=F(base15,64),
        rank15_base_saving_per_component=F(base_gap,64),rank15_cross_savings=[dict(partner=k,saving=F(v,64))for k,v in sorted(pair_gaps.items())if v],
        all_other_later_positive_families_disjoint=True,complete_new_positive_family_upper=max(maximum,guarded_D17_overlap),
        coarse_maximizers=[dict(counts=r['counts'],k5=r['k5'],upper=r['upper'])for r in rows if r['upper']==coarse_maximum])


def run():
    li=local();ii=ideal();ff=finite();prior=V.load(Q.TAG)['finite'];previous=F(prior['complete_final_finite_maximum'])
    witness=list(map(tuple,prior['attaining_forest']));assert value(witness)==64*previous
    sources={Path(__file__),Path(Q.__file__),Path(C.__file__),ROOT/f'results/{Q.TAG}.json',ROOT/f'results/{LOCAL}.json',
        ROOT/f'scripts/{LOCAL}.py',ROOT/f'results/{LOCAL}.md'}
    return V.O.enc(dict(status='PASS_ORDER16_D_QUARTER_EPSILON_IDEAL57_SCOUT',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(sources)},
        local=li,ideal=ii,finite=ff,previous_finite_maximum=previous,unchanged_formal_attainer=witness,
        conditional_complete_finite_upper=max(previous,ff['complete_new_positive_family_upper']),
        universal_composition_complete=False,independent_review_complete=False,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);Q.N.A.R.M.R.resources();U.memory_guard();p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();Q.N.A.R.M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['ideal']['new_negative_count'],out['ideal']['normalized_minimum'],out['finite']['classes'],
        out['finite']['base_order17_and_quarterD10_D15_family_upper'],out['conditional_complete_finite_upper'],flush=True)
