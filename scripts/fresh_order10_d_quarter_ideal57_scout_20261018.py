#!/usr/bin/env python3
"""Separate quarter-D10 ideal/profile and finite-family scout."""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
import json
from pathlib import Path
import signal
import fresh_a16_d16_spectator_extension_scout_20261018 as N

ROOT=N.ROOT;C=N.C;V=N.V;U=N.U
TAG='fresh_order10_d_quarter_ideal57_scout_20261018'
LOCAL='fresh_order10_d_opposite_fork_fractional_transfer_20261018'


def increment(fs):
    co=Counter(fs);return co['D',9]*(12*co['A',1]+24*sum(co['A',r]for r in range(2,9)))


def value(fs):return N.value(fs)+increment(fs)


def ideal():
    old=V.load(V.O.TAG)['ideal'];assert old['complete_profiles']==2270 and old['new_negative_count']==34
    rows=[]
    for r in old['rows']:
        saving=F(3*r['ranks'].count(9),80);margin=F(r['normalized_margin'])+saving
        rows.append(dict(ranks=r['ranks'],old_margin=r['normalized_margin'],normalized_saving=saving,normalized_margin=margin))
    bad=[r for r in rows if r['normalized_margin']<0];zero=[r for r in rows if not r['normalized_margin']]
    assert len(bad)==18 and not zero
    assert min(r['normalized_margin']for r in rows)==-F(17719,255552)
    return dict(complete_profiles=2270,rows=rows,negative_profiles=bad,zero_profiles=zero,
        old_negative_count=34,new_negative_count=18,new_order10_family_totals=dict(A='85/2',B='44',C='44',D='44'),
        old_rank9_envelope='71/16',new_rank9_envelope='22/5',normalized_saving_per_rank9_factor='3/80',
        worst_profiles=[r['ranks']for r in rows if r['normalized_margin']==-F(17719,255552)],
        normalized_minimum='-17719/255552',no_negative_surplus_extension=True)


def finite():
    rep=C.representative;upper_cross=C.C.upper_cross;extra={};pairs=0
    for a in C.TYPES:
        for b in C.TYPES:
            inc=increment([a,b]);need=max(0,C.actual_cross(a,b)+inc-upper_cross(rep(a),rep(b)))
            key=tuple(sorted((rep(a),rep(b))))
            if need:extra[key]=max(extra.get(key,0),need)
            pairs+=1
    assert extra=={tuple(sorted((t,('B',9)))):24 for t in(('A',5),('E',6),('A',7),('A',8))}
    rows=[];checks=preds=ordinary_options=quarter_options=0
    for cs,used,key in C.fixed_domain():
        if not cs[0]:continue
        for k5 in range((57-used)//5+1):
            fixed=C.fixed_components(cs,key)+[('A',5)]*k5;rem=57-used-5*k5
            fc=C.upper_value(fixed)+24*cs[0]*k5
            ws={r:C.item(r,key,k5)+(24*cs[0]if r in(6,7,8)else 0)for r in C.ITEMS}
            tabs=[]
            for no2 in(False,True):
                ranks=[r for r in C.ITEMS if r<=rem and(not no2 or r!=2)]
                ts=[rep(('A'if r<=8 else'D',r))for r in ranks]
                tab,a,b=C.additive(ts,[ws[r]for r in ranks],rem);checks+=a;preds+=b;tab['exclude2']=no2;tabs.append(tab)
            opts=[(tabs[0]['values'][rem],0,0,rem,False)]
            for a1 in(1,2):
                oldguard=a1==2 and bool(cs[1]);newguard=bool(cs[0]and k5)
                if rem<a1 or not(oldguard or newguard)or(a1==2 and not oldguard):continue
                for a2 in range((rem-a1)//2+1):
                    left=rem-a1-2*a2;v=tabs[1]['values'][left]
                    if v is not None:opts.append((a1*ws[1]+a2*ws[2]+v+(36*a2*(a2-1)if oldguard else 0)+(64*a2*k5 if newguard else 0),a1,a2,left,True))
            ordinary_options+=len(opts);best=max(z for z in opts if z[0]is not None)
            order17=16*k5*key[9]+(8*k5*rem if cs[2]else 0)
            ordinary=fc+best[0]+order17;quarter=None
            if cs[1]and rem>=7:
                for a1 in range(1,rem-5):
                    for a4 in range(1,(rem-a1-2)//4+1):
                        for a2 in range(1,(rem-a1-4*a4)//2+1):
                            left=rem-a1-4*a4-2*a2;v=tabs[1]['values'][left]
                            if v is None:continue
                            z=fc+a1*ws[1]+a2*ws[2]+a4*ws[4]+v+order17+(36*a2*(a2-1)if a1>=2 else 0)+(64*a2*k5 if cs[0]else 0)+64*a2*a4
                            quarter=z if quarter is None else max(quarter,z);quarter_options+=1
            rows.append(dict(counts=cs,k5=k5,a12=key[9],key=key,remaining_rank=rem,tables=tabs,
                ordinary_upper=F(ordinary,64),quarter_upper=None if quarter is None else F(quarter,64),
                upper=F(max(ordinary,quarter if quarter is not None else -1),64)))
    maximum=max(r['upper']for r in rows)
    assert len(rows)==1181 and maximum==F(729,2)
    current=V.load(N.TAG);families=current['complete_coupon_families']
    assert max(F(families[k])for k in('old_A16_family_maximum','new_D16_family_maximum'))==F(741,2)
    # D9 cannot fit the D16/R17/A17 new-positive families. Its only
    # coupled-A16 overlap has A1+A13+R16+D18+D9 (one extra A1/D9 pair).
    a16_overlap=F(741,2)+F(3,16)
    # Guarded-D17 positives contain D>=18 and A13 or A14. If x D9s
    # occur, at most26-9x total ranks can be charged A1..A8 factors.
    d17_extra=max(F(3*x*(26-9*x),16)for x in(1,2))
    d17_overlap=F(365)+d17_extra
    increasing=max(maximum,a16_overlap,d17_overlap);assert increasing==F(5931,16)<F(1483,4)
    witness=list(map(tuple,current['finite']['attaining_forest']));assert value(witness)==64*F(1483,4)
    return dict(ordered_actual_type_pairs=pairs,representative_pair_increments=[dict(types=k,increment=F(v,64))for k,v in sorted(extra.items())],
        classes=len(rows),rows=rows,Bellman_inequalities=checks,attaining_predecessors=preds,
        ordinary_context_options=ordinary_options,quarter_context_options=quarter_options,
        base_and_order17_D9_family_upper=maximum,
        quarter_classes=sum(r['quarter_upper']is not None for r in rows),
        maximizers=[dict(counts=r['counts'],k5=r['k5'],a12=r['a12'],upper=r['upper'])for r in rows if r['upper']==maximum],
        coupled_A16_later_overlap_upper=a16_overlap,guarded_D17_later_overlap_extra=d17_extra,
        guarded_D17_later_overlap_upper=d17_overlap,complete_new_positive_family_upper=increasing,
        complete_final_finite_maximum='1483/4',attaining_forest=witness,
        later_positive_family_overlap_proof_in_written_report=True)


def local():
    old=V.load(LOCAL)
    assert F(old['new_cost'])==F(355,8)
    assert F(old['delta'])==F(5,32)and 1-F(1,4)>=F(1,2)
    assert F(C.actual_base(('A',17)),64)==F(old['A17_new_cost'])-F(5,32)+F(1,2)
    weights={tuple(r['pair']):F(r['weight'])for r in old['weighted_pairs']}
    changed={(6,8),(6,9),(7,8),(7,9)}
    assert all(weights[p]==F(27,32)for p in changed)
    assert sum(F(3,4)if p in changed else w for p,w in weights.items())==44
    return dict(delta='1/4',old_delta='5/32',changed_ideal_pairs=[[6,8],[6,9],[7,8],[7,9]],
        ideal_old_weight='27/32',ideal_new_weight='3/4',affine_D_total='44',
        finite_increment_per_mark_pair='3/32',finite_A1_D9_increment='3/16',finite_A2_through_A8_D9_increment='3/8',
        required_finite_A17_endpoint='1/4',actual_retained_A17_endpoint='1/2',extra_A17_cost='0',
        inherited_geometry_parameter_free=True,all_weight_sum_inequalities_hold_for_delta_at_most_half=True)


def run():
    li=local();ii=ideal();ff=finite()
    paths={Path(__file__),Path(N.__file__),Path(C.__file__),ROOT/f'results/{N.TAG}.json',ROOT/f'results/{V.O.TAG}.json',ROOT/f'results/{LOCAL}.json',
        ROOT/f'scripts/{LOCAL}.py',ROOT/f'results/{LOCAL}.md',
        ROOT/'results/fresh_order10_low8_ambient58_frontier_20261018.json'}
    return V.O.enc(dict(status='PASS_ORDER10_D_QUARTER_IDEAL57_SCOUT',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},
        local=li,ideal=ii,finite=ff,universal_composition_complete=False,independent_review_complete=False,
        complete_final_finite_bound_proved=True,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);N.A.R.M.R.resources();U.memory_guard();p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();N.A.R.M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['ideal']['new_negative_count'],out['ideal']['normalized_minimum'],
        out['finite']['classes'],out['finite']['base_and_order17_D9_family_upper'],
        len(out['finite']['maximizers']),out['finite']['maximizers'][:1],flush=True)
