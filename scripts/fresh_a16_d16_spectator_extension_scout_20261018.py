#!/usr/bin/env python3
"""Separate coupled-A16 extension to an actual whole D16 spectator."""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
from itertools import combinations_with_replacement
import json
from pathlib import Path
import signal
import fresh_a17_inner_unit_rank16_rank17_guard_scout_20261018 as A

ROOT=A.ROOT;C=A.C;V=A.V;U=A.U
TAG='fresh_a16_d16_spectator_extension_scout_20261018'


def credit(fs):
    co=Counter(fs)
    return(64*co['A',1]*co['A',13]+(64*co['A',1]*(co['A',13]+co['D',13])if co['A',2]else 0)
        +(160*co['A',4]*co['A',12]+128*co['A',2]*co['A',14]if co['A',1]else 0)
        +128*co['A',4]*(co['A',13]+co['D',13])+64*co['A',15])


def value(fs):
    co=Counter(fs)
    if not V.has_guard(fs)or not co['D',16]:return A.value(fs)
    debit=-co['A',16]*(192+(80 if co['A',1]else 0))if co['A',16]<2 else 0
    return A.value(fs)+debit+(credit(fs)if not co['A',16]else 0)


def families():
    # The union of old-A16 and new-D16 positive families is evaluated under
    # the FINAL rule, including all later coupons and endpoint restorations.
    partners=[('A',12),('A',13),('D',13),('A',14),('A',15)]
    small=[t for t in C.TYPES if t[1]<=8 and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in small for b in small)
    rows=[];checks=preds=pair_checks=options=0
    for witness in(('A',16),('D',16)):
        for partner in partners:
            for s in range(18,42-partner[1]):
                rem=41-partner[1]-s
                extras=[None]+[t for t in C.TYPES if 9<=t[1]<=rem]
                for extra in extras:
                    fixed=[witness,partner,('D',s)]+([extra]if extra else[]);left=rem-(extra[1]if extra else 0);tabs=[];opts=[]
                    for a1 in range(left+1):
                        for a2 in range((left-a1)//2+1):
                            if partner==('A',12)and not a1 or partner==('A',14)and not(a1 and a2):continue
                            need4=partner==('A',12)or partner==('A',13)and not a1 or partner==('D',13)and not(a1 and a2)
                            cap=left-a1-2*a2-4*need4
                            if cap<0:continue
                            anchor=fixed+[('A',1)]*a1+[('A',2)]*a2+([('A',4)]if need4 else[]);bc=value(anchor)
                            assert credit(anchor)>0
                            ts=[t for t in small if t[1]<=cap];ws=[value(anchor+[t])-bc for t in ts]
                            tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b
                            for i,x in enumerate(ts):
                                for j in range(i,len(ts)):
                                    y=ts[j]
                                    if x[1]+y[1]<=cap:assert value(anchor+[x,y])==bc+ws[i]+ws[j];pair_checks+=1
                            tabs.append(dict(a1=a1,a2=a2,reserved_A4=bool(need4),table=tab))
                            if tab['values'][cap]is not None:opts.append((bc+tab['values'][cap],len(tabs)-1,cap));options+=1
                    if not opts:continue
                    z,ti,left=max(opts);t=tabs[ti];forest=fixed+[('A',1)]*t['a1']+[('A',2)]*t['a2']+([('A',4)]if t['reserved_A4']else[]);tab=t['table']
                    while left:x=tab['types'][tab['last'][left]];forest.append(x);left-=x[1]
                    assert value(forest)==z and sum(t[1]for t in forest)==57 and credit(forest)>0
                    rows.append(dict(witness=witness,partner=partner,reserved_D=s,extra_large=extra,tables=tabs,
                        maximum=F(z,64),attaining_forest=forest))
    return dict(rows=rows,Bellman_inequalities=checks,attaining_predecessors=preds,
        feasible_pair_additivity_checks=pair_checks,explicit_options=options,
        old_A16_family_maximum=max(r['maximum']for r in rows if r['witness']==('A',16)),
        new_D16_family_maximum=max(r['maximum']for r in rows if r['witness']==('D',16)))


def finite():
    old=V.load(A.R.TAG)['finite']['original_compact_class_maxima'];high=[r for r in old if F(r['maximum'])>F(1483,4)]
    assert len(high)==6 and all(sum(r['counts'][:7])==r['k5']==r['a12']==0 for r in high)
    original=[r for r in high if not r['counts'][8]];assert{r['counts'][7]for r in original}=={0,1,2}
    before=V.load(A.TAG)['finite']['classes'];assert max(F(r['maximum'])for r in before)<=F(1483,4)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS};res=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    assert all(C.actual_cross(a,b)==0 for a in res+[('A',1)]for b in res+[('A',1)])
    fibers=[('A',16),('B',16),('D',16)];out=[];checks=preds=options=0
    for k16 in(0,1,2):
        rows=[]
        for fixed0 in combinations_with_replacement(fibers,k16):
            fixed=list(fixed0);rem=57-16*k16;tabs=[];opts=[]
            for a1 in range(rem+1):
                for guard in(False,True):
                    if guard and rem-a1<18:continue
                    anchor=fixed+[('A',1)]*a1+([('D',18)]if guard else[]);cap=rem-a1-(18 if guard else 0);bc=value(anchor)
                    ts=[t for t in res if t[1]<=cap and(guard or not(t[0]=='D'and t[1]>=18))];ws=[value(anchor+[t])-bc for t in ts]
                    tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b;tabs.append(dict(a1=a1,guard=guard,table=tab))
                    for s in(range(18,rem-a1+1)if guard else[None]):
                        left=rem-a1-(s if guard else 0);z=tab['values'][left]
                        if z is None:continue
                        fs=fixed+[('A',1)]*a1+([('D',s)]if guard else[]);base=value(fs)
                        for t,w in zip(ts,ws):
                            if t[1]<=left:assert value(fs+[t])-base==w
                        opts.append((base+z,a1,guard,s,left,len(tabs)-1));options+=1
            z,a1,guard,s,left,ti=max(opts);forest=fixed+[('A',1)]*a1+([('D',s)]if guard else[]);tab=tabs[ti]['table']
            while left:t=tab['types'][tab['last'][left]];forest.append(t);left-=t[1]
            assert value(forest)==z and sum(t[1]for t in forest)==57
            rows.append(dict(fixed=fixed,tables=tabs,maximum=F(z,64),attaining_forest=forest))
        out.append(dict(fixed_rank16_count=k16,rows=rows,maximum=max(r['maximum']for r in rows)))
    witness=[('A',16),('D',17),('D',24)];assert value(witness)==64*F(1483,4)
    return dict(original_above_new_target_classes=high,classes=out,
        inherited_A17_class_uppers=[r['maximum']for r in before],Bellman_inequalities=checks,
        attaining_predecessors=preds,explicit_options=options,attaining_forest=witness,
        conditional_complete_new_maximum='1483/4')


def retention(physical):
    old=physical['current_later_interfaces'];literal=physical['current_literal_interface']
    assert old['later_corners']==564 and len(old['target_hits'])==20
    assert len(physical['local']['selected_and_later_seals'])==17
    assert all(c['A16_outputs']==0 for c in old['chosen_D12'])
    latest=V.load(A.TAG)['geometry']['consumer'];positives=0
    for c in latest['corners']:
        if not F(c['payment']):continue
        assert c['restore']and c['canonical_images']==[[0,14]]
        positives+=1
    d16=V.load(A.L.N.TAG)['geometry'];new_positive=0
    for p in d16['full']+d16['proper']:
        groups=[p.get('corners',[])]+[mp.get('corners',[])for mp in p.get('maps',[])]
        for cs in groups:
            for c in cs:
                if not F(c['payment']):continue
                assert c['coupon']and not c['selected_D16_zero']
                a,b=[next(b for b in c['components']if m in b['nodes'])for m in c['marks']]
                assert a is not b and{a['type'],b['type']}=={'A2','A13'}
                new_positive+=1
    assert positives==1 and new_positive>0
    return dict(inherited_later_corners=564,inherited_affected_A16_roles=20,
        inherited_sealed_selected_and_consumer_contexts=17,
        inherited_chosen_D12_A16_outputs=0,new_A17_payment_A15_endpoint=True,
        new_D16_and_R17_positive_outputs_are_A2_A13=new_positive,
        guard_extension_uses_the_same_full_core_E=True)


def run():
    fam=families();top=finite();covers=list(map(F,V.load(A.TAG)['sequential_increasing_family_uppers']))
    assert covers[2]==371;covers[2]=fam['old_A16_family_maximum'];covers.append(fam['new_D16_family_maximum'])
    assert max(covers)<=F(1483,4)and max(r['maximum']for r in top['classes'])<=F(1483,4)
    physical=V.load(A.S.TAG);ret=retention(physical)
    sources={Path(__file__),Path(A.__file__),Path(A.R.__file__),Path(A.S.__file__),
        ROOT/f'results/{A.TAG}.json',ROOT/f'results/{A.S.TAG}.json',ROOT/f'results/{A.L.TAG}.json'}
    return V.O.enc(dict(status='PASS_A16_D16_SPECTATOR_EXTENSION_SCOUT',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(sources)},
        rule='Extend the unmarked whole A16 witness in every coupled-A16 debit/coupon to an unmarked whole A16 orD16; retain D>=18.',
        complete_coupon_families=fam,finite=top,sequential_increasing_family_uppers=covers,
        unchanged_physical_seals=len(physical['local']['selected_and_later_seals']),
        current_retention=ret,geometric_proof_in_separate_report=True,
        independent_review_complete=False,universal_composition_complete=False,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);A.R.M.R.resources();U.memory_guard();p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();A.R.M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_coupon_families']['old_A16_family_maximum'],out['complete_coupon_families']['new_D16_family_maximum'],
        out['finite']['conditional_complete_new_maximum'],[(r['fixed_rank16_count'],r['maximum'])for r in out['finite']['classes']],flush=True)
