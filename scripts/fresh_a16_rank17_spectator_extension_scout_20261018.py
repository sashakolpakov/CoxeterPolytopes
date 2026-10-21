#!/usr/bin/env python3
"""Whole A17/D17 spectator extension of the sealed coupled A16 rule."""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations_with_replacement
import json
from pathlib import Path
import signal
import fresh_order17_d13_shortleaf_seven_thirtyseconds_scout_20261018 as H

ROOT=H.ROOT;C=H.C;V=H.V;U=H.U;N=H.E.Q.N;A=N.A
TAG='fresh_a16_rank17_spectator_extension_scout_20261018'
TARGET=F(11861,32)


def delta(fs):
    co=Counter(fs);r16=co['A',16]+co['D',16];r17=co['A',17]+co['D',17]
    if not V.has_guard(fs)or not r17:return 0
    debit=-co['A',16]*(192+80*bool(co['A',1]))if co['A',16]==1 and not co['D',16]else 0
    return debit+(N.credit(fs)if not r16 else 0)


def value(fs):return H.value(fs)+delta(fs)


def geometry():
    old=V.load(A.S.TAG);rows=[];hits=0
    for p in old['local']['selected_and_later_seals']:
        found=[]
        for variant in p['variants']:
            for c in variant['corners']:
                if F(c['delta'])>=0:continue
                a,b=[next(b for b in c['components']if m in b['nodes'])for m in c['marks']]
                assert a is b and a['type']=='A16'and len(a['nodes'])==16
                outside=set(p['white'])|set(c['marks']);outside-=set(a['nodes'])
                assert len(outside)<=1
                found.append(dict(external_A1=variant['external_A1'],marks=c['marks'],
                    displayed_outside_A16=sorted(outside),payment_before=c['current_before']))
        assert found;hits+=len(found)
        rows.append(dict(index=p['index'],seal_method=p['seal']['method'],A16_debits=found))
    assert len(rows)==17
    ret=N.retention(old)
    return dict(seals=rows,sealed_contexts=17,actual_debit_variants=hits,inherited_retention=ret,
        spectator_types=['A16','D16','A17','D17'],large_D_unchanged=True,
        proof='At each A16 debit corner the displayed remainder has at most one root. The full inherited core seal makes the extra-white forest orthogonal. Any distinct whole rank17 witness and D>=18 therefore lie wholly in the same extra forest at all replacement corners.')


def families():
    partners=[('A',12),('A',13),('D',13),('A',14),('A',15)]
    small=[t for t in C.TYPES if t[1]<=8 and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in small for b in small)
    rows=[];checks=preds=pairs=options=0
    for witness in(('A',17),('D',17)):
        for partner in partners:
            for s in range(18,41-partner[1]):
                left=40-partner[1]-s;tabs=[];opts=[]
                # A rank9 extra cannot coexist with any required small role.
                for a1 in range(left+1):
                    for a2 in range((left-a1)//2+1):
                        if partner==('A',12)and not a1 or partner==('A',14)and not(a1 and a2):continue
                        need4=partner==('A',12)or partner==('A',13)and not a1 or partner==('D',13)and not(a1 and a2)
                        cap=left-a1-2*a2-4*need4
                        if cap<0:continue
                        assert cap<=8
                        anchor=[witness,partner,('D',s)]+[('A',1)]*a1+[('A',2)]*a2+([('A',4)]if need4 else[])
                        assert N.credit(anchor)>0
                        bc=value(anchor);ts=[t for t in small if t[1]<=cap];ws=[value(anchor+[t])-bc for t in ts]
                        for i,a in enumerate(ts):
                            for j in range(i,len(ts)):
                                b=ts[j]
                                if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
                        tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b
                        tabs.append(dict(a1=a1,a2=a2,reserved_A4=bool(need4),table=tab))
                        if tab['values'][cap]is not None:opts.append((bc+tab['values'][cap],len(tabs)-1,anchor));options+=1
                if not opts:continue
                z,ti,att=max(opts);att=list(att);tab=tabs[ti]['table'];cap=len(tab['values'])-1
                while cap:t=tab['types'][tab['last'][cap]];att.append(t);cap-=t[1]
                assert value(att)==z and sum(t[1]for t in att)==57
                rows.append(dict(witness=witness,partner=partner,D_rank=s,tables=tabs,maximum=F(z,64),attaining_forest=sorted(att)))
    return dict(rows=rows,complete_cases=len(rows),Bellman_inequalities=checks,attaining_predecessors=preds,
        feasible_pair_additivity_checks=pairs,explicit_options=options,
        A17_witness_maximum=max(r['maximum']for r in rows if r['witness']==('A',17)),
        D17_witness_maximum=max(r['maximum']for r in rows if r['witness']==('D',17)),
        complete_positive_family_maximum=max(r['maximum']for r in rows))


def top():
    old=V.load(A.R.TAG)['finite']['original_compact_class_maxima']
    high=[r for r in old if F(r['maximum'])>TARGET]
    assert len(old)==38 and len(high)==6
    assert all(sum(r['counts'][:7])==r['k5']==r['a12']==0 for r in high)
    previous=V.load(A.TAG)['finite']['classes'];assert max(F(r['maximum'])for r in previous)<=F(741,2)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    res=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    assert ('A',17)not in res and ('D',17)in res
    assert all(C.actual_cross(a,b)==0 for a in res+[('A',1)]for b in res+[('A',1)])
    out=[];checks=preds=options=pairs=0
    for k16 in(0,1,2):
        rows=[]
        for fixed0 in combinations_with_replacement([('A',16),('B',16),('D',16)],k16):
            fixed=list(fixed0);rem=57-16*k16;tabs=[];opts=[]
            for a1 in range(rem+1):
                for guard in(False,True):
                    if guard and rem-a1<18:continue
                    anchor=fixed+[('A',1)]*a1+([('D',18)]if guard else[]);cap=rem-a1-(18 if guard else 0);bc=value(anchor)
                    ts=[t for t in res if t[1]<=cap and(guard or not(t[0]=='D'and t[1]>=18))];ws=[value(anchor+[t])-bc for t in ts]
                    if a1<=1:
                        for i,a in enumerate(ts):
                            for j in range(i,len(ts)):
                                b=ts[j]
                                if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
                    tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b;tabs.append(dict(a1=a1,guard=guard,table=tab))
                    for s in(range(18,rem-a1+1)if guard else[None]):
                        left=rem-a1-(s if guard else 0);z=tab['values'][left]
                        if z is None:continue
                        fs=fixed+[('A',1)]*a1+([('D',s)]if guard else[]);base=value(fs)
                        for t,w in zip(ts,ws):
                            if t[1]<=left:assert value(fs+[t])-base==w
                        opts.append((base+z,a1,guard,s,left,len(tabs)-1));options+=1
            z,a1,guard,s,left,ti=max(opts);att=fixed+[('A',1)]*a1+([('D',s)]if guard else[]);tab=tabs[ti]['table']
            while left:t=tab['types'][tab['last'][left]];att.append(t);left-=t[1]
            assert value(att)==z and sum(t[1]for t in att)==57
            rows.append(dict(fixed=fixed,tables=tabs,maximum=F(z,64),attaining_forest=att))
        out.append(dict(fixed_rank16_count=k16,rows=rows,maximum=max(r['maximum']for r in rows)))
    return dict(original_high_classes=high,classes=out,inherited_A17_class_uppers=[r['maximum']for r in previous],
        Bellman_inequalities=checks,attaining_predecessors=preds,explicit_options=options,
        feasible_pair_additivity_checks=pairs,maximum=max(r['maximum']for r in out))


def complete():
    t=top();f=families();old=V.load(N.TAG);covers=list(map(F,old['sequential_increasing_family_uppers']))
    assert max(covers)==F(741,2)
    d10=V.load(H.E.Q.TAG)['finite'];examples=[]
    for w in(('A',16),('D',16)):
        fs=[('A',1),('A',13),w,('D',18),('D',9)]
        assert sum(t[1]for t in fs)==57
        examples.append(dict(forest=fs,maximum=F(value(fs),64)))
    # These are the ONLY coupled-A16/D16 overlap forests for D10: their
    # minimum ranks already sum57, and all other partners need more ranks.
    d10upper=max(F(d10['base_and_order17_D9_family_upper']),F(d10['guarded_D17_later_overlap_upper']),max(r['maximum']for r in examples))
    eps=F(V.load(H.E.TAG)['finite']['complete_new_positive_family_upper'])
    current=V.load(H.TAG);last=max(F(current['complete_base_increase_family']['complete_base_family_upper']),F(current['later_positive_overlaps']['maximum']))
    allbounds=covers+[t['maximum'],max(map(F,t['inherited_A17_class_uppers'])),f['complete_positive_family_maximum'],d10upper,eps,last]
    assert max(allbounds)==TARGET
    att=[('D',13)]*3+[('A',4)]*4+[('A',2)]
    assert sum(t[1]for t in att)==57 and value(att)==64*TARGET
    return dict(top=t,new_positive_families=f,D10_exact_coupled_overlap_forests=examples,D10_positive_family_upper=d10upper,
        inherited_sequential_family_uppers=covers,epsilon_positive_family_upper=eps,D17_distal_positive_family_upper=last,
        complete_finite_maximum=TARGET,attaining_forest=att,
        former_attainer=dict(forest=[('A',16),('D',17),('D',24)],old=F(1483,4),new=F(value([('A',16),('D',17),('D',24)]),64)))


def run():
    g=geometry();f=complete()
    sources={Path(__file__),Path(H.__file__),Path(N.__file__),ROOT/f'results/{H.TAG}.json',ROOT/f'results/{N.TAG}.json',
        ROOT/f'results/{A.S.TAG}.json',ROOT/f'results/{H.E.TAG}.json',ROOT/f'results/{H.E.Q.TAG}.json'}
    return V.O.enc(dict(status='PASS_A16_RANK17_SPECTATOR_EXTENSION_SCOUT',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(sources)},
        geometry=g,finite=f,ideal_weights_unchanged=True,geometric_proof_in_separate_report=True,
        universal_composition_complete=False,independent_review_complete=False,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);N.A.R.M.R.resources();U.memory_guard()
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();N.A.R.M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['finite']['complete_finite_maximum'],out['finite']['new_positive_families']['complete_positive_family_maximum'],flush=True)
