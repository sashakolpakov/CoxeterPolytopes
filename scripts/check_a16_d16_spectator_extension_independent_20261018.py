#!/usr/bin/env python3
"""Independent final-cost/domain checks for the D16-spectator extension."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations_with_replacement
import signal
import check_a17_inner_unit_guard_independent_20261018 as Y

ROOT=Y.ROOT;W=Y.W;C=Y.C;OLD=Y.OLD;K=Y.X.K
TAG='fresh_a16_d16_spectator_extension_scout_20261018'
PHYSICAL='fresh_a16_double_large_d_physical_scout_20261018'


def credit(co):
    return(64*co['A',1]*co['A',13]+(64*co['A',1]*(co['A',13]+co['D',13])if co['A',2]else 0)
        +(160*co['A',4]*co['A',12]+128*co['A',2]*co['A',14]if co['A',1]else 0)
        +128*co['A',4]*(co['A',13]+co['D',13])+64*co['A',15])


def value(fs):
    fs=list(map(tuple,fs));co=Counter(fs);a,d=co['A',16],co['D',16]
    g=any(k=='D'and n>=18 for k,n in fs)
    debit=a*(192+80*bool(co['A',1]))if g and d and a<2 else 0
    inc=credit(co)if g and d and not a else 0
    return Y.value(fs)-debit+inc


def pair_check(anchor,ts,ws,cap):
    base=value(anchor);n=0
    for i,a in enumerate(ts):
        for j in range(i,len(ts)):
            b=ts[j]
            if a[1]+b[1]<=cap:assert value(anchor+[a,b])==base+ws[i]+ws[j];n+=1
    return n


def families(doc):
    fam=doc['complete_coupon_families'];partners={('A',12):5,('A',13):1,('D',13):3,('A',14):3,('A',15):0}
    small=[t for t in C.TYPES if t[1]<=8 and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in small for b in small)
    expected=set()
    for w in(('A',16),('D',16)):
        for p,m in partners.items():
            for s in range(18,42-p[1]):
                rem=41-p[1]-s
                for e in[None]+[t for t in C.TYPES if 9<=t[1]<=rem]:
                    if rem-(e[1]if e else 0)>=m:expected.add((w,p,s,e))
    actual={(tuple(r['witness']),tuple(r['partner']),r['reserved_D'],tuple(r['extra_large'])if r['extra_large']else None)for r in fam['rows']}
    assert actual==expected and len(actual)==len(fam['rows'])
    checks=preds=pairs=options=0
    for r in fam['rows']:
        w,p=tuple(r['witness']),tuple(r['partner']);s=r['reserved_D'];e=tuple(r['extra_large'])if r['extra_large']else None
        fixed=[w,p,('D',s)]+([e]if e else[]);left=41-p[1]-s-(e[1]if e else 0)
        expected_tables=set()
        for a1 in range(left+1):
            for a2 in range((left-a1)//2+1):
                if p==('A',12)and not a1 or p==('A',14)and not(a1 and a2):continue
                need4=p==('A',12)or p==('A',13)and not a1 or p==('D',13)and not(a1 and a2)
                if left-a1-2*a2-4*need4>=0:expected_tables.add((a1,a2,need4))
        assert{(t['a1'],t['a2'],t['reserved_A4'])for t in r['tables']}==expected_tables
        best=None
        for t in r['tables']:
            a1,a2,need4=t['a1'],t['a2'],t['reserved_A4'];anchor=fixed+[('A',1)]*a1+[('A',2)]*a2+([('A',4)]if need4 else[])
            cap=left-a1-2*a2-4*need4;bc=value(anchor);assert credit(Counter(anchor))>0
            ts=[x for x in small if x[1]<=cap];ws=[value(anchor+[x])-bc for x in ts];tab=t['table']
            assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
            pairs+=pair_check(anchor,ts,ws,cap);a,b=K.bellman(tab);checks+=a;preds+=b
            if tab['values'][cap]is not None:z=bc+tab['values'][cap];best=z if best is None else max(best,z);options+=1
        assert F(best,64)==F(r['maximum'])and value(r['attaining_forest'])==best
    assert(checks,preds,pairs,options)==(fam['Bellman_inequalities'],fam['attaining_predecessors'],fam['feasible_pair_additivity_checks'],fam['explicit_options'])
    assert max(F(r['maximum'])for r in fam['rows']if r['witness']==['A',16])==F(1481,4)==F(fam['old_A16_family_maximum'])
    assert max(F(r['maximum'])for r in fam['rows']if r['witness']==['D',16])==F(741,2)==F(fam['new_D16_family_maximum'])
    return len(actual),checks,preds,pairs,options


def finite(doc):
    top=doc['finite'];old=W.load(Y.X.TAG)['finite']['original_compact_class_maxima'];high=[r for r in old if F(r['maximum'])>F(1483,4)]
    assert high==top['original_above_new_target_classes']and len(high)==6
    before=W.load(Y.TAG)['finite']['classes'];assert list(map(F,top['inherited_A17_class_uppers']))==[F(r['maximum'])for r in before]and max(F(r['maximum'])for r in before)<=F(1483,4)
    assert all(sum(r['counts'][:7])==r['k5']==r['a12']==0 for r in high)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS};res=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    assert all(C.actual_cross(a,b)==0 for a in res+[('A',1)]for b in res+[('A',1)])
    fibers=[t for t in C.TYPES if C.representative(t)==('D',16)];assert fibers==[('A',16),('B',16),('D',16)]
    checks=preds=options=pairs=0
    for cl in top['classes']:
        k=cl['fixed_rank16_count'];assert{tuple(map(tuple,r['fixed']))for r in cl['rows']}==set(combinations_with_replacement(fibers,k))
        for r in cl['rows']:
            fixed=list(map(tuple,r['fixed']));rem=57-16*k;best=None
            assert{(t['a1'],t['guard'])for t in r['tables']}=={(a,g)for a in range(rem+1)for g in(False,True)if not g or rem-a>=18}
            for t in r['tables']:
                a1=t['a1'];g=t['guard'];anchor=fixed+[('A',1)]*a1+([('D',18)]if g else[]);cap=rem-a1-(18 if g else 0)
                ts=[x for x in res if x[1]<=cap and(g or not(x[0]=='D'and x[1]>=18))];ws=[value(anchor+[x])-value(anchor)for x in ts];tab=t['table']
                assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
                # A1 multiplicity contributes only linear terms after its
                # presence test is fixed; test all pairs at both presence anchors.
                if a1<=1:pairs+=pair_check(anchor,ts,ws,cap)
                a,b=K.bellman(tab);checks+=a;preds+=b
                for s in(range(18,rem-a1+1)if g else[None]):
                    n=rem-a1-(s if g else 0);z=tab['values'][n]
                    if z is None:continue
                    fs=fixed+[('A',1)]*a1+([('D',s)]if g else[]);bc=value(fs)
                    for x,w in zip(ts,ws):
                        if x[1]<=n:assert value(fs+[x])-bc==w
                    z+=bc;best=z if best is None else max(best,z);options+=1
            assert F(best,64)==F(r['maximum'])and value(r['attaining_forest'])==best
        assert max(F(r['maximum'])for r in cl['rows'])==F(cl['maximum'])
    assert(checks,preds,options)==(top['Bellman_inequalities'],top['attaining_predecessors'],top['explicit_options'])
    assert[F(r['maximum'])for r in top['classes']]==[F(370),F(1483,4),F(370)]
    assert value(top['attaining_forest'])==64*F(1483,4)and top['conditional_complete_new_maximum']=='1483/4'
    return checks,preds,options,pairs


def main():
    doc=W.load(TAG);hashes=0
    for p,h in doc['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    phys=W.load(PHYSICAL);pd=W.geometry(phys);ret=doc['current_retention']
    assert ret['inherited_later_corners']==564 and ret['inherited_affected_A16_roles']==20 and ret['inherited_sealed_selected_and_consumer_contexts']==17
    assert ret['new_A17_payment_A15_endpoint']and ret['new_D16_and_R17_positive_outputs_are_A2_A13']>0
    a=families(doc);b=finite(doc)
    assert list(map(F,doc['sequential_increasing_family_uppers']))==[F(739,2),F(365),F(1481,4),F(11513,32),F(22917,64),F(1481,4),F(741,2)]
    OLD.U.memory_guard();print('PASS_INDEPENDENT_A16_D16_SPECTATOR_EXTENSION',hashes,pd,*a,*b,'1483/4',flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
