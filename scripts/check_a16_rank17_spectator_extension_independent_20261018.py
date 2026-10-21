#!/usr/bin/env python3
"""Independent arithmetic and full inherited seals for R17 A16 extension."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations_with_replacement
import signal
import check_order17_d13_shortleaf_seven_thirtyseconds_independent_20261018 as I
import check_a16_double_large_d_independent_20261018 as G

ROOT=I.ROOT;W=I.W;C=I.C;K=I.K;OLD=I.OLD
TAG='fresh_a16_rank17_spectator_extension_scout_20261018'
PHYSICAL='fresh_a16_double_large_d_physical_scout_20261018'
PREVIOUS='fresh_a16_d16_spectator_extension_scout_20261018'
R17='fresh_d16_rank17_witness_extension_scout_20261018'
A17='fresh_a17_inner_unit_rank16_rank17_guard_scout_20261018'
H='fresh_order17_d13_shortleaf_seven_thirtyseconds_scout_20261018'
Q='fresh_order10_d_quarter_ideal57_scout_20261018'
EPS='fresh_order16_d_quarter_epsilon_ideal57_scout_20261018'
TARGET=F(11861,32)


def credit(co):
    return(64*co['A',1]*co['A',13]+(64*co['A',1]*(co['A',13]+co['D',13])if co['A',2]else 0)
        +(160*co['A',4]*co['A',12]+128*co['A',2]*co['A',14]if co['A',1]else 0)
        +128*co['A',4]*(co['A',13]+co['D',13])+64*co['A',15])


def value(fs):
    fs=list(map(tuple,fs));co=Counter(fs);r16=co['A',16]+co['D',16];r17=co['A',17]+co['D',17]
    d=0
    if r17 and any(k=='D'and n>=18 for k,n in fs):
        if co['A',16]==1 and not co['D',16]:d-=192+80*bool(co['A',1])
        if not r16:d+=credit(co)
    return I.value(fs)+d


def geometry(doc):
    old=W.load(PHYSICAL);pd=G.geometry(old);rows=[];hits=0
    for p in old['local']['selected_and_later_seals']:
        found=[];g=[list(map(F,row))for row in p['gram']]
        for var in p['variants']:
            for c in var['corners']:
                if F(c['delta'])>=0:continue
                bs=G.blocks(g,sorted(set(p['white'])|set(c['marks'])))
                a,b=[next(b for b in bs if m in b['nodes'])for m in c['marks']]
                assert a is b and a['type']=='A16'and len(a['nodes'])==16
                outside=(set(p['white'])|set(c['marks']))-set(a['nodes']);assert len(outside)<=1
                found.append(dict(external_A1=var['external_A1'],marks=c['marks'],displayed_outside_A16=sorted(outside),payment_before=c['current_before']))
        assert found;hits+=len(found);rows.append(dict(index=p['index'],seal_method=p['seal']['method'],A16_debits=found))
    assert rows==doc['geometry']['seals']and len(rows)==17==doc['geometry']['sealed_contexts']
    assert hits==doc['geometry']['actual_debit_variants']
    assert doc['geometry']['inherited_retention']==W.load(PREVIOUS)['current_retention']
    return pd,len(rows),hits


def families(doc):
    ff=doc['finite']['new_positive_families'];partners={('A',12):5,('A',13):1,('D',13):3,('A',14):3,('A',15):0}
    expected={(w,p,s)for w in(('A',17),('D',17))for p,m in partners.items()for s in range(18,41-p[1]-m)}
    assert expected=={(tuple(r['witness']),tuple(r['partner']),r['D_rank'])for r in ff['rows']}
    assert len(expected)==ff['complete_cases']
    small=[t for t in C.TYPES if t[1]<=8 and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in small for b in small)
    checks=preds=pairs=options=0
    for r in ff['rows']:
        w,p=tuple(r['witness']),tuple(r['partner']);s=r['D_rank'];left=40-p[1]-s;expected_tables=set();best=None
        for a1 in range(left+1):
            for a2 in range((left-a1)//2+1):
                if p==('A',12)and not a1 or p==('A',14)and not(a1 and a2):continue
                need4=p==('A',12)or p==('A',13)and not a1 or p==('D',13)and not(a1 and a2)
                if left-a1-2*a2-4*need4>=0:expected_tables.add((a1,a2,bool(need4)))
        assert expected_tables=={(t['a1'],t['a2'],t['reserved_A4'])for t in r['tables']}
        for t in r['tables']:
            a1,a2,need4=t['a1'],t['a2'],t['reserved_A4'];cap=left-a1-2*a2-4*need4;assert cap<=8
            anchor=[w,p,('D',s)]+[('A',1)]*a1+[('A',2)]*a2+([('A',4)]if need4 else[])
            assert credit(Counter(anchor))>0;bc=value(anchor);ts=[x for x in small if x[1]<=cap];ws=[value(anchor+[x])-bc for x in ts];tab=t['table']
            assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
            for i,a in enumerate(ts):
                for j in range(i,len(ts)):
                    b=ts[j]
                    if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
            a,b=K.bellman(tab);checks+=a;preds+=b;z=tab['values'][cap]
            if z is not None:options+=1;best=bc+z if best is None else max(best,bc+z)
        assert F(best,64)==F(r['maximum'])and value(r['attaining_forest'])==best
    assert(checks,preds,pairs,options)==(ff['Bellman_inequalities'],ff['attaining_predecessors'],ff['feasible_pair_additivity_checks'],ff['explicit_options'])
    assert max(F(r['maximum'])for r in ff['rows']if r['witness']==['A',17])==F(ff['A17_witness_maximum'])
    assert max(F(r['maximum'])for r in ff['rows']if r['witness']==['D',17])==F(ff['D17_witness_maximum'])
    assert max(F(r['maximum'])for r in ff['rows'])==F(ff['complete_positive_family_maximum'])
    return len(expected),checks,preds,pairs,options


def paircheck(anchor,ts,ws,cap):
    n=0;bc=value(anchor)
    for i,a in enumerate(ts):
        for j in range(i,len(ts)):
            b=ts[j]
            if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j];n+=1
    return n


def top(doc):
    ff=doc['finite']['top'];old=W.load(R17)['finite']['original_compact_class_maxima']
    high=[r for r in old if F(r['maximum'])>TARGET];assert len(old)==38 and len(high)==6 and high==ff['original_high_classes']
    assert all(sum(r['counts'][:7])==r['k5']==r['a12']==0 for r in high)
    before=W.load(A17)['finite']['classes'];assert list(map(F,ff['inherited_A17_class_uppers']))==[F(r['maximum'])for r in before]
    assert max(map(F,ff['inherited_A17_class_uppers']))<=F(741,2)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    res=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    assert ('A',17)not in res and('D',17)in res
    assert all(C.actual_cross(a,b)==0 for a in res+[('A',1)]for b in res+[('A',1)])
    checks=preds=options=pairs=0
    assert{r['fixed_rank16_count']for r in ff['classes']}=={0,1,2}
    for cl in ff['classes']:
        k16=cl['fixed_rank16_count'];domain=set(combinations_with_replacement([('A',16),('B',16),('D',16)],k16))
        assert domain=={tuple(map(tuple,r['fixed']))for r in cl['rows']}
        for r in cl['rows']:
            fixed=list(map(tuple,r['fixed']));rem=57-16*k16;best=None
            assert{(t['a1'],t['guard'])for t in r['tables']}=={(a,g)for a in range(rem+1)for g in(False,True)if not g or rem-a>=18}
            for t in r['tables']:
                a1,g=t['a1'],t['guard'];anchor=fixed+[('A',1)]*a1+([('D',18)]if g else[]);cap=rem-a1-(18 if g else 0)
                ts=[x for x in res if x[1]<=cap and(g or not(x[0]=='D'and x[1]>=18))];ws=[value(anchor+[x])-value(anchor)for x in ts];tab=t['table']
                assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
                if a1<=1:pairs+=paircheck(anchor,ts,ws,cap)
                a,b=K.bellman(tab);checks+=a;preds+=b
                for s in(range(18,rem-a1+1)if g else[None]):
                    left=rem-a1-(s if g else 0);z=tab['values'][left]
                    if z is None:continue
                    fs=fixed+[('A',1)]*a1+([('D',s)]if g else[]);bc=value(fs)
                    for x,w in zip(ts,ws):
                        if x[1]<=left:assert value(fs+[x])-bc==w
                    options+=1;best=bc+z if best is None else max(best,bc+z)
            assert F(best,64)==F(r['maximum'])and value(r['attaining_forest'])==best
        assert max(F(r['maximum'])for r in cl['rows'])==F(cl['maximum'])
    assert(checks,preds,options,pairs)==(ff['Bellman_inequalities'],ff['attaining_predecessors'],ff['explicit_options'],ff['feasible_pair_additivity_checks'])
    assert max(F(r['maximum'])for r in ff['classes'])==F(ff['maximum'])
    return checks,preds,options,pairs


def complete(d):
    f=d['finite'];old=W.load(PREVIOUS);covers=list(map(F,old['sequential_increasing_family_uppers']))
    assert max(covers)==F(741,2)and covers==list(map(F,f['inherited_sequential_family_uppers']))
    assert{tuple(map(tuple,r['forest']))for r in f['D10_exact_coupled_overlap_forests']}=={(('A',1),('A',13),w,('D',18),('D',9))for w in(('A',16),('D',16))}
    for r in f['D10_exact_coupled_overlap_forests']:assert F(value(r['forest']),64)==F(r['maximum'])
    q=W.load(Q)['finite'];qupper=max(F(q['base_and_order17_D9_family_upper']),F(q['guarded_D17_later_overlap_upper']),max(F(r['maximum'])for r in f['D10_exact_coupled_overlap_forests']))
    assert qupper==F(f['D10_positive_family_upper'])
    eps=F(W.load(EPS)['finite']['complete_new_positive_family_upper']);assert eps==F(f['epsilon_positive_family_upper'])
    h=W.load(H);last=max(F(h['complete_base_increase_family']['complete_base_family_upper']),F(h['later_positive_overlaps']['maximum']))
    assert last==F(f['D17_distal_positive_family_upper'])
    upper=max(covers+[F(f['top']['maximum']),max(map(F,f['top']['inherited_A17_class_uppers'])),F(f['new_positive_families']['complete_positive_family_maximum']),qupper,eps,last])
    assert upper==TARGET==F(f['complete_finite_maximum'])and value(f['attaining_forest'])==64*TARGET
    assert F(value(f['former_attainer']['forest']),64)==F(f['former_attainer']['new'])<F(f['former_attainer']['old'])


def main():
    d=W.load(TAG);hashes=0
    for p,h in d['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    a=geometry(d);b=families(d);c=top(d);complete(d);OLD.U.memory_guard()
    print('PASS_INDEPENDENT_A16_RANK17_SPECTATOR_EXTENSION',hashes,*a,*b,*c,TARGET,flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
