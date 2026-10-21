#!/usr/bin/env python3
"""Check R17 certificates without importing the new R17 producer."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations_with_replacement
import signal
import check_d16_node3_guard_independent_20261018 as K

ROOT=K.ROOT;C=K.C;W=K.W;OLD=K.OLD
TAG='fresh_d16_rank17_witness_extension_scout_20261018'


def value(fs):
    fs=list(map(tuple,fs));co=Counter(fs)
    r16=co['A',16]+co['D',16];r17=co['A',17]+co['D',17]
    g=any(k=='D'and n>=18 for k,n in fs)
    debit=64*co['D',16]if g and(r16>=2 or r16 and r17)else 0
    credit=256*co['A',2]*co['A',13]if g and co['A',1]and(r16 or r17)else 0
    return OLD.value(fs)-debit+credit


def coupon(doc):
    up=doc['coupon_upper'];cs=(0,0,0,0,1,0,0,0,1)
    key=next(key for counts,used,key in C.fixed_domain()if counts==cs and key[9]==0)
    assert C.actual_base(('D',17))-C.upper_value([('A',17)])==16
    assert all(C.actual_cross(t,('D',17))==0 for t in C.TYPES if t[1]<=6 or t==('A',13))
    assert{(r['k5'],r['reserved_D'])for r in up['rows']}=={(k,s)for k in(0,1)for s in range(18,25-5*k)}
    checks=preds=0
    for row in up['tables']:
        k=row['k5'];tab=row['table'];ws={r:C.item(r,key,k)+(256 if r==2 else 128 if r==4 else 0)for r in C.ITEMS}
        assert dict(tab['weights'])==ws
        v=tab['values'];last=tab['last'];assert len(v)==len(last)==58 and v[0]==0 and last[0]is None
        for n in range(1,58):
            opts=[r for r in ws if r<=n and v[n-r]is not None]
            for r in opts:assert v[n]>=v[n-r]+ws[r];checks+=1
            if opts:r=last[n];assert r in opts and v[n]==v[n-r]+ws[r];preds+=1
            else:assert v[n]is None and last[n]is None
        fixed=C.fixed_components(cs,key)+[('A',5)]*k;fc=C.upper_value(fixed)+16
        for r in[x for x in up['rows']if x['k5']==k]:
            s=r['reserved_D'];n=27-5*k-s;assert n==r['residual']
            opts=[]
            for a1 in range(1,n-1):
                for a2 in range(1,(n-a1)//2+1):
                    rest=n-a1-2*a2
                    if v[rest]is not None:opts.append(fc+C.actual_base(('D',s))+a1*ws[1]+a2*ws[2]+v[rest])
            assert F(max(opts),64)==F(r['upper'])
    assert checks==up['Bellman_inequalities']and preds==up['attaining_predecessors']
    assert max(F(r['upper'])for r in up['rows'])==F(22917,64)
    return checks,preds


def finite(doc):
    top=doc['finite'];old=top['original_compact_class_maxima'];assert len(old)==38
    high=[r for r in old if F(r['maximum'])>F(743,2)]
    assert len(high)==5 and all(r['k5']==r['a12']==0 for r in high)
    assert{tuple(r['counts'])for r in high}=={tuple([0]*7+[k,0])for k in(0,1,2)}|{tuple([0]*7+[k,1])for k in(0,1)}
    assert all(r['counts'][4]==0 for r in high)
    previous=W.load(K.J);assert max(F(r['maximum'])for r in previous['top_classes']['classes'])==371
    assert top['newly_refined_old_classes']==[r for r in high if r['counts'][8]]
    assert C.representative(('A',17))==('A',17)and C.representative(('D',17))==('D',17)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    res=[t for t in C.TYPES if C.representative(t)in reps and t!=('A',1)]
    fibers=[t for t in C.TYPES if C.representative(t)==('D',16)]
    assert fibers==[('A',16),('B',16),('D',16)]and top['actual_rank16_fiber']==list(map(list,fibers))
    assert ('A',17)not in res and ('D',17)in res
    assert all(C.actual_cross(a,b)==0 for a in res+[('A',1)]for b in res+[('A',1)])
    checks=preds=options=pair_checks=0
    for cl in top['classes']:
        k=cl['fixed_rank16_count']
        assert{tuple(map(tuple,r['fixed']))for r in cl['rows']}=={(('A',17),)+f for f in combinations_with_replacement(fibers,k)}
        for row in cl['rows']:
            fixed=list(map(tuple,row['fixed']));rem=40-16*k;best=None
            expected={(a,g)for a in range(rem+1)for g in(False,True)if not g or rem-a>=18}
            assert{(r['a1'],r['guard'])for r in row['tables']}==expected
            for tt in row['tables']:
                a1=tt['a1'];g=tt['guard'];anchor=fixed+[('A',1)]*a1+([('D',18)]if g else[])
                cap=rem-a1-(18 if g else 0);bc=value(anchor)
                ts=[t for t in res if t[1]<=cap and(g or not(t[0]=='D'and t[1]>=18))]
                ws=[value(anchor+[t])-bc for t in ts];tab=tt['table']
                assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
                # The written fixed-test polynomial argument proves arbitrary
                # multiplicities; these independently check every feasible pair.
                for i,a in enumerate(ts):
                    for j in range(i,len(ts)):
                        b=ts[j]
                        if a[1]+b[1]<=cap:
                            assert value(anchor+[a,b])==bc+ws[i]+ws[j];pair_checks+=1
                a,b=K.bellman(tab);checks+=a;preds+=b
                for s in(range(18,rem-a1+1)if g else[None]):
                    n=rem-a1-(s if g else 0);v=tab['values'][n]
                    if v is None:continue
                    chosen=fixed+[('A',1)]*a1+([('D',s)]if g else[]);z=value(chosen)+v
                    for x,w in zip(ts,ws):
                        if x[1]<=n:assert value(chosen+[x])-value(chosen)==w
                    best=z if best is None else max(best,z);options+=1
            assert F(best,64)==F(row['maximum'])and value(row['attaining_forest'])==best
        assert max(F(r['maximum'])for r in cl['rows'])==F(cl['maximum'])
    assert(checks,preds,options)==(top['Bellman_inequalities'],top['attaining_predecessors'],top['explicit_options'])
    assert [F(c['maximum'])for c in top['classes']]==[F(1483,4),F(743,2)]
    assert top['conditional_complete_new_maximum']=='743/2'and value(top['attaining_forest'])==64*F(743,2)
    return checks,preds,options,pair_checks


def main():
    doc=W.load(TAG);hashes=0
    for p,h in doc['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    a=coupon(doc);b=finite(doc)
    covers=list(map(F,doc['sequential_increasing_family_uppers']))
    assert covers==[F(739,2),F(365),F(371),F(11513,32),F(22917,64)]and max(covers)<F(743,2)
    OLD.U.memory_guard();print('PASS_INDEPENDENT_D16_RANK17_WITNESS',hashes,*a,*b,'743/2',flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
