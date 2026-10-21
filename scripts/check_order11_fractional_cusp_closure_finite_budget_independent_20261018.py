#!/usr/bin/env python3
"""Independent exact family cover for the fractional order11 finite rule."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations_with_replacement, product
from math import prod
import signal
import fresh_a16_rank17_spectator_extension_scout_20261018 as R

ROOT=R.ROOT;C=R.C;H=R.H;V=R.V;U=R.U
TAG='fresh_order11_fractional_cusp_closure_finite_budget_20261018'
TARGET=F(11861,32)


def extra(fs):
    a=Counter(fs)
    return 22*a['A',1]*a['D',10]+44*sum(a['A',r]for r in(2,3,4))*a['D',10]+22*a['A',2]*a['B',10]+44*(2*a['A',2]+a['A',3])*a['B',8]


@lru_cache(maxsize=8192)
def ordered_cost(counts):
    fs=[t for t,n in counts for _ in range(n)]
    twice=sum(H.O.actual_cross(a,b)*(na*nb-(na if a==b else 0))for a,na in counts for b,nb in counts)
    assert twice%2==0
    ordinary=sum(H.O.actual_base(t)*n for t,n in counts)+twice//2+H.O.context(fs)
    return 8*(ordinary+H.E.Q.increment(fs)+H.E.increment(fs)+H.increment(fs))+extra(fs)


def old_upper(fs,guard=None):
    z=ordered_cost(tuple(sorted(Counter(fs).items())))
    return z+(8*V.context_delta(fs)if guard is not None else 0)


def dp(ts,weights,cap):
    values=[None]*(cap+1);values[0]=0;checks=preds=0
    for rank in range(1,cap+1):
        cand=[values[rank-t[1]]+w for t,w in zip(ts,weights)if t[1]<=rank and values[rank-t[1]]is not None]
        checks+=len(cand)
        if cand:values[rank]=max(cand);preds+=1
    return values,checks,preds


def coarse(doc):
    rows=doc['coarse']['rows'];position=0;high=[];counts=[0,0,0,0];lowmax=F()
    for a in C.TYPES:
        for b in C.TYPES:
            rp={C.representative(a),C.representative(b)};pair={a,b};padding=0
            if ('D',10)in rp:padding+=22 if ('A',1)in rp else 44 if rp&{('A',2),('A',3),('A',4)}else 0
            if rp=={('A',2),('A',8)}:padding+=88
            if rp=={('A',3),('A',8)}:padding+=44
            if rp=={('A',1),('D',15)}:padding+=96
            if rp in({('A',4),('D',13)},{('A',2),('D',14)}):padding+=1024
            if ('B',9)in rp and rp&{('A',5),('E',6),('A',7),('A',8)}:padding+=192
            charge=8*C.actual_cross(a,b)+8*(H.E.Q.increment([a,b])+H.E.increment([a,b])+H.increment([a,b]))+extra([a,b])
            if pair in({('A',4),('A',13)},{('A',2),('A',14)}):charge+=1024
            assert charge<=8*C.C.upper_cross(C.representative(a),C.representative(b))+padding
    for cs,used,key in C.fixed_domain():
        for k5 in range((57-used)//5+1):
            for k8 in range((57-used-5*k5)//8+1):
                if cs[1]==0 and k8==0:continue
                row=rows[position];position+=1
                rem=57-used-5*k5-8*k8
                assert row['counts']==list(cs)and row['key']==list(key)and(row['k5'],row['k8'],row['remaining_rank'],row['a12'])==(k5,k8,rem,key[9])
                fixed=C.fixed_components(cs,key)+[('A',5)]*k5+[('A',8)]*k8
                constant=8*(C.upper_value(fixed)+24*cs[0]*(k5+k8))
                weights={r:8*C.item(r,key,k5)for r in C.ITEMS if r!=8}
                weights[1]+=22*cs[1]+96*cs[6]
                weights[2]+=44*cs[1]+88*k8+1024*cs[5]
                weights[3]+=44*cs[1]+44*k8
                weights[4]+=44*cs[1]+1024*cs[4]
                for r in(6,7):weights[r]+=192*cs[0]
                tabs=[]
                for forbid2 in(False,True):
                    ranks=[r for r in weights if r<=rem and(not forbid2 or r!=2)]
                    ts=[C.representative(('A'if r<=8 else'D',r))for r in ranks]
                    tab,a,b=dp(ts,[weights[r]for r in ranks],rem);tabs.append(tab);counts[0]+=a;counts[1]+=b
                candidates=[tabs[0][rem]]
                for a1 in(1,2):
                    quadratic=a1==2 and bool(cs[1]);b9=bool(cs[0]and k5)
                    if rem<a1 or not(quadratic or b9)or(a1==2 and not quadratic):continue
                    for a2 in range((rem-a1)//2+1):
                        z=tabs[1][rem-a1-2*a2]
                        if z is not None:candidates.append(a1*weights[1]+a2*weights[2]+z+(288*a2*(a2-1)if quadratic else 0)+(512*a2*k5 if b9 else 0))
                counts[2]+=len(candidates);o17=128*k5*key[9]+(64*k5*rem if cs[2]else 0)
                ordinary=constant+max(v for v in candidates if v is not None)+o17;quarter=None
                if cs[1]and rem>=7:
                    for a1 in range(1,rem-5):
                        for a4 in range(1,(rem-a1-2)//4+1):
                            for a2 in range(1,(rem-a1-4*a4)//2+1):
                                z=tabs[1][rem-a1-4*a4-2*a2]
                                if z is None:continue
                                q=constant+a1*weights[1]+a2*weights[2]+a4*weights[4]+z+o17+(288*a2*(a2-1)if a1>=2 else 0)+(512*a2*k5 if cs[0]else 0)+512*a2*a4
                                quarter=q if quarter is None else max(quarter,q);counts[3]+=1
                upper=F(max(ordinary,quarter if quarter is not None else -1),512)
                assert F(row['ordinary_upper'])==F(ordinary,512)and(row['quarter_upper']is None if quarter is None else F(row['quarter_upper'])==F(quarter,512))and F(row['upper'])==upper
                if upper>TARGET:high.append(row)
                else:lowmax=max(lowmax,upper)
    assert position==len(rows)==doc['coarse']['complete_classes']
    assert counts==[doc['coarse'][k]for k in('Bellman_inequalities','attaining_predecessors','ordinary_context_options','quarter_context_options')]
    return high,lowmax


def actual(doc,high):
    fibers={}
    for t in C.TYPES:fibers.setdefault(C.representative(t),[]).append(t)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS if r!=8}
    residual=[t for t in C.TYPES if C.representative(t)in reps and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in residual+[('A',1),('A',2)]for b in residual+[('A',1),('A',2)])
    classes=doc['actual']['classes'];position=0;fibercount=0;globalmax=F();checks=preds=options=singles=reserved=0
    for coarse_row in high:
        U.memory_guard();fixed=C.fixed_components(coarse_row['counts'],coarse_row['key'])+[('A',5)]*coarse_row['k5']+[('A',8)]*coarse_row['k8']
        groups=sorted(Counter(fixed).items());rem=coarse_row['remaining_rank']
        domain=sorted({tuple(sorted(t for part in parts for t in part))for parts in product(*(combinations_with_replacement(fibers[t],n)for t,n in groups))})
        domain=[fs for fs in domain if set(fs)&{('B',8),('B',10),('D',10)}]
        if not domain:continue
        entry=classes[position];position+=1
        assert all(entry[k]==coarse_row[k]for k in('counts','k5','k8','a12'))and len(entry['rows'])==len(domain)
        classmax=F()
        for fs,row in zip(domain,entry['rows']):
            assert [list(t)for t in fs]==row['fixed_actual'];co=Counter(fs);cache={};best=None;fibercount+=1
            assert old_upper(list(fs))==8*H.basevalue(list(fs))+extra(list(fs))
            raw={t:old_upper(list(fs)+[t])-old_upper(list(fs))-8*(C.context(list(fs)+[t])-C.context(list(fs)))for t in residual}
            for a1 in range(rem+1):
                for a2 in range((rem-a1)//2+1):
                    branches=[(False,0)]if a1 or a2 or not co['D',14]else[(True,0),(False,1)]
                    for forbid3,reserve3 in branches:
                        left=rem-a1-2*a2-3*reserve3
                        if left<0:continue
                        anchor=list(fs)+[('A',1)]*a1+[('A',2)]*a2+[('A',3)]*reserve3
                        for guard in((False,True)if co['A',13]or co['A',14]else[None]):
                            if guard and left<18:continue
                            cap=left-(18 if guard else 0);prefix=anchor+([('D',18)]if guard else[])
                            ts=[t for t in residual if t[1]<=cap and not(forbid3 and t==('A',3))and(guard is not False or not(t[0]=='D'and t[1]>=18))]
                            weights=[raw[t]+(8*(32*co['A',12]+64*a2*bool(co['A',10]))if a1 and t==('A',4)else 0)+(1024*co['A',13]if guard and t==('A',4)else -512 if guard and t==('D',17)else 0)for t in ts]
                            prefixbase=old_upper(prefix,guard)
                            for t,w in zip(ts,weights):assert old_upper(prefix+[t],guard)-prefixbase==w
                            singles+=len(ts);key=tuple(ts),tuple(weights),cap
                            if key not in cache:
                                tab,a,b=dp(ts,weights,cap);cache[key]=tab;checks+=a;preds+=b
                            tab=cache[key]
                            for ds in(range(18,left+1)if guard else[None]):
                                n=left-(ds if guard else 0)
                                if tab[n]is None:continue
                                aa=anchor+([('D',ds)]if guard else[]);b0=old_upper(aa,guard)
                                for t,w in zip(ts,weights):
                                    if t[1]<=n:assert old_upper(aa+[t],guard)-b0==w;reserved+=1
                                options+=1;z=b0+tab[n];best=z if best is None else max(best,z)
            assert best is not None and F(best,512)==F(row['maximum'])
            att=[tuple(t)for t in row['attaining_forest']]
            assert old_upper(att)==8*H.basevalue(att)+extra(att)
            assert sum(t[1]for t in att)==57 and F(extra(att),512)==F(row['added_increment'])and F(8*R.value(att)+extra(att),512)==F(row['actual_current'])
            assert F(row['actual_current'])<=F(row['maximum'])
            classmax=max(classmax,F(best,512))
        assert classmax==F(entry['maximum']);globalmax=max(globalmax,classmax)
    assert position==len(classes)and fibercount==doc['actual']['actual_fixed_fibers']
    assert [checks,preds,options,singles,reserved]==[doc['actual'][k]for k in('Bellman_inequalities','attaining_predecessors','explicit_options','singleton_context_checks','reserved_D_singleton_checks')]
    assert globalmax==F(doc['actual']['maximum'])
    return globalmax,fibercount,checks,preds,options,singles,reserved


def main():
    doc=V.load(TAG)
    for p,h in doc['sources'].items():assert U.digest(ROOT/p)==h
    assert F(doc['delta'])==F(11,512)and doc['scale']==512
    high,low=coarse(doc);a=actual(doc,high);bound=max(low,a[0])
    assert bound==F(doc['new_positive_family_upper'])and max(TARGET,bound)==F(doc['complete_finite_upper'])
    att=[tuple(t)for t in doc['unchanged_attaining_forest']]
    assert extra(att)==0 and R.value(att)==64*TARGET
    assert doc['complete_finite_maximum']==(str(TARGET)if bound<=TARGET else None)
    local=V.load('fresh_order11_fractional_cusp_closure_local_20261018');absolute=[]
    for row in local['ideal']['rows']:
        volume=prod(r+1 for r in row['ranks']);absolute.append((volume*F(row['normalized_margin']),row['ranks'],volume))
    minimum=min(x[0]for x in absolute)
    assert minimum>0 and minimum==F(doc['ideal_absolute_terminal_minimum'])
    expected=[dict(ranks=r,vertices=v,surplus=str(m))for m,r,v in absolute if m==minimum]
    assert expected==doc['ideal_absolute_minimizers']
    print('PASS_INDEPENDENT_ORDER11_FRACTIONAL_FINITE',len(doc['sources']),len(doc['coarse']['rows']),len(high),*a,bound,minimum,flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);R.N.A.R.M.R.resources();U.memory_guard();main();U.memory_guard();R.N.A.R.M.R.resources(start=False)
