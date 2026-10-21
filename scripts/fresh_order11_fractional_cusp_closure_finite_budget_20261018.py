#!/usr/bin/env python3
"""Complete conditional finite-family budget for the fractional order11 rule."""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement,product
import json
from math import prod
from pathlib import Path
import signal
import fresh_a16_rank17_spectator_extension_scout_20261018 as R

ROOT=R.ROOT;C=R.C;H=R.H;V=R.V;U=R.U
TAG='fresh_order11_fractional_cusp_closure_finite_budget_20261018'
TARGET=F(11861,32);SCALE=512


def increment(fs):
    co=Counter(fs);a1,a2,a3,a4=(co['A',r]for r in range(1,5))
    return(22*a1+44*(a2+a3+a4))*co['D',10]+22*a2*co['B',10]+(88*a2+44*a3)*co['B',8]


def current(fs):return 8*R.value(fs)+increment(fs)


@lru_cache(maxsize=8192)
def multiplicity_base(items):
    # This is exactly H.basevalue, grouping its ordinary unordered
    # component pairs by type. Equal-type pairs have binomial multiplicity.
    fs=[t for t,n in items for _ in range(n)]
    z=sum(H.O.actual_base(t)*n for t,n in items)
    for i,(a,na)in enumerate(items):
        z+=H.O.actual_cross(a,a)*na*(na-1)//2
        z+=sum(H.O.actual_cross(a,b)*na*nb for b,nb in items[i+1:])
    return 8*(z+H.O.context(fs)+H.E.Q.increment(fs)+H.E.increment(fs)+H.increment(fs))+increment(fs)


def base(fs):return multiplicity_base(tuple(sorted(Counter(fs).items())))


def guarded(fs):return base(fs)+8*V.context_delta(fs)


def tabhash(tab):return sha256(json.dumps(C.enc(tab),separators=(',',':')).encode()).hexdigest()


def coarse():
    rows=[];checks=preds=ordinary_options=quarter_options=typed=0
    for a in C.TYPES:
        for b in C.TYPES:
            pair={a,b};rp={C.representative(a),C.representative(b)};extra=0
            if ('D',10)in rp:
                extra+=22 if ('A',1)in rp else 44 if rp&{('A',2),('A',3),('A',4)}else 0
            if rp=={('A',2),('A',8)}:extra+=88
            if rp=={('A',3),('A',8)}:extra+=44
            if rp=={('A',1),('D',15)}:extra+=96
            if rp in({('A',4),('D',13)},{('A',2),('D',14)}):extra+=1024
            if ('B',9)in rp and rp&{('A',5),('E',6),('A',7),('A',8)}:extra+=192
            actual=8*C.actual_cross(a,b)+8*(H.E.Q.increment([a,b])+H.E.increment([a,b])+H.increment([a,b]))+increment([a,b])
            if pair in({('A',4),('A',13)},{('A',2),('A',14)}):actual+=1024
            assert actual<=8*C.C.upper_cross(C.representative(a),C.representative(b))+extra
            typed+=1
    for cs,used,key in C.fixed_domain():
        for k5 in range((57-used)//5+1):
            for k8 in range((57-used-5*k5)//8+1):
                if not(cs[1]or k8):continue
                rem=57-used-5*k5-8*k8
                fixed=C.fixed_components(cs,key)+[('A',5)]*k5+[('A',8)]*k8
                fc=8*(C.upper_value(fixed)+24*cs[0]*(k5+k8))
                ws={r:8*C.item(r,key,k5)+(22*cs[1]+96*cs[6]if r==1 else
                    44*cs[1]+88*k8+1024*cs[5]if r==2 else
                    44*cs[1]+44*k8 if r==3 else
                    44*cs[1]+1024*cs[4]if r==4 else
                    192*cs[0]if r in(6,7)else 0)for r in C.ITEMS if r!=8}
                tabs=[]
                for no2 in(False,True):
                    ranks=[r for r in ws if r<=rem and(not no2 or r!=2)]
                    ts=[C.representative(('A'if r<=8 else'D',r))for r in ranks]
                    tab,a,b=C.additive(ts,[ws[r]for r in ranks],rem);checks+=a;preds+=b;tabs.append(tab)
                opts=[tabs[0]['values'][rem]]
                for a1 in(1,2):
                    og=a1==2 and bool(cs[1]);ng=bool(cs[0]and k5)
                    if rem<a1 or not(og or ng)or(a1==2 and not og):continue
                    for a2 in range((rem-a1)//2+1):
                        v=tabs[1]['values'][rem-a1-2*a2]
                        if v is not None:opts.append(a1*ws[1]+a2*ws[2]+v+(288*a2*(a2-1)if og else 0)+(512*a2*k5 if ng else 0))
                ordinary_options+=len(opts);order17=128*k5*key[9]+(64*k5*rem if cs[2]else 0)
                ordinary=fc+max(v for v in opts if v is not None)+order17;quarter=None
                if cs[1]and rem>=7:
                    for a1 in range(1,rem-5):
                        for a4 in range(1,(rem-a1-2)//4+1):
                            for a2 in range(1,(rem-a1-4*a4)//2+1):
                                v=tabs[1]['values'][rem-a1-4*a4-2*a2]
                                if v is None:continue
                                z=fc+a1*ws[1]+a2*ws[2]+a4*ws[4]+v+order17+(288*a2*(a2-1)if a1>=2 else 0)+(512*a2*k5 if cs[0]else 0)+512*a2*a4
                                quarter=z if quarter is None else max(quarter,z);quarter_options+=1
                rows.append(dict(counts=cs,k5=k5,k8=k8,a12=key[9],key=key,remaining_rank=rem,
                    table_hashes=[tabhash(t)for t in tabs],ordinary_upper=F(ordinary,512),
                    quarter_upper=None if quarter is None else F(quarter,512),upper=F(max(ordinary,quarter if quarter is not None else -1),512)))
    return dict(rows=rows,complete_classes=len(rows),actual_type_pair_domination_checks=typed,
        Bellman_inequalities=checks,attaining_predecessors=preds,
        ordinary_context_options=ordinary_options,quarter_context_options=quarter_options,
        maximum=max(r['upper']for r in rows))


def refinement(coarse):
    fibers={}
    for t in C.TYPES:fibers.setdefault(C.representative(t),[]).append(t)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS if r!=8}
    residual=[t for t in C.TYPES if C.representative(t)in reps and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in residual+[('A',1),('A',2)]for b in residual+[('A',1),('A',2)])
    out=[];checks=preds=options=singletons=reserve_checks=0
    for old in coarse['rows']:
        if old['upper']<=TARGET:continue
        U.memory_guard();fixed=C.fixed_components(old['counts'],old['key'])+[('A',5)]*old['k5']+[('A',8)]*old['k8']
        groups=sorted(Counter(fixed).items());rem=old['remaining_rank']
        domain=sorted({tuple(sorted(t for group in choice for t in group))for choice in product(*(
            combinations_with_replacement(fibers[t],n)for t,n in groups))})
        domain=[fs for fs in domain if any(t in fs for t in[('B',8),('B',10),('D',10)])];actual=[]
        for fs in domain:
            assert base(list(fs))==8*H.basevalue(list(fs))+increment(list(fs))
            co=Counter(fs);positive_D17=bool(co['A',13]or co['A',14]);cache={};best=None;hashes=[]
            raw={t:base(list(fs)+[t])-base(list(fs))-8*(C.context(list(fs)+[t])-C.context(list(fs)))for t in residual}
            for a1 in range(rem+1):
                for a2 in range((rem-a1)//2+1):
                    modes=[(False,0)]if a1 or a2 or not co['D',14]else[(True,0),(False,1)]
                    for no3,res3 in modes:
                        left=rem-a1-2*a2-3*res3
                        if left<0:continue
                        anchor=list(fs)+[('A',1)]*a1+[('A',2)]*a2+[('A',3)]*res3
                        for guard in((False,True)if positive_D17 else[None]):
                            if guard and left<18:continue
                            cap=left-(18 if guard else 0);fun=guarded if guard is not None else base
                            anch=anchor+([('D',18)]if guard else[]);bc=fun(anch)
                            ts=[t for t in residual if t[1]<=cap and not(no3 and t==('A',3))and(guard is not False or not(t[0]=='D'and t[1]>=18))]
                            # Actual residual ordinary cross charges vanish.
                            # Once presences are fixed, C.context only adds
                            # the standard A4 term; the guarded D17 terms are
                            # linear in A4/A2 and in the number of D17 factors.
                            ws=[raw[t]+(8*(32*co['A',12]+64*a2*bool(co['A',10]))if a1 and t==('A',4)else 0)
                                +(1024*co['A',13]if guard and t==('A',4)else -512 if guard and t==('D',17)else 0)for t in ts]
                            for t,w in zip(ts,ws):assert fun(anch+[t])-bc==w
                            singletons+=len(ts);key=tuple(ts),tuple(ws),cap
                            if key not in cache:
                                tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b;cache[key]=tab;hashes.append(tabhash(tab))
                            tab=cache[key]
                            for ds in(range(18,left+1)if guard else[None]):
                                n=left-(ds if guard else 0);z=tab['values'][n]
                                if z is None:continue
                                aa=anchor+([('D',ds)]if guard else[]);b0=fun(aa)
                                for t,w in zip(ts,ws):
                                    if t[1]<=n:assert fun(aa+[t])-b0==w;reserve_checks+=1
                                options+=1
                                if best is None or b0+z>best[0]:best=(b0+z,aa,tab,n,guard)
            assert best is not None
            z,att,tab,n,guard=best;att=list(att);fun=guarded if guard is not None else base
            while n:t=tab['types'][tab['last'][n]];att.append(t);n-=t[1]
            assert sum(t[1]for t in att)==57 and fun(att)==z
            actual.append(dict(fixed_actual=fs,table_hashes=hashes,maximum=F(z,512),attaining_forest=sorted(att),
                added_increment=F(increment(att),512),actual_current=F(current(att),512)))
        if actual:out.append(dict(counts=old['counts'],k5=old['k5'],k8=old['k8'],a12=old['a12'],rows=actual,maximum=max(r['maximum']for r in actual)))
    return dict(classes=out,actual_fixed_fibers=sum(len(c['rows'])for c in out),Bellman_inequalities=checks,
        attaining_predecessors=preds,explicit_options=options,singleton_context_checks=singletons,
        reserved_D_singleton_checks=reserve_checks,maximum=max((c['maximum']for c in out),default=F()))


def run():
    local=R.V.load('fresh_order11_fractional_cusp_closure_local_20261018');assert F(local['delta'])==F(11,512)
    absolute=[dict(ranks=r['ranks'],vertices=prod(j+1 for j in r['ranks']),
        surplus=prod(j+1 for j in r['ranks'])*F(r['normalized_margin']))for r in local['ideal']['rows']]
    amin=min(r['surplus']for r in absolute)
    print('ORDER11_IDEAL_ABSOLUTE',amin,flush=True)
    coarse_doc=coarse()
    print('ORDER11_FINITE_COARSE',coarse_doc['complete_classes'],sum(r['upper']>TARGET for r in coarse_doc['rows']),flush=True)
    actual=refinement(coarse_doc)
    bound=max(max(r['upper']for r in coarse_doc['rows']if r['upper']<=TARGET),actual['maximum'])
    att=[('D',13)]*3+[('A',4)]*4+[('A',2)];assert current(att)==512*TARGET and increment(att)==0
    paths={Path(__file__),Path(R.__file__),ROOT/f'results/{R.TAG}.json',ROOT/'results/order11_fractional_cusp_closure_research_20261018.md',
        ROOT/'results/fresh_order11_fractional_cusp_closure_local_20261018.json'}
    return C.enc(dict(status='PASS_ORDER11_FRACTIONAL_CUSP_CLOSURE_FINITE_BUDGET',sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},
        delta='11/512',scale=512,coarse=coarse_doc,actual=actual,new_positive_family_upper=bound,
        complete_finite_upper=max(TARGET,bound),complete_finite_maximum=TARGET if bound<=TARGET else None,
        unchanged_attaining_forest=att,ideal_absolute_terminal_minimum=amin,
        ideal_absolute_minimizers=[r for r in absolute if r['surplus']==amin],
        new_dimension_bound=None,mathematical_proof_in_separate_report=True))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);R.N.A.R.M.R.resources();U.memory_guard()
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();R.N.A.R.M.R.resources(start=False)
    if args.verify:assert R.V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['coarse']['complete_classes'],out['actual']['actual_fixed_fibers'],out['new_positive_family_upper'],out['complete_finite_maximum'],flush=True)
