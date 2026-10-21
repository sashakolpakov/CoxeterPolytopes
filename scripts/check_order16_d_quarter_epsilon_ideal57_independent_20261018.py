#!/usr/bin/env python3
"""Independent sparse-D epsilon arithmetic and fresh rank56 convolutions."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations_with_replacement,product
from math import prod
import signal
import check_order10_d_quarter_ideal57_independent_20261018 as X

ROOT=X.ROOT;W=X.W;C=X.C;K=X.K;OLD=X.OLD
TAG='fresh_order16_d_quarter_epsilon_ideal57_scout_20261018'
LOCAL='fresh_order16_cd_sparse_sixteenth_transfer_20261018'


def increment(fs):
    co=Counter(map(tuple,fs));return 12*co['A',1]*co['D',15]


def value(fs):return X.Z.value(fs)+X.increment(fs)+increment(fs)


def refinement(ff,fibers):
    coarse={(tuple(r['counts']),r['k5'],r['a12']):r for r in ff['rows']if F(r['upper'])>F(1483,4)}
    cert=ff['actual_refinement'];assert len(coarse)==4==len(cert['classes'])
    assert set(coarse)=={(tuple(c['counts']),c['k5'],c['a12'])for c in cert['classes']}
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    residual=[t for t in C.TYPES if C.representative(t)in reps and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in residual for b in residual)
    checks=preds=options=contexts=pairs=fixed_count=0
    for cl in cert['classes']:
        old=coarse[tuple(cl['counts']),cl['k5'],cl['a12']]
        assert old['counts'][6]==old['counts'][1]==1 and old['k5']==0
        rem=old['remaining_rank'];assert rem==cl['remaining_rank']
        groups=sorted(Counter(C.fixed_components(old['counts'],old['key'])).items())
        expected={tuple(sorted(t for group in choice for t in group))for choice in product(*(
            combinations_with_replacement(fibers[t],n)for t,n in groups))}
        expected={fs for fs in expected if('D',15)in fs}
        assert expected=={tuple(map(tuple,r['fixed_actual']))for r in cl['rows']}
        fixed_count+=len(expected)
        for r in cl['rows']:
            fs=list(map(tuple,r['fixed_actual']));best=None
            assert{(t['a1'],t['a2'])for t in r['tables']}=={(a,b)for a in range(1,rem+1)for b in range((rem-a)//2+1)}
            for tab in r['tables']:
                a1,a2=tab['a1'],tab['a2'];anchor=fs+[('A',1)]*a1+[('A',2)]*a2
                cap=rem-a1-2*a2;bc=value(anchor);assert bc==C.value(anchor)+increment(anchor)
                ts=[t for t in residual if t[1]<=cap];ws=[value(anchor+[t])-bc for t in ts]
                assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
                for t,w in zip(ts,ws):assert value(anchor+[t])==C.value(anchor+[t])+increment(anchor+[t])
                contexts+=len(ts)
                if a1==1:
                    for i,a in enumerate(ts):
                        for j in range(i,len(ts)):
                            b=ts[j]
                            if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
                a,b=K.bellman(tab);checks+=a;preds+=b
                z=tab['values'][cap]
                if z is not None:options+=1;best=bc+z if best is None else max(best,bc+z)
            assert F(best,64)==F(r['maximum'])and value(r['attaining_forest'])==best
        assert max(F(r['maximum'])for r in cl['rows'])==F(cl['maximum'])
    assert(checks,preds,options,contexts,pairs,fixed_count)==(cert['Bellman_inequalities'],cert['attaining_predecessors'],cert['explicit_options'],cert['singleton_context_checks'],cert['feasible_pair_additivity_checks'],cert['actual_fixed_fibers'])
    maximum=max(F(c['maximum'])for c in cert['classes']);assert maximum==F(cert['actual_refined_maximum'])
    return checks,preds,options,contexts,pairs,fixed_count,maximum


def local(doc):
    old=W.load(LOCAL);d=next(r for r in old['rules']if r['kind']=='D')
    changed={(1,14),(1,15),(10,12),(10,13)}
    assert set(map(tuple,d['changed_pairs']))==changed
    weights={tuple(r['pair']):F(r['weight'])for r in d['weighted_pairs']}
    assert sum(weights.values())==F(411,4)and all(weights[p]==F(7,16)for p in changed)
    new={p:F(1,4)if p in changed else w for p,w in weights.items()}
    assert sum(new.values())==102
    assert all(new[tuple(p)]==weights[tuple(p)]for p in d['retained_basic_pairs'])
    assert all(new[tuple(p)]==1 for p in d['retained_triangle_pairs'])
    outcomes=Counter();payments=[]
    for c in d['cases']:
        for r in c['completions']:outcomes[r['proof']['outcome']]+=1
        for r in c['current_payments']:
            corners=[]
            for z in r['corners']:
                if z['kind']=='affineD15_order16':p=F(1,4)
                else:p=sum((F(t['payment'])for t in z['existing_coupon']),F())+F(len(z['new_coupon']),4)
                corners.append(dict(marks=z['marks'],kind=z['kind'],payment=str(p)))
            total=sum(F(z['payment'])for z in corners)
            assert total==(1 if c['case']in(2368,2512)else F(7,4))
            payments.append(dict(case=c['case'],corners=corners,total=str(total)))
    assert outcomes==dict(excluded_projection_above_norm=32,excluded_lorentz_index=16,sealed_completion=32)
    assert payments==doc['local']['whole_map_payments']
    assert doc['local']['new_family_totals']==dict(A='80',B='102',C='101',D='102')
    assert F(doc['local']['new_rank15_envelope'])==F(51,8)
    return sum(outcomes.values()),len(payments)


def profiles(doc):
    env={int(k):F(v)for k,v in X.object_at(ROOT/'results/fresh_order10_low8_ambient58_frontier_20261018.json',b'"current_envelopes": {').items()}
    assert env[9]==F(71,16)and env[16]==F(227,34)and env[15]==F(411,64)
    env[9]=F(22,5);env[16]=F(113,17);previous=dict(env);env[15]=F(51,8)
    domain=sorted((lo,)+p if lo else p for lo in range(8)for p in X.parts(56-lo))
    ii=doc['ideal'];assert len(domain)==2270==ii['complete_profiles']
    assert domain==[tuple(r['ranks'])for r in ii['rows']]
    bad=[];zeros=[];checks=0
    for rs,row in zip(domain,ii['rows']):
        h=[1]
        for r in rs:
            out=[0]*(len(h)+r)
            for i,a in enumerate(h):
                for j in range(r+1):out[i+j]+=a;checks+=1
            h=out
        v=prod(r+1 for r in rs);assert len(h)==57 and h==h[::-1]and sum(h)==v
        g=sum(min(i,56-i)*a for i,a in enumerate(h));curv=F(5*g-56*v,2*v)
        m=curv-sum(env[r]for r in rs);before=curv-sum(previous[r]for r in rs)
        assert m==F(row['normalized_margin'])and before==F(row['old_margin'])
        assert m-before==F(3*rs.count(15),64)==F(row['normalized_saving'])
        if m<0:bad.append(row)
        if not m:zeros.append(row)
    assert bad==ii['negative_profiles']and zeros==ii['zero_profiles']
    assert len(bad)==8 and not zeros
    assert min(F(r['normalized_margin'])for r in ii['rows'])==-F(1553,22627)==F(ii['normalized_minimum'])
    return checks,len(bad)


def finite(doc):
    ff=doc['finite'];rep=C.representative;fibers={}
    for t in C.TYPES:fibers.setdefault(rep(t),[]).append(t)
    f15=fibers['D',15];assert f15==[('A',15),('B',15),('D',15)]
    base=max(C.actual_base(t)+(32 if t in(('A',15),('D',15))else 0)for t in f15)
    gap=C.C.upper_base(('D',15))-base;assert gap>=0
    assert F(base,64)==F(ff['rank15_typed_base_with_context_padding'])and F(gap,64)==F(ff['rank15_base_saving_per_component'])
    assert list(map(tuple,ff['rank15_actual_fiber']))==f15
    cross={}
    for b,ts in fibers.items():
        bound=max(C.actual_cross(a,t)+X.increment([a,t])+increment([a,t])for a in f15 for t in ts)
        cross[b]=C.C.upper_cross(('D',15),b)+(12 if b==('A',1)else 0)-bound
        assert cross[b]>=0
    assert {(tuple(r['partner']),F(r['saving']))for r in ff['rank15_cross_savings']}=={(b,F(n,64))for b,n in cross.items()if n}
    extra={}
    for a in C.TYPES:
        for b in C.TYPES:
            need=max(0,C.actual_cross(a,b)+X.increment([a,b])+increment([a,b])-C.C.upper_cross(rep(a),rep(b)))
            if need:key=tuple(sorted((rep(a),rep(b))));extra[key]=max(extra.get(key,0),need)
    expected={tuple(sorted((t,('B',9)))):24 for t in(('A',5),('E',6),('A',7),('A',8))}
    expected[tuple(sorted((('A',1),('D',15))))]=12;assert extra==expected
    domain={(cs,k,key[9]):key for cs,u,key in C.fixed_domain()if cs[6]for k in range((57-u)//5+1)}
    assert len(domain)==606==ff['classes']and{(tuple(r['counts']),r['k5'],r['a12'])for r in ff['rows']}==set(domain)
    checks=preds=ordinary_options=quarter_options=0
    for r in ff['rows']:
        cs=tuple(r['counts']);k=r['k5'];key=domain[cs,k,r['a12']]
        rem=57-sum(a*b for a,b in zip(range(9,18),cs))-5*k
        assert rem==r['remaining_rank']and list(key)==r['key']
        fixed=C.fixed_components(cs,key)+[('A',5)]*k
        saving=sum(cross[b]if a==('D',15)else cross[a]if b==('D',15)else 0 for i,a in enumerate(fixed)for b in fixed[i+1:])
        fc=C.upper_value(fixed)+24*cs[0]*k-gap*cs[6]-saving
        ws={j:C.item(j,key,k)+(12*cs[6]if j==1 else 24*cs[0]if j in(6,7,8)else 0)-cs[6]*cross[rep(('A'if j<=8 else'D',j))]for j in C.ITEMS}
        tabs={}
        for t in r['tables']:
            no2=t['exclude2'];js=[j for j in C.ITEMS if j<=rem and(not no2 or j!=2)]
            assert list(map(tuple,t['types']))==[rep(('A'if j<=8 else'D',j))for j in js]
            assert t['weights']==[ws[j]for j in js]and len(t['values'])==rem+1
            a,b=K.bellman(t);checks+=a;preds+=b;tabs[no2]=t
        assert set(tabs)=={False,True}
        opts=[tabs[False]['values'][rem]]
        for a1 in(1,2):
            og=a1==2 and bool(cs[1]);ng=bool(cs[0]and k)
            if rem<a1 or not(og or ng)or(a1==2 and not og):continue
            for a2 in range((rem-a1)//2+1):
                v=tabs[True]['values'][rem-a1-2*a2]
                if v is not None:opts.append(a1*ws[1]+a2*ws[2]+v+(36*a2*(a2-1)if og else 0)+(64*a2*k if ng else 0))
        ordinary_options+=len(opts);order17=16*k*key[9]+(8*k*rem if cs[2]else 0)
        ordinary=fc+max(v for v in opts if v is not None)+order17;quarter=None
        if cs[1]and rem>=7:
            for a1 in range(1,rem-5):
                for a4 in range(1,(rem-a1-2)//4+1):
                    for a2 in range(1,(rem-a1-4*a4)//2+1):
                        v=tabs[True]['values'][rem-a1-4*a4-2*a2]
                        if v is None:continue
                        z=fc+a1*ws[1]+a2*ws[2]+a4*ws[4]+v+order17+(36*a2*(a2-1)if a1>=2 else 0)+(64*a2*k if cs[0]else 0)+64*a2*a4
                        quarter=z if quarter is None else max(quarter,z);quarter_options+=1
        assert F(ordinary,64)==F(r['ordinary_upper'])
        assert(None if quarter is None else F(quarter,64))==(None if r['quarter_upper']is None else F(r['quarter_upper']))
        assert F(max(ordinary,quarter if quarter is not None else -1),64)==F(r['upper'])
    assert(checks,preds,ordinary_options,quarter_options)==(ff['Bellman_inequalities'],ff['attaining_predecessors'],ff['ordinary_context_options'],ff['quarter_context_options'])
    assert max(F(r['upper'])for r in ff['rows'])==F(ff['coarse_upper_before_actual_refinement'])
    refined=refinement(ff,fibers)
    maximum=max(max(F(r['upper'])for r in ff['rows']if F(r['upper'])<=F(1483,4)),refined[-1])
    assert maximum==F(ff['base_order17_and_quarterD10_D15_family_upper'])
    overlap=365+F(3*(57-18-13-15),16);assert overlap==F(5873,16)==F(ff['guarded_D17_overlap_upper'])
    assert max(maximum,overlap)==F(ff['complete_new_positive_family_upper'])
    assert value(doc['unchanged_formal_attainer'])==64*F(1483,4)==64*F(doc['previous_finite_maximum'])
    assert max(F(1483,4),maximum,overlap)==F(doc['conditional_complete_finite_upper'])
    return checks,preds,ordinary_options,quarter_options,*refined


def main():
    doc=W.load(TAG);hashes=0
    for p,h in doc['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    a=local(doc);b=profiles(doc);c=finite(doc);OLD.U.memory_guard()
    print('PASS_INDEPENDENT_ORDER16_D_QUARTER_EPSILON_IDEAL57',hashes,*a,*b,*c,doc['conditional_complete_finite_upper'],flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
