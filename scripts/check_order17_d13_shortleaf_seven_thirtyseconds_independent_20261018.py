#!/usr/bin/env python3
"""Independent D17 distal-coupon finite arithmetic and ideal convolutions."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations_with_replacement,product
from math import prod
import signal
import check_order16_d_quarter_epsilon_ideal57_independent_20261018 as I
import fresh_d17_partial_weight_transfer_20261017 as P

ROOT=I.ROOT;W=I.W;C=I.C;K=I.K;OLD=I.OLD
TAG='fresh_order17_d13_shortleaf_seven_thirtyseconds_scout_20261018'
ORIGINAL='fresh_order17_d5_alternate_half_zero_transfer_20261018'
SUPPORT='fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018'
CHANGE={(1,15),(1,16),(11,13),(11,14)}
CASES=[1974,2367,2504,2511]


def increment(fs):
    co=Counter(map(tuple,fs));return 72*co['A',4]*co['D',13]


def value(fs):return I.value(fs)+increment(fs)


def basevalue(fs):
    fs=list(map(tuple,fs));co=Counter(fs)
    old17=16*co['D',5]*co['A',12]
    if co['A',11]and not(co['A',9]and co['A',2]):old17+=8*co['A',1]*co['D',5]
    return C.value(fs)+old17+I.X.increment(fs)+I.increment(fs)+increment(fs)


def local(doc):
    old=W.load(ORIGINAL)['local_proof'];weights={tuple(r['pair']):F(r['weight'])for r in old['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values())==113 and all(weights[p]==F(1,2)for p in CHANGE)
    new={p:F(7,32)if p in CHANGE else w for p,w in weights.items()}
    assert sum(new.values())==F(895,8)
    for k in('retained_basic_pairs','retained_short8','retained_triangle_pairs'):assert all(new[tuple(p)]==1 for p in old[k])
    s=W.load(SUPPORT);support=s['reconstructed_support']
    assert s['complete_eligible_patterns']==140 and len(support)==71 and sum(len(r['embeddings'])for r in support)==488
    hits=[r for r in support if set(map(tuple,r['marked_images']))&CHANGE]
    assert[r['index']for r in hits]==CASES
    assert all(len(r['embeddings'])==8 and set(map(tuple,r['marked_images']))<=CHANGE for r in hits)
    cases={r['index']:r for r in W.load(P.TAG)['cases']}
    assert[r['index']for r in doc['local']['proofs']]==CASES
    for r in doc['local']['proofs']:
        old=cases[r['index']];g=[[F(x)for x in row]for row in old['gram']];white=old['white']
        assert r['gram']==old['gram']and r['white']==white and len(g)==19
        assert P.Q.inertia(g)==[17,1,1]and P.Q.inertia(P.L.sub(g,white))==[15,0,0]
        black=sorted(set(range(19))-set(white))
        assert all(P.Q.inertia(P.L.sub(g,white+[b]))==[16,0,0]for b in black)
        z=list(map(F,old['radical']));ns=old['affine_nodes'];kernel=list(map(F,old['positive_affine_kernel']))
        assert all(sum(g[i][j]*z[j]for j in range(19))==0 for i in range(19))and z[3]*z[13]<0
        assert ns==sorted(set(range(19))-{3,13})and all(k>0 for k in kernel)
        assert all(sum(g[i][j]*kernel[t]for t,j in enumerate(ns))==0 for i in ns)
        assert r['radical']==old['radical']and r['affine_nodes']==ns and r['positive_affine_kernel']==old['positive_affine_kernel']
        selected=next(c for c in old['corners']if c['kind']=='affineD16_order17')
        assert selected['marks']==r['selected_marks']and F(r['selected_payment'])==F(7,32)
        marks=r['finite_coupon']['marks'];cc=P.L.blocks(g,sorted(set(white)|set(marks)),[1]*19)
        aa,bb=[next(c for c in cc if m in c['nodes'])for m in marks]
        if aa['type']=='D13':aa,bb=bb,aa;marks=marks[::-1]
        assert aa['type']=='A4'and bb['type']=='D13'
        assert aa['roles'][str(marks[0])]=='endpoint'and bb['roles'][str(marks[1])]=='short_fork_leaf'
        assert F(r['finite_coupon']['old'])==F(1,2)and F(r['finite_coupon']['new'])==F(25,32)
        assert r['finite_coupon']['marks']!=r['selected_marks']and F(r['total'])==F(7,32)+F(25,32)==1
    assert F(doc['local']['new_total'])==F(895,8)and F(doc['local']['new_envelope'])==F(895,136)
    return 140,488,32,4


def profiles(doc):
    env={int(k):F(v)for k,v in I.X.object_at(ROOT/'results/fresh_order10_low8_ambient58_frontier_20261018.json',b'"current_envelopes": {').items()}
    env[9]=F(22,5);env[15]=F(51,8);env[16]=F(113,17);before=dict(env);env[16]=F(895,136)
    domain=sorted((lo,)+p if lo else p for lo in range(8)for p in I.X.parts(56-lo))
    ii=doc['ideal'];assert len(domain)==2270==ii['complete_profiles']and domain==[tuple(r['ranks'])for r in ii['rows']]
    bad=[];zero=[];checks=0
    for rs,row in zip(domain,ii['rows']):
        h=[1]
        for r in rs:
            out=[0]*(len(h)+r)
            for i,a in enumerate(h):
                for j in range(r+1):out[i+j]+=a;checks+=1
            h=out
        v=prod(r+1 for r in rs);assert len(h)==57 and h==h[::-1]and sum(h)==v
        g=sum(min(i,56-i)*a for i,a in enumerate(h));curv=F(5*g-56*v,2*v)
        m=curv-sum(env[r]for r in rs);old=curv-sum(before[r]for r in rs)
        assert m==F(row['normalized_margin'])and old==F(row['old_margin'])
        assert m-old==F(9*rs.count(16),136)==F(row['normalized_saving'])
        if m<0:bad.append(row)
        if not m:zero.append(row)
    assert len(bad)==5 and not zero and bad==ii['negative_profiles']and zero==ii['zero_profiles']
    assert min(F(r['normalized_margin'])for r in ii['rows'])==F(ii['normalized_minimum'])
    return checks,len(bad)


def refinement(ff):
    cert=ff['actual_refinement'];target=F(1483,4)
    costly={(tuple(r['counts']),r['k5'],r['a12']):r for r in ff['rows']if F(r['upper'])>target}
    assert set(costly)=={(tuple(c['counts']),c['k5'],c['a12'])for c in cert['classes']}
    fibers={}
    for t in C.TYPES:fibers.setdefault(C.representative(t),[]).append(t)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    residual=[t for t in C.TYPES if C.representative(t)in reps and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in residual for b in residual)
    checks=preds=options=pairs=singletons=fixed_count=0
    for cl in cert['classes']:
        old=costly[tuple(cl['counts']),cl['k5'],cl['a12']];rem=old['remaining_rank']
        groups=sorted(Counter(C.fixed_components(old['counts'],old['key'])+[('A',5)]*old['k5']).items())
        domain={tuple(sorted(t for group in choice for t in group))for choice in product(*(
            combinations_with_replacement(fibers[t],n)for t,n in groups))}
        domain={fs for fs in domain if('D',13)in fs}
        assert domain=={tuple(map(tuple,r['fixed_actual']))for r in cl['rows']};fixed_count+=len(domain)
        for r in cl['rows']:
            fs=list(map(tuple,r['fixed_actual']));co=Counter(fs);tabs={};visited=set();best=None
            for t in r['tables']:
                key=tuple(map(tuple,t['types'])),tuple(t['weights']),len(t['values'])-1
                assert key not in tabs;tabs[key]=t
            for a1 in range(rem-3):
                for a2 in range((rem-a1-4)//2+1):
                    modes=[(False,0)]if a1 or a2 or not co['D',14]else[(True,0),(False,1)]
                    for no3,res3 in modes:
                        cap=rem-a1-2*a2-4-3*res3
                        if cap<0:continue
                        anchor=fs+[('A',1)]*a1+[('A',2)]*a2+[('A',4)]+[('A',3)]*res3;bc=basevalue(anchor)
                        ts=[t for t in residual if t[1]<=cap and not(no3 and t==('A',3))]
                        ws=[basevalue(anchor+[t])-bc for t in ts];key=tuple(ts),tuple(ws),cap
                        assert key in tabs;tab=tabs[key];singletons+=len(ts)
                        if key not in visited:
                            visited.add(key);a,b=K.bellman(tab);checks+=a;preds+=b
                            for i,a in enumerate(ts):
                                for j in range(i,len(ts)):
                                    b=ts[j]
                                    if a[1]+b[1]<=cap:assert basevalue(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
                        z=tab['values'][cap]
                        if z is not None:options+=1;best=bc+z if best is None else max(best,bc+z)
            assert visited==set(tabs)and F(best,64)==F(r['maximum'])and basevalue(r['attaining_forest'])==best
        assert max(F(r['maximum'])for r in cl['rows'])==F(cl['maximum'])
    assert(checks,preds,options,pairs,singletons,fixed_count)==(cert['Bellman_inequalities'],cert['attaining_predecessors'],cert['explicit_options'],cert['feasible_pair_additivity_checks'],cert['singleton_context_checks'],cert['actual_fixed_fibers'])
    maximum=max((F(c['maximum'])for c in cert['classes']),default=F());assert maximum==F(cert['maximum'])
    return checks,preds,options,pairs,singletons,fixed_count,maximum


def basefamily(doc):
    ff=doc['complete_base_increase_family'];rep=C.representative;typed=0
    for a in C.TYPES:
        for b in C.TYPES:
            rs={rep(a),rep(b)}
            allowance=(72 if rs=={('A',4),('D',13)}else 12 if rs=={('A',1),('D',15)}else
                       24 if ('B',9)in rs and bool(rs&{('A',5),('E',6),('A',7),('A',8)})else 0)
            assert C.actual_cross(a,b)+I.X.increment([a,b])+I.increment([a,b])+increment([a,b])<=C.C.upper_cross(rep(a),rep(b))+allowance
            typed+=1
    assert typed==ff['actual_type_pair_domination_checks']
    domain={(cs,k,key[9]):key for cs,u,key in C.fixed_domain()if cs[4]for k in range((57-u)//5+1)if 57-u-5*k>=4}
    assert len(domain)==ff['complete_classes']and set(domain)=={(tuple(r['counts']),r['k5'],r['a12'])for r in ff['rows']}
    checks=preds=ordinary_options=quarter_options=0
    for r in ff['rows']:
        cs=tuple(r['counts']);k=r['k5'];key=domain[cs,k,r['a12']]
        rem=57-sum(a*b for a,b in zip(range(9,18),cs))-5*k
        assert rem==r['remaining_rank']and list(key)==r['key']
        ws={j:C.item(j,key,k)+(12*cs[6]if j==1 else 24*cs[0]if j in(6,7,8)else 72*cs[4]if j==4 else 0)for j in C.ITEMS};tabs={}
        for t in r['tables']:
            no2=t['exclude2'];js=[j for j in C.ITEMS if j<=rem and(not no2 or j!=2)]
            assert list(map(tuple,t['types']))==[rep(('A'if j<=8 else'D',j))for j in js]and t['weights']==[ws[j]for j in js]
            assert len(t['values'])==rem+1;a,b=K.bellman(t);checks+=a;preds+=b;tabs[no2]=t
        assert set(tabs)=={False,True}
        fixed=C.fixed_components(cs,key)+[('A',5)]*k;fc=C.upper_value(fixed)+24*cs[0]*k
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
        assert F(max(ordinary,quarter if quarter is not None else-1),64)==F(r['upper'])
    assert(checks,preds,ordinary_options,quarter_options)==(ff['Bellman_inequalities'],ff['attaining_predecessors'],ff['ordinary_context_options'],ff['quarter_context_options'])
    assert max(F(r['upper'])for r in ff['rows'])==F(ff['coarse_maximum'])
    a=refinement(ff);maximum=max(max(F(r['upper'])for r in ff['rows']if F(r['upper'])<=F(1483,4)),a[-1])
    assert maximum==F(ff['complete_base_family_upper'])
    return len(domain),typed,checks,preds,ordinary_options,quarter_options,*a


def overlaps(doc):
    ff=doc['later_positive_overlaps'];small=[t for t in C.TYPES if t[1]<=8 and t not in[('A',1),('A',2)]]
    witnesses={'coupled_A16':('A',16),'coupled_D16':('D',16),'guarded_D17_A13':('A',13),'guarded_D17_A14':('A',14)}
    expected=set()
    for kind,w in witnesses.items():
        for s in range(18,58-w[1]-13-4):
            rem=57-w[1]-13-s
            for e in[None]+[t for t in C.TYPES if 9<=t[1]<=rem-4]:
                left=rem-(e[1]if e else 0)
                if left>=(7 if kind=='guarded_D17_A14'else 4):expected.add((kind,s,e))
    assert expected=={(r['kind'],r['D_rank'],tuple(r['extra'])if r['extra']else None)for r in ff['rows']}
    assert len(expected)==ff['complete_cases']
    checks=preds=options=pairs=0
    for r in ff['rows']:
        w=witnesses[r['kind']];s=r['D_rank'];e=tuple(r['extra'])if r['extra']else None
        assert tuple(r['witness'])==w
        fixed=[w,('D',13),('D',s)]+([e]if e else[]);left=57-s-w[1]-13-(e[1]if e else 0)
        domain={(a,b)for a in range(left+1)for b in range((left-a)//2+1)if left-a-2*b>=4 and(r['kind']!='guarded_D17_A14'or(a and b))}
        assert domain=={(t['a1'],t['a2'])for t in r['tables']};best=None
        for t in r['tables']:
            a1,a2=t['a1'],t['a2'];anchor=fixed+[('A',1)]*a1+[('A',2)]*a2+[('A',4)];cap=left-a1-2*a2-4
            bc=value(anchor);ts=[x for x in small if x[1]<=cap];ws=[value(anchor+[x])-bc for x in ts];tab=t['table']
            assert list(map(tuple,tab['types']))==ts and tab['weights']==ws and len(tab['values'])==cap+1
            for i,a in enumerate(ts):
                for j in range(i,len(ts)):
                    b=ts[j]
                    if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
            a,b=K.bellman(tab);checks+=a;preds+=b
            z=tab['values'][cap]
            if z is not None:options+=1;best=bc+z if best is None else max(best,bc+z)
        assert F(best,64)==F(r['maximum'])and value(r['attaining_forest'])==best
    assert(checks,preds,options,pairs)==(ff['Bellman_inequalities'],ff['attaining_predecessors'],ff['explicit_options'],ff['feasible_pair_additivity_checks'])
    for k in witnesses:assert max(F(r['maximum'])for r in ff['rows']if r['kind']==k)==F(ff['family_maxima'][k])
    assert max(map(F,ff['family_maxima'].values()))==F(ff['maximum'])
    return len(expected),checks,preds,options,pairs


def main():
    d=W.load(TAG);hashes=0
    for p,h in d['sources'].items():assert W.digest(ROOT/p)==h;hashes+=1
    a=local(d);b=profiles(d);c=basefamily(d);e=overlaps(d)
    upper=max(F(1483,4),F(d['complete_base_increase_family']['complete_base_family_upper']),F(d['later_positive_overlaps']['maximum']))
    assert upper==F(1483,4)==F(d['complete_finite_upper'])==F(d['complete_finite_maximum'])
    assert value(d['unchanged_formal_attainer'])==64*upper;OLD.U.memory_guard()
    print('PASS_INDEPENDENT_D17_D13_SHORTLEAF_SEVEN_THIRTYSECONDS',hashes,*a,*b,*c,*e,upper,flush=True)


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);OLD.R.resources();main();OLD.R.resources(start=False)
