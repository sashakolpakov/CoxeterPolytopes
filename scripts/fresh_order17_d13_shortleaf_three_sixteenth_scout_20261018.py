#!/usr/bin/env python3
"""Separate four-half D17 repair with an explicit shared finite budget."""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as F
from itertools import combinations_with_replacement,product
import json
from pathlib import Path
import signal
import fresh_order16_d_quarter_epsilon_ideal57_scout_20261018 as E

ROOT=E.ROOT;C=E.C;V=E.V;U=E.U;O=V.O
TAG='fresh_order17_d13_shortleaf_three_sixteenth_scout_20261018'
PARTIAL='fresh_d17_partial_weight_transfer_20261017'
CHANGE={(1,15),(1,16),(11,13),(11,14)}
CASES=[1974,2367,2504,2511]


def increment(fs):
    co=Counter(fs);return 80*co['A',4]*co['D',13]


def value(fs):return E.value(fs)+increment(fs)


def basevalue(fs):return O.value(fs)+E.Q.increment(fs)+E.increment(fs)+increment(fs)


def local():
    old=V.load(O.TAG)['local_proof'];sup=V.load(O.M.TAG)
    weights={tuple(r['pair']):F(r['weight'])for r in old['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values())==113
    assert all(weights[p]==F(1,2)for p in CHANGE)
    new={p:F(3,16)if p in CHANGE else w for p,w in weights.items()}
    assert sum(new.values())==F(447,4)
    for k in('retained_basic_pairs','retained_short8','retained_triangle_pairs'):
        assert all(new[tuple(p)]==1 for p in old[k])
    target=old['target'];_,dist,_=O.M.P.distance_data(target)
    assert all(dist[a][b]==12 for a,b in CHANGE)
    for a in O.M.P.embeddings(target,target):
        assert{tuple(sorted(a[i]for i in p))for p in CHANGE}==CHANGE
    support=sup['reconstructed_support'];assert sup['complete_eligible_patterns']==140
    assert len(support)==71 and sum(len(r['embeddings'])for r in support)==488
    hits=[r for r in support if set(map(tuple,r['marked_images']))&CHANGE]
    assert[r['index']for r in hits]==CASES
    assert all(len(r['embeddings'])==8 and set(map(tuple,r['marked_images']))<=CHANGE for r in hits)
    prior={r['index']:r for r in V.load(PARTIAL)['cases']};rows=[]
    for index in CASES:
        p=prior[index];g=V.A.matrix(p['gram']);white=p['white'];z=list(map(F,p['radical']))
        assert len(g)==19 and V.Q.inertia(g)==[17,1,1]
        assert V.Q.inertia(V.A.sub(g,white))==[15,0,0]
        black=sorted(set(range(19))-set(white))
        assert all(V.Q.inertia(V.A.sub(g,white+[b]))==[16,0,0]for b in black)
        assert all(sum(g[i][j]*z[j]for j in range(19))==0 for i in range(19))and z[3]*z[13]<0
        ns=p['affine_nodes'];kernel=list(map(F,p['positive_affine_kernel']))
        assert ns==sorted(set(range(19))-{3,13})and all(k>0 for k in kernel)
        assert all(sum(g[i][j]*kernel[t]for t,j in enumerate(ns))==0 for i in ns)
        selected=[];paying=[]
        for q in p['corners']:
            if q['kind']=='affineD16_order17':selected.append(q['marks']);continue
            if not q['coupon_roles']:continue
            blocks=V.D.finite_blocks(g,sorted(set(white)|set(q['marks'])))
            a,b=[next(b for b in blocks if m in b['nodes'])for m in q['marks']]
            if a['type']=='D13':a,b=b,a;marks=q['marks'][::-1]
            else:marks=q['marks']
            assert a['type']=='A4'and b['type']=='D13'
            ai=V.T.mark_pairs(g,a['nodes'],[marks[0]],4,'A')
            bi=V.T.mark_pairs(g,b['nodes'],[marks[1]],13,'D')
            assert ai==[(0,),(3,)]and bi==[(11,),(12,)]
            paying.append(dict(marks=q['marks'],A4_images=ai,D13_images=bi,old='1/2',new='13/16'))
        assert len(selected)==len(paying)==1 and selected[0]!=paying[0]['marks']
        rows.append(dict(index=index,gram=p['gram'],white=white,radical=p['radical'],
            affine_nodes=ns,positive_affine_kernel=p['positive_affine_kernel'],
            selected_marks=selected[0],selected_payment='3/16',finite_coupon=paying[0],total='1'))
    return dict(complete_patterns=140,successful_contexts=71,complete_maps=488,affected_maps=32,
        changed_pairs=sorted(CHANGE),old_weight='1/2',new_weight='3/16',old_total='113',new_total='447/4',
        family_totals=dict(A='85',B='111',C='108',D='447/4'),new_envelope='447/68',
        finite_coupon_old='1/2',finite_coupon_new='13/16',component_pair_increment='5/4',
        retained_basic_short8_triangle_units=True,proofs=rows)


def ideal():
    old=V.load(E.TAG)['ideal'];assert old['new_negative_count']==8;rows=[]
    for r in old['rows']:
        saving=F(5*r['ranks'].count(16),68);margin=F(r['normalized_margin'])+saving
        rows.append(dict(ranks=r['ranks'],old_margin=r['normalized_margin'],normalized_saving=saving,normalized_margin=margin))
    bad=[r for r in rows if r['normalized_margin']<0];zero=[r for r in rows if not r['normalized_margin']]
    assert len(bad)==4 and not zero
    return dict(complete_profiles=len(rows),rows=rows,old_negative_count=8,new_negative_count=4,
        negative_profiles=bad,zero_profiles=zero,normalized_minimum=min(r['normalized_margin']for r in rows),
        normalized_saving_per_rank16='5/68',no_negative_surplus_extension=True)


def overlaps():
    small=[t for t in C.TYPES if t[1]<=8 and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in small for b in small)
    rows=[];checks=preds=options=pairs=0
    for kind,witness in(('coupled_A16',('A',16)),('coupled_D16',('D',16)),
                         ('guarded_D17_A13',('A',13)),('guarded_D17_A14',('A',14))):
        for s in range(18,58-witness[1]-13-4):
            rem=57-witness[1]-13-s
            for extra in[None]+[t for t in C.TYPES if 9<=t[1]<=rem-4]:
                fixed=[witness,('D',13),('D',s)]+([extra]if extra else[])
                left=rem-(extra[1]if extra else 0);tabs=[];opts=[]
                for a1 in range(left+1):
                    for a2 in range((left-a1)//2+1):
                        if kind=='guarded_D17_A14'and not(a1 and a2):continue
                        cap=left-a1-2*a2-4
                        if cap<0:continue
                        anchor=fixed+[('A',1)]*a1+[('A',2)]*a2+[('A',4)];bc=value(anchor)
                        ts=[t for t in small if t[1]<=cap];ws=[value(anchor+[t])-bc for t in ts]
                        for i,a in enumerate(ts):
                            for j in range(i,len(ts)):
                                b=ts[j]
                                if a[1]+b[1]<=cap:assert value(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
                        tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b
                        tabs.append(dict(a1=a1,a2=a2,table=tab))
                        if tab['values'][cap]is not None:opts.append((bc+tab['values'][cap],len(tabs)-1,anchor));options+=1
                if not opts:continue
                z,ti,att=max(opts);tab=tabs[ti]['table'];att=list(att);left=len(tab['values'])-1
                while left:t=tab['types'][tab['last'][left]];att.append(t);left-=t[1]
                assert value(att)==z and sum(t[1]for t in att)==57
                rows.append(dict(kind=kind,witness=witness,D_rank=s,extra=extra,tables=tabs,
                    maximum=F(z,64),attaining_forest=sorted(att)))
    return dict(rows=rows,complete_cases=len(rows),Bellman_inequalities=checks,attaining_predecessors=preds,
        explicit_options=options,feasible_pair_additivity_checks=pairs,
        family_maxima={k:max(r['maximum']for r in rows if r['kind']==k)for k in sorted({r['kind']for r in rows})},
        maximum=max(r['maximum']for r in rows))


def refinement(rows,target):
    fibers={}
    for t in C.TYPES:fibers.setdefault(C.representative(t),[]).append(t)
    reps={C.representative(('A'if r<=8 else'D',r))for r in C.ITEMS}
    residual=[t for t in C.TYPES if C.representative(t)in reps and t not in[('A',1),('A',2)]]
    assert all(C.actual_cross(a,b)==0 for a in residual for b in residual)
    out=[];checks=preds=options=pairs=singletons=0
    for old in rows:
        if old['upper']<=target:continue
        U.memory_guard();fixed=C.fixed_components(old['counts'],old['key'])+[('A',5)]*old['k5']
        groups=sorted(Counter(fixed).items());rem=old['remaining_rank']
        domain=sorted({tuple(sorted(t for group in choice for t in group))for choice in product(*(
            combinations_with_replacement(fibers[t],n)for t,n in groups))})
        domain=[fs for fs in domain if('D',13)in fs];actual=[]
        for fs in domain:
            co=Counter(fs);certs=[];cache={};best=None
            for a1 in range(rem-3):
                for a2 in range((rem-a1-4)//2+1):
                    modes=[(False,0)]if a1 or a2 or not co['D',14]else[(True,0),(False,1)]
                    for no3,res3 in modes:
                        cap=rem-a1-2*a2-4-3*res3
                        if cap<0:continue
                        anchor=list(fs)+[('A',1)]*a1+[('A',2)]*a2+[('A',4)]+[('A',3)]*res3
                        bc=basevalue(anchor);ts=[t for t in residual if t[1]<=cap and not(no3 and t==('A',3))]
                        ws=[basevalue(anchor+[t])-bc for t in ts];key=(tuple(ts),tuple(ws),cap)
                        if key not in cache:
                            tab,a,b=C.additive(ts,ws,cap);checks+=a;preds+=b
                            cache[key]=len(certs);certs.append(tab)
                            for i,a in enumerate(ts):
                                for j in range(i,len(ts)):
                                    b=ts[j]
                                    if a[1]+b[1]<=cap:assert basevalue(anchor+[a,b])==bc+ws[i]+ws[j];pairs+=1
                        ti=cache[key];tab=certs[ti];singletons+=len(ts)
                        v=tab['values'][cap]
                        if v is not None:
                            options+=1
                            if best is None or bc+v>best[0]:best=(bc+v,anchor,ti)
            assert best is not None
            z,att,ti=best;att=list(att);tab=certs[ti];left=len(tab['values'])-1
            while left:t=tab['types'][tab['last'][left]];att.append(t);left-=t[1]
            assert sum(t[1]for t in att)==57 and basevalue(att)==z
            actual.append(dict(fixed_actual=fs,tables=certs,maximum=F(z,64),attaining_forest=sorted(att)))
        assert actual
        out.append(dict(counts=old['counts'],k5=old['k5'],a12=old['a12'],rows=actual,
            maximum=max(r['maximum']for r in actual)))
    return dict(classes=out,actual_fixed_fibers=sum(len(c['rows'])for c in out),Bellman_inequalities=checks,
        attaining_predecessors=preds,explicit_options=options,feasible_pair_additivity_checks=pairs,
        singleton_context_checks=singletons,maximum=max((c['maximum']for c in out),default=F()))


def basefamily():
    rows=[];checks=preds=ordinary_options=quarter_options=typed_pairs=0
    for a in C.TYPES:
        for b in C.TYPES:
            rs={C.representative(a),C.representative(b)}
            allow=(80 if rs=={('A',4),('D',13)}else
                   12 if rs=={('A',1),('D',15)}else
                   24 if ('B',9)in rs and bool(rs&{('A',5),('E',6),('A',7),('A',8)})else 0)
            assert C.actual_cross(a,b)+E.Q.increment([a,b])+E.increment([a,b])+increment([a,b])<=C.C.upper_cross(C.representative(a),C.representative(b))+allow
            typed_pairs+=1
    for cs,used,key in C.fixed_domain():
        if not cs[4]:continue
        for k5 in range((57-used)//5+1):
            fixed=C.fixed_components(cs,key)+[('A',5)]*k5;rem=57-used-5*k5
            if rem<4:continue
            fc=C.upper_value(fixed)+24*cs[0]*k5
            ws={r:C.item(r,key,k5)+(12*cs[6]if r==1 else 24*cs[0]if r in(6,7,8)else 80*cs[4]if r==4 else 0)for r in C.ITEMS};tabs=[]
            for no2 in(False,True):
                ranks=[r for r in C.ITEMS if r<=rem and(not no2 or r!=2)]
                ts=[C.representative(('A'if r<=8 else'D',r))for r in ranks]
                tab,a,b=C.additive(ts,[ws[r]for r in ranks],rem);checks+=a;preds+=b;tab['exclude2']=no2;tabs.append(tab)
            opts=[tabs[0]['values'][rem]]
            for a1 in(1,2):
                og=a1==2 and bool(cs[1]);ng=bool(cs[0]and k5)
                if rem<a1 or not(og or ng)or(a1==2 and not og):continue
                for a2 in range((rem-a1)//2+1):
                    v=tabs[1]['values'][rem-a1-2*a2]
                    if v is not None:opts.append(a1*ws[1]+a2*ws[2]+v+(36*a2*(a2-1)if og else 0)+(64*a2*k5 if ng else 0))
            ordinary_options+=len(opts);order17=16*k5*key[9]+(8*k5*rem if cs[2]else 0)
            ordinary=fc+max(v for v in opts if v is not None)+order17;quarter=None
            if cs[1]and rem>=7:
                for a1 in range(1,rem-5):
                    for a4 in range(1,(rem-a1-2)//4+1):
                        for a2 in range(1,(rem-a1-4*a4)//2+1):
                            v=tabs[1]['values'][rem-a1-4*a4-2*a2]
                            if v is None:continue
                            z=fc+a1*ws[1]+a2*ws[2]+a4*ws[4]+v+order17+(36*a2*(a2-1)if a1>=2 else 0)+(64*a2*k5 if cs[0]else 0)+64*a2*a4
                            quarter=z if quarter is None else max(quarter,z);quarter_options+=1
            rows.append(dict(counts=cs,k5=k5,a12=key[9],key=key,remaining_rank=rem,tables=tabs,
                ordinary_upper=F(ordinary,64),quarter_upper=None if quarter is None else F(quarter,64),
                upper=F(max(ordinary,quarter if quarter is not None else -1),64)))
    target=F(1483,4);refined=refinement(rows,target)
    maximum=max(max(r['upper']for r in rows if r['upper']<=target),refined['maximum'])
    return dict(rows=rows,complete_classes=len(rows),actual_type_pair_domination_checks=typed_pairs,
        Bellman_inequalities=checks,attaining_predecessors=preds,
        ordinary_context_options=ordinary_options,quarter_context_options=quarter_options,
        coarse_maximum=max(r['upper']for r in rows),actual_refinement=refined,complete_base_family_upper=maximum)


def run():
    a=local();b=ideal();c=overlaps();d=basefamily();previous=F(V.load(E.TAG)['conditional_complete_finite_upper'])
    upper=max(previous,c['maximum'],d['complete_base_family_upper'])
    att=[('A',16),('D',17),('D',24)];assert value(att)==64*previous
    paths={Path(__file__),Path(E.__file__),ROOT/f'results/{E.TAG}.json',ROOT/f'results/{O.TAG}.json',
        ROOT/f'results/{O.M.TAG}.json',ROOT/f'results/{PARTIAL}.json',ROOT/f'results/{PARTIAL}.md'}
    return O.enc(dict(status='PASS_ORDER17_D13_SHORTLEAF_THREE_SIXTEENTH_SCOUT',
        sources={str(p.relative_to(ROOT)):U.digest(p)for p in sorted(paths)},local=a,ideal=b,later_positive_overlaps=c,
        complete_base_increase_family=d,complete_base_increase_family_checked=True,
        complete_finite_upper=upper,complete_finite_maximum=upper if upper==previous else None,
        unchanged_formal_attainer=att,
        universal_composition_complete=False,new_dimension_bound=None))


if __name__=='__main__':
    lock=(ROOT/'formal/.proof-batches/exclusive.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    signal.alarm(59);E.Q.N.A.R.M.R.resources();U.memory_guard()
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists();out=run();U.memory_guard();E.Q.N.A.R.M.R.resources(start=False)
    if args.verify:assert V.load(TAG)==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['ideal']['new_negative_count'],out['complete_base_increase_family']['complete_base_family_upper'],
        out['later_positive_overlaps']['maximum'],out['complete_finite_maximum'],flush=True)
