#!/usr/bin/env python3
"""Replace part of four sealed D5/D12 payments by distinct whole corners.

Exact local reconstruction, all current literal interfaces, and a complete
rank59 finite certificate with rank5 counts fixed before residual counting.
No new ambient dimension or geometric realization claim.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_rank13_contextual_coupon_rebalance_20261018 as R
import fresh_small_d16_coupon_and_ideal17_rebalance_20261018 as S
import fresh_order14_d_near_fork_half_transfer_20261018 as H

ROOT=R.ROOT;T=R.T;A=R.A;Q=R.Q;D=R.D;P=R.P
TAG='fresh_d5_d12_alternate_coupon_rebalance_20261018'
N=59;SCALE=16;DELTA=F(3,8)
STEMS=tuple(dict.fromkeys((*R.STEMS,R.TAG,*S.STEMS,S.TAG,*H.STEMS,H.TAG)))
enc=R.enc


def changes(g,white,marks):
    ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids))
    if sig!=[len(ids),0,0]:return dict(delta=F(),roles=[],components=[])
    blocks=D.finite_blocks(g,ids)
    aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
    delta=F();roles=[]
    if aa is not bb:
        for small,other,m,n in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
            if small['type']!='D5':continue
            sm=T.mark_pairs(g,small['nodes'],[m],5,'D')
            if sm!=[(0,)]:continue
            images=[];amount=F();role=None
            if other['type']=='D12':
                images=T.mark_pairs(g,other['nodes'],[n],12,'D')
                if set(images)=={(10,),(11,)}:amount=-DELTA;role='D5_long_D12_short_debit'
            elif other['type']=='A12':
                images=T.mark_pairs(g,other['nodes'],[n],12,'A')
                if set(images)=={(0,),(11,)}:amount=DELTA;role='D5_long_A12_endpoint_credit'
            elif other['type']=='A1':
                images=[(0,)];amount=DELTA;role='D5_long_A1_credit'
            if amount:
                delta+=amount;roles.append(dict(role=role,small=small,small_mark=m,small_images=sm,
                    other=other,other_mark=n,other_images=images,delta=amount))
    assert len(roles)<=1
    return dict(delta=delta,roles=roles,components=blocks)


def local_seals(docs):
    # Reconstruct the four full radical seals using the inherited exact
    # native Gram routine, and compare with the frozen seven-eighth rule.
    prior=S.seals(docs);assert enc(prior)==docs[S.TAG]['local_proofs']
    proofs=[]
    for old in prior['proofs']:
        core=old['old_seal'];g=A.matrix(core['gram']);white=core['white'];corners=[]
        assert Q.inertia(g)==[17,1,1]
        assert Q.inertia(A.sub(g,white))==[len(white),0,0]
        black=sorted(set(range(len(g)))-set(white))
        assert len(black)==4
        assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0] for b in black)
        for previous in old['corners']:
            ch=changes(g,white,previous['marks'])
            before=previous['ideal_payment']+previous['finite_coupon']
            after=before+ch['delta'];assert after>=0
            corners.append(dict(marks=previous['marks'],previous=before,**ch,payment=after))
        assert sorted(c['delta'] for c in corners)==[-DELTA,F(),F(),DELTA]
        assert sorted(c['payment'] for c in corners)==[F(),F(1,8),DELTA,F(1,2)]
        assert sum(c['payment'] for c in corners)==1
        positive=next(c for c in corners if c['delta']>0)
        wanted='A12' if old['index'] in (1978,2535) else 'A1'
        assert positive['roles'][0]['other']['type']==wanted
        proofs.append(dict(index=old['index'],old_seal=core,corners=corners,total='1'))
    assert [p['index'] for p in proofs]==[1978,2371,2535,2539]
    role=prior['disjoint_coupon_roles']
    assert role['this_packet']==[[0,9]] and role['preceding_D114']==[[0,10],[0,11]]
    assert role['finite_D17']==[[3,9],[3,10],[3,11],[4,9],[4,10],[4,11]]
    return dict(proofs=proofs,complete_inherited_selected_patterns=140,complete_inherited_maps=488,
        affected_maps=32,historical_map_search_repeated=False,disjoint_coupon_roles=role,
        ideal_weight_unchanged='1/8',finite_D5_long_D12_short_new_weight='1/2',
        new_D5_long_A12_endpoint_increment='3/8',new_D5_long_A1_increment='3/8',
        complete_contact_domain=[3,13],nonzero_contact_alternative='length2 path between the adjacent black sides of actual corner(3,13)',
        otherwise_all_extra_whites_orthogonal_to_full_core=True)


def literal(docs):
    old={R.key(r['source'],r['case']):r for r in docs[R.TAG]['combined_literal_interface']['rows']}
    h=docs[H.TAG]['current_literal_interface'];assert h['configurations']==216 and h['changed']==[]
    hrows={R.key(r['source'],r['case']):r for r in h['rows']};assert set(hrows)==set(old)
    rows=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=old[R.key(stem,index)];other=hrows[R.key(stem,index)];corners=[]
        assert len(prev['corners'])==len(other['corners'])==len(cs)
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks']==other['corners'][j]['marks']
            assert F(other['corners'][j]['debit'])==F(other['corners'][j]['increment'])==0
            ch=changes(g,white,marks);before=F(prev['corners'][j]['retained_payment_lower_bound'])
            after=before+ch['delta'];assert after>=0
            hits.update(x['role'] for x in ch['roles'])
            corners.append(dict(marks=marks,previous_lower=before,**ch,retained_payment_lower_bound=after))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216 and len(old)==len(hrows)==216
    assert hits==dict(D5_long_D12_short_debit=4,D5_long_A12_endpoint_credit=2,D5_long_A1_credit=4),dict(hits)
    return dict(configurations=216,changed_roles=dict(hits),rows=rows,
        changed=[r for r in rows if any(c['delta'] for c in r['corners'])],
        minimum_total=min(r['total_lower_bound'] for r in rows))


def actual_base(t):return R.actual_base(t)
def upper_base(t):return R.upper_base(t)


def actual_cross(a,b):
    value=R.actual_cross(a,b,True);pair={a,b}
    if pair=={('D',5),('D',12)}:value-=12
    if pair=={('D',5),('A',12)}:value+=12
    if pair=={('A',1),('D',5)}:value+=6
    return value


def upper_cross(a,b):
    value=R.upper_cross(a,b,True);pair={a,b}
    if pair=={('A',5),('D',12)}:value-=12
    if pair=={('A',1),('A',5)}:value+=6
    return value


def value(fs,upper=False):
    base,cross=(upper_base,upper_cross) if upper else (actual_base,actual_cross)
    return sum(base(t) for t in fs)+sum(cross(a,b) for i,a in enumerate(fs) for b in fs[i+1:])+R.context_cost(fs,upper)


ITEMS=tuple(r for r in P.ITEMS if r!=5)


def item(r,key,k5):
    assert r in ITEMS
    return R.item(r,key,True)+(6*k5 if r==1 else 0)


def table(key,k5,exclude2):
    ws={r:item(r,key,k5) for r in ITEMS if not exclude2 or r!=2}
    v=[0]+[-10**9]*N
    for r,w in ws.items():
        for n in range(r,N+1):v[n]=max(v[n],v[n-r]+w)
    last=[0]+[next(r for r,w in ws.items() if r<=n and v[n]==v[n-r]+w) for n in range(1,N+1)]
    return dict(key=list(key),fixed_A5_count=k5,exclude2=exclude2,values=v,last=last)


def check_table(t,key,k5,exclude2):
    assert t['key']==list(key) and t['fixed_A5_count']==k5 and t['exclude2']==exclude2
    v,last=t['values'],t['last'];assert len(v)==len(last)==N+1 and v[0]==last[0]==0
    ws={r:item(r,key,k5) for r in ITEMS if not exclude2 or r!=2};checks=0
    for n in range(1,N+1):
        for r,w in ws.items():
            if r<=n:assert v[n]>=v[n-r]+w;checks+=1
        r=last[n];assert r in ws and r<=n and v[n]==v[n-r]+ws[r]
    return checks


def finite(docs,saved=None):
    rep=P.Z.B.B.representative;types=[t for t in P.Z.B.B.types() if t[1]<=N]
    for a in types:
        assert upper_base(rep(a))>=actual_base(a)>=0
        for b in types:assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)==actual_cross(b,a)>=0
    dom=P.fixed_domain();assert len(dom)==961 and len({k for cs,u,k in dom})==960
    kd=sorted({(key,k5) for cs,used,key in dom for k5 in range((N-used)//5+1)})
    tabkeys={(key,k5,e) for key,k5 in kd for e in (False,True)}
    tabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in saved['tables']} if saved else {
        (k,k5,e):table(k,k5,e) for k,k5,e in sorted(tabkeys)}
    assert set(tabs)==tabkeys and (saved is None or len(saved['tables'])==len(tabkeys))
    checks=sum(check_table(tabs[k,k5,e],k,k5,e) for k,k5,e in sorted(tabkeys))
    residual=[rep(('A' if r<=8 else 'D',r)) for r in ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];expanded=[];items=choices=0
    previous={tuple(row['counts']):row for row in docs[R.TAG]['conditional_neighbor_preflight']['alternatives']}
    assert len(previous)==961
    for cs,used,key in dom:
        large=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for _,r in large)==used and key[2]==cs[4] and key[6]==cs[1]
        options5=[]
        for k5 in range((N-used)//5+1):
            fixed=large+[('A',5)]*k5;rem=N-used-5*k5
            for r,t in zip(ITEMS,residual):
                assert item(r,key,k5)==upper_base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
            normal=tabs[key,k5,False];no2=tabs[key,k5,True]
            options=[(normal['values'][rem],0,0,False)]
            if cs[1] and rem>=2:
                options += [(2*item(1,key,k5)+j*item(2,key,k5)+no2['values'][rem-2-2*j]+4*j*(j-1),j,2,True)
                    for j in range((rem-2)//2+1)]
            choices+=len(options);best,j,a,e=max(options);tab=tabs[key,k5,e]
            fs=fixed+[('A',2)]*j+[('A',1)]*a;rest=rem-2*j-a
            while rest:
                r=tab['last'][rest];fs.append(rep(('A' if r<=8 else 'D',r)));rest-=r
            fixed_cost=sum(upper_base(t) for t in fixed)+sum(upper_cross(a,b) for i,a in enumerate(fixed) for b in fixed[i+1:])
            upper=fixed_cost+best
            assert sum(r for _,r in fs)==N and value(fs,True)==upper and value(fs)<=upper
            row=dict(counts=list(cs),key=list(key),fixed_A5_count=k5,forest=sorted(k+str(r) for k,r in fs),
                sixteenths_upper=upper,sixteenths_actual=value(fs),triggered=e,
                fixed_A2_count=j,reserved_A1_count=a)
            expanded.append(row);options5.append(row)
        best=max(options5,key=lambda x:(x['sixteenths_upper'],x['fixed_A5_count']))
        rows.append(best)
    maximum=max(r['sixteenths_upper'] for r in rows)
    assert maximum==max(r['sixteenths_actual'] for r in rows)==6360
    attainers=[['D16','D17','D26'],['A17','D16','D26'],['D16','D16','D27']]
    assert all(sum(int(s[1:]) for s in fs)==N and value([(s[0],int(s[1:])) for s in fs])==maximum for fs in attainers)
    witnesses=[]
    for n1,n5,cost in ((3,4,F(391)),(8,3,F(783,2))):
        fs=[('A',1)]*n1+[('D',5)]*n5+[('D',12)]*3
        assert sum(r for _,r in fs)==N and F(value(fs),SCALE)==cost
        witnesses.append(dict(forest=sorted(k+str(r) for k,r in fs),rank=N,
            old=str(F(R.value(fs,False,True),SCALE)),new=str(cost)))
    old_above=sum(r['sixteenths_upper']>6264 for r in previous.values())
    above=[r for r in rows if r['sixteenths_upper']>6264]
    assert old_above==14 and len(above)==12
    return dict(exact_maximum=str(F(maximum,SCALE)),curvature='783/2',fixed_large_vectors=961,
        original_interaction_keys=960,fixed_large_and_A5_vectors=len(expanded),interaction_keys=len(kd),
        tables=[tabs[k,k5,e] for k,k5,e in sorted(tabkeys)],simultaneous_type_pairs=len(types)**2,
        residual_item_checks=items,Bellman_upper_inequalities=checks,attaining_predecessors=len(tabkeys)*N,
        context_options=choices,expanded_alternatives=expanded,alternatives=rows,
        previous_above_curvature=old_above,above_curvature=above,actual_attainers=attainers,witnesses=witnesses,
        component_costs={name:str(F(actual_cross(*pair),SCALE)) for name,pair in {
            'D5_D12':(('D',5),('D',12)),'D5_A12':(('D',5),('A',12)),
            'A1_D5':(('A',1),('D',5)),'A1_D13':(('A',1),('D',13))}.items()},
        A5_counts_exhaustive=True,residual_excludes_rank5=True,
        contextual_A2_pair_payment='1/8',context_trigger='whole unmarked A10 and at least two whole A1',
        neighbor_half_is_proved_by_bound_D14_packet=True,
        formal_forest_attainment_not_geometric_realization=True,replay_optimizes=False)


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[R.TAG]['conditional_neighbor_preflight']['exact_maximum']=='795/2'
    assert docs[H.TAG]['new_cost']=='80' and docs[H.TAG]['new_coupon']['component_pair_increment']=='1/2'
    paths={Path(__file__),Path(R.__file__),Path(S.__file__),Path(S.M.__file__),Path(S.M.O.__file__),
        Path(H.__file__),Path(T.__file__),Path(D.__file__),Path(A.__file__),Path(Q.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    local=local_seals(docs);lit=literal(docs)
    return dict(status='PASS_D5_D12_ALTERNATE_COUPON_REBALANCE',sources=sources,delta=str(DELTA),
        local_proofs=local,combined_literal_interface=lit,finite=finite(docs,saved['finite'] if saved else None),
        triangle_minimum=2,quadrilateral_minimum=1,ideal_weights_unchanged=True,
        current_order14_family_costs=docs[H.TAG]['current_order14_family_costs'],
        same_component_finite_weights_unchanged=True,all_old_broad_halves_preserved=True,
        completed_unrestricted_bound=59,new_dimension_bound=None,new_finite_volume_example=False,
        scope='Four exact full-core seals and all216 current literal interfaces; inherited complete selected/nonselected geometric scope and current coupon-output inventory. Complete rank59 finite certificate fixes every rank5 count and separately treats the contextual A2 interaction. No new geometry census, ideal convolution, dimension bound, or geometric forest realization.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    if not args.verify:assert not path.exists(),'Certificate already exists; use --verify.'
    out=enc(run(saved))
    if args.verify:assert out==saved
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'above',len(out['finite']['above_curvature']),
        'literal',out['combined_literal_interface']['changed_roles'],flush=True)
