#!/usr/bin/env python3
"""Complete simultaneous finite/ideal frontier; all geometry is bound locally.

New residual weights require fresh exact Bellman certificates. Fixed A17,
residual D17, actual D16 and its A16-dominating allowance stay distinct.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import prod
from pathlib import Path
import fresh_d5_d12_alternate_coupon_rebalance_20261018 as X
import fresh_order15_bd_basic_fractional_transfer_20261018 as I
import fresh_d17_existing_three_corner_transfer_20261018 as V
import fresh_finite_d17_leaf_a16_half_transfer_20261018 as H
import fresh_order14_contextual_ambient59_frontier_20261018 as U

ROOT=X.ROOT;T=X.T;A=X.A;Q=X.Q;D=X.D;P=X.P;R=X.R
TAG='fresh_order15_d16_d17_ambient59_frontier_20261018'
D16='fresh_finite_d16_alternate_half_zero_transfer_20261018'
A17='fresh_a17_five_pair_alternate_transfer_20261018'
N=59;SCALE=16;ITEMS=X.ITEMS
LOCAL_MDS={I.TAG:'a8e10b9a35bd5485eb26593f8c1fedbb3c6375809246067cc59b35a157dda4ed',
 V.TAG:'03606f81300acb76cac429d80045de042b560cf216001059d1352666e56b837a',
 H.TAG:'56ccbce32eb877f3de0d3b1b4b474bd850d0cf48b1046ee67051eb4b17769a76',
 D16:'dd940e8d23642709a1af782dcca9457286abb9e8b7477875a24045ac5242e11c',
 A17:'7a05ac84e36a33c738f816899828a4c45a39219f9fe27840bb7aab5f03d205ae'}
STEMS=tuple(dict.fromkeys((*X.STEMS,X.TAG,*I.STEMS,I.TAG,V.TAG,H.TAG,D16,A17,U.TAG,U.FRONT,U.OLD)))
enc=R.enc
A17_DEBITS={(2,12):F(1),(4,14):F(1),(2,13):F(1,2),(3,14):F(1,2),(2,14):F(1)}


def actual_base(t):return X.actual_base(t)-(16 if t==('D',16) else 30 if t==('D',17) else 64 if t==('A',17) else 0)
def upper_base(t):return X.upper_base(t)-(8 if t==('D',16) else 30 if t==('D',17) else 64 if t==('A',17) else 0)


def actual_cross(a,b):
    pair={a,b}
    return X.actual_cross(a,b)+(16 if pair in ({('A',1),('A',15)},{('A',1),('A',16)})
        else 1 if pair=={('A',1),('D',14)} else 32 if pair=={('A',1),('A',14)}
        else 64 if pair=={('D',5),('A',12)} else 0)


def upper_cross(a,b):return X.upper_cross(a,b)+(1 if {a,b}=={('A',1),('D',14)} else 0)


def value(fs,upper=False):
    base,cross=(upper_base,upper_cross) if upper else (actual_base,actual_cross)
    return sum(base(t) for t in fs)+sum(cross(a,b) for i,a in enumerate(fs) for b in fs[i+1:])+R.context_cost(fs,upper)


def item(r,key,k5):return X.item(r,key,k5)-(30 if r==17 else 0)+(key[4] if r==1 else 0)


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
    assert actual_base(('D',16))==1756 and upper_base(('D',16))==actual_base(('A',16))==1764
    assert actual_base(('D',17))==upper_base(('D',17))==1838
    dom=P.fixed_domain();assert len(dom)==961 and len({k for cs,u,k in dom})==960
    kd=sorted({(key,k5) for cs,used,key in dom for k5 in range((N-used)//5+1)})
    tabkeys={(key,k5,e) for key,k5 in kd for e in (False,True)}
    tabs={(tuple(t['key']),t['fixed_A5_count'],t['exclude2']):t for t in saved['tables']} if saved else {
        (k,k5,e):table(k,k5,e) for k,k5,e in sorted(tabkeys)}
    assert set(tabs)==tabkeys and (saved is None or len(saved['tables'])==len(tabkeys))
    checks=sum(check_table(tabs[k,k5,e],k,k5,e) for k,k5,e in sorted(tabkeys))
    residual=[rep(('A' if r<=8 else 'D',r)) for r in ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    prior={tuple(r['counts']):r for r in docs[X.TAG]['finite']['alternatives']}
    rows=[];expanded=[];identity_checks=choices=0
    for cs,used,key in dom:
        large=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for _,r in large)==used and key[2]==cs[4] and key[4]==cs[5] and key[6]==cs[1]
        options5=[]
        for k5 in range((N-used)//5+1):
            fixed=large+[('A',5)]*k5;rem=N-used-5*k5
            for r,t in zip(ITEMS,residual):
                assert item(r,key,k5)==upper_base(t)+sum(upper_cross(t,b) for b in fixed);identity_checks+=1
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
                sixteenths_upper=upper,sixteenths_actual=value(fs),triggered=e,fixed_A2_count=j,reserved_A1_count=a)
            expanded.append(row);options5.append(row)
        best=max(options5,key=lambda x:(x['sixteenths_upper'],x['fixed_A5_count']))
        best['previous_sixteenths_upper']=prior[cs]['sixteenths_upper'];rows.append(best)
    maximum=max(r['sixteenths_upper'] for r in rows)
    witness=[('A',16),('A',16),('D',27)]
    assert sum(r for k,r in witness)==N and value(witness)==maximum
    above=[r for r in rows if r['sixteenths_upper']>6264]
    return dict(exact_maximum=F(maximum,16),curvature=F(783,2),actual_attainer=['A16','A16','D27'],
        actual_D16=F(actual_base(('D',16)),16),representative_rank16=F(upper_base(('D',16)),16),actual_A16=F(actual_base(('A',16)),16),
        actual_D17=F(actual_base(('D',17)),16),fixed_A17=F(actual_base(('A',17)),16),
        fixed_large_vectors=len(dom),original_interaction_keys=960,fixed_large_and_A5_vectors=len(expanded),interaction_keys=len(kd),
        tables=[tabs[k,k5,e] for k,k5,e in sorted(tabkeys)],simultaneous_type_pairs=len(types)**2,
        residual_item_checks=identity_checks,Bellman_upper_inequalities=checks,attaining_predecessors=len(tabkeys)*N,context_options=choices,
        alternatives=rows,expanded_alternatives=expanded,above_curvature=above,
        component_costs={name:F(actual_cross(*pair),16) for name,pair in {
            'A1_A14':(('A',1),('A',14)),'A1_A15':(('A',1),('A',15)),'A1_A16':(('A',1),('A',16)),
            'A1_D14':(('A',1),('D',14)),'D5_D12':(('D',5),('D',12)),
            'D5_A12':(('D',5),('A',12)),'A1_D5':(('A',1),('D',5)),
            'A2_D13':(('A',2),('D',13)),'A1_D13':(('A',1),('D',13))}.items()},
        fixed_rank17='A17',residual_rank17='D17',rank5_multiplicity_exhaustive=True,
        context_rule='k*(k-1)/4 only with a whole unmarked A10 and at least two whole unmarked A1 components; any rank10 is an upper witness',
        geometric_attainment_claimed=False,
        certificate_method='Fresh item-first generation; separate verifier checks every last-item inequality and attaining predecessor')


def ideal(docs):
    old=U.ideal(docs[U.FRONT]['ideal']);assert enc(old)==docs[U.TAG]['ideal']
    assert old['current_envelopes'][14]==F(92,15)
    env={**old['current_envelopes'],14:F(367,60)};rows=[]
    for r in old['rows']:
        rs=r['ranks'];saving=F(rs.count(14),60);cost=sum(env[x] for x in rs)
        assert r['normalized_cost']-cost==saving
        margin=r['curvature']/r['link_vertices']-cost
        assert margin==r['normalized_margin']+saving
        rows.append(dict(ranks=rs,link_vertices=r['link_vertices'],G=r['G'],curvature=r['curvature'],
            previous_normalized_cost=r['normalized_cost'],normalized_cost=cost,order15_saving=saving,
            previous_normalized_margin=r['normalized_margin'],normalized_margin=margin,absolute_margin=margin*r['link_vertices']))
    summary=U.summary(rows);assert summary['negative_count']==6 and summary['zero_count']==0
    assert next(r['normalized_margin'] for r in rows if r['ranks']==(8,9,9,9,9,14))==F(8,1125)
    plans=[{12:F(1),15:F(1,2),16:F(3,4)},{12:F(1),16:F(5,4)}]
    scenarios=[U.conditional(rows,p) for p in plans]
    assert all(p['summary']['negative_count']==p['summary']['zero_count']==0 for p in scenarios)
    return dict(complete_profiles=1890,complete_low9_domain=True,saved_polynomials_checked=True,fresh_convolution_performed=False,
        current_order15_family_costs=docs[I.TAG]['new_family_costs'],current_envelopes=env,rows=rows,summary=summary,
        previous_negative_count=7,conditional_scenarios=scenarios,
        coordinated_thresholds=[U.threshold(rows,r,{k:v for k,v in plans[0].items() if k!=r}) for r in plans[0]],
        negative_absolute_margin_not_extended_to_arbitrary_profiles=True)


def a15_coupon(g,blocks,marks):
    a,b=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if a is b:return []
    out=[]
    for sm,la,m,n in ((a,b,*marks),(b,a,marks[1],marks[0])):
        if sm['type']=='A1' and la['type']=='A15' and T.mark_pairs(g,la['nodes'],[n],15,'A')==[(2,),(12,)]:
            out.append(dict(small=sm,large=la,marks=[m,n],increment='1/2'))
    assert len(out)<=1;return out


def a17_coupons(g,blocks,marks):
    a,b=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if a is b:return []
    out=[]
    for sm,la,m,n in ((a,b,*marks),(b,a,marks[1],marks[0])):
        if sm['type']=='A1' and la['type']=='A14':
            images=T.mark_pairs(g,la['nodes'],[n],14,'A')
            if images in ([(2,),(11,)],[(3,),(10,)]):
                out.append(dict(rule='A1_A14',amount=F(1,2),small=sm,large=la,marks=[m,n],images=images))
        if sm['type']=='D5' and la['type']=='A12':
            si=T.mark_pairs(g,sm['nodes'],[m],5,'D');li=T.mark_pairs(g,la['nodes'],[n],12,'A')
            if si==[(3,),(4,)] and li==[(2,),(9,)]:
                out.append(dict(rule='D5_A12',amount=F(1),small=sm,large=la,marks=[m,n],images=[si,li]))
    assert len(out)<=1;return out


def literal(docs):
    before={R.key(r['source'],r['case']):r for r in docs[X.TAG]['combined_literal_interface']['rows']}
    local={s:{R.key(r['source'],r['case']):r for r in docs[s][field]['rows']}
        for s,field in ((I.TAG,'current_literal_interface'),(V.TAG,'combined_literal_interface'),
            (H.TAG,'current_literal_interface'),(D16,'combined_literal_interface'),(A17,'current_literal_interface'))}
    assert len(before)==216 and all(set(m)==set(before) and len(m)==216 for m in local.values())
    targets={r['kind']:r['target'] for r in docs[I.TAG]['rules']}
    rows=[];hits=Counter();local_comparisons=0
    for stem,index,g,white,cs in T.literal_inputs(docs):
        rowkey=R.key(stem,index);old=before[rowkey];assert len(cs)==len(old['corners']);corners=[]
        assert all(len(m[rowkey]['corners'])==len(cs) for m in local.values())
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks'];ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids))
            debits=[];credits=[];images={}
            if sig==[len(ids),0,0]:
                bs=D.finite_blocks(g,ids);b=next((b for b in bs if set(marks)<=set(b['nodes'])),None)
                if b and b['type'] in ('D16','D17','A17'):
                    n=int(b['type'][1:]);ps=T.mark_pairs(g,b['nodes'],marks,n,b['type'][0]);images[b['type']]=ps
                    rules=A17_DEBITS if b['type']=='A17' else {(2,14):F(1,2),(2,15):F(1,2)} if n==16 else {**V.DEBITS,**{p:F(1,2) for p in H.CHANGE}}
                    ds={rules.get(p,F()) for p in ps};assert len(ds)==1;debit=ds.pop()
                    if debit:debits.append(dict(rule=b['type'],amount=debit));hits[b['type']+'_debits']+=1
                for rule,cp in [('A1_A15',a15_coupon(g,bs,marks)),('A1_A16',H.coupon(g,bs,marks)),('A1_D14',I.coupon(g,ids,marks))]:
                    if cp:
                        amount=F(1,16) if rule=='A1_D14' else F(1,2)
                        assert len(cp)==1;credits.append(dict(rule=rule,amount=amount,roles=cp));hits[rule+'_credits']+=1
                for cp in a17_coupons(g,bs,marks):
                    credits.append(dict(rule=cp['rule'],amount=cp['amount'],roles=[cp]));hits[cp['rule']+'_credits']+=1
            elif sig[1]==0:
                block=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if block and len(block)==15:
                    labels=T.T.labels(g,block)
                    if labels is not None:
                        for kind,target in targets.items():
                            ps=sorted({tuple(sorted(mp[block.index(m)] for m in marks)) for mp in I.P.embeddings(labels,target)})
                            if ps:
                                images['ideal_'+kind]=ps;changed=I.BDROP if kind=='B' else I.DCHANGES
                                if set(ps)&changed:
                                    assert set(ps)<=changed;debits.append(dict(rule='ideal_'+kind,amount=F(1) if kind=='B' else F(1,16)));hits['ideal15_debits']+=1
            payment=F(old['corners'][j]['retained_payment_lower_bound'])-sum((z['amount'] for z in debits),F())+sum((z['amount'] for z in credits),F())
            assert payment>=0
            debit_by_rule={z['rule']:z['amount'] for z in debits}
            credit_by_rule={z['rule']:z['amount'] for z in credits}
            assert len(debit_by_rule)==len(debits) and len(credit_by_rule)==len(credits)
            native_D17=images.get('D17',[])
            node2_values={V.DEBITS.get(p,F()) for p in native_D17} or {F()}
            leaf_values={F(1,2) if p in H.CHANGE else F() for p in native_D17} or {F()}
            assert len(node2_values)==len(leaf_values)==1
            node2_debit=node2_values.pop();leaf_debit=leaf_values.pop()
            expected={
                I.TAG:(sum(v for k,v in debit_by_rule.items() if k.startswith('ideal_')),credit_by_rule.get('A1_D14',F())),
                V.TAG:(node2_debit,F()),H.TAG:(leaf_debit,credit_by_rule.get('A1_A16',F())),
                D16:(debit_by_rule.get('D16',F()),credit_by_rule.get('A1_A15',F())),
                A17:(debit_by_rule.get('A17',F()),credit_by_rule.get('A1_A14',F())+credit_by_rule.get('D5_A12',F()))}
            assert debit_by_rule.get('D17',F())==node2_debit+leaf_debit
            for s,(debit,increment) in expected.items():
                c=local[s][rowkey]['corners'][j];assert c['marks']==marks
                field='coupon_increment' if s==D16 else 'increment'
                assert F(c['debit'])==debit and F(c.get(field,'0'))==increment,(s,stem,index,j)
                local_comparisons+=1
            # The four other local deltas are independently zero on this domain.
            # Thus this is also a direct check of the leaf certificate's baseline.
            assert payment==F(local[H.TAG][rowkey]['corners'][j]['retained_payment_lower_bound'])
            corners.append(dict(marks=marks,previous_lower=old['corners'][j]['retained_payment_lower_bound'],images=images,
                debits=debits,credits=credits,retained_payment_lower_bound=payment))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216 and +hits==Counter(D17_debits=16,A1_A16_credits=8)
    return dict(configurations=216,changes=dict(+hits),rows=rows,minimum_total=min(r['total_lower_bound'] for r in rows),
        bound_local_corner_delta_comparisons=local_comparisons,local_keysets_and_marks_identical=True)


def run(saved=None):
    assert all(h!='PENDING' for h in LOCAL_MDS.values()),'Wait for all local reports to freeze.'
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    for s,h in LOCAL_MDS.items():assert sha256((ROOT/f'results/{s}.md').read_bytes()).hexdigest()==h,s
    assert docs[U.OLD]['dimension_bound']==59 and docs[U.TAG]['strict_global_constant']=='4350'
    assert docs[I.TAG]['new_coupon']['component_pair_increment']=='1/16'
    assert docs[H.TAG]['new_cross_component_coupons'][0]['component_pair_increment']=='1'
    assert docs[V.TAG]['component_saving']=='7/8' and docs[H.TAG]['total_saving']=='1'
    assert docs[D16]['new_cost']=='439/4' and docs[D16]['changed_pairs']==[[2,14],[2,15]]
    assert docs[D16]['typed_costs']['new_rank16_representative_upper']=='441/4'
    assert docs[D16]['new_coupon']['large_marks']==[2,12] and docs[D16]['new_coupon']['component_pair_increment']=='1'
    assert docs[A17]['new_cost']=='451/4' and docs[A17]['component_saving']=='4'
    assert docs[A17]['all_cross_upper_costs_unchanged'] and docs[A17]['all_ideal_weights_unchanged']
    assert docs[D16]['all_ideal_weights_unchanged'] and docs[H.TAG]['all_ideal_weights_unchanged']
    assert docs[V.TAG]['ideal_weights_unchanged']
    assert not set(V.CHANGES)&H.CHANGE
    weights={tuple(r['pair']):F(r['weight']) for r in docs[V.TAG]['weighted_pairs']}
    for p in H.CHANGE:assert weights[p]==F(1,2);weights[p]=F()
    assert sum(weights.values())==F(919,8)
    paths={Path(__file__),Path(X.__file__),Path(I.__file__),Path(V.__file__),Path(H.__file__),Path(U.__file__),
        Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__),Path(R.__file__)}
    paths.update(ROOT/f'{f}/{s}.{e}' for s in STEMS for f,e in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_ORDER15_D16_D17_AMBIENT59_FRONTIER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        ambient_dimension=59,link_dimension=58,completed_unrestricted_bound=59,new_dimension_bound=None,
        finite=finite(docs,saved['finite'] if saved else None),ideal=ideal(docs),combined_literal_interface=literal(docs),
        combined_D17_weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        combined_D16_weighted_pairs=docs[D16]['D16_weighted_pairs'],combined_A17_weighted_pairs=docs[A17]['weighted_pairs'],
        triangle_minimum=2,quadrilateral_minimum=1,finite_curvature='783/2',strict_global_constant='4350',constant_counted_once=True,
        scope='Compatible proved local transfers, complete fresh rank59 forest certificates and saved-polynomial1890 ideal update. Actual D16 versus representative A16 allowance and fixedA17/residualD17 are separate. No negative-margin extension, new ambient bound, geometric forest realization or global arithmetic optimization.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    out=run(saved)
    if args.verify:
        assert out==saved
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'above',len(out['finite']['above_curvature']),
        'ideal negative',out['ideal']['summary']['negative_count'],flush=True)
