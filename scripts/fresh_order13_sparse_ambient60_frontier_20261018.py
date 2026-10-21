#!/usr/bin/env python3
"""Count the full order13 family improvement with its new finite A4/B9 cost."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_balanced_rank16_order14_ambient60_frontier_20261018 as B
import fresh_order13_bc_sparse_terminal_transfer_20261018 as K
import fresh_order13_d_center_leaf_half_transfer_20261018 as H

ROOT=B.ROOT; Z=B.Z; DIM=60
TAG='fresh_order13_sparse_ambient60_frontier_20261018'
STEMS=tuple(dict.fromkeys((B.TAG,K.TAG,H.TAG,*B.STEMS)))


def base(t):return B.base(t)
def actual_cross(a,b):return B.actual_cross(a,b)+(8 if {a,b}=={('A',4),('B',9)} else 0)
def upper_cross(a,b):return B.upper_cross(a,b)+(8 if {a,b}=={('A',4),('B',9)} else 0)
def item(r,key):return B.item(r,key)+(8*key[5] if r==4 else 0)


def table(key):
    values=[0]+[-10**9]*DIM
    for r in Z.Z.ITEMS:
        w=item(r,key)
        for n in range(r,DIM+1):values[n]=max(values[n],values[n-r]+w)
    last=[0]+[next(r for r in Z.Z.ITEMS if r<=n and values[n]==values[n-r]+item(r,key)) for n in range(1,DIM+1)]
    return dict(key=list(key),values=values,last=last)


def check_table(t,key):
    assert t['key']==list(key)
    v,last=t['values'],t['last'];assert len(v)==len(last)==61 and v[0]==last[0]==0
    count=0
    for n in range(1,61):
        for r in Z.Z.ITEMS:
            if r<=n:assert v[n]>=v[n-r]+item(r,key);count+=1
        r=last[n];assert r in Z.Z.ITEMS and r<=n and v[n]==v[n-r]+item(r,key)
    return count


def finite(saved=None):
    rep=Z.Z.B.B.representative;types=[t for t in Z.Z.B.B.types() if t[1]<=DIM]
    for a in types:
        assert base(rep(a))>=base(a)
        for b in types:
            assert actual_cross(a,b)==actual_cross(b,a)
            assert upper_cross(rep(a),rep(b))>=actual_cross(a,b)
    assert actual_cross(('A',4),('B',9))==upper_cross(('A',4),('B',9))==8
    assert actual_cross(('D',5),('A',11))==16 and upper_cross(('A',5),('D',11))==76
    domain=Z.Z.large_domain();keys=sorted({key for cs,used,key in domain})
    assert len(domain)==1044 and len(keys)==792
    tables={tuple(t['key']):t for t in saved['tables']} if saved else {key:table(key) for key in keys}
    assert set(tables)==set(keys) and (saved is None or len(saved['tables'])==len(keys))
    checks=sum(check_table(tables[key],key) for key in keys)
    residual=[rep(('A' if r<=8 else 'D',r)) for r in Z.Z.ITEMS]
    assert all(upper_cross(a,b)==0 for a in residual for b in residual)
    rows=[];items=0
    for cs,used,key in domain:
        fixed=[('B' if r==9 else 'A' if r==17 else 'D',r) for r,c in zip(range(9,18),cs) for _ in range(c)]
        assert sum(r for k,r in fixed)==used and key[5]==cs[0]
        assert all(upper_cross(a,b)==0 for a in fixed for b in fixed)
        for r,t in zip(Z.Z.ITEMS,residual):
            assert item(r,key)==base(t)+sum(upper_cross(t,b) for b in fixed);items+=1
        rem=DIM-used;t=tables[key];value=sum(base(x) for x in fixed)+t['values'][rem]
        forest=fixed[:]
        while rem:
            r=t['last'][rem];forest.append(rep(('A' if r<=8 else 'D',r)));rem-=r
        assert sum(r for k,r in forest)==DIM
        upper=sum(base(t) for t in forest)+sum(upper_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        actual=sum(base(t) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        assert upper==value and actual<=upper
        rows.append(dict(counts=list(cs),key=list(key),forest=sorted(k+str(r) for k,r in forest),
            octupled_upper=value,octupled_actual=actual,deep7_deletions_included=sum(B.C.deep(t) for t in forest)))
    assert max(r['octupled_upper'] for r in rows)==3276
    witnesses=[]
    for ranks in ((16,16,28),(16,17,27),(17,17,26)):
        forest=[('D',r) for r in ranks]
        value=sum(base(t) for t in forest)+sum(actual_cross(a,b) for i,a in enumerate(forest) for b in forest[i+1:])
        assert value==3276;witnesses.append(dict(forest=['D'+str(r) for r in ranks],cost='819/2'))
    assert sum(r['octupled_upper']>3240 for r in rows)==15
    return dict(exact_maximum='819/2',previous_maximum='819/2',curvature='405',finite_deficit='9/2',
        complete_fixed_vectors=1044,interaction_keys=792,simultaneous_type_pairs=len(types)**2,
        residual_item_checks=items,Bellman_upper_inequalities=checks,attaining_predecessors=60*len(keys),
        tables=[tables[k] for k in keys],alternatives=rows,above_curvature=15,actual_attainers=witnesses,
        exceptional_costs=dict(A16='110',D16='441/4',A17='117',D17='469/4'),
        new_A4_B9_actual_cost='1',new_A4_B9_upper_allowance='1',finite_cost_unit='1/8',
        deep7_deletions_counted_in_upper=True,new_tables_generated=True,replay_optimizes=False)


def ideal(docs):
    previous=B.ideal(docs);assert previous==docs[B.TAG]['ideal'];rows=[]
    for old in previous['rows']:
        row=dict(old);saving=F(row['ranks'].count(12),26)
        row['normalized_cost']=str(F(row['normalized_cost'])-saving)
        row['normalized_margin']=str(F(row['normalized_margin'])+saving)
        row['absolute_margin']=str(row['link_vertices']*F(row['normalized_margin']))
        row['order13_normalized_saving']=str(saving);rows.append(row)
    summary=Z.Z.summarize(rows)
    assert len(rows)==2113 and summary['negative_count']==10 and summary['zero_count']==0
    assert summary['normalized_minimum']=='-1117/7680'
    conditional=Z.Z.conditional(rows,{11:F(1,2),15:F(5,4)})
    assert conditional['summary']['negative_count']==conditional['summary']['zero_count']==0
    assert conditional['summary']['normalized_minimum']=='73/104000'
    return dict(complete_profiles=2113,previous_negative_profiles=17,resolved_negative_profiles=7,
        rows=rows,summary=summary,remaining_conditional_target=conditional,low9_concentration_only=True,
        negative_terminal_absolute_margin_is_universal_bound=False)


def literal(docs):
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in B.C.L.STEMS}
    previous={(r['source'],B.C.L.N.normalize(r['case'])):r for r in docs[B.TAG]['combined_literal_interface']['rows']}
    targets={r['kind']:r['target'] for r in docs[K.TAG]['records']};targets['D']=docs[H.TAG]['target']
    changes={kind:{pair} for kind,pair in K.PAIRS.items()};changes['D']=H.CHANGE
    rows=[];changed=[];debits=credits=0
    for stem,index,g,white,cs in K.T.literal_inputs(inputs):
        key=(stem,B.C.L.N.normalize(index));before=previous[key];corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks))
            images={};debit=F();credit=F();coupons=[];sig=K.Q.inertia(K.A.sub(g,ids))
            if sig==[len(ids),0,0]:coupons=K.coupon(g,ids,marks);credit=F(len(coupons),2)
            elif sig[1]==0:
                shared=next((p for p in K.A.components(g,ids) if set(marks)<=set(p)),None)
                if shared and len(shared)==13:
                    graph=K.T.T.labels(g,shared)
                    if graph is not None:
                        for kind,target in targets.items():
                            pairs={tuple(sorted(mp[shared.index(m)] for m in marks)) for mp in K.P.embeddings(graph,target)}
                            images[kind]=sorted(pairs)
                            if pairs&changes[kind]:
                                assert pairs<=changes[kind];debit=F(1) if kind=='C' else F(1,2)
            value=F(before['corners'][j]['retained_payment_lower_bound'])+credit-debit;assert value>=0
            debits+=debit>0;credits+=credit>0
            corners.append(dict(marks=marks,canonical_images=images,order13_debit=debit,A4_B9_increment=credit,
                coupons=coupons,retained_payment_lower_bound=value))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        row=dict(source=stem,case=index,corners=corners,total_lower_bound=total);rows.append(row)
        if any(c['order13_debit'] or c['A4_B9_increment'] for c in corners):changed.append(row)
    assert len(rows)==216 and debits==credits==0
    return dict(configurations=216,debited_corners=debits,credited_corners=credits,changed=changed,rows=rows,
        inherited_D16_eighth_credits=docs[B.TAG]['combined_literal_interface']['inherited_D16_eighth_credits'])


def run(saved=None):
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[K.TAG]['simultaneous_order13_family_costs']==dict(A='65',B='145/2',C='72',D='71')
    assert docs[K.TAG]['new_coupon']['increment']=='1' and docs[H.TAG]['new_cost']=='71'
    assert docs[K.TAG]['triangle_minimum']==docs[H.TAG]['triangle_minimum']==2
    assert docs[K.TAG]['quadrilateral_minimum']==docs[H.TAG]['quadrilateral_minimum']==1
    assert docs[B.TAG]['completed_unrestricted_dimension_bound']==60
    assert F(5*sum(min(j,59-j) for j in range(60))-59*60,2)==405
    assert sum(2*j-59 for j in range(60) if 2*j>59)==900
    paths={Path(__file__),Path(B.__file__),Path(K.__file__),Path(H.__file__),Path(Z.__file__),
           Path(Z.Z.__file__),Path(Z.Z.B.__file__),Path(Z.Z.B.B.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    if saved:assert sources==saved['sources']
    return K.enc(dict(status='PASS_ORDER13_SPARSE_AMBIENT60_FRONTIER',sources=sources,
        completed_unrestricted_dimension_bound=60,ambient_dimension=60,link_dimension=59,
        finite=finite(saved['finite'] if saved else None),ideal=ideal(docs),combined_literal_interface=literal(docs),
        current_ideal_order13_costs=docs[K.TAG]['simultaneous_order13_family_costs'],
        current_ideal_order14_costs=docs[B.TAG]['current_ideal_order14_costs'],
        current_ideal_order17_costs=docs[B.TAG]['current_ideal_order17_costs'],
        strict_global_constant='4500',constant_counted_once=True,new_dimension_bound=None,new_finite_volume_example=False,
        scope='Complete simultaneous finite count includes the added A4/B9 cost1. All1044 fixed vectors/792 new Bellman tables,2113 ideal convolutions and216 current literal interfaces are reconstructed or checked. Conditional order12/order16 targets are not theorems. No new geometric census or dimension59 conclusion.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';saved=json.loads(path.read_text()) if args.verify else None
    out=run(saved)
    if args.verify:assert out==saved
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'finite',out['finite']['exact_maximum'],'negative',out['ideal']['summary']['negative_count'],flush=True)
