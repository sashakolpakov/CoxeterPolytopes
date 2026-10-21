#!/usr/bin/env python3
"""Finite D17 cost119 ->235/2, with an explicit +3 coupon.

Validate the complete saved selected maps, reconstruct six radical proofs
and all216 current literal interfaces. No forest optimization here.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_rank15_sparse_payment_exchange_20261018 as R
import fresh_finite_d17_five_unit_reduction_20261018 as N

ROOT=R.ROOT;T=N.T;A=N.A;Q=N.Q;D=N.D;P=N.P
TAG='fresh_finite_d17_d5_d12_unit_transfer_20261018'
CHANGE={(4,14),(4,15),(4,16)}
IDEAL15='fresh_order15_bd_terminal_coupon_transfer_20261018'
CURRENT='fresh_rank15_order15_ambient60_frontier_20261018'
STEMS=tuple(dict.fromkeys((*R.STEMS,R.TAG,N.TAG,N.COUPON,IDEAL15,CURRENT)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def coupon(g,blocks,marks):
    aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
    if aa is bb:return []
    out=[]
    for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        if sm['type']!='D5' or la['type']!='D12':continue
        si=T.mark_pairs(g,sm['nodes'],[m],5)
        li=T.mark_pairs(g,la['nodes'],[lm],12)
        if si==[(3,),(4,)] and li in ([(9,)],[(10,),(11,)]):
            out.append(dict(small_type='whole D5',small_nodes=sm['nodes'],small_mark=m,
                small_images=si,large_type='whole D12',large_nodes=la['nodes'],large_mark=lm,
                large_images=li,old_weight='1/2',new_weight='1'))
    assert len(out)<=1
    return out


def selected_support(saved,over,target):
    expected=[c['index'] for c in over['isolation']['cases'] if c.get('component_order_bound')==17]
    assert [r['index'] for r in saved['rows']]==expected and len(expected)==140
    maps=0;affected=[]
    for row in saved['rows']:
        c=over['isolation']['cases'][row['index']]
        assert row['base_nodes']==c['marked_component'] and row['marks']==c['preferred_corner']
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        base=c['marked_component'];graph=T.T.labels(g,base);seen=set();pairs=set()
        for m in row['maps']:
            mp=m['images'];assert len(mp)==len(base)==len(set(mp)) and tuple(mp) not in seen
            assert all(0<=i<17 for i in mp);seen.add(tuple(mp));maps+=1
            assert graph is not None
            assert all(graph[i][j]==target[mp[i]][mp[j]] for i in range(len(base)) for j in range(len(base)))
            pair=tuple(sorted(mp[base.index(i)] for i in c['preferred_corner']))
            assert list(pair)==m['marked_images'];pairs.add(pair)
        assert sorted(pairs)==list(map(tuple,row['marked_images']))
        if pairs&CHANGE:
            assert pairs<=CHANGE and len(base)==17
            affected.append(dict(index=c['index'],maps=row['maps'],marked_images=row['marked_images']))
    assert maps==saved['complete_maps']==168
    assert [r['index'] for r in affected]==N.RADICAL_IDS
    assert sum(len(r['maps']) for r in affected)==12
    return dict(complete_patterns=140,validated_saved_maps=maps,affected=affected,
        induced_map_search_repeated=False)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[CURRENT]['completed_unrestricted_dimension_bound']==60
    assert docs[CURRENT]['current_exceptional_finite_costs']==dict(A15='101',B15='97',D15='102',A17='116')
    assert docs[IDEAL15]['new_family_costs']==dict(A=75,B=92,C=91,D=92)
    assert not docs[IDEAL15]['uses_finite_same_component_units']
    old={tuple(r['pair']):F(r['weight']) for r in docs[R.J.TAG]['D17_weighted_pairs']}
    assert len(old)==136 and sum(old.values())==119 and all(old[p]==F(1,2) for p in CHANGE)
    new={p:(F() if p in CHANGE else w) for p,w in old.items()}
    assert sum(new.values())==F(235,2)
    target=docs[N.TAG]['target'];adj,dist,_=P.distance_data(target)
    assert all(dist[i][j]>=10 for i,j in CHANGE)
    autos=docs[N.TAG]['automorphisms']
    assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGE}==CHANGE for mp in autos)
    assert all(new[tuple(p)]==1 for p in docs[N.TAG]['retained_basic_pairs'])
    assert new[1,15]==new[1,16]==F(1,2)
    over=docs[A.OVER];support=selected_support(docs[N.TAG]['complete_selected_support'],over,target)
    proofs=[]
    for i in N.RADICAL_IDS:
        c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
        proof=N.seal(c,core)
        assert enc(proof)==next(p for p in docs[N.TAG]['radical_transfers'] if p['index']==i)
        g=A.matrix(proof['gram']);corners=[]
        for r in proof['corners']:
            matches=coupon(g,r['components'],r['marks'])
            assert bool(matches)==bool(r['existing_coupons'])
            payment=F(1) if matches else F()
            corners.append(dict(marks=r['marks'],components=r['components'],
                selected=bool(F(r['selected_half'])),upgraded_coupons=matches,payment=payment))
        assert sum(r['payment'] for r in corners)==1
        assert sum(bool(r['payment']) for r in corners)==1
        assert all(not r['payment'] for r in corners if r['selected'])
        proof.update(corners=corners,new_coupon_cost='3',total='1')
        proofs.append(proof)
    # The upgrade is a subset of the existing whole-small<=5/D12 half rule.
    oldcoupon=docs[N.COUPON]
    assert [r['rank'] for r in oldcoupon['small_components']]==[1,2,3,4,5]
    assert next(r for r in oldcoupon['old_large_components'] if r['type']=='D12')['eligible_marks']==3
    assert 'one has rank<=5' in oldcoupon['rule'] and 'fork center or short leaf' in oldcoupon['rule']
    previous={(r['source'],N.normalize(r['case'])):r for r in docs[R.TAG]['literal_rows']}
    assert len(previous)==216
    literals=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=previous[stem,N.normalize(index)];assert len(cs)==len(prev['corners']);corners=[]
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks']
            ids=sorted(set(white)|set(marks));debit=F();inc=F();images=[];coupons=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                blocks=D.finite_blocks(g,ids)
                shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
                if shared and shared['type']=='D17':
                    images=T.mark_pairs(g,shared['nodes'],marks,17)
                    if set(images)&CHANGE:
                        assert set(images)<=CHANGE;debit=F(1,2);hits['D17_debits']+=1
                coupons=coupon(g,blocks,marks)
                if coupons:inc=F(1,2);hits['coupon_increments']+=1
            val=F(prev['corners'][j]['retained_payment_lower_bound'])-debit+inc
            assert val>=0
            corners.append(dict(marks=marks,whole_D17_images=images,debit=debit,
                coupon_increment=inc,upgraded_coupons=coupons,retained_payment_lower_bound=val))
        total=sum(c['retained_payment_lower_bound'] for c in corners);assert total>=1
        literals.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(literals)==216 and not hits
    # Explicitly retain the positive same-component outputs of later finite
    # exchanges. The six old selected halves above are the only replacements.
    interfaces=[]
    for stem in (R.R.TAG,R.R.E.TAG,R.U.E.TAG,R.J.TAG,R.A16,R.A15,R.TAG):
        for proof in docs[stem].get('selected_proofs',[]):
            if 'corners' not in proof:continue
            for corner in proof['corners']:
                if not F(corner.get('payment','0')):continue
                for b in corner['components']:
                    if b['type']=='D17' and set(corner['marks'])<=set(b['nodes']):
                        images=T.mark_pairs(A.matrix(proof['gram']),b['nodes'],corner['marks'],17)
                        assert not set(images)&CHANGE
                        interfaces.append(dict(source=stem,case=proof['index'],marks=corner['marks'],images=images))
    paths={Path(__file__),Path(R.__file__),Path(N.__file__),Path(T.__file__),Path(T.T.__file__),
        Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_FINITE_D17_D5_D12_UNIT_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        old_cost='119',new_cost='235/2',saving='3/2',zero_pairs=sorted(CHANGE),target=target,
        automorphisms=autos,weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(new.items())],
        retained_basic_pairs=docs[N.TAG]['retained_basic_pairs'],retained_D16_alternative_pairs=[[1,15],[1,16]],
        complete_selected_support=support,selected_proofs=proofs,
        upgraded_coupon=dict(small='whole D5',small_marks=[3,4],large='whole D12',large_marks=[9,10,11],
            old_weight='1/2',new_weight='1',changed_marked_pairs=6,component_pair_cost_increase='3',
            old_total_component_pair_cost='15/2',new_total_component_pair_cost='21/2',
            D5_long_arm_endpoint_unchanged=True),
        literal_configurations=216,literal_changed_corners=dict(hits),literal_rows=literals,
        retained_positive_interfaces=interfaces,triangle_minimum=2,quadrilateral_minimum=1,
        all_ideal_weights_unchanged=True,compatible_finite7_feature_depth=10,
        current_order15_maximum='92',new_dimension_bound=None,new_finite_volume_example=False,
        finite_forest_bound_proved=False,
        scope='Complete saved140-pattern/168-map validation, six reconstructed arbitrary-support radical seals, current216 literal and later finite-positive interfaces. D17three halves becomezero; six D5shortleaf/D12fork-role halves becomeunits, cost+3. Written universal priorities and nonselected reductions remain inherited inputs. No new map census, forest budget or dimension theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'D17',out['new_cost'],'coupon+3','literal',out['literal_changed_corners'],flush=True)
