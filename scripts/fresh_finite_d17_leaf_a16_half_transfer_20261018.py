#!/usr/bin/env python3
"""Delete two D17 leaf halves using retained A17 units and one typed coupon.

Exact inherited selected-map validation, four full seals, current216
literal interfaces and pointwise coupon absorption. No forest optimization.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_finite_d17_a4_d13_half_transfer_20261018 as L
import fresh_d5_d12_alternate_coupon_rebalance_20261018 as C
ROOT=L.ROOT;T=L.T;A=L.A;Q=L.Q;D=L.D;P=L.P
TAG='fresh_finite_d17_leaf_a16_half_transfer_20261018'
BALANCE='fresh_center_reversion_rank16_balance_ambient60_frontier_20261018'
CHANGE={(3,15),(3,16)};CASES=[1944,2032,2041,2337]
STEMS=tuple(dict.fromkeys((*C.STEMS,C.TAG,L.TAG,A.OVER,BALANCE)))
enc=L.enc


def coupon(g,blocks,marks):
    a,b=[next(b for b in blocks if m in b['nodes'])for m in marks];out=[]
    if a is b:return out
    for aa,bb,m,n in((a,b,*marks),(b,a,marks[1],marks[0])):
        if aa['type']!='A1'or bb['type']!='A16':continue
        images=T.mark_pairs(g,bb['nodes'],[n],16,'A')
        if images==[(3,),(12,)]:out.append(dict(small=aa,small_mark=m,large=bb,large_mark=n,
            large_images=images,old_weight='0',new_weight='1/2'))
    assert len(out)<=1
    return out


def selected(saved,over,target):
    expected=[c['index']for c in over['isolation']['cases']if c.get('component_order_bound')==17]
    assert len(expected)==140 and [r['index']for r in saved['rows']]==expected
    total=0;hits=[]
    for row in saved['rows']:
        c=over['isolation']['cases'][row['index']];base=c['marked_component']
        assert row['base_nodes']==base and row['marks']==c['preferred_corner']
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram']);labels=T.T.labels(g,base);seen=set();images=set()
        for m in row['maps']:
            mp=m['images'];assert len(mp)==len(base)==len(set(mp))and tuple(mp)not in seen
            assert all(0<=i<17 for i in mp);seen.add(tuple(mp));total+=1
            assert labels is not None and all(labels[i][j]==target[mp[i]][mp[j]]for i in range(len(base))for j in range(len(base)))
            pair=tuple(sorted(mp[base.index(i)]for i in c['preferred_corner']));assert list(pair)==m['marked_images'];images.add(pair)
        assert sorted(images)==list(map(tuple,row['marked_images']))
        if images&CHANGE:
            assert images<=CHANGE and len(base)==17
            hits.append(dict(index=c['index'],maps=row['maps'],marked_images=row['marked_images']))
    assert total==168 and [r['index']for r in hits]==CASES and sum(len(r['maps'])for r in hits)==8
    return dict(complete_patterns=140,validated_saved_maps=168,affected_maps=8,affected=hits,
        historical_map_search_repeated=False)


def seals(docs,weightsA):
    over=docs[A.OVER];out=[]
    for i in CASES:
        c=over['isolation']['cases'][i];core=over['isolation']['cores'][c['core_key']]
        proof=L.seal(c,core);assert enc(proof)==next(p for p in docs[L.TAG]['selected_proofs']if p['index']==i)
        g=A.matrix(proof['gram']);corners=[]
        for r in proof['corners']:
            matches=coupon(g,r['components'],r['marks']);block=next((b for b in r['components']if b['type']=='A17'and set(r['marks'])<=set(b['nodes'])),None)
            images=[];unit=F()
            if block:
                images=T.mark_pairs(g,block['nodes'],r['marks'],17,'A')
                assert set(images)=={(0,13),(3,16)}and all(weightsA[p]==1 for p in images);unit=F(1)
            oldhalf=F(1,2)if r['existing_coupons']else F()
            value=unit if i in(1944,2032)else oldhalf+sum((F(m['new_weight'])for m in matches),F())
            assert value>=0 and not(r['selected']and value)
            corners.append(dict(marks=r['marks'],components=r['components'],selected=r['selected'],
                existing_A4_D13_halves=r['existing_coupons'],new_A1_A16_coupons=matches,
                retained_whole_A17_images=images,retained_whole_A17_weight=str(unit),payment=value))
        assert sum(r['payment']for r in corners)==1
        assert sum(bool(r['new_A1_A16_coupons'])for r in corners)==(0 if i in(1944,2032)else 1)
        proof.update(corners=corners,total='1');out.append(proof)
    return out


def literals(docs):
    previous={C.R.key(r['source'],r['case']):r for r in docs[C.TAG]['combined_literal_interface']['rows']}
    assert len(previous)==216;rows=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        before=previous[C.R.key(stem,index)];assert len(before['corners'])==len(cs);corners=[]
        for j,(marks,*_)in enumerate(cs):
            assert marks==before['corners'][j]['marks'];ids=sorted(set(white)|set(marks));debit=F();images=[];coupons=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                blocks=D.finite_blocks(g,ids);b=next((b for b in blocks if b['type']=='D17'and set(marks)<=set(b['nodes'])),None)
                if b:
                    images=T.mark_pairs(g,b['nodes'],marks,17)
                    if set(images)&CHANGE:assert set(images)<=CHANGE;debit=F(1,2);hits['D17_debits']+=1
                coupons=coupon(g,blocks,marks)
            inc=sum((F(z['new_weight'])for z in coupons),F());hits['A1_A16_credits']+=bool(inc)
            val=F(before['corners'][j]['retained_payment_lower_bound'])-debit+inc;assert val>=0
            corners.append(dict(marks=marks,previous_lower=before['corners'][j]['retained_payment_lower_bound'],
                whole_D17_images=images,debit=debit,new_coupons=coupons,increment=inc,retained_payment_lower_bound=val))
        total=sum(r['retained_payment_lower_bound']for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216 and hits==dict(D17_debits=16,A1_A16_credits=8)
    return dict(configurations=216,changed_roles=dict(hits),rows=rows,minimum_total=min(r['total_lower_bound']for r in rows))


def positive_interfaces(docs):
    examined=[];hits=[]
    for stem,d in docs.items():
        for name in('selected_proofs','radical_transfers','whole_lorentz_cases','local_proofs'):
            ps=d.get(name,[])
            if isinstance(ps,dict):ps=ps.get('proofs',[])
            if not isinstance(ps,list):continue
            for proof in ps:
                if 'gram'not in proof:continue
                g=A.matrix(proof['gram'])
                for corner in proof.get('corners',[]):
                    value=corner.get('payment',corner.get('used_new_payment',corner.get('retained_payment_lower_bound','0')))
                    try:positive=F(value)>0
                    except(TypeError,ValueError):continue
                    if not positive or 'components'not in corner:continue
                    marks=corner['marks'];b=next((b for b in corner['components']if b['type']=='D17'and set(marks)<=set(b['nodes'])),None)
                    if not b:continue
                    images=T.mark_pairs(g,b['nodes'],marks,17)
                    r=dict(source=stem,array=name,index=proof.get('index'),marks=marks,images=images,payment=str(value));examined.append(r)
                    if set(images)&CHANGE:hits.append(r)
    assert [(r['source'],r['index'])for r in hits]==[(L.TAG,i)for i in CASES]
    return dict(examined_positive_whole_D17_corners=examined,only_changed_selected_proofs=hits,
        scope='Direct positive-corner arrays in bound local packets with stored Grams. This supplements, not replaces, the written complete universal/selected/later-output inventory.')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    old=docs[BALANCE]['local_proofs']['finite_weighted_pairs'];wa={tuple(r['pair']):F(r['weight'])for r in old['A17']};wd={tuple(r['pair']):F(r['weight'])for r in old['D17']}
    assert wa[0,13]==wa[3,16]==1 and all(wd[p]==F(1,2)for p in CHANGE)
    assert sum(wd.values())==F(C.actual_base(('D',17)),16)==F(467,4)
    new={p:(F()if p in CHANGE else w)for p,w in wd.items()};assert sum(new.values())==F(463,4)
    target=docs[L.TAG]['target'];adj,dist,_=P.distance_data(target)
    assert all(dist[i][j]>=10 for i,j in CHANGE)
    autos=docs[L.TAG]['automorphisms'];assert all({tuple(sorted(mp[i]for i in p))for p in CHANGE}==CHANGE for mp in autos)
    assert all(new[p]==w for p,w in wd.items()if p not in CHANGE)
    support=selected(docs[L.TAG]['complete_selected_support'],docs[A.OVER],target)
    proofs=seals(docs,wa);lit=literals(docs);interfaces=positive_interfaces(docs)
    rep=C.P.Z.B.B.representative;types=[t for t in C.P.Z.B.B.types()if t[1]<=59]
    for a in types:
        assert C.upper_base(rep(a))-(16 if rep(a)==('D',17)else 0)>=C.actual_base(a)-(16 if a==('D',17)else 0)
        for b in types:
            add=16 if {a,b}=={('A',1),('A',16)}else 0
            assert C.actual_cross(a,b)+add<=C.upper_cross(rep(a),rep(b))
    assert C.actual_cross(('A',1),('A',16))==0 and C.upper_cross(('A',1),('D',16))==24
    paths={Path(__file__),Path(L.__file__),Path(C.__file__),Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{f}/{s}.{e}'for s in STEMS for f,e in(('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_FINITE_D17_LEAF_A16_HALF_TRANSFER',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sorted(paths)},
        zero_pairs=sorted(CHANGE),target=target,automorphisms=autos,old_cost='467/4',new_cost='463/4',total_saving='1',
        weighted_pairs=[dict(pair=p,weight=w)for p,w in sorted(new.items())],retained_whole_A17_unit_pairs=[[0,13],[3,16]],
        new_cross_component_coupons=[dict(small='whole A1',small_nodes=[0],large='whole A16',large_nodes=[3,12],
            weight='1/2',component_pair_increment='1',old_actual_total='0',new_actual_total='1',existing_representative_allowance='3/2')],
        selected_support=support,selected_proofs=proofs,current_literal_interface=lit,positive_interface_audit=interfaces,
        unchanged_upper_cross_domination=dict(simultaneous_actual_type_pairs=len(types)**2,rank_limit=59,all_pass=True),
        all_ideal_weights_unchanged=True,triangle_minimum=2,quadrilateral_minimum=1,
        finite_forest_count_performed=False,new_dimension_bound=None,
        scope='Four inherited arbitrary-support full radical seals, complete saved140/168 selected-map validation, current216 literals, direct later positive-corner checks and unchanged upper coupon domination. Universal/protected/current-priority proofs remain bound inputs. No forest optimum, historical recensus or ambient theorem.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert json.loads(p.read_text())==out
    else:assert not p.exists();p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_cost'],out['current_literal_interface']['changed_roles'],flush=True)
