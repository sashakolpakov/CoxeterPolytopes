#!/usr/bin/env python3
"""Finite B15/D15 cost103, using only existing payments.

The complete saved map domain and current literal interfaces are inputs;
forced-completion and whole-component arguments are recorded explicitly.
No new forest calculation or ambient bound is asserted.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_rank15_existing_half_transfer_20261018 as R
import fresh_finite_d17_a4_d13_half_transfer_20261018 as E

ROOT=R.ROOT;T=R.T;A=R.A;Q=R.Q;D=R.D;P=R.P;X=R.X
TAG='fresh_finite_b15_d15_existing_payment_transfer_20261018'
DEPTH='fresh_depth8_finite8_zero_transfer_20261018'
CHANGES={'B':{(0,13)},'D':{(1,13),(1,14)}}
PROPER=[1646,1666,1933,2326]
WHOLE=[1162,1175,1910,2303]
STEMS=tuple(dict.fromkeys((*E.STEMS,E.TAG,DEPTH)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def validate_support(rec,over):
    assert rec['complete_eligible_patterns']==536
    target=rec['target'];kind=rec['kind'];total=0;affected=[]
    for row in rec['selected_patterns']:
        c=over['isolation']['cases'][row['index']]
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        nodes=c['marked_component'];graph=T.T.labels(g,nodes)
        pairs=set();seen=set();hit=[]
        for m in row['maps']:
            mp=m['images'];assert len(mp)==len(nodes)==len(set(mp))
            assert all(0<=i<15 for i in mp) and tuple(mp) not in seen
            seen.add(tuple(mp));total+=1
            assert all(graph[i][j]==target[mp[i]][mp[j]] for i in range(len(mp)) for j in range(len(mp)))
            pair=tuple(sorted(mp[nodes.index(i)] for i in c['preferred_corner']))
            assert list(pair)==m['marked_images'];pairs.add(pair)
            if pair in CHANGES[kind]:hit.append(mp)
        assert sorted(pairs)==list(map(tuple,row['marked_images']))
        if hit:
            assert pairs<=CHANGES[kind]
            affected.append(dict(index=c['index'],maps=hit))
    assert total=={'B':240,'D':870}[kind]
    assert [r['index'] for r in affected]==([] if kind=='B' else sorted(PROPER+WHOLE))
    return dict(kind=kind,patterns=536,maps=total,affected=affected)


def proof(c,core,maps,target):
    g,J,du,dv,end=X.constraints(c,core)
    n=len(g);u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    white=sorted(set(range(n))-set(u+v));nodes=c['marked_component']
    assert Q.inertia(g)==[n-1,1,0]
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0] for b in u+v)
    if c['index'] in PROPER:
        assert n==17 and len(nodes)==14
        witnesses=[]
        for mp in maps:
            image={mp[j]:i for j,i in enumerate(nodes)}
            missing=sorted(set(range(15))-set(image));assert missing==[0]
            contacts=[(i,image[i],target[0][i]) for i in image if target[0][i]]
            assert len(contacts)==1
            ti,ci,label=contacts[0]
            assert ti==1 and label==3 and ci in c['possible_contacts']
            assert J[ci][ci]==F(7,2)>2
            assert ci in c['preferred_corner']
            witnesses.append(dict(images=mp,missing=missing,forced_old_target_node=ti,
                core_contact=ci,label=label,projection_lower_bound=J[ci][ci]))
        return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],
            white=white,method='impossible_proper_completion',end_projection=end,
            maps=witnesses,old_other_corner_nodes=end['other_old_corner_nodes'],
            conclusion='A required new supporting white has norm2 but projection at least7/2 unless its contacts cross both ends and give short8.')
    assert c['index'] in WHOLE and n==18 and len(nodes)==15
    remaining=sorted(set(c['possible_contacts'])-set(nodes))
    assert sorted(J[i][i] for i in remaining)==[4,4,14]
    selected=T.mark_pairs(g,nodes,c['preferred_corner'],15,'D')
    assert selected==[(1,13),(1,14)]
    corners=[];coupons=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids)
        aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
        value=F();reason=None;roles=[]
        if marks==c['preferred_corner']:
            value=F(1,2);reason='selected_D15_half'
        elif aa is not bb:
            for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if sm['type']=='A2' and la['type']=='D13':
                    si=T.mark_pairs(g,sm['nodes'],[m],2,'A')
                    li=T.mark_pairs(g,la['nodes'],[lm],13,'D')
                    assert si==[(0,),(1,)] and li==[(11,),(12,)]
                    assert not value
                    value=F(1,2);reason='existing_whole_A2_D13_short_fork_leaf_half'
                    roles=[R.S.marked_role(g,sm,m),R.S.marked_role(g,la,lm)]
                    coupons.append(dict(small=sm['type'],small_images=si,large=la['type'],large_images=li,weight='1/2'))
        corners.append(dict(marks=marks,components=blocks,payment=value,reason=reason,roles=roles))
    assert len(coupons)==1 and sum(r['payment'] for r in corners)==1
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],white=white,
        method='whole_D15_sealing',end_projection=end,selected_nodes=nodes,
        remaining_contact_nodes=remaining,remaining_inverse_diagonals=[[i,J[i][i]] for i in remaining],
        selected_images=selected,corners=corners,coupons=coupons,total='1')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[E.TAG]['new_cost']=='118'
    assert docs[DEPTH]['feature_depth']==8 and docs[DEPTH]['quadrilateral_minimum']==1
    assert [r['rank'] for r in docs[E.N.COUPON]['small_components']]==[1,2,3,4,5]
    assert next(r for r in docs[E.N.COUPON]['old_large_components'] if r['type']=='D13')['eligible_marks']==3
    oldrules={r['kind']:r for r in docs[R.TAG]['rules']}
    scouts={r['kind']:r for r in docs[R.R.TAG]['records']};over=docs[A.OVER]
    rules=[];supports=[];proofs=[]
    for kind in ('B','D'):
        rec=oldrules[kind];old={tuple(r['pair']):F(r['weight']) for r in rec['weighted_pairs']}
        assert len(old)==105 and sum(old.values())==104 and all(old[p]==1 for p in CHANGES[kind])
        new={p:(F() if kind=='B' else F(1,2)) if p in CHANGES[kind] else w for p,w in old.items()}
        assert sum(new.values())==103
        target=rec['target'];_,dist,_=P.distance_data(target)
        assert all(dist[i][j]>=10 for i,j in CHANGES[kind])
        assert all(new[tuple(p)]==1 for p in rec['retained_basic_pairs'])
        assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES[kind]}==CHANGES[kind] for mp in rec['automorphisms'])
        sup=validate_support(scouts[kind],over);supports.append(sup)
        for r in sup['affected']:
            c=over['isolation']['cases'][r['index']]
            proofs.append(proof(c,over['isolation']['cores'][c['core_key']],r['maps'],target))
        rules.append(dict(kind=kind,rank=15,target=target,automorphisms=rec['automorphisms'],
            changed_pairs=sorted(CHANGES[kind]),old_cost='104',new_cost='103',
            weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(new.items())],retained_basic_pairs=rec['retained_basic_pairs']))
    previous={(r['source'],E.N.normalize(r['case'])):r for r in docs[E.TAG]['literal_rows']}
    assert len(previous)==216
    literal=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(docs):
        prev=previous[stem,E.N.normalize(index)];rows=[]
        assert len(cs)==len(prev['corners'])
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks'];ids=sorted(set(white)|set(marks))
            debit=F();typ=None;images=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                block=next((b for b in D.finite_blocks(g,ids) if set(marks)<=set(b['nodes'])),None)
                if block and block['type'] in ('B15','D15'):
                    typ=block['type'];images=T.mark_pairs(g,block['nodes'],marks,15,typ[0])
                    if set(images)&CHANGES[typ[0]]:
                        assert set(images)<=CHANGES[typ[0]]
                        debit=F(1) if typ=='B15' else F(1,2);hits[typ]+=1
            value=F(prev['corners'][j]['retained_payment'])-debit;assert value>=0
            rows.append(dict(marks=marks,type=typ,images=images,debit=debit,retained_payment=value))
        total=sum(r['retained_payment'] for r in rows);assert total>=1
        literal.append(dict(source=stem,case=index,corners=rows,total=total))
    assert len(literal)==216
    # Later positive finite alternatives never spend one of these long pairs.
    for p in docs[R.TAG]['selected_proofs']+docs[R.E.TAG]['selected_proofs']+docs[E.TAG]['selected_proofs']:
        for c in p['corners']:
            if not F(c['payment']):continue
            for b in c['components']:
                if b['type'] in ('B15','D15') and set(c['marks'])<=set(b['nodes']):
                    images=T.mark_pairs(A.matrix(p['gram']),b['nodes'],c['marks'],15,b['type'][0])
                    assert not set(images)&CHANGES[b['type'][0]]
    routing=docs[R.TAG]['complete_high_overlap_routing']
    assert routing==dict(total=2614,bound_at_most14=1691,protected_or_adaptive=83,bound_at_least15=840,displayed_larger_than15=304,eligible_rank15=536)
    paths=[Path(__file__),Path(R.__file__),Path(E.__file__),Path(X.__file__),Path(T.__file__),Path(A.__file__),Path(Q.__file__),Path(D.__file__),Path(P.__file__)]
    paths +=[ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return enc(dict(status='PASS_FINITE_B15_D15_EXISTING_PAYMENT_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        new_costs=dict(B15='103',D15='103'),rules=rules,complete_selected_support=supports,
        selected_proofs=proofs,proper_completion_cases=len(PROPER),whole_component_cases=len(WHOLE),
        complete_high_overlap_routing=routing,literal_rows=literal,literal_configurations=216,
        literal_changed_corners=dict(hits),new_cross_component_coupons=[],
        compatible_D17_cost='118',compatible_feature_depth=8,
        all_other_finite_and_ideal_weights_unchanged=True,triangle_minimum=2,quadrilateral_minimum=1,
        finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Complete saved536-pattern B/D15 maps; four forced new-white projection exclusions, four whole-D15 arbitrary-support seals, and current216 simultaneous literal debits. Existing wholeA2/D13 short-leaf halves only. Universal and later routes are separate written inputs; no forest or ambient claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:
        assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['new_costs'],'literal',result['literal_changed_corners'],flush=True)
