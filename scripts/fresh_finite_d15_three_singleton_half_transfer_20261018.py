#!/usr/bin/env python3
"""Three whole D15 singleton units become halves using existing payments.

Re-enumerate the entire selected native-map domain, prove eleven arbitrary-
support cases, and reconstruct the current216 literal interfaces. No forest
optimization is performed here. In particular A15=101 still exceeds D15.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_rank15_sparse_payment_exchange_20261018 as S
import fresh_order12_order16_depth8_ambient60_frontier_20261018 as C
import fresh_d16_d17_eighth_balance_20261018 as B

R=S.R; ROOT=R.ROOT; T=R.T; A=R.A; Q=R.Q; P=R.P; D=R.D; X=R.X
TAG='fresh_finite_d15_three_singleton_half_transfer_20261018'
CHANGES={(0,11),(0,12),(1,12)}
WHOLE={1147,1269,1283,1622,1682,2501}
PROPER={1030,1284,1613,1623}
SATURATED=1641
EXPECTED=sorted(WHOLE|PROPER|{SATURATED})
INTERFACES=(R.TAG,R.E.TAG,S.U.E.TAG,S.J.TAG,S.A16,S.A15,S.U.TAG,S.TAG)
STEMS=tuple(dict.fromkeys((C.TAG,B.TAG,S.TAG,R.R.TAG,*C.STEMS,*S.STEMS)))


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x


def selected_support(rec,over):
    patterns=[c for c in over['isolation']['cases']
              if (c.get('component_order_bound') or 0)>=15 and len(c['marked_component'])<=15]
    assert len(patterns)==rec['complete_eligible_patterns']==536
    assert [c['index'] for c in patterns]==[r['index'] for r in rec['selected_patterns']]
    count=0;affected=[]
    for c,row in zip(patterns,rec['selected_patterns']):
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        nodes=c['marked_component'];graph=T.T.labels(g,nodes)
        maps=P.embeddings(graph,rec['target']) if graph is not None else []
        assert sorted(maps)==sorted(m['images'] for m in row['maps'])
        pairs={tuple(sorted(mp[nodes.index(i)] for i in c['preferred_corner'])) for mp in maps}
        assert sorted(pairs)==list(map(tuple,row['marked_images']))
        for m in row['maps']:
            assert m['marked_images']==sorted(m['images'][nodes.index(i)] for i in c['preferred_corner'])
        count+=len(maps)
        hits=[m for m in row['maps'] if tuple(m['marked_images']) in CHANGES]
        if hits:
            # A proper displayed component can have other target completions;
            # only these maps use a newly decreased actual marked pair.
            affected.append(dict(index=c['index'],maps=hits))
    assert count==870 and [r['index'] for r in affected]==EXPECTED
    assert sum(len(r['maps']) for r in affected)==24
    return dict(patterns=536,maps=870,native_search_rerun=True,affected=affected)


def forced_contacts(image,target,t):
    contacts=[(j,image[j],target[t][j]) for j in sorted(image) if target[t][j]]
    assert all(label==3 for _,_,label in contacts)
    return contacts


def proof(c,core,maps,target,d16):
    g,J,du,dv,end=X.constraints(c,core);n=len(g)
    left=c['left'][1];right=[c['right_copy_map'][i] for i in c['right'][1]]
    white=sorted(set(range(n))-set(left+right));nodes=c['marked_component']
    assert Q.inertia(g)==[n-1,1,0]
    assert Q.inertia(A.sub(g,white))==[len(white),0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[len(white)+1,0,0] for b in left+right)
    common=dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],white=white,
                selected_nodes=nodes,end_projection=end,core_inertia=[n-1,1,0])
    if c['index'] in WHOLE:
        assert n==18 and len(nodes)==15
        remaining=sorted(set(c['possible_contacts'])-set(nodes))
        assert sorted(J[i][i] for i in remaining)==[4,4,14]
        selected=T.mark_pairs(g,nodes,c['preferred_corner'],15,'D')
        assert selected and set(selected)<=CHANGES
        corners=[];alternatives=0
        for old in c['all_corners']:
            marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
            assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
            blocks=D.finite_blocks(g,ids)
            aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
            payment=F();reason=None;images=[];roles=[]
            if marks==c['preferred_corner']:
                payment=F(1,2);reason='selected_D15_half';images=selected
            elif aa is bb and aa['type']=='D16':
                images=T.mark_pairs(g,aa['nodes'],marks,16,'D')
                wanted=(1,13) if c['index'] in (1283,1622) else (1,12)
                assert images==[wanted] and d16[wanted]==F(1,2)
                payment=F(1,2);reason='existing_whole_D16_half';alternatives+=1
            elif aa is not bb and c['index'] in (1147,1269):
                for sm,la,m,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                    if sm['type']=='A2' and la['type']=='D13':
                        si=T.mark_pairs(g,sm['nodes'],[m],2,'A')
                        li=T.mark_pairs(g,la['nodes'],[lm],13,'D')
                        assert si==[(0,),(1,)] and li==[(10,)] and not payment
                        payment=F(1,2);reason='existing_whole_A2_D13_fork_center_half'
                        roles=[dict(type='A2',images=si),dict(type='D13',images=li)];alternatives+=1
            corners.append(dict(marks=marks,components=blocks,payment=payment,reason=reason,canonical_images=images,roles=roles))
        assert alternatives==1 and sum(r['payment'] for r in corners)==1
        return dict(**common,method='whole_D15_seals_all_displayed_factors',
            remaining_contacts=remaining,remaining_inverse_diagonals=[[i,J[i][i]] for i in remaining],
            selected_images=selected,corners=corners,total='1')
    witnesses=[]
    for m in maps:
        image={j:k for k,j in zip(nodes,m['images'])};missing=sorted(set(range(15))-set(image))
        contacts={t:forced_contacts(image,target,t) for t in missing}
        assert all(i in c['possible_contacts'] for cs in contacts.values() for _,i,_ in cs)
        if c['index'] in PROPER:
            assert n==17 and len(nodes)==14 and len(missing)==1
            t=missing[0];cs=contacts[t];assert len(cs)==1
            tj,i,label=cs[0];expected=F(6) if c['index'] in (1030,1613) else F(7,2)
            assert J[i][i]==expected>2
            assert (t in (13,14) and tj==12) if c['index'] in (1030,1613) else (t==0 and tj==1)
            witnesses.append(dict(images=m['images'],marked_images=m['marked_images'],missing_nodes=missing,
                required_contacts=cs,projection_lower_bound=expected,outcome='forced_new_white_projection_exceeds_norm2'))
            continue
        assert c['index']==SATURATED and n==16 and len(nodes)==12
        assert end['end_sets']==[[0,1,2,3],[10,11,12,13]] and end['both_ends_shortcut_bound']==4
        pair=tuple(m['marked_images'])
        if pair==(1,12):
            assert missing==[0,13,14]
            assert len(contacts[13])==len(contacts[14])==1 and contacts[13]==contacts[14]
            tj,i,label=contacts[13][0];assert tj==12 and i in (1,11) and J[i][i]==2
            # Both new leaves have sole unit contact i and no residual vector.
            # Their full product is therefore2, contradicting distinct acute walls.
            witnesses.append(dict(images=m['images'],marked_images=m['marked_images'],missing_nodes=missing,
                required_contacts=contacts,contact_node=i,projection_square='2',residual_square='0',
                forced_leaf_product='2',outcome='two_saturated_short_leaves_coincide'))
        else:
            assert pair==(0,11) and missing==[12,13,14]
            assert len(contacts[12])==1 and not contacts[13] and not contacts[14]
            tj,i,label=contacts[12][0];assert tj==11 and i in (1,11) and J[i][i]==2
            same=next(e for e in end['end_sets'] if i in e)
            opposite=next(e for e in end['end_sets'] if i not in e)
            assert all(J[i][j]>0 for j in same)
            assert [J[i][j] for j in opposite]==[-6,-3,-4,-2]
            assert max(J[i][j] for j in opposite)==-2
            assert target[12][13]==target[12][14]==3
            # The saturated center's sole old-root contact is mapped target11.
            # No other old root can supply either missing adjacent leaf.
            # Any new leaf has zero, positive, or <=-2 product with the center
            # under the one-end contact alternative, never the required -1.
            witnesses.append(dict(images=m['images'],marked_images=m['marked_images'],missing_nodes=missing,
                required_contacts=contacts,contact_node=i,projection_square='2',residual_square='0',
                opposite_end=opposite,opposite_inverse_row=[J[i][j] for j in opposite],
                required_center_leaf_product='-1',maximum_negative_product='-2',
                outcome='saturated_fork_center_has_no_possible_simple_leaf'))
    assert len(witnesses)==len(maps)
    return dict(**common,method='impossible_proper_completion' if c['index'] in PROPER else 'saturated_completion_exclusion',maps=witnesses)


def literal(docs):
    inputs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in C.R.B.C.L.STEMS}
    previous={(r['source'],json.dumps(r['case'],sort_keys=True)):r for r in docs[C.TAG]['combined_literal_interface']['rows']}
    rows=[];hits=Counter()
    for stem,index,g,white,cs in T.literal_inputs(inputs):
        prev=previous[stem,json.dumps(index,sort_keys=True)];corners=[]
        assert len(cs)==len(prev['corners'])
        for j,(marks,*_) in enumerate(cs):
            assert marks==prev['corners'][j]['marks'];ids=sorted(set(white)|set(marks))
            debit=F();images=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                block=next((b for b in D.finite_blocks(g,ids) if b['type']=='D15' and set(marks)<=set(b['nodes'])),None)
                if block:
                    images=T.mark_pairs(g,block['nodes'],marks,15,'D')
                    if set(images)&CHANGES:
                        assert set(images)<=CHANGES;debit=F(1,2);hits[str(images[0])]+=1
            value=F(prev['corners'][j]['retained_payment_lower_bound'])-debit;assert value>=0
            corners.append(dict(marks=marks,D15_images=images,previous_lower=prev['corners'][j]['retained_payment_lower_bound'],
                debit=debit,retained_payment_lower_bound=value))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==216 and hits=={'(0, 12)':8}
    return dict(configurations=216,changed_corners=dict(hits),rows=rows,all_totals_at_least_one=True)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[C.TAG]['finite']['exact_maximum']=='409'
    rule=next(r for r in docs[S.TAG]['rules'] if r['kind']=='D')
    weights={tuple(r['pair']):F(r['weight']) for r in rule['weighted_pairs']}
    assert sum(weights.values())==102 and all(weights[p]==1 for p in CHANGES)
    for p in CHANGES:weights[p]=F(1,2)
    assert sum(weights.values())==F(201,2)
    rec=next(r for r in docs[R.R.TAG]['records'] if r['kind']=='D');target=rec['target']
    assert target==rule['target']
    adj,dist,_=P.distance_data(target)
    assert sorted(dist[i][j] for i,j in CHANGES)==[11,11,12]
    assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES}==CHANGES for mp in P.embeddings(target,target))
    assert all(weights[tuple(p)]==1 for p in rec['retained_basic_pairs'])
    assert len(weights)==105 and all(weights[p]>=0 for p in combinations(range(15),2))
    d16={tuple(r['pair']):F(r['weight']) for r in docs[B.TAG]['D16_weighted_pairs']}
    assert d16[1,12]==d16[1,13]==F(1,2)
    coupon=docs[S.U.E.N.COUPON]
    assert [r['rank'] for r in coupon['small_components']]==[1,2,3,4,5]
    assert next(r for r in coupon['old_large_components'] if r['type']=='D13')['eligible_marks']==3
    assert docs[S.TAG]['upgraded_coupon']['large_marks']==[11,12] # center10 unchanged
    over=docs[A.OVER];support=selected_support(rec,over);proofs=[]
    for row in support['affected']:
        c=over['isolation']['cases'][row['index']]
        proofs.append(proof(c,over['isolation']['cores'][c['core_key']],row['maps'],target,d16))
    later=[]
    for stem in INTERFACES:
        for p in docs[stem].get('selected_proofs',[]):
            for corner in p.get('corners',[]):
                if not F(corner.get('payment',0)):continue
                for b in corner.get('components',[]):
                    if b['type']=='D15' and set(corner['marks'])<=set(b['nodes']):
                        im=T.mark_pairs(A.matrix(p['gram']),b['nodes'],corner['marks'],15,'D')
                        assert not set(im)&CHANGES,(stem,p['index'],im)
                        later.append(dict(source=stem,index=p['index'],images=im,payment=corner['payment']))
    cases=over['isolation']['cases']
    counts=[sum(0<(c.get('component_order_bound') or 0)<=14 for c in cases),
            sum(c.get('component_order_bound') is None for c in cases),
            sum((c.get('component_order_bound') or 0)>=15 and len(c['marked_component'])>15 for c in cases),536]
    assert len(cases)==2614 and counts==[1691,83,304,536] and sum(counts)==2614
    paths={Path(__file__),Path(R.__file__),Path(S.__file__),Path(C.__file__),Path(B.__file__),
           Path(T.__file__),Path(T.T.__file__),Path(A.__file__),Path(Q.__file__),Path(P.__file__),Path(D.__file__),Path(X.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_FINITE_D15_THREE_SINGLETON_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        target=target,changed_pairs=sorted(CHANGES),old_weight='1',new_weight='1/2',old_cost='102',new_cost='201/2',
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        complete_selected_support=support,selected_proofs=proofs,whole_cases=6,proper_projection_cases=4,saturated_case=1641,
        complete_high_overlap_routing=dict(total=2614,bound_at_most14=1691,protected_or_adaptive=83,displayed_larger_than15=304,eligible_rank15=536),
        later_positive_D15_interfaces=later,combined_literal_interface=literal(docs),
        existing_D16_half_pairs=[[1,12],[1,13]],existing_A2_D13_role='fork_center10',
        no_new_coupons=True,triangle_minimum=2,quadrilateral_minimum=1,all_other_weights_unchanged=True,
        rank15_actual_costs=dict(A15='101',B15='97',D15='201/2'),required_rank15_upper_base='101',
        finite_forest_recomputed=False,new_dimension_bound=None,
        scope='Fresh536-pattern/870-map native enumeration; six whole-factor seals, four strict projection exclusions and four saturated maps in one context; current216 literal corners and later positive interfaces. Universal actual-face routing is an explicit inherited input. No forest optimization or new dimension theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    result=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'D15',result['new_cost'],'contexts',len(result['selected_proofs']),'maps',result['complete_selected_support']['maps'],flush=True)
