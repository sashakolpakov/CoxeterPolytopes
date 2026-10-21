#!/usr/bin/env python3
"""Finite D17 cost125 ->120, using only existing whole-component coupons.

Complete selected and literal interfaces; arbitrary-support scope is in
the written proof. No forest optimization or ambient theorem here.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018 as S

ROOT=S.ROOT;T=S.T;A=S.A;P=S.P;Q=S.Q;D=S.D
TAG='fresh_finite_d17_five_unit_reduction_20261018'
ZERO={(0,12),(0,13)}
HALF={(0,14),(0,15),(0,16),(4,14),(4,15),(4,16)}
RADICAL_IDS=[1297,1950,2063,2072,2078,2343]
CURRENT='fresh_sparse_coupon_order11_13_dimension62_20261018'
FRONT='fresh_ambient62_three_transfer_frontier_20261018'
D17='fresh_order17_e6_unit_coupon_zero_transfer_20261018'
COUPON='fresh_finite8_a17_d11_coupon_budget_20261017'
STEMS=tuple(dict.fromkeys((*S.STEMS,S.TAG,CURRENT,FRONT,D17,COUPON,
    'fresh_order13_d_center_pair_zero_transfer_20261018',
    'fresh_order16_affine_b_singleton_projection_transfer_20261018')))


def normalize(x):
    return tuple(x) if isinstance(x,list) else x


def seal(c,core):
    g=A.matrix(core['gram']);n=len(g)
    assert c['method']=='radical' and c['core_key']=='(25, 25, ((9, 9),))'
    assert n==19 and core['gauge']==[1]*n and Q.inertia(g)==[17,1,1]
    z=list(map(F,core['radical'][0]))
    assert all(sum(g[i][j]*z[j] for j in range(n))==0 for i in range(n))
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    black=set(u+v);white=sorted(set(range(n))-black)
    assert len(white)==15 and Q.inertia(A.sub(g,white))==[15,0,0]
    assert all(Q.inertia(A.sub(g,sorted(set(white)|{b})))==[16,0,0] for b in black)
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for k in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[k]==c[k]
    assert len(c['marked_component'])==17
    images=T.mark_pairs(g,c['marked_component'],c['preferred_corner'],17)
    assert images and set(images)<={(4,14),(4,15),(4,16)}
    allowed=c['possible_contacts'];assert all(z[i] for i in allowed)
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v);shortcuts=[]
    for i,j in combinations(allowed,2):
        if z[i]*z[j]<0:
            length=min(du[i]+dv[j],du[j]+dv[i])+2
            assert length<=8
            shortcuts.append(dict(nodes=[i,j],length=length))
    assert shortcuts
    corners=[];coupon_count=0
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[17,0,0]
        blocks=D.finite_blocks(g,ids)
        selected=marks==c['preferred_corner'];coupons=[]
        if not selected:
            aa,bb=[next(b for b in blocks if m in b['nodes']) for m in marks]
            for small,large,sm,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
                if small['type']!='D5' or large['type']!='D12':continue
                assert small is not large
                role=large['roles'][str(lm)]
                assert role in ('fork_center','short_fork_leaf')
                assert S.marked_role(g,small,sm)['canonical_images']==[3,4]
                canonical=T.mark_pairs(g,large['nodes'],[lm],12)
                assert canonical==([(9,)] if role=='fork_center' else [(10,),(11,)])
                coupons.append(dict(small_type='whole D5',small_nodes=small['nodes'],small_mark=sm,
                    small_role='short fork leaf',large_type='whole D12',large_nodes=large['nodes'],
                    large_mark=lm,large_role=role,large_canonical_images=canonical,weight='1/2'))
        assert len(coupons)<=1
        coupon_count+=len(coupons)
        value=F(1,2) if selected or coupons else F()
        corners.append(dict(marks=marks,components=blocks,selected_half='1/2' if selected else '0',
                            existing_coupons=coupons,payment=str(value)))
    assert coupon_count==1 and sum(F(r['payment']) for r in corners)==1
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],radical=list(map(str,z)),
                white=white,selected_images=images,allowed_contacts=allowed,
                mixed_sign_shortcuts=shortcuts,corners=corners,total='1',new_coupon_cost='0')


def literal_check(data,old_weights,new_weights):
    old={(r['source'],normalize(r['case'])):r for r in data[S.TAG]['literal_interfaces']['rows']}
    assert len(old)==216
    out=[];affected=[];d17_halves=[]
    ideal_target=data[D17]['target']
    ideal_weights={tuple(r['pair']):F(r['weight']) for r in data[D17]['weighted_pairs']}
    for stem,index,g,white,cs in T.literal_inputs(data):
        previous=old[(stem,normalize(index))];corners=[]
        assert len(cs)==len(previous['corners']) and len(cs) in (3,4)
        for j,(marks,unit,_) in enumerate(cs):
            baseline=previous['corners'][j]
            assert marks==baseline['marks']
            ids=sorted(set(white)|set(marks));sig=Q.inertia(A.sub(g,ids))
            images=[];debit=F()
            if sig==[len(ids),0,0]:
                parts=A.components(g,ids)
                shared=next((p for p in parts if set(marks)<=set(p)),None)
                if shared and len(shared)==17:
                    images=T.mark_pairs(g,shared,marks,17)
                    if images:
                        values={old_weights[p]-new_weights[p] for p in images}
                        assert len(values)==1;debit=values.pop()
            value=F(baseline['retained_weight'])-debit
            assert value>=0
            corners.append(dict(marks=marks,old_payment=baseline['retained_weight'],
                                whole_D17_images=images,debit=str(debit),retained_payment=str(value)))
        total=sum(F(c['retained_payment']) for c in corners)
        assert total>=1,(stem,index,total)
        row=dict(source=stem,case=index,corners=corners,total=str(total));out.append(row)
        if any(F(c['debit']) for c in corners):
            affected.append(row)
            if stem=='fresh_d17_partial_weight_transfer_20261017':
                # The older partial-weight literal uses one ideal half. Prove
                # it is not one of the four newly zeroed E6-coupon orbit pairs.
                hits=[]
                for marks,unit,_ in cs:
                    ids=sorted(set(white)|set(marks));p=T.T.labels(g,ids)
                    maps=[] if p is None else P.embeddings(p,ideal_target)
                    if maps:
                        pairids=[ids.index(m) for m in marks]
                        pairs=sorted({tuple(sorted(mp[i] for i in pairids)) for mp in maps})
                        assert all(ideal_weights[pair]==F(1,2) for pair in pairs)
                        hits.append(dict(marks=marks,canonical_pairs=pairs,retained_current_weight='1/2'))
                assert len(hits)==1
                d17_halves.append(dict(case=index,ideal_corner=hits[0]))
    assert len(out)==216 and len(affected)==9 and len(d17_halves)==4
    assert sorted(r['case'] for r in affected if r['source']==T.T.A.OLD17)==[393,647,662,675,828]
    assert sorted(r['case'] for r in d17_halves)==[1974,1978,2504,2535]
    assert min(F(r['total']) for r in out)==1
    return dict(complete_configurations=216,affected_configurations=9,rows=out,
                affected=affected,current_D17_ideal_halves=d17_halves,minimum_total='1')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,doc in docs.items():
        assert doc['status'].startswith('PASS'),stem
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    assert docs[CURRENT]['dimension_bound']==62 and docs[FRONT]['inherited_completed_ceiling']==62
    assert [r['rank'] for r in docs[COUPON]['small_components']]==[1,2,3,4,5]
    assert next(r for r in docs[COUPON]['old_large_components'] if r['type']=='D12')['eligible_marks']==3
    assert 'one has rank<=5' in docs[COUPON]['rule']
    assert 'fork center or short leaf' in docs[COUPON]['rule']
    old=next(r for r in docs[S.TAG]['rules'] if (r['kind'],r['rank'])==('D',17))
    weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values())==125
    new=dict(weights)
    assert all(weights[p]==1 for p in ZERO|HALF)
    for p in ZERO:new[p]=F()
    for p in HALF:new[p]=F(1,2)
    assert sum(new.values())==120
    graph=T.target(17);adj,dist,_=P.distance_data(graph);autos=P.embeddings(graph,graph)
    for group in (ZERO,HALF):
        assert all({tuple(sorted(mp[i] for i in p)) for p in group}==group for mp in autos)
    assert all(dist[i][j]>=10 for i,j in ZERO|HALF)
    basic={p for p in weights if dist[p[0]][p[1]]<=8 or
           (dist[p[0]][p[1]]==9 and any(len(adj[i])==1 for i in p))}
    assert all(new[p]==1 for p in basic) and not basic&(ZERO|HALF)
    over=docs[A.OVER];support=T.selected(over,17)
    assert support['complete_eligible_patterns']==140 and support['complete_maps']==168
    hits=[r for r in support['rows'] if set(map(tuple,r['marked_images']))&(ZERO|HALF)]
    assert [r['index'] for r in hits]==RADICAL_IDS
    assert all(set(map(tuple,r['marked_images']))<={(4,14),(4,15),(4,16)} for r in hits)
    assert sum(len(r['maps']) for r in hits)==12
    proofs=[seal(over['isolation']['cases'][i],
                 over['isolation']['cores'][over['isolation']['cases'][i]['core_key']]) for i in RADICAL_IDS]
    literals=literal_check(docs,weights,new)
    # The reused whole D5/D12 coupon is an unchanged member of this complete
    # small<=5 / D11..16 fork-center-or-short-leaf rule.
    paths=[Path(__file__),Path(S.__file__),Path(T.__file__),Path(A.__file__),
           Path(P.__file__),Path(Q.__file__),Path(D.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_FINITE_D17_FIVE_UNIT_COUPON_FREE_REDUCTION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        target=graph,automorphisms=autos,zero_pairs=sorted(ZERO),half_pairs=sorted(HALF),
        old_cost='125',new_cost='120',total_saving='5',
        weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(new.items())],
        retained_basic_pairs=sorted(basic),complete_selected_support=support,
        affected_selected_contexts=RADICAL_IDS,affected_maps=12,radical_transfers=proofs,
        literal_interfaces=literals,existing_coupon_rule=docs[COUPON]['rule'],
        new_finite_coupons=[],all_other_pair_weights_unchanged=True,
        baseline_ideal_rule=CURRENT,compatible_three_transfer_rule=FRONT,
        triangle_minimum=2,quadrilateral_minimum=1,finite_budget_recomputed=False,
        new_dimension_bound=None,
        scope='Complete140-pattern/168-map selected domain, six sealed full radical contexts and216 current literal configurations. Existing whole D5/D12 half coupons suffice. Written arbitrary-support and universal/priority composition inherited and checked; no finite-forest maximum or ambient theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists(),'Existing output; use --verify.'
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'cost',out['new_cost'],'maps',out['affected_maps'],
          'literal',out['literal_interfaces']['affected_configurations'],flush=True)
