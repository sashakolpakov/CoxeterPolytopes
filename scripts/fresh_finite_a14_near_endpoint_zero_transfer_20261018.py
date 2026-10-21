#!/usr/bin/env python3
"""Delete the A14 near-endpoint orbit using an existing whole A15 unit.

Fresh complete selected maps, arbitrary-contact seals and current native
literal interfaces. No coupon or ideal weight changes; no new dimension bound.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_finite_rank14_existing_payment_transfer_20261018 as L
import fresh_depth6_finite6_order13_order17_ambient59_frontier_20261018 as C

ROOT=L.ROOT; R=L.R; T=L.T; A=L.A; P=L.P; Q=L.Q; D=L.D; X=L.X
TAG='fresh_finite_a14_near_endpoint_zero_transfer_20261018'
A15='fresh_rank15_sparse_payment_exchange_20261018'
CHANGES={(0,12),(1,13)}
STEMS=tuple(dict.fromkeys((*L.STEMS,L.TAG,*C.STEMS,C.TAG,A15)))
enc=L.enc


def support(rec,over):
    patterns=[c for c in over['isolation']['cases']
              if (c.get('component_order_bound') or 0)>=14 and len(c['marked_component'])<=14]
    assert len(patterns)==rec['complete_eligible_patterns']==600
    assert [c['index'] for c in patterns]==[r['index'] for r in rec['selected_patterns']]
    rows=[];hits=[];total=0
    for c,old in zip(patterns,rec['selected_patterns']):
        g=A.matrix(over['isolation']['cores'][c['core_key']]['gram'])
        maps=T.T.mapped(g,c['marked_component'],c['preferred_corner'],rec['target'])
        assert maps==old['maps'];total+=len(maps)
        changed=[m for m in maps if tuple(m['marked_images']) in CHANGES]
        rows.append(dict(index=c['index'],maps=maps))
        if changed:
            assert len(changed)==len(maps)
            hits.append(dict(index=c['index'],maps=changed))
    assert total==552 and [r['index'] for r in hits]==[1157,1617,1640]
    assert [len(r['maps']) for r in hits]==[2,2,4]
    return dict(complete_patterns=600,complete_maps=552,fresh_native_search=True,
                rows=rows,affected=hits,affected_maps=8)


def whole(c,core,a15):
    g,white,left,right=L.basic(c,core);gg,J,du,dv,end=X.constraints(c,core)
    assert gg==g and c['index'] in (1157,1617)
    assert len(g)==18 and len(c['marked_component'])==14
    assert Q.inertia(g)==[17,1,0]
    remaining=sorted(set(c['possible_contacts'])-set(c['marked_component']))
    assert remaining==[2,3,12,13]
    assert [J[i][i] for i in remaining]==[14,4,14,4]
    corners=[]
    for old in c['all_corners']:
        marks=old['preferred_corner'];ids=sorted(set(white)|set(marks))
        assert Q.inertia(A.sub(g,ids))==[len(ids),0,0]
        blocks=D.finite_blocks(g,ids)
        shared=next((b for b in blocks if set(marks)<=set(b['nodes'])),None)
        payment=F();images=[];reason=None
        if marks==c['preferred_corner']:
            assert shared and shared['type']=='A14'
            images=T.mark_pairs(g,shared['nodes'],marks,14,'A')
            assert set(images)==CHANGES;reason='new_A14_zero'
        elif shared and shared['type']=='A15':
            images=T.mark_pairs(g,shared['nodes'],marks,15,'A')
            assert images==[(1,13)] and a15[1,13]==1
            payment=F(1);reason='retained_whole_A15_unit'
        corners.append(dict(marks=marks,components=blocks,canonical_images=images,
                            payment=payment,reason=reason))
    assert sum(r['payment'] for r in corners)==1
    assert sum(r['reason']=='retained_whole_A15_unit' for r in corners)==1
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],white=white,
        signature=[17,1,0],sealing=dict(method='whole_A14_positive_end_projection',
        whole_selected_nodes=c['marked_component'],end_projection=end,
        remaining_contacts=remaining,remaining_diagonals=[[i,J[i][i]] for i in remaining],
        minimum_nonzero_projection='4',maximum_further_core_contact_whites=0),
        corners=corners,total='1')


def literal(docs):
    prior={C.C.R.key(r['source'],r['case']):r
           for r in docs[C.TAG]['combined_literal_interface']['rows']}
    rows=[];changed=0
    for stem,index,g,white,cs in T.literal_inputs(docs):
        old=prior[C.C.R.key(stem,index)];corners=[]
        assert len(cs)==len(old['corners'])
        for j,(marks,*_) in enumerate(cs):
            assert marks==old['corners'][j]['marks']
            ids=sorted(set(white)|set(marks));images=[]
            if Q.inertia(A.sub(g,ids))==[len(ids),0,0]:
                block=next((b for b in A.components(g,ids) if set(marks)<=set(b)),None)
                if block and len(block)==14:
                    images=T.mark_pairs(g,block,marks,14,'A')
            debit=F(bool(set(images)&CHANGES))
            if debit:assert set(images)<=CHANGES;changed+=1
            pay=F(old['corners'][j]['retained_payment_lower_bound'])-debit
            assert pay>=0
            corners.append(dict(marks=marks,A14_images=images,debit=debit,
                retained_payment_lower_bound=pay))
        total=sum(r['retained_payment_lower_bound'] for r in corners);assert total>=1
        rows.append(dict(source=stem,case=index,corners=corners,total_lower_bound=total))
    assert len(rows)==len(prior)==216 and changed==0
    return dict(configurations=216,changed_corners=changed,native_Grams_reconstructed=True,
                rows=rows,minimum_total=min(r['total_lower_bound'] for r in rows))


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,d in docs.items():
        assert d['status'].startswith('PASS'),stem
        for p,h in d['sources'].items():
            assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(stem,p)
    arule=next(r for r in docs[A15]['rules'] if r['kind']=='A')
    a15={tuple(r['pair']):F(r['weight']) for r in arule['weighted_pairs']}
    assert a15[1,13]==1 and sum(a15.values())==101
    assert C.actual_base(('A',15))==16*101 and C.deleted(('A',15))==0
    assert docs[L.TAG]['new_costs']==dict(A14='91',B14='87',D14='181/2')
    rec=next(r for r in docs[R.TAG]['records'] if r['kind']=='A')
    over=docs[A.OVER];sup=support(rec,over);proofs=[]
    for hit in sup['affected']:
        c=over['isolation']['cases'][hit['index']];core=over['isolation']['cores'][c['core_key']]
        if len(c['marked_component'])==14:
            pr=whole(c,core,a15)
        else:
            assert c['index']==1640 and len(c['marked_component'])==13
            pr=L.proper(c,core,hit['maps'],rec['target'],'A')
            assert len(pr['maps'])==4
            assert all(r['projection_lower_bound']==F(7,2) for r in pr['maps'])
            assert {r['forced_node'] for r in pr['maps']}=={1,11}
        pr['affected_maps']=hit['maps'];proofs.append(pr)
    weights={p:F(p not in CHANGES) for p in combinations(range(14),2)}
    assert sum(weights.values())==89
    _,dist,features=P.distance_data(rec['target'])
    assert all(dist[i][j]==12 for i,j in CHANGES)
    assert all(weights[tuple(p)]==1 for p in rec['retained_basic_pairs'])
    assert all({tuple(sorted(mp[i] for i in p)) for p in CHANGES}==CHANGES
               for mp in rec['automorphisms'])
    assert all(weights[p]==1 for p in weights if dist[p[0]][p[1]]<=9)
    types=[t for t in C.C.P.Z.B.B.types() if t[1]<=59];rep=C.C.P.Z.B.B.representative
    def abase(t):return C.actual_base(t)-(32 if t==('A',14) else 0)
    def ubase(t):return C.upper_base(t)-(8 if t==('D',14) else 0)
    for a in types:
        assert ubase(rep(a))>=abase(a)>=0
        for b in types:assert C.upper_cross(rep(a),rep(b))>=C.actual_cross(a,b)>=0
    assert [F(abase((k,14)),16) for k in ('A','B','D')]==[89,87,F(181,2)]
    assert ubase(('D',14))==16*F(181,2)
    # Preserve the earlier rejected A2/A13 role shortcut explicitly. This
    # theorem uses a different corner, whose two marks are in whole A15.
    rejected=L.rejected_A14_shortcut(over)
    assert enc(rejected)==docs[L.TAG]['rejected_A14_shortcut']
    paths={Path(__file__),Path(L.__file__),Path(C.__file__),Path(R.__file__),Path(T.__file__),
           Path(T.T.__file__),Path(A.__file__),Path(P.__file__),Path(Q.__file__),Path(D.__file__),Path(X.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS
                 for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return enc(dict(status='PASS_FINITE_A14_NEAR_ENDPOINT_ZERO_EXISTING_A15_UNIT',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        target=rec['target'],automorphisms=rec['automorphisms'],changed_pairs=sorted(CHANGES),
        weighted_pairs=[dict(pair=p,weight=w) for p,w in sorted(weights.items())],
        new_costs=dict(A14='89',B14='87',D14='181/2'),rank14_representative_upper='181/2',
        complete_selected_support=sup,selected_proofs=proofs,
        complete_high_overlap_partition=docs[L.TAG]['complete_high_overlap_partition'],
        retained_A15_pair=[1,13],retained_A15_weight='1',new_cross_component_coupons=[],
        previous_rejected_A2_A13_shortcut=rejected,
        current_literal_interface=literal(docs),simultaneous_type_pairs=len(types)**2,
        all_ideal_weights_unchanged=True,all_other_finite_weights_unchanged=True,
        triangle_minimum=2,quadrilateral_minimum=1,finite_forest_recomputed=False,
        completed_unrestricted_bound=59,new_dimension_bound=None,
        scope='Complete600-pattern/552-map A14 search, two whole Lorentz seals and four forced proper-completion exclusions, current216 native literal configurations. Universal routing and independent short8 escape as explained in the report. No new coupon, forest certificate or dimension bound.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json'
    if not args.verify:assert not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['new_costs'],out['current_literal_interface']['changed_corners'],flush=True)
