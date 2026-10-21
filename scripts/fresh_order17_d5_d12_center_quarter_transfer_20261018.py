#!/usr/bin/env python3
"""Two affine D16 halves become quarters using a whole D5/D12 3/4 coupon.

Complete selected native-map reconstruction and two full-core seals.
The added finite cost is exactly1/4; a global count remains separate.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018 as M

ROOT=M.ROOT;O=M.O;D=M.D;P=M.P;Q=M.Q;L=M.L
TAG='fresh_order17_d5_d12_center_quarter_transfer_20261018'
FINITE='fresh_finite_d17_d5_d12_unit_transfer_20261018'
CHANGE={(0,10),(2,12)}
EXPECTED=[1314,2529]
STEMS=(M.TAG,O.TAG,O.SMALL,O.OVER,FINITE)


def center_role(g,b,mark):
    ns=b['nodes'];assert b['type']=='D12' and len(ns)==12
    adj={i:[j for j in ns if i!=j and g[i][j]] for i in ns}
    assert all(g[i][i]==2 for i in ns)
    assert all(g[i][j]==-1 for i in ns for j in adj[i])
    assert sum(map(len,adj.values()))==22
    centers=[i for i in ns if len(adj[i])==3];assert centers==[mark]
    ds=Q.graph_dist(L.sub(g,ns),[ns.index(mark)])
    assert all(x is not None for x in ds)
    assert sorted(ds[ns.index(i)] for i in ns if len(adj[i])==1)==[1,1,9]
    return dict(type='D12',nodes=ns,actual_mark=mark,role='fork_center',
                eligible_marks=[mark],leaf_arm_lengths=[1,1,9])


def inspect(c,core,old):
    inherited=O.inspect(c,core);assert inherited==old
    g=[[F(x) for x in row] for row in core['gram']];n=len(g);assert n==19
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[key]==c[key]
    assert c['possible_contacts']==[3,13] and set(c['marked_component'])==set(range(n))-{3,13}
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v)
    assert min(du[3]+dv[13],du[13]+dv[3])+2==2
    corners=[];selected=None;other=None;roles=None
    for row in inherited['corners']:
        if row['kind']=='affineD16_order17':
            selected=row['marks']
            corners.append(dict(marks=selected,kind=row['kind'],ideal_payment='1/4',finite_coupon='0'))
            continue
        assert row['signature']==[17,0,0]
        blocks=D.finite_blocks(g,sorted(set(inherited['white'])|set(row['marks'])))
        assert [(b['nodes'],b['type']) for b in blocks]==[(b['nodes'],b['type']) for b in row['components']]
        matches=[]
        for r in row['coupon_roles']:
            assert r['small_type']=='D5' and r['large_type']=='D12' and r['large_role']=='fork_center'
            sm=next(b for b in blocks if r['small_mark'] in b['nodes'])
            la=next(b for b in blocks if r['large_mark'] in b['nodes'])
            assert sm is not la and all(not g[i][j] for i in sm['nodes'] for j in la['nodes'])
            matches.append(dict(small=M.role(g,sm,r['small_mark']),large=center_role(g,la,r['large_mark']),whole_components=True))
        assert len(matches)<=1
        if matches:
            assert other is None and F(row['finite_half'])==F(1,2)
            other=row['marks'];roles=matches[0]
        corners.append(dict(marks=row['marks'],kind='finite',components=blocks,ideal_payment='0',
            finite_coupon='3/4' if matches else '0',same_component_payment_used=False,coupon_roles=matches))
    assert selected is not None and other is not None and selected!=other
    assert sum(F(r['ideal_payment']) for r in corners)==F(1,4)
    assert sum(F(r['finite_coupon']) for r in corners)==F(3,4)
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],signature=[17,1,1],
        radical=inherited['radical'],white=inherited['white'],affine_nodes=inherited['affine_nodes'],
        positive_affine_kernel=inherited['positive_affine_kernel'],complete_contact_domain=[3,13],
        opposite_sign_contact_shortcut=2,original_ideal_marks=selected,distinct_coupon_marks=other,
        whole_coupon_role=roles,corners=corners,total='1',uses_finite_same_component_payments=False)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    old=docs[M.TAG]
    assert old['new_order17_costs']==dict(A=85,B=111,C=108,D=114)
    assert old['upgraded_coupon']['component_pair_cost_increase']=='1'
    assert old['upgraded_coupon']['small_marks']==[0] and old['upgraded_coupon']['large_marks']==[10,11]
    ft=docs[FINITE];fc=ft['upgraded_coupon']
    assert ft['new_cost']=='235/2' and fc['small_marks']==[3,4] and fc['large_marks']==[9,10,11]
    assert fc['component_pair_cost_increase']=='3'
    new_coupon_marks={(0,9)};preceding_coupon_marks={(0,10),(0,11)}
    finite_coupon_marks={(i,j) for i in (3,4) for j in (9,10,11)}
    assert not new_coupon_marks&preceding_coupon_marks and not (new_coupon_marks|preceding_coupon_marks)&finite_coupon_marks
    weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values())==114
    for p in CHANGE:assert weights[p]==F(1,2);weights[p]=F(1,4)
    assert sum(weights.values())==F(227,2)
    target=old['target'];autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert all({tuple(sorted(mp[i] for i in pair)) for pair in CHANGE}==CHANGE for mp in autos)
    assert all(dist[i][j]==10 for i,j in CHANGE)
    assert all(weights[tuple(p)]==1 for p in old['retained_basic_pairs'])
    assert all(weights[tuple(p)]==1 for p in old['retained_short8'])
    assert all(weights[tuple(p)]==1 for p in old['retained_triangle_pairs'])
    small=docs[O.SMALL];over=docs[O.OVER]
    patterns=[p for p in small['patterns'] if p['order_bound']>=17 and len(p['nodes'])<=17]
    assert len(patterns)==140
    rec=next(r for r in small['records'] if r['order']==17 and r['kind']=='D')
    assert rec['target']==target
    saved={r['case']:r for r in rec['case_embeddings']};support=[];hits=[]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings'] and pairs==list(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        changed=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in CHANGE]
        if changed:
            assert len(changed)==8 and len(p['nodes'])==17
            hits.append(dict(index=p['index'],embeddings=changed,marked_images=pairs))
    assert json.loads(json.dumps(support))==old['reconstructed_support']
    assert len(support)==len(saved)==71 and sum(len(r['embeddings']) for r in support)==488
    assert [r['index'] for r in hits]==EXPECTED and sum(len(r['embeddings']) for r in hits)==16
    oldcases={r['index']:r for r in docs[O.TAG]['cases']};proofs=[]
    for hit in hits:
        c=over['isolation']['cases'][hit['index']];core=over['isolation']['cores'][c['core_key']]
        proof=inspect(c,core,oldcases[hit['index']]);proof['maps']=hit['embeddings'];proofs.append(proof)
    paths={Path(__file__),Path(M.__file__),Path(O.__file__),Path(D.__file__),Path(P.__file__),Path(Q.__file__),Path(L.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in
        (('scripts','py'),('results','json'),('results','md')))
    return dict(status='PASS_ORDER17_D227_OVER2_D5_D12_CENTER_QUARTER_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        target=target,automorphisms=autos,changed_pairs=sorted(CHANGE),old_weight='1/2',new_weight='1/4',
        old_cost='114',new_cost='227/2',saving='1/2',
        new_order17_costs=dict(A='85',B='111',C='108',D='227/2'),new_order17_envelope='227/34',
        weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(weights.items())],
        retained_basic_pairs=old['retained_basic_pairs'],retained_short8=old['retained_short8'],
        retained_triangle_pairs=old['retained_triangle_pairs'],complete_eligible_patterns=140,
        reconstructed_support=support,affected_contexts=EXPECTED,affected_maps=16,selected_proofs=proofs,
        inherited_complete_high_overlap_routing=old['complete_high_overlap_routing'],
        upgraded_coupon=dict(small='whole D5',small_role='unique_long_arm_endpoint',small_marks=[0],
            large='whole D12',large_role='fork_center',large_marks=[9],old_weight='1/2',new_weight='3/4',
            changed_marked_pairs=1,component_pair_cost_increase='1/4',
            preceding_D114_total_component_pair_cost='17/2',new_total_without_finite_transfer='35/4',
            simultaneous_finite_transfer_total_component_pair_cost='47/4'),
        disjoint_coupon_roles=dict(this_packet=sorted(new_coupon_marks),preceding_D114=sorted(preceding_coupon_marks),
            finite_D17=sorted(finite_coupon_marks)),uses_same_component_finite_payments=False,
        all_other_finite_weights_unchanged=True,compatible_finite_D17_cost='235/2',
        triangle_minimum=2,quadrilateral_minimum=1,finite_budget_proved=False,new_dimension_bound=None,
        scope='Complete140 native searches and16 affected maps in two rechecked whole-affine radical seals. Each ideal quarter plus a distinct whole D5-long-end/D12-center three-quarter coupon totals1; added finite cost1/4. Disjoint simultaneous finite/ideal coupon roles checked. No new census, forest or dimension claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists(),'Frozen certificate exists; use --verify.'
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'patterns',out['complete_eligible_patterns'],'maps',out['affected_maps'],
        'ideal',out['new_cost'],'finite increment',out['upgraded_coupon']['component_pair_cost_increase'],flush=True)
