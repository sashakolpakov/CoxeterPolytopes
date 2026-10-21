#!/usr/bin/env python3
"""Four affine D16 halves vanish using a precisely marked whole D5/D12 unit.

Reconstruct the complete140-pattern selected map domain and four whole-core
radical/finite-role certificates. The new +1 component-pair finite cost is
explicit; this file does not prove a finite forest or new ambient bound.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_d17_partial_weight_transfer_20261017 as O
import fresh_d16_partial_weight_transfer_20261017 as D
import fresh_rank15_order15_ambient60_frontier_20261018 as C

ROOT=O.ROOT;P=O.P;Q=O.Q;L=O.L
TAG='fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018'
CURRENT='fresh_rank15_order15_ambient60_frontier_20261018'
OLD='fresh_order17_e6_unit_coupon_zero_transfer_20261018'
DEPTH='fresh_depth10_finite7_zero_transfer_20261018'
CHANGE={(2,15),(2,16),(10,13),(10,14)}
EXPECTED=[1978,2371,2535,2539]
STEMS=(CURRENT,OLD,O.TAG,D.TAG,O.OVER,O.SMALL,DEPTH,
       'fresh_order15_bd_terminal_coupon_transfer_20261018',
       'fresh_rank15_sparse_payment_exchange_20261018',
       'fresh_order11_order16_dimension60_20261018',
       'fresh_finite8_a17_d11_coupon_budget_20261017',
       'fresh_finite7_ideal5_triangle_payments_20261017')


def role(g,component,mark):
    ns=component['nodes'];typ=component['type'];assert typ in ('D5','D12')
    assert len(ns)==int(typ[1:])
    adj={i:[j for j in ns if i!=j and g[i][j]] for i in ns}
    assert all(g[i][i]==2 for i in ns)
    assert all(g[i][j]==-1 for i in ns for j in adj[i])
    assert sum(map(len,adj.values()))==2*(len(ns)-1)
    forks=[i for i in ns if len(adj[i])==3];assert len(forks)==1
    center=forks[0];dist=Q.graph_dist(L.sub(g,ns),[ns.index(center)])
    assert len(dist)==len(ns) and all(v is not None for v in dist)
    leaves=[i for i in ns if len(adj[i])==1]
    lengths=sorted(dist[ns.index(i)] for i in leaves)
    assert lengths==([1,1,2] if typ=='D5' else [1,1,9])
    wanted=2 if typ=='D5' else 1
    eligible=sorted(i for i in leaves if dist[ns.index(i)]==wanted)
    assert len(eligible)==(1 if typ=='D5' else 2) and mark in eligible
    return dict(type=typ,nodes=ns,actual_mark=mark,fork_center=center,
        leaf_arm_lengths=lengths,eligible_marks=eligible,
        role='unique_long_arm_endpoint' if typ=='D5' else 'short_fork_leaf')


def inspect(c,core,old):
    inherited=O.inspect(c,core);assert inherited==old
    g=[[F(x) for x in row] for row in core['gram']];n=len(g);assert n==19
    u=c['left'][1];v=[c['right_copy_map'][i] for i in c['right'][1]]
    actual=Q.inspect_corner(g,[1]*n,u+v,*c['preferred_corner'])
    for key in ('blocks','marked_component','omitted_black','possible_contacts','preferred_corner'):
        assert actual[key]==c[key]
    assert c['possible_contacts']==[3,13] and set(c['marked_component'])==set(range(n))-{3,13}
    du=Q.graph_dist(g,u);dv=Q.graph_dist(g,v)
    shortcut=min(du[3]+dv[13],du[13]+dv[3])+2;assert shortcut==2
    corners=[];selected=None;other=None;roles=None
    for corner in inherited['corners']:
        if corner['kind']=='affineD16_order17':
            selected=corner['marks']
            corners.append(dict(marks=selected,kind='affineD16_order17',ideal_payment='0',finite_coupon='0'))
            continue
        assert corner['signature']==[17,0,0]
        blocks=D.finite_blocks(g,sorted(set(inherited['white'])|set(corner['marks'])))
        assert [(b['nodes'],b['type']) for b in blocks]==[(b['nodes'],b['type']) for b in corner['components']]
        matches=[]
        for r in corner['coupon_roles']:
            assert r['small_type']=='D5' and r['large_type']=='D12' and r['large_role']=='short_fork_leaf'
            sm=next(b for b in blocks if r['small_mark'] in b['nodes'])
            la=next(b for b in blocks if r['large_mark'] in b['nodes'])
            assert sm is not la and all(not g[i][j] for i in sm['nodes'] for j in la['nodes'])
            matches.append(dict(small=role(g,sm,r['small_mark']),large=role(g,la,r['large_mark']),whole_components=True))
        assert len(matches)<=1
        if matches:
            assert other is None and F(corner['finite_half'])==F(1,2)
            other=corner['marks'];roles=matches[0]
        corners.append(dict(marks=corner['marks'],kind='finite',components=blocks,
            ideal_payment='0',finite_coupon=str(int(bool(matches))),
            same_component_payment_used=False,coupon_roles=matches))
    assert selected is not None and other is not None and selected!=other
    assert sum(F(x['finite_coupon']) for x in corners)==1
    return dict(index=c['index'],core_key=c['core_key'],gram=core['gram'],signature=[17,1,1],
        radical=inherited['radical'],white=inherited['white'],affine_nodes=inherited['affine_nodes'],
        positive_affine_kernel=inherited['positive_affine_kernel'],complete_contact_domain=[3,13],
        opposite_sign_contact_shortcut=shortcut,original_ideal_marks=selected,distinct_unit_marks=other,
        whole_coupon_role=roles,corners=corners,total='1',uses_finite_same_component_payments=False)


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[CURRENT]['completed_unrestricted_dimension_bound']==60
    assert docs[CURRENT]['current_ideal_order15_costs']==dict(A=75,B=92,C=91,D=92)
    coupon_base=docs['fresh_finite8_a17_d11_coupon_budget_20261017']
    assert 5 in [r['rank'] for r in coupon_base['small_components']]
    assert next(r for r in coupon_base['old_large_components'] if r['type']=='D12')['eligible_marks']==3
    assert C.actual_cross(('D',5),('D',12))==15  # The frozen counting code uses doubled costs.
    assert C.actual_cross(('D',12),('D',5))==15
    assert docs[DEPTH]['feature_depth']==10 and docs[DEPTH]['ideal_weights_unchanged']
    old=docs[OLD];assert old['new_order17_costs']==dict(A=85,B=111,C=108,D=116)
    weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
    assert len(weights)==136 and sum(weights.values())==116
    for p in CHANGE:assert weights[p]==F(1,2);weights[p]=F()
    assert sum(weights.values())==114
    small=docs[O.SMALL];over=docs[O.OVER]
    rec=next(r for r in small['records'] if r['order']==17 and r['kind']=='D')
    target=rec['target'];assert target==old['target']==docs[O.TAG]['target']
    autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert all({tuple(sorted(mp[i] for i in pair)) for pair in CHANGE}==CHANGE for mp in autos)
    assert all(dist[i][j]==11 for i,j in CHANGE)
    basic=L.S.basic_pairs(target);assert len(basic)==110 and all(weights[p]==1 for p in basic)
    triangles=docs[O.TAG]['triangle_pairs_retained_at_weight_one']
    assert all(weights[tuple(p)]==1 for p in triangles)
    short8=[p for p in weights if dist[p[0]][p[1]]<=8];assert all(weights[p]==1 for p in short8)
    patterns=[p for p in small['patterns'] if p['order_bound']>=17 and len(p['nodes'])<=17]
    assert len(patterns)==140
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
    assert len(support)==len(saved)==71 and sum(len(r['embeddings']) for r in support)==488
    assert [r['index'] for r in hits]==EXPECTED and sum(len(r['embeddings']) for r in hits)==32
    oldcases={r['index']:r for r in docs[O.TAG]['cases']};proofs=[]
    for hit in hits:
        c=over['isolation']['cases'][hit['index']];core=over['isolation']['cores'][c['core_key']]
        proof=inspect(c,core,oldcases[hit['index']]);proof['maps']=hit['embeddings'];proofs.append(proof)
    cs=over['isolation']['cases'];assert len(cs)==2614
    low=[c for c in cs if 0<(c.get('component_order_bound') or 0)<=16]
    protected=[c for c in cs if c.get('component_order_bound') is None]
    large=[c for c in cs if (c.get('component_order_bound') or 0)>=17 and len(c['marked_component'])>17]
    eligible=[c for c in cs if (c.get('component_order_bound') or 0)>=17 and len(c['marked_component'])<=17]
    assert [len(x) for x in (low,protected,large,eligible)]==[2391,83,0,140]
    assert Counter(c['method'] for c in protected)==dict(finite_leaf9_or_ideal_feature5=74,nine_context_adaptive_finite_payment=9)
    paths={Path(__file__),Path(O.__file__),Path(D.__file__),Path(P.__file__),Path(Q.__file__),Path(L.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in
        (('scripts','py'),('results','json'),('results','md')))
    return dict(status='PASS_ORDER17_D114_WHOLE_D5_D12_UNIT_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        target=target,automorphisms=autos,changed_pairs=sorted(CHANGE),old_weight='1/2',new_weight='0',
        old_cost='116',new_cost='114',new_order17_costs=dict(A=85,B=111,C=108,D=114),new_order17_envelope='114/17',
        weighted_pairs=[dict(pair=p,weight=str(w)) for p,w in sorted(weights.items())],
        retained_basic_pairs=sorted(basic),retained_short8=short8,retained_triangle_pairs=triangles,
        complete_eligible_patterns=140,reconstructed_support=support,affected_contexts=EXPECTED,affected_maps=32,
        complete_high_overlap_routing=dict(total=2614,bound_at_most16=2391,protected_or_adaptive=83,
            displayed_larger_than17=0,eligible_order17=140),selected_proofs=proofs,
        upgraded_coupon=dict(small='whole D5',small_role='unique_long_arm_endpoint',small_marks=[0],
            large='whole D12',large_role='short_fork_leaf',large_marks=[10,11],
            old_weight='1/2',new_weight='1',changed_marked_pairs=2,component_pair_cost_increase='1',
            old_total_component_pair_cost='15/2',new_total_component_pair_cost='17/2'),
        uses_same_component_finite_payments=False,all_other_finite_weights_unchanged=True,
        unchanged_order15_costs=docs[CURRENT]['current_ideal_order15_costs'],compatible_finite7_zero_feature_depth=10,
        triangle_minimum=2,quadrilateral_minimum=1,finite_budget_proved=False,new_dimension_bound=None,
        scope='Complete140-pattern native reconstruction and four rechecked arbitrary-support radical seals. Each changed ideal half is replaced by a distinct whole D5-long-end/D12-shortleaf unit. Exactly+1 finite component-pair cost requires separate counting; all basic/short8/triangle and inherited nonselected/later priorities are retained. No census or ambient-bound claim.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists(),'Frozen certificate exists; use --verify.'
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'patterns',out['complete_eligible_patterns'],'maps',out['affected_maps'],
          'ideal',out['new_cost'],'finite increment',out['upgraded_coupon']['component_pair_cost_increase'],flush=True)
