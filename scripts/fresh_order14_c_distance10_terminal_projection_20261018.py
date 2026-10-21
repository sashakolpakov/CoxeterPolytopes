#!/usr/bin/env python3
"""C13/order14 cost83 ->81 with two coupon-free terminal exclusions.

Reconstruct every600 eligible native search, including failures, and all
28 affected maps. Complete positive-end contact bounds use physical norm2.
No finite counting or ambient-dimension theorem is asserted.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_order14_c_inner_endpoint_projection_transfer_20261018 as C

X=C.X;ROOT=C.ROOT;P=C.P
TAG='fresh_order14_c_distance10_terminal_projection_20261018'
CURRENT='fresh_d17_fractional_ambient60_frontier_20261018'
FAMILIES='fresh_order14_ambient62_ideal_margin_closure_20261018'
DROP={(1,11),(2,12)}
EXPECTED=[366,380,392,531,614,1115,1159,1503,1505,1508,1541,1619,1642]
STEMS=(CURRENT,FAMILIES,C.TAG,X.SMALL,X.OVER,X.COMBINED,
       'fresh_depth10_finite7_zero_transfer_20261018',
       'fresh_order15_bd_terminal_coupon_transfer_20261018',
       'fresh_order17_d5_d12_center_quarter_transfer_20261018')


def inspect(c,core,maps,target):
    g,inv,du,dv,ends=X.constraints(c,core)
    assert c['method'] in ('lorentz_residual','relevant_component')
    base=c['marked_component'];rows=[]
    for mp in maps:
        image={t:base[i] for i,t in enumerate(mp)}
        assert set(image) in (set(range(1,13)),set(range(1,12)),set(range(2,13)))
        assert len(image)==len(base) in (11,12)
        entries=[]
        for new in sorted(set(range(14))-set(image)):
            for old,node in sorted(image.items()):
                bond=target[new][old]
                if not bond:continue
                assert bond in (3,4) and node in c['possible_contacts']
                assert new not in image
                q2=F(1 if bond==3 else 2);lower=q2*inv[node][node]
                entries.append(dict(new_target_node=new,old_target_node=old,core_node=node,
                    native_label=bond,physical_contact_square=str(q2),
                    inverse_diagonal=str(inv[node][node]),projection_lower_bound=str(lower),
                    old_node_is_marked=node in c['preferred_corner']))
        assert len(entries)==2 and len({e['new_target_node'] for e in entries})==2
        bad=[e for e in entries if e['native_label']==4 and F(e['projection_lower_bound'])>2]
        assert bad and min(F(e['projection_lower_bound']) for e in bad)>=F(5,2)
        assert all(e['new_target_node'] in (0,13) for e in bad)
        if len(base)==12:assert len(bad)==2
        else:assert c['index']==1642 and len(bad)==1
        rows.append(dict(map=mp,displayed_target_nodes=sorted(image),forced_entries=entries,
            obstruction=bad[0],outcome='genuinely_new_norm2_doubled_terminal_has_projection_above2'))
    return dict(index=c['index'],method=c['method'],core_key=c['core_key'],core_gram=core['gram'],core_inverse=core['inverse'],
        base_nodes=base,preferred_corner=c['preferred_corner'],end_projection_proof=ends,maps=rows)


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert data[CURRENT]['completed_unrestricted_dimension_bound']==60
    assert data[CURRENT]['finite']['exact_maximum']=='411' and data[CURRENT]['ideal']['summary']['negative_count']==28
    assert data[FAMILIES]['order14_family_costs']==dict(A=70,B=82,C=83,D=82)
    previous=data[C.TAG];assert previous['new_cost']==83
    prior=set(map(tuple,previous['retained_unit_pairs']));assert len(prior)==83 and DROP<=prior
    retained=prior-DROP;assert len(retained)==81
    small=data[X.SMALL];over=data[X.OVER]
    rec=next(r for r in small['records'] if (r['order'],r['kind'])==(14,'C'))
    target=rec['target'];assert target==previous['target']
    autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert all({tuple(sorted(a[i] for i in pair)) for pair in DROP}==DROP for a in autos)
    assert all(dist[i][j]==10 for i,j in DROP)
    basic=X.L.S.basic_pairs(target);assert len(basic)==81 and basic==retained
    short8={p for p in combinations(range(14),2) if dist[p[0]][p[1]]<=8};assert short8<=retained
    tri=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if(r['order'],r['kind'])==(14,'C'))
    assert set(map(tuple,tri['required_pairs']))<=retained
    patterns=[p for p in small['patterns'] if p['order_bound']>=14 and len(p['nodes'])<=14]
    assert len(patterns)==600
    saved={r['case']:r for r in rec['case_embeddings']};support=[];cases=[];outcomes=Counter();hist=Counter()
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings'] and pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        selected=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in DROP]
        if not selected:continue
        c=over['isolation']['cases'][p['index']];core=over['isolation']['cores'][c['core_key']]
        proof=inspect(c,core,selected,target);cases.append(proof)
        for row in proof['maps']:
            outcomes[len(row['displayed_target_nodes'])]+=1
            for e in row['forced_entries']:
                if e['native_label']==4:hist[e['projection_lower_bound']]+=1
    assert len(support)==len(saved)==43 and sum(len(r['embeddings']) for r in support)==92
    assert json.loads(json.dumps(support))==previous['reconstructed_support']
    assert [r['index'] for r in cases]==EXPECTED
    assert sum(len(r['maps']) for r in cases)==28 and outcomes=={12:24,11:4}
    assert all(F(x)>2 for x in hist)
    cs=over['isolation']['cases'];assert len(cs)==2614
    low=[c for c in cs if 0<(c.get('component_order_bound') or 0)<=13]
    protected=[c for c in cs if c.get('component_order_bound') is None]
    large=[c for c in cs if (c.get('component_order_bound') or 0)>=14 and len(c['marked_component'])>14]
    eligible=[c for c in cs if (c.get('component_order_bound') or 0)>=14 and len(c['marked_component'])<=14]
    assert [len(x) for x in (low,protected,large,eligible)]==[1378,83,553,600]
    assert Counter(c['method'] for c in protected)==dict(finite_leaf9_or_ideal_feature5=74,nine_context_adaptive_finite_payment=9)
    paths={Path(__file__),Path(C.__file__),Path(X.__file__),Path(X.L.__file__),Path(X.Q.__file__),Path(P.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in
        (('scripts','py'),('results','json'),('results','md')))
    return dict(status='PASS_ORDER14_C81_DISTANCE10_TERMINAL_PROJECTION_DELETION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        affine_type='C13',order=14,target=target,automorphisms=autos,deleted_pairs=sorted(DROP),
        old_cost='83',new_cost='81',new_order14_family_costs=dict(A=70,B=82,C=81,D=82),new_order14_maximum='82',
        old_order14_envelope='83/14',new_order14_envelope='41/7',
        retained_unit_pairs=sorted(retained),weighted_pairs=[dict(pair=p,weight=str(int(p in retained))) for p in combinations(range(14),2)],
        basic_pairs_retained=sorted(basic),short8_pairs_retained=sorted(short8),triangle_pairs_retained=tri['required_pairs'],
        complete_eligible_patterns=600,successful_patterns=43,complete_successful_maps=92,reconstructed_support=support,
        affected_contexts=13,affected_maps=28,cases=cases,displayed_order_histogram=dict(outcomes),
        forced_doubled_terminal_bound_histogram=dict(hist),
        complete_high_overlap_routing=dict(total=2614,bound_at_most13=1378,protected_or_adaptive=83,
            displayed_larger_than14=553,eligible_order14=600),
        finite_weights_unchanged=True,added_finite_coupons=[],triangle_minimum=2,quadrilateral_minimum=1,
        compatible_finite7_feature_depth=10,compatible_order17_cost='227/2',
        completed_dimension_bound_unchanged=60,new_dimension_bound=None,
        scope='Complete600 native searches, including557 failures, and all28 affected maps. Genuinely new physical norm2 doubled terminals have forced projection>2 under complete arbitrary-contact end bounds. Short/basic/triangle and inherited nonselected/current-priority interfaces retained. No new finite payment, forest calculation or ambient bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:
        assert not path.exists(),'Frozen certificate exists; use --verify.'
        path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'patterns',out['complete_eligible_patterns'],'maps',out['affected_maps'],
        'C',out['new_cost'],'maximum',out['new_order14_maximum'],out['forced_doubled_terminal_bound_histogram'],flush=True)
