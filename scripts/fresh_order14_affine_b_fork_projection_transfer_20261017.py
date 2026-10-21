#!/usr/bin/env python3
"""Two further order14 affine-B pairs vanish by forced terminal projections.

All52 affected maps are checked with complete necessary contact domains.
No finite coupon, finite weight or other ideal factor is changed.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=X.ROOT;P=X.P;Q=X.Q;L=X.L
TAG='fresh_order14_affine_b_fork_projection_transfer_20261017'
DROP={(10,12),(10,13)}
EXPECTED=[624,629,637,674,713,719,736,778,842,1555,1559,1567,1572,1584,
          1641,1647,1667,1848,1855,1934,2241,2248,2327]


def run():
    stems=(X.TAG,X.SMALL,X.OVER,X.COMBINED,
        'fresh_ideal5_weighted_dimension66_20261017',
        'fresh_cycle12_ideal5_transfer_20261017',
        'fresh_affine_d11_leaf_pair_transfer_20261017',
        'fresh_order12_affine_b_half_transfer_20261017',
        'fresh_order12_affine_c_endpoint_reduction_20261017')
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in stems}
    for stem,d in data.items():
        assert d['status'].startswith('PASS'),stem
        for path,h in d['sources'].items():
            assert sha256((ROOT/path).read_bytes()).hexdigest()==h,(stem,path)
    previous=next(r for r in data[X.TAG]['records'] if r['order']==14)
    rec=next(r for r in data[X.SMALL]['records'] if r['order']==14 and r['kind']=='B')
    target=rec['target'];assert target==previous['target']
    paid=set(map(tuple,previous['new_paid_pairs']));assert len(paid)==86 and DROP<=paid
    reduced=paid-DROP;assert len(reduced)==84
    autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert all({tuple(sorted(a[i] for i in p)) for p in DROP}==DROP for a in autos)
    assert all(dist[i][j]==11 for i,j in DROP)
    short8={p for p in combinations(range(14),2) if dist[p[0]][p[1]]<=8}
    basic=L.S.basic_pairs(target)
    assert short8<=reduced and basic<=reduced and len(basic)==82
    triangle=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if r['order']==14 and r['kind']=='B')
    assert set(map(tuple,triangle['required_pairs']))<=reduced
    affected=[r for r in rec['case_embeddings'] if DROP&set(map(tuple,r['marked_images']))]
    assert [r['case'] for r in affected]==EXPECTED
    over=data[X.OVER];rows=[];histogram=Counter();map_count=0;cores=set()
    patterns={r['index']:r for r in data[X.SMALL]['patterns']}
    for saved in affected:
        c=over['isolation']['cases'][saved['case']]
        assert c['method']=='lorentz_residual'
        core=over['isolation']['cores'][c['core_key']]
        g,inv,du,dv,endproof=X.constraints(c,core)
        cores.add(c['core_key']);ids=c['marked_component']
        native=[[int(2-g[i][j]) if i!=j and g[i][j] else 0 for j in ids] for i in ids]
        assert native==patterns[c['index']]['labels']
        all_maps=P.embeddings(native,target)
        assert all_maps==saved['embeddings']
        maps=[]
        for mp in all_maps:
            image=tuple(sorted(mp[ids.index(i)] for i in c['preferred_corner']))
            if image not in DROP:
                continue
            row=X.inspect_map(c,g,inv,du,dv,target,mp)
            assert row['witness']['new_target_node']==11 and row['witness']['old_target_node']==10
            assert row['witness']['coxeter_label']==4
            maps.append(row);map_count+=1
            histogram[row['witness']['forced_projection_lower_bound']]+=1
        assert maps
        rows.append(dict(index=c['index'],core_key=c['core_key'],method=c['method'],
            left=c['left'],right=c['right'],match=c['match'],core_gram=core['gram'],
            core_inverse=core['inverse'],base_nodes=ids,preferred_corner=c['preferred_corner'],
            contact_domain=c['possible_contacts'],end_sets=c['end_sets'],
            end_projection_proof=endproof,complete_pattern_maps=all_maps,affected_maps=maps))
    assert map_count==52 and histogram==Counter({'4':12,'7':12,'16':20,'52':8})
    costs=dict(data['fresh_ideal5_weighted_dimension66_20261017']['small_order_costs']['14'])
    assert costs==dict(A=70,B=86,C=84,D=82)
    costs['B']=84;assert max(costs.values())==84
    paths=[Path(__file__),Path(X.__file__),Path(L.__file__),Path(Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in stems for folder,ext in (
        ('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER14_AFFINE_B_TWO_FORK_PAIR_PROJECTION_EXCLUSIONS',
        affine_type='B13',diagram_order=14,target=target,automorphisms=autos,
        inherited_cost=86,new_cost=84,removed_pairs=sorted(DROP),retained_pairs=sorted(reduced),
        affected_cases=EXPECTED,affected_contexts=len(rows),affected_map_count=map_count,
        distinct_core_keys=len(cores),projection_histogram=dict(histogram),cases=rows,
        triangle_pairs_retained=triangle['required_pairs'],short8_pairs_retained=sorted(short8),
        basic_pairs_retained=sorted(basic),finite_weights_unchanged=True,added_finite_coupons=[],
        new_order14_costs=costs,new_rank13_normalized_envelope='6',new_ambient_bound=None,
        domain_scope='Complete saved high support; low whole-component bound12; disjoint exclusion; protected/basic ideal paths<=9 retained; universal immediate affine outputs<=13; post69 cycle/D11/B12/C12 outputs retained short paths, finite payments, or smaller affine factors. No finite or other ideal payment changes.',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        scope='All52 affected induced maps in23 complete selected contexts force a genuinely new norm2 doubled terminal with projection square>=4. Arbitrary extra contacts in the same end increase the bound; opposite-end contacts give a retained short path. No unspecified contacts are set zero, no new finite coupon and no ambient theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert out==json.loads(p.read_text())
    else:
        p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['affected_contexts'],out['affected_map_count'],out['projection_histogram'],flush=True)
