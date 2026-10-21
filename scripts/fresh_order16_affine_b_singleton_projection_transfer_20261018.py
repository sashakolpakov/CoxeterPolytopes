#!/usr/bin/env python3
"""Delete two order16 affine-B singleton units by forced norm-two entries.

Complete selected maps are reconstructed. Written composition retains all
nonselected and current62 interfaces; no finite coupon is introduced.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=X.ROOT;P=X.P
TAG='fresh_order16_affine_b_singleton_projection_transfer_20261018'
CURRENT='fresh_sparse_coupon_order11_13_dimension62_20261018'
TABLES='fresh_deep_finite8_dimension64_20261018'
FINITE='fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018'
D16='fresh_d16_partial_weight_transfer_20261017'
DROP={(1,11),(2,12)}
EXPECTED={(1,11):[1031,1614],(2,12):[504,1031,1537,1614]}
STEMS=(CURRENT,TABLES,FINITE,D16,X.SMALL,X.OVER,X.COMBINED,
       'fresh_order11_bcd_half_pair_transfer_20261018',
       'fresh_cycle11_half_pair_transfer_20261018',
       'fresh_order17_e6_half_coupon_transfer_20261018',
       'fresh_depth11_finite8_zero_transfer_20261018',
       'fresh_order11_c_endpoint_zero_transfer_20261018',
       'fresh_order12_b_center_terminal_transfer_20261018',
       'fresh_order13_d_center_pair_half_transfer_20261018')


def run():
    data={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,doc in data.items():
        assert doc['status'].startswith('PASS'),tag
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    assert data[CURRENT]['dimension_bound']==62
    small=data[X.SMALL];old=next(r for r in small['records'] if (r['order'],r['kind'])==(16,'B'))
    target=old['target'];paid=set(map(tuple,old['retained_paid_pairs']));assert len(paid)==107 and DROP<=paid
    current=next(r for r in data[TABLES]['short7_small_factor_interfaces'] if (r['order'],r['kind'])==(16,'B'))
    assert target==current['target']
    assert {tuple(r['pair']) for r in current['weighted_pairs'] if r['weight']=='1'}==paid
    assert all(r['weight']=='1' for r in current['weighted_pairs'])
    retained=paid-DROP;assert len(retained)==105
    autos=P.embeddings(target,target);adj,ds,_=P.distance_data(target)
    assert all({tuple(sorted(a[i] for i in p)) for p in DROP}==DROP for a in autos)
    assert all(ds[i][j]==10 for i,j in DROP)
    basic=X.L.S.basic_pairs(target);assert len(basic)==100 and basic<=retained
    short8={p for p in combinations(range(16),2) if ds[p[0]][p[1]]<=8};assert short8<=retained
    triangle=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if (r['order'],r['kind'])==(16,'B'))
    assert set(map(tuple,triangle['required_pairs']))<=retained
    patterns=[p for p in small['patterns'] if p['order_bound']>=16 and len(p['nodes'])<=16]
    saved={r['case']:r for r in old['case_embeddings']};support=[];affected={p:[] for p in DROP}
    over=data[X.OVER];proofs=[];hist=Counter()
    for pattern in patterns:
        maps=P.embeddings(pattern['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in pattern['marks'])) for mp in maps})
        if not maps:assert pattern['index'] not in saved;continue
        assert maps==saved[pattern['index']]['embeddings']
        assert pairs==sorted(map(tuple,saved[pattern['index']]['marked_images']))
        support.append(dict(index=pattern['index'],embeddings=maps,marked_images=pairs))
        for pair in sorted(DROP):
            selected=[mp for mp in maps if tuple(sorted(mp[i] for i in pattern['marks']))==pair]
            if not selected:continue
            c=over['isolation']['cases'][pattern['index']];core=over['isolation']['cores'][c['core_key']]
            assert c['method']=='lorentz_residual'
            g,inv,du,dv,ends=X.constraints(c,core);base=c['marked_component'];rows=[]
            for mp in selected:
                image={mp[t]:i for t,i in enumerate(base)};entries=[]
                for new in sorted(set(range(16))-set(image)):
                    for old_node,node in sorted(image.items()):
                        bond=target[new][old_node]
                        if bond:
                            assert bond in (3,4) and node in c['possible_contacts']
                            q2=F(2 if bond==4 else 1)
                            entries.append(dict(new_target_node=new,old_target_node=old_node,old_core_node=node,
                                                native_label=bond,physical_contact_square=str(q2),
                                                inverse_diagonal=str(inv[node][node]),
                                                projection_lower_bound=str(q2*inv[node][node])))
                bad=next((e for e in entries if F(e['projection_lower_bound'])>2),None);assert bad
                if pair==(1,11):
                    assert bad['new_target_node'] in (14,15) and bad['old_target_node']==0
                    assert bad['native_label']==3 and bad['projection_lower_bound']=='10/3'
                else:
                    assert bad['new_target_node']==13 and bad['old_target_node']==12
                    assert bad['native_label']==4 and bad['projection_lower_bound'] in ('4','16')
                hist[bad['projection_lower_bound']]+=1
                rows.append(dict(map=mp,forced_entries=entries,witness=bad,outcome='projection_exceeds_physical_norm_two'))
            affected[pair].append(c['index'])
            proofs.append(dict(index=c['index'],changed_pair=pair,core_key=c['core_key'],core_gram=core['gram'],
                               core_inverse=core['inverse'],base_nodes=base,preferred_corner=c['preferred_corner'],
                               end_projection_proof=ends,maps=rows))
    assert len(patterns)==394 and len(support)==len(saved)==125
    assert sum(len(r['embeddings']) for r in support)==296
    assert affected==EXPECTED and sum(len(r['maps']) for r in proofs)==10
    assert dict(hist)=={'10/3':4,'16':4,'4':2}
    costs={r['kind']:F(r['total_weight']) for r in data[TABLES]['short7_small_factor_interfaces'] if r['order']==16}
    assert costs==dict(A=F(80),B=F(107),C=F(103),D=F(106))
    costs['B']=F(105);assert max(costs.values())==106
    paths=[Path(__file__),Path(X.__file__),Path(X.L.__file__),Path(X.Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER16_B_TWO_SINGLETON_PROJECTION_DELETIONS',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                affine_type='B15',order=16,target=target,automorphisms=autos,deleted_pairs=sorted(DROP),
                old_cost=107,new_cost=105,retained_unit_pairs=sorted(retained),basic_pairs_retained=sorted(basic),
                short8_pairs_retained=sorted(short8),triangle_pairs_retained=triangle['required_pairs'],
                complete_eligible_patterns=394,successful_patterns=125,complete_successful_maps=296,
                reconstructed_support=support,affected_contexts=4,affected_target_contexts=6,affected_maps=10,
                cases=proofs,projection_histogram=dict(hist),new_order16_family_costs={k:str(v) for k,v in costs.items()},
                new_order16_maximum='106',new_rank15_envelope='53/8',
                finite_weights_unchanged=True,added_finite_coupons=[],
                triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
                scope='Complete394-pattern induced-map reconstruction and all10 physical norm-two forced-entry exclusions. Written full-domain/current62 composition retains all short/basic/triangle, finite and later ideal interfaces. No finite charge or ambient bound is added.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_eligible_patterns'],out['affected_maps'],out['projection_histogram'],flush=True)
