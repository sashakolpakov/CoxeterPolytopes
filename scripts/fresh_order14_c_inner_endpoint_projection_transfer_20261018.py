#!/usr/bin/env python3
"""Delete the C13 pair at both inner neighbors of the doubled terminals.

All eligible selected patterns, including failed maps, are reconstructed.
Norm-two first-entry projection excludes every changed occurrence.
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
TAG='fresh_order14_c_inner_endpoint_projection_transfer_20261018'
CURRENT='fresh_sparse_coupon_order11_13_dimension62_20261018'
TABLES='fresh_deep_finite8_dimension64_20261018'
DROP=(1,12)
STEMS=(CURRENT,TABLES,X.SMALL,X.OVER,X.COMBINED,
       'fresh_order14_affine_b_fork_projection_transfer_20261017',
       'fresh_order16_affine_b_singleton_projection_transfer_20261018',
       'fresh_order13_d_center_pair_zero_transfer_20261018',
       'fresh_order11_order17_dimension63_20261018',
       'fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018')


def run():
    data={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,doc in data.items():
        assert doc['status'].startswith('PASS'),tag
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(tag,p)
    assert data[CURRENT]['dimension_bound']==62
    small=data[X.SMALL]
    rec=next(r for r in small['records'] if (r['order'],r['kind'])==(14,'C'))
    current=next(r for r in data[TABLES]['short7_small_factor_interfaces'] if (r['order'],r['kind'])==(14,'C'))
    target=rec['target'];assert target==current['target']
    paid=set(map(tuple,rec['retained_paid_pairs']));assert len(paid)==84 and DROP in paid
    assert {tuple(r['pair']) for r in current['weighted_pairs'] if r['weight']=='1'}==paid
    assert all(r['weight']=='1' for r in current['weighted_pairs'])
    retained=paid-{DROP};autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert len(retained)==83 and dist[1][12]==11
    assert all(tuple(sorted(a[i] for i in DROP))==DROP for a in autos)
    basic=X.L.S.basic_pairs(target);assert len(basic)==81 and basic<=retained
    short8={p for p in combinations(range(14),2) if dist[p[0]][p[1]]<=8};assert short8<=retained
    triangle=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if (r['order'],r['kind'])==(14,'C'))
    assert set(map(tuple,triangle['required_pairs']))<=retained
    patterns=[p for p in small['patterns'] if p['order_bound']>=14 and len(p['nodes'])<=14]
    saved={r['case']:r for r in rec['case_embeddings']};support=[];proofs=[];hist=Counter()
    over=data[X.OVER]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings']
        assert pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        selected=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks']))==DROP]
        if not selected:continue
        c=over['isolation']['cases'][p['index']];core=over['isolation']['cores'][c['core_key']]
        assert c['method']=='lorentz_residual'
        g,inv,du,dv,ends=X.constraints(c,core);base=c['marked_component'];rows=[]
        for mp in selected:
            image={mp[t]:i for t,i in enumerate(base)}
            assert set(image)==set(range(1,13))
            entries=[]
            for new,old in [(0,1),(13,12)]:
                node=image[old];assert target[new][old]==4 and node in c['possible_contacts']
                assert node in c['preferred_corner']
                lower=2*inv[node][node];assert lower>2
                entries.append(dict(new_target_node=new,old_target_node=old,core_node=node,
                                    native_label=4,physical_contact_square='2',
                                    inverse_diagonal=str(inv[node][node]),projection_lower_bound=str(lower)))
                hist[str(lower)]+=1
            rows.append(dict(map=mp,forced_terminal_entries=entries,
                             outcome='either_missing_terminal_projection_exceeds_physical_norm_two'))
        proofs.append(dict(index=c['index'],core_key=c['core_key'],core_gram=core['gram'],
                           core_inverse=core['inverse'],base_nodes=base,preferred_corner=c['preferred_corner'],
                           end_projection_proof=ends,maps=rows))
    assert len(patterns)==600 and len(support)==len(saved)==43
    assert sum(len(r['embeddings']) for r in support)==92
    assert [r['index'] for r in proofs]==[624,1555,1641]
    assert sum(len(r['maps']) for r in proofs)==6 and hist==Counter({'4':4,'7':4,'16':4})
    costs={r['kind']:F(r['total_weight']) for r in data[TABLES]['short7_small_factor_interfaces'] if r['order']==14}
    assert costs==dict(A=F(70),B=F(84),C=F(84),D=F(82))
    costs['C']=F(83);assert max(costs.values())==84
    paths=[Path(__file__),Path(X.__file__),Path(X.L.__file__),Path(X.Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER14_C_INNER_ENDPOINT_PAIR_PROJECTION_DELETION',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        affine_type='C13',order=14,target=target,automorphisms=autos,deleted_pair=DROP,
        old_cost=84,new_cost=83,retained_unit_pairs=sorted(retained),basic_pairs_retained=sorted(basic),
        short8_pairs_retained=sorted(short8),triangle_pairs_retained=triangle['required_pairs'],
        complete_eligible_patterns=600,successful_patterns=43,complete_successful_maps=92,
        reconstructed_support=support,affected_contexts=3,affected_maps=6,
        cases=proofs,forced_terminal_projection_histogram=dict(hist),
        new_order14_family_costs={k:str(v) for k,v in costs.items()},new_order14_maximum='84',
        finite_weights_unchanged=True,added_finite_coupons=[],triangle_minimum=2,quadrilateral_minimum=1,
        compatible_with_D13_zero=True,compatible_with_B16_singleton_deletions=True,
        new_dimension_bound=None,
        scope='Complete600-pattern induced-map reconstruction, including557 failed patterns, and all6 affected maps. Both missing physical norm-two doubled terminals have forced projections>2. Arbitrary same-end contacts only increase the bound; opposite-end contacts give the inherited protected short8 route. Full nonselected/current-rule interfaces are inherited explicitly in the report; no finite or global gain is asserted.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_eligible_patterns'],out['affected_maps'],out['forced_terminal_projection_histogram'],flush=True)
