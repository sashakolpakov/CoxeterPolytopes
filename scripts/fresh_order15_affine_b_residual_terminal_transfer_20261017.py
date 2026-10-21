#!/usr/bin/env python3
"""Order15 B96 to94 by forced entries and a residual-zero terminal exclusion."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=X.ROOT;P=X.P
TAG='fresh_order15_affine_b_residual_terminal_transfer_20261017'
BASE='fresh_ideal5_weighted_dimension65_20261017'
DROP={(10,13),(10,14)}
EXPECTED=[629,637,674,778,1163,1176,1559,1572,1641,1647,1667,1848,
          1911,1934,2241,2304,2327]
SATURATED={1647,1667,1934,2327}
STEMS=(BASE,X.TAG,X.SMALL,X.OVER,X.COMBINED,
       'fresh_order15_affine_b_pair_reduction_20261017',
       'fresh_order12_bc_unit_coupon_transfer_20261017',
       'fresh_order13_affine_b_fork_projection_transfer_20261017',
       'fresh_order14_affine_b_fork_projection_transfer_20261017',
       'fresh_large_finite_d_leaf_half_transfer_20261017')


def inspect(c,core,maps,target):
    g,inv,du,dv,ends=X.constraints(c,core);base=c['marked_component'];rows=[]
    assert c['method']=='lorentz_residual'
    for mp in maps:
        image={t:base[i] for i,t in enumerate(mp)}
        entries=[]
        for new in sorted(set(range(15))-set(mp)):
            for old,node in sorted(image.items()):
                bond=target[new][old]
                if bond:
                    assert bond in (3,4) and node in c['possible_contacts']
                    square=F(1 if bond==3 else 2);bound=square*inv[node][node]
                    entries.append(dict(new_target_node=new,old_target_node=old,old_core_node=node,
                        native_label=bond,norm_two_contact_square=str(square),
                        inverse_diagonal=str(inv[node][node]),projection_lower_bound=str(bound)))
        assert entries
        bad=next((z for z in entries if F(z['projection_lower_bound'])>2),None)
        if bad:
            assert c['index'] not in SATURATED
            rows.append(dict(map=mp,forced_entries=entries,outcome='projection_exclusion',witness=bad))
            continue
        assert c['index'] in SATURATED and len(entries)==1
        e=entries[0];i=e['old_core_node']
        assert (e['new_target_node'],e['old_target_node'],e['native_label'])==(11,10,3)
        assert inv[i][i]==2 and i in c['preferred_corner']
        assert 11 not in image and 12 not in image and target[11][12]==4
        same=next(end for end in c['end_sets'] if i in end)
        opposite=next(end for end in c['end_sets'] if i not in end)
        assert len(opposite)==1 and all(inv[i][j]>0 for j in same)
        j=opposite[0];assert inv[i][j]==-2
        # Saturation forces only the unit i-contact and zero residual. Every
        # other new white is orthogonal to H, on the same positive end, or on
        # the singleton opposite end. Its actual product with the forced
        # root is respectively0, positive, or<=-2. No required label4 bond
        # is possible; the last alternative violates strict white PD.
        rows.append(dict(map=mp,forced_entries=entries,outcome='residual_zero_no_terminal',
            residual_square='0',same_end=same,
            same_end_projection_products={str(k):str(inv[i][k]) for k in same},
            opposite_end=opposite,opposite_unit_projection_product='-2',
            required_next_target_node=12,required_next_native_label=4,
            core_orthogonal_product='0'))
    return dict(core_gram=core['gram'],core_inverse=core['inverse'],base_nodes=base,
        preferred_corner=c['preferred_corner'],end_projection_proof=ends,maps=rows)


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    previous=next(r for r in data[X.TAG]['records'] if r['order']==15)
    small=data[X.SMALL];rec=next(r for r in small['records'] if r['order']==15 and r['kind']=='B')
    target=rec['target'];assert target==previous['target']
    prior=set(map(tuple,previous['new_paid_pairs']));assert len(prior)==96 and DROP<=prior
    paid=prior-DROP;assert len(paid)==94
    autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert all({tuple(sorted(a[i] for i in p)) for p in DROP}==DROP for a in autos)
    assert all(dist[i][j]==11 for i,j in DROP)
    basic=X.L.S.basic_pairs(target);assert basic<=paid and len(basic)==91
    short8={p for p in combinations(range(15),2) if dist[p[0]][p[1]]<=8};assert short8<=paid
    triangle=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if r['kind']=='B' and r['order']==15)
    assert set(map(tuple,triangle['required_pairs']))<=paid
    patterns=[p for p in small['patterns'] if p['order_bound']>=15 and len(p['nodes'])<=15]
    saved={r['case']:r for r in rec['case_embeddings']};all_support=[];affected=[]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:
            assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings']
        assert pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
        all_support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        selected=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in DROP]
        if selected:affected.append((p['index'],selected))
    assert len(all_support)==len(saved) and [i for i,mp in affected]==EXPECTED
    over=data[X.OVER];rows=[];outcomes=Counter();bounds=Counter()
    for index,maps in affected:
        c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']]
        result=inspect(c,core,maps,target)
        for r in result['maps']:
            outcomes[r['outcome']]+=1
            if r['outcome']=='projection_exclusion':bounds[r['witness']['projection_lower_bound']]+=1
        rows.append(dict(index=index,core_key=c['core_key'],proof=result))
    assert outcomes==Counter(projection_exclusion=28,residual_zero_no_terminal=8)
    costs=dict(data[BASE]['small_order_costs']['15']);assert costs==dict(A=75,B=96,C=95,D=96)
    costs['B']=94;assert max(costs.values())==96
    paths=[Path(__file__),Path(X.__file__),Path(X.L.__file__),Path(X.Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER15_B94_RESIDUAL_TERMINAL_EXCLUSIONS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        affine_type='B14',diagram_order=15,target=target,automorphisms=autos,
        old_cost=96,new_cost=94,removed_pairs=sorted(DROP),retained_unit_pairs=sorted(paid),
        basic_pairs_retained=sorted(basic),short8_pairs_retained=sorted(short8),
        triangle_pairs_retained=triangle['required_pairs'],complete_eligible_patterns=len(patterns),
        reconstructed_high_support=all_support,affected_contexts=17,affected_maps=36,
        cases=rows,outcomes=dict(outcomes),direct_projection_histogram=dict(bounds),
        new_order15_costs=costs,order15_envelope_unchanged='32/5',
        finite_weights_unchanged=True,added_finite_coupons=[],
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        scope='Complete eligible selected maps; forced native norm2 entry exclusion or saturated first root with no possible terminal among any supporting whites. All short/basic/triangle and inherited low/immediate/later interfaces retained. Original geometric census is inherited. C95/D96 remain order15 bottlenecks; no new ideal envelope or ambient theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_eligible_patterns'],out['outcomes'],out['direct_projection_histogram'],flush=True)
