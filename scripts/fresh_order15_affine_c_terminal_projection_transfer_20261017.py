#!/usr/bin/env python3
"""Four order15 affine-C payments vanish by forced norm-two terminal entries."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=X.ROOT; P=X.P
TAG='fresh_order15_affine_c_terminal_projection_transfer_20261017'
BASE='fresh_ideal5_weighted_dimension65_20261017'
B15='fresh_order15_affine_b_residual_terminal_transfer_20261017'
DROP={(1,11),(3,13),(1,12),(2,13)}
EXPECTED=[366,611,986,1031,1159,1503,1614,1619,1641,1642]
STEMS=(BASE,B15,X.TAG,X.SMALL,X.OVER,X.COMBINED,
       'fresh_order12_bc_unit_coupon_transfer_20261017',
       'fresh_order13_affine_b_fork_projection_transfer_20261017',
       'fresh_order14_affine_b_fork_projection_transfer_20261017',
       'fresh_large_finite_d_leaf_half_transfer_20261017')


def inspect(c,core,maps,target):
    g,J,du,dv,ends=X.constraints(c,core);base=c['marked_component'];rows=[]
    for mp in maps:
        image={t:base[i] for i,t in enumerate(mp)}
        entries=[]
        for new in sorted(set(range(15))-set(mp)):
            for old,node in sorted(image.items()):
                bond=target[new][old]
                if not bond:continue
                assert bond in (3,4) and node in c['possible_contacts']
                square=F(1 if bond==3 else 2)
                entries.append(dict(new_target_node=new,old_target_node=old,old_core_node=node,
                    native_label=bond,norm_two_contact_square=str(square),
                    inverse_diagonal=str(J[node][node]),projection_lower_bound=str(square*J[node][node])))
        assert entries
        bad=next(z for z in entries if F(z['projection_lower_bound'])>2)
        assert bad['native_label']==4 and bad['norm_two_contact_square']=='2'
        rows.append(dict(map=mp,forced_entries=entries,outcome='projection_exclusion',witness=bad))
    return dict(core_gram=core['gram'],core_inverse=core['inverse'],base_nodes=base,
        preferred_corner=c['preferred_corner'],end_projection_proof=ends,maps=rows)


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    small=data[X.SMALL]
    rec=next(r for r in small['records'] if r['order']==15 and r['kind']=='C')
    target=rec['target'];prior=set(map(tuple,rec['retained_paid_pairs']))
    assert len(prior)==95 and DROP<=prior
    paid=prior-DROP;assert len(paid)==91
    autos=P.embeddings(target,target);_,dist,_=P.distance_data(target)
    assert all({tuple(sorted(a[i] for i in p)) for p in DROP}==DROP for a in autos)
    assert sorted(dist[i][j] for i,j in DROP)==[10,10,11,11]
    basic=X.L.S.basic_pairs(target);assert basic<=paid and len(basic)==90
    assert paid-basic=={(2,12)}
    short8={p for p in combinations(range(15),2) if dist[p[0]][p[1]]<=8};assert short8<=paid
    triangle=next(r for r in data[X.COMBINED]['triangle_pair_retention'] if r['kind']=='C' and r['order']==15)
    assert set(map(tuple,triangle['required_pairs']))<=paid
    patterns=[p for p in small['patterns'] if p['order_bound']>=15 and len(p['nodes'])<=15]
    assert len(patterns)==536
    saved={r['case']:r for r in rec['case_embeddings']};support=[];affected=[]
    for p in patterns:
        maps=P.embeddings(p['labels'],target)
        pairs=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if not maps:
            assert p['index'] not in saved;continue
        assert maps==saved[p['index']]['embeddings']
        assert pairs==sorted(map(tuple,saved[p['index']]['marked_images']))
        support.append(dict(index=p['index'],embeddings=maps,marked_images=pairs))
        selected=[mp for mp in maps if tuple(sorted(mp[i] for i in p['marks'])) in DROP]
        if selected:affected.append((p['index'],selected))
    assert len(support)==len(saved) and [i for i,mp in affected]==EXPECTED
    over=data[X.OVER];rows=[];bounds=Counter();methods=Counter();maps_count=0
    for index,maps in affected:
        c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']]
        result=inspect(c,core,maps,target);methods[c['method']]+=1
        for z in result['maps']:bounds[z['witness']['projection_lower_bound']]+=1;maps_count+=1
        rows.append(dict(index=index,core_key=c['core_key'],method=c['method'],proof=result))
    assert maps_count==24
    costs=dict(data[B15]['new_order15_costs']);assert costs==dict(A=75,B=94,C=95,D=96)
    costs['C']=91;assert max(costs.values())==96
    paths=[Path(__file__),Path(X.__file__),Path(X.L.__file__),Path(X.Q.__file__),Path(P.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER15_C91_TERMINAL_PROJECTION_EXCLUSIONS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        affine_type='C14',diagram_order=15,target=target,automorphisms=autos,
        old_cost=95,new_cost=91,removed_pairs=sorted(DROP),retained_unit_pairs=sorted(paid),
        basic_pairs_retained=sorted(basic),short8_pairs_retained=sorted(short8),
        triangle_pairs_retained=triangle['required_pairs'],complete_eligible_patterns=len(patterns),
        reconstructed_high_support=support,affected_contexts=len(rows),affected_maps=maps_count,
        cases=rows,projection_histogram=dict(bounds),affected_methods=dict(methods),
        simultaneous_order15_costs=costs,order15_envelope_unchanged='32/5',
        finite_weights_unchanged=True,added_finite_coupons=[],
        triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        scope='Complete selected-map reconstruction and forced norm2 native label4 projection exclusions, including arbitrary additional contacts and supporting whites. All basic/short/triangle and inherited low/immediate/later interfaces retained. Original all-label geometric census inherited. D96 still obstructs the order15 maximum94 target; no ambient theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['complete_eligible_patterns'],out['affected_maps'],out['projection_histogram'],flush=True)
