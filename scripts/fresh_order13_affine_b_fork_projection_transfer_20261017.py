#!/usr/bin/env python3
"""Delete the last two nonbasic pairs of order13 affine B without new coupons.

The four radical completions use the already proved whole A2/D11 unit.
All proper completions fail a forced norm2 doubled-terminal projection.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_order13_affine_b_singleton_transfer_20261017 as L
import fresh_b14_b15_singleton_completion_transfer_20261017 as X

ROOT=L.ROOT;P=L.P
TAG='fresh_order13_affine_b_fork_projection_transfer_20261017'
DROP={(9,11),(9,12)}
EXPECTED=[149,151,625,630,638,715,721,1356,1364,1556,1642,1648,1668,
          1742,1849,1857,1922,1935,2135,2242,2250,2315,2328]
BASE='fresh_ideal5_weighted_dimension65_20261017'
BC='fresh_order12_bc_unit_coupon_transfer_20261017'
STEMS=(BASE,BC,L.TAG,L.DIAG,L.OVER,L.D14,
       'fresh_affine_d12_leaf_pair_transfer_20261017',
       'fresh_finite8_ideal5_combined_payments_20261017',
       'fresh_order14_affine_b_fork_projection_transfer_20261017',
       'fresh_large_finite_d_leaf_half_transfer_20261017')


def radical_payment(c,core):
    row=L.inspect_case(c,core);payments=[]
    for corner in row['corners']:
        if corner['kind']!='finite':
            assert corner['kind']=='affineB12_order13'
            payments.append(dict(marks=corner['marks'],kind=corner['kind'],weight='0'))
            continue
        cc=corner['components'];marks=corner['marks']
        a,b=[next(z for z in cc if m in z['nodes']) for m in marks]
        unit=int(a is b);assert unit==corner['finite_unit']
        if unit:assert a['type']=='B13' and len(a['nodes'])==13
        coupon=0;roles=None
        for small,large,lm in ((a,b,marks[1]),(b,a,marks[0])):
            if small['type']=='A2' and large['type']=='D11' and large['roles'].get(str(lm))=='short_fork_leaf':
                assert a is not b
                coupon=1;roles=dict(small_nodes=small['nodes'],large_nodes=large['nodes'],
                                    large_mark=lm,large_role='short_fork_leaf')
        payments.append(dict(marks=marks,kind='finite',components=cc,
            finite_unit=unit,existing_A2_D11_unit=coupon,coupon_roles=roles,weight=str(unit+coupon)))
    total=sum(F(p['weight']) for p in payments)
    assert total==(2 if c['index'] in (149,1742) else 1)
    assert sum(p.get('existing_A2_D11_unit',0) for p in payments)==1
    return dict(geometry=row,payments=payments,finite_total=str(total),ideal_weight='0')


def proper_maps(c,core,maps,target):
    g,inv,du,dv,ends=X.constraints(c,core)
    base=c['marked_component'];corner=sorted(set(range(len(g)))-set(c['omitted_black']))
    other=sorted(set(corner)-set(base))
    assert all(g[i][j]==0 for i in base for j in other)
    records=[]
    for mp in maps:
        assert 9 in mp and 10 not in mp
        i=base[mp.index(9)]
        assert i in c['preferred_corner'] and i in c['possible_contacts']
        assert target[9][10]==4
        value=2*inv[i][i];assert value>2
        records.append(dict(map=mp,new_target_node=10,old_target_node=9,
            old_core_node=i,native_label=4,norm_two_contact_square='2',
            inverse_diagonal=str(inv[i][i]),forced_projection_lower_bound=str(value)))
    return dict(core_gram=core['gram'],core_inverse=core['inverse'],
        base_nodes=base,preferred_corner=c['preferred_corner'],
        end_projection_proof=ends,maps=records)


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,d in data.items():
        assert d['status'].startswith('PASS'),s
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert data[BASE]['dimension_bound']==65
    assert data[BC]['coupon_upgrade'].startswith('Whole A2')
    assert data[BC]['new_coupon_increment_per_A2_D11_component_pair']=='2'
    diag=data[L.DIAG];over=data[L.OVER]
    rec=next(r for r in diag['records'] if r['kind']=='B');target=rec['target']
    assert diag['eligible_patterns']==len(diag['patterns'])==801
    prior=set(map(tuple,rec['required_union']))-{(0,10),(10,11),(10,12)}
    assert len(prior)==75 and DROP<=prior
    paid=prior-DROP;assert paid==set(map(tuple,rec['basic_pairs'])) and len(paid)==73
    autos=P.embeddings(target,target);_,dist,features=P.distance_data(target)
    assert all({tuple(sorted(a[i] for i in p)) for p in DROP}==DROP for a in autos)
    assert all(dist[i][j]==10 for i,j in DROP)
    short8={p for p in combinations(range(13),2) if dist[p[0]][p[1]]<=8}
    assert short8<=paid
    triangle={p for p in combinations(range(13),2) if dist[p[0]][p[1]]<=5 or
              dist[p[0]][p[1]]<=8 and any(dist[i][f]<=5 for i in p for f in features)}
    assert triangle==set(map(tuple,data[L.TAG]['triangle_pairs_retained'])) and triangle<=paid
    saved={r['case']:r for r in rec['case_embeddings']};all_support=[];affected=[]
    for p in diag['patterns']:
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
    rows=[];hist=Counter();methods=Counter();mapcount=0
    for index,maps in affected:
        c=over['isolation']['cases'][index];core=over['isolation']['cores'][c['core_key']]
        assert c['index']==index
        if c['method']=='radical':
            assert index in (149,151,1742,2135) and all(len(mp)==13 for mp in maps)
            result=radical_payment(c,core);outcome='sealed_A2_D11_unit'
        else:
            result=proper_maps(c,core,maps,target);outcome='forced_doubled_terminal_projection'
            hist.update(r['forced_projection_lower_bound'] for r in result['maps'])
        rows.append(dict(index=index,core_key=c['core_key'],method=c['method'],
            embeddings=maps,outcome=outcome,proof=result))
        methods[outcome]+=1;mapcount+=len(maps)
    assert mapcount==52 and methods==Counter(sealed_A2_D11_unit=4,forced_doubled_terminal_projection=19)
    assert sum(hist.values())==44
    assert not(DROP&set(map(tuple,rec['immediate_pairs'])))
    costs=dict(data[BASE]['small_order_costs']['13'])
    assert costs==dict(A=65,B=75,C=73,D=74)
    costs['B']=73;assert max(costs.values())==74
    paths=[Path(__file__),Path(L.__file__),Path(X.__file__),Path(P.__file__),Path(L.Q.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS
              for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_ORDER13_B_TWO_FORK_PAIR_REMOVALS_WITH_EXISTING_A2_D11',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        affine_type='B12',diagram_order=13,target=target,automorphisms=autos,
        old_cost=75,new_cost=73,removed_pairs=sorted(DROP),retained_unit_pairs=sorted(paid),
        short8_pairs_retained=sorted(short8),triangle_pairs_retained=sorted(triangle),
        complete_high_patterns=801,reconstructed_high_support=all_support,
        affected_contexts=23,affected_maps=52,cases=rows,outcomes=dict(methods),
        projection_histogram=dict(hist),immediate_D14_affected_pairs=[],
        new_order13_costs=costs,new_rank12_envelope='74/13',
        finite_weights_unchanged=True,added_finite_coupons=[],
        existing_coupon_dependency=BC,triangle_minimum=2,quadrilateral_minimum=1,new_dimension_bound=None,
        scope='All801 eligible high patterns re-enumerated;19 proper-completion contexts excluded by complete positive-end norm2 projection, four full-affine radical contexts sealed and paid by existing whole A2/D11 units. Low whole-order bound12, retained short/basic units and all universal/later outputs covered in the report. Original all-label geometric census is inherited. No new finite cost or ambient theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['outcomes'],out['projection_histogram'],out['new_order13_costs'],flush=True)
