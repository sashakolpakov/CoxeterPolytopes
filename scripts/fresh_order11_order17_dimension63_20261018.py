#!/usr/bin/env python3
"""Compose order11/order17 transfers, finite464 and positive ambient64 cusps.

Written actual-corner case splitting supplies the simultaneous geometric proof.
Exact factor, coupon and scalar interfaces are checked here.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order11_order17_dimension63_20261018'
ORDER11='fresh_order11_bcd_half_pair_transfer_20261018'
CYCLE='fresh_cycle11_half_pair_transfer_20261018'
ORDER17='fresh_order17_e6_half_coupon_transfer_20261018'
BUDGET='fresh_e6_b10_half_coupon_forest_budget_20261018'
IDEAL='fresh_ambient64_ideal_margin_closure_20261018'
BASE='fresh_ambient64_finite_closure_20261018'
ZERO='fresh_depth11_finite8_zero_transfer_20261018'
EXCHANGE='fresh_finite_d16_d17_two_unit_exchange_20261018'
TRIANGLE='fresh_finite7_ideal5_triangle_payments_20261017'
TAIL='fresh_deep_finite8_dimension64_20261018'
STEMS=(ORDER11,CYCLE,ORDER17,BUDGET,IDEAL,BASE,ZERO,EXCHANGE,TRIANGLE,TAIL)


def distances(graph):
    out=[]
    for start in range(len(graph)):
        row=[None]*len(graph);row[start]=0;queue=deque([start])
        while queue:
            i=queue.popleft()
            for j,label in enumerate(graph[i]):
                if label and row[j] is None:row[j]=row[i]+1;queue.append(j)
        assert all(x is not None for x in row);out.append(row)
    return out


def order11_graph(kind):
    if kind=='A':edges=[(i,(i+1)%11,3) for i in range(11)]
    elif kind=='B':edges=[(i,i+1,4 if i==7 else 3) for i in range(8)]+[(0,9,3),(0,10,3)]
    elif kind=='C':edges=[(i,i+1,4 if i in (0,9) else 3) for i in range(10)]
    else:
        assert kind=='D';edges=[(i,i+1,3) for i in range(6)]+[(0,7,3),(0,8,3),(6,9,3),(6,10,3)]
    graph=[[0]*11 for _ in range(11)]
    for i,j,m in edges:graph[i][j]=graph[j][i]=m
    return graph


def run():
    data={stem:json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem,doc in data.items():
        assert doc['status'].startswith('PASS'),stem
        for name,h in doc['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    eleven=data[ORDER11];cycle=data[CYCLE];seventeen=data[ORDER17]
    factor_rows=[]
    for kind in 'ABCD':
        rows=cycle['weighted_pairs'] if kind=='A' else next(r for r in eleven['factor_rules'] if r['kind']==kind)['weighted_pairs']
        weights={tuple(r['pair']):F(r['weight']) for r in rows}
        assert set(weights)==set(combinations(range(11),2))
        graph=order11_graph(kind);dist=distances(graph)
        halved=[p for p,w in weights.items() if w==F(1,2)]
        assert len(halved)==dict(A=11,B=2,C=1,D=4)[kind]
        assert all(dist[i][j]==dict(A=5,B=9,C=10,D=8)[kind] for i,j in halved)
        assert all(w in (F(1,2),F(1)) for w in weights.values())
        short7=[p for p in weights if dist[p[0]][p[1]]<=7]
        missing=[p for p in short7 if weights[p]!=1]
        assert missing==(halved if kind=='A' else [])
        total=sum(weights.values(),F());assert str(total)==eleven['new_order11_costs'][kind]
        factor_rows.append(dict(order=11,kind=kind,total_cost=str(total),halved_pairs=halved,
                                changed_short7_pairs=missing,short7_interface='cycle11 independent two-corner theorem' if kind=='A' else 'all retained units'))
    assert max(F(r['total_cost']) for r in factor_rows)==F(109,2)
    old17=next(r for r in data[TAIL]['short7_small_factor_interfaces'] if (r['order'],r['kind'])==(17,'D'))
    assert seventeen['target']==old17['target']
    weights={tuple(r['pair']):F(r['weight']) for r in seventeen['weighted_pairs']};dist=distances(seventeen['target'])
    assert set(weights)==set(combinations(range(17),2)) and sum(weights.values(),F())==118
    assert all(weights[i,j]==1 for i,j in weights if dist[i][j]<=8)
    factor_rows.append(dict(order=17,kind='D',total_cost='118',halved_new_pairs=seventeen['changed_pairs'],all_short8_units_retained=True))
    assert seventeen['uses_finite_same_component_units'] is False
    # The order17 proper-core branch can yield distance8, not only7.
    # Every inherited small tree of orders12..17 retains these units.
    short8_trees=[]
    for row in data[TAIL]['short7_small_factor_interfaces']:
        if row['kind']=='A':continue
        ds=distances(row['target']);ws={tuple(r['pair']):F(r['weight']) for r in row['weighted_pairs']}
        pairs=[p for p in combinations(range(row['order']),2) if ds[p[0]][p[1]]<=8]
        assert all(ws.get(p,F())==1 for p in pairs)
        short8_trees.append(dict(order=row['order'],kind=row['kind'],short8_unit_pairs=len(pairs)))
    assert len(short8_trees)==18
    # Current local outputs are protected finite units or typed half coupons.
    small=eleven['selected_whole_small_corners'];assert len(small)==12
    assert all(F(r['weight'])>=F(1,2) for r in small)
    assert all(r['attached_spectator_A_rank']<=3 for r in small)
    mixed=eleven['sealed_far_spectator_corners'];assert len(mixed)==12
    for row in mixed:
        if row['outcome']=='excluded by actual edge/corner':continue
        assert row['outcome']=='retained finite unit'
        paid=[c for c in row['corners'] if c['retained_finite_unit']];assert paid
        for corner in paid:
            part=next(p for p in corner['components'] if set(corner['marks'])<=set(p))
            assert len(part)<=18 and corner['marked_distance']<=9
    budget=data[BUDGET]
    assert budget['finite_maximum']==budget['finite_curvature']=='464' and budget['finite_slack']=='0'
    assert budget['attaining_forest']==['D16']*4
    assert [(r['types'],r['new_component_pair_cost']) for r in budget['added_coupons']]==[(['E6','D11'],'2'),(['A1','B10'],'1/2')]
    assert seventeen['finite_coupon']['new_total']=='2'
    assert data[BASE]['finite_budget_hypotheses_discharged'] and data[ZERO]['changed_weight']=='0'
    assert data[EXCHANGE]['new_finite_component_costs']['D16']==116
    assert data[EXCHANGE]['new_finite_component_costs']['D17']==127
    triangle=data[TRIANGLE]
    assert triangle['finite_distance_at_most']==7 and triangle['whole_ideal_factor_order_at_most']==9
    assert triangle['triangle_minimum']==2
    ideal=data[IDEAL]
    assert ideal['complete_terminal_profiles']==2419 and ideal['all_profiles_positive']
    assert ideal['minimum_normalized_margin']['normalized_margin']=='2385/146432'
    assert ideal['minimum_absolute_margin']['absolute_margin']=='240/23'
    assert ideal['global_additive_constant']=='5120'
    assert data[TAIL]['dimension_bound']==64 and data[TAIL]['unrestricted_finite_volume']
    n=64;G=sum(min(i,n-1-i) for i in range(n));finite_curvature=F(5*G-(n-1)*n,2)
    assert finite_curvature==464 and 5*(n*n//4)==5120
    assert F(ideal['minimum_absolute_margin']['absolute_margin'])>0
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
                            for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md',
              ROOT/'results/fresh_cusp_link_shape_curvature_20261017.md',
              ROOT/'results/fresh_ambient65_ideal_margin_closure_20261017.md']
    return dict(status='PASS_ORDER11_ORDER17_UNRESTRICTED_DIMENSION63',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                dimension_bound=63,excluded_dimension=64,inherited_higher_tail_bound=64,
                unrestricted_finite_volume=True,arithmetic_hypothesis=False,
                triangle_minimum=2,quadrilateral_minimum=1,new_factor_interfaces=factor_rows,
                retained_small_tree_short8_interfaces=short8_trees,
                finite_maximum='464',finite_curvature='464',finite_surplus='0',
                all_cusp_absolute_surplus='240/23',strict_global_constant='5120',
                impossible_inequality='0 >= (240/23)*v_i + 5120',
                quadrilateral_case_split=['Any changed cycle11 corner: its independent half theorem supplies total>=1.',
                    'Otherwise any changed order11 tree corner: its independent coupled theorem supplies two halves or a protected finite alternative.',
                    'Otherwise any changed order17D corner: retained short8 replacement or sealed wholeE6/D11 half plus original half.',
                    'No changed new ideal corner: baseline finite464 face coverage, with only nonnegative coupon increments.'],
                unchanged_compact_arithmetic_bound=21,unchanged_arbitrary_compact_bound=29,
                scope='Complete written simultaneous local-to-global composition, exact enlarged finite forest certificate and independently reconstructed2419 ideal profiles. The actual-corner split avoids circular reuse of changed units. Inherited geometric and all-rank classification/concentration inputs retain their stated scopes. No new example, external refereeing or proof-assistant verification.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'bound',out['dimension_bound'],'interfaces',len(out['new_factor_interfaces']),flush=True)
