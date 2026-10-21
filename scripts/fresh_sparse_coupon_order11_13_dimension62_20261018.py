#!/usr/bin/env python3
"""Compose the sparse finite rule and three ideal refinements at ambient63.

The accompanying proof supplies the exhaustive actual-corner priority split.
This binder checks exact frozen interfaces; it does not rerun old censuses,
the 2204-profile calculation or the finite-forest dynamic programs.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sparse_coupon_order11_13_dimension62_20261018'
CURRENT='fresh_order11_order17_dimension63_20261018'
C11='fresh_order11_c_endpoint_zero_transfer_20261018'
B12='fresh_order12_b_center_terminal_transfer_20261018'
D13='fresh_order13_d_center_pair_half_transfer_20261018'
FINITE='fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018'
BUDGET='fresh_sparse_coupon_rank63_forest_budget_20261018'
IDEAL='fresh_ambient63_ideal_margin_closure_20261018'
OLD_IDEAL='fresh_ambient64_ideal_margin_closure_20261018'
TABLES='fresh_deep_finite8_dimension64_20261018'
BCD='fresh_order11_bcd_half_pair_transfer_20261018'
CYCLE='fresh_cycle11_half_pair_transfer_20261018'
ZERO='fresh_depth11_finite8_zero_transfer_20261018'
TRIANGLE='fresh_finite7_ideal5_triangle_payments_20261017'
STEMS=(CURRENT,C11,B12,D13,FINITE,BUDGET,IDEAL,OLD_IDEAL,TABLES,BCD,CYCLE,ZERO,TRIANGLE)


def graph(kind,n):
    if kind=='A': edges=[(i,(i+1)%n,3) for i in range(n)]
    elif kind=='B': edges=[(i,i+1,4 if i==n-4 else 3) for i in range(n-3)]+[(0,n-2,3),(0,n-1,3)]
    elif kind=='C': edges=[(i,i+1,4 if i in (0,n-2) else 3) for i in range(n-1)]
    else:
        assert kind=='D'
        edges=[(i,i+1,3) for i in range(n-5)]+[(0,n-4,3),(0,n-3,3),(n-5,n-2,3),(n-5,n-1,3)]
    g=[[0]*n for _ in range(n)]
    for i,j,m in edges:g[i][j]=g[j][i]=m
    return g


def distances(g):
    ds=[]
    for start in range(len(g)):
        d=[None]*len(g);d[start]=0;q=[start]
        for i in q:
            for j,m in enumerate(g[i]):
                if m and d[j] is None:d[j]=d[i]+1;q.append(j)
        assert all(x is not None for x in d);ds.append(d)
    return ds


def factor_interfaces(data):
    entries={}
    for row in data[TABLES]['short7_small_factor_interfaces']:
        entries[row['order'],row['kind']]={tuple(r['pair']):F(r['weight']) for r in row['weighted_pairs']}
    for kind in 'ABCD':
        rows=data[CYCLE]['weighted_pairs'] if kind=='A' else next(r for r in data[BCD]['factor_rules'] if r['kind']==kind)['weighted_pairs']
        entries[11,kind]={tuple(r['pair']):F(r['weight']) for r in rows}
    old_entries={key:dict(value) for key,value in entries.items()}
    entries[11,'C']={tuple(r['pair']):F(r['weight']) for r in data[C11]['weighted_pairs']}
    entries[12,'B']={tuple(p):F(1) for p in data[B12]['new_paid_pairs']}
    entries[13,'D']={tuple(r['pair']):F(r['weight']) for r in data[D13]['weighted_pairs']}
    expected={11:['99/2','54','54','53'],12:['60','63','63','62'],13:['65','73','73','147/2']}
    changed_expected={(11,'C'):[([0,10],'1/2','0')],(12,'B'):[([0,9],'1','0')],(13,'D'):[([0,8],'1','1/2')]}
    rows=[]
    for n,kind in product_order():
        g=graph(kind,n);ds=distances(g);ws=entries[n,kind];old=old_entries[n,kind]
        total=sum(ws.values(),F());assert str(total)==expected[n]['ABCD'.index(kind)]
        changes=[(list(p),str(old.get(p,F())),str(ws.get(p,F()))) for p in combinations(range(n),2)
                 if old.get(p,F())!=ws.get(p,F())]
        assert changes==changed_expected.get((n,kind),[])
        short4=[p for p in combinations(range(n),2) if ds[p[0]][p[1]]<=4]
        assert all(ws.get(p,F())==1 for p in short4)
        rows.append(dict(order=n,kind=kind,total_cost=str(total),changed_pairs=changes,
                         retained_short4_units=len(short4)))
    envelopes={str(n):str(max(F(r['total_cost']) for r in rows if r['order']==n)) for n in (11,12,13)}
    assert envelopes=={'11':'54','12':'63','13':'147/2'}
    return rows,envelopes,entries


def product_order():
    return [(n,k) for n in (11,12,13) for k in 'ABCD']


def run():
    data={tag:json.loads((ROOT/f'results/{tag}.json').read_text()) for tag in STEMS}
    for tag,doc in data.items():
        assert doc['status'].startswith('PASS'),tag
        for path,digest in doc['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==digest,(tag,path)
    factors,envelopes,weights=factor_interfaces(data)
    c=data[C11];assert c['new_cost']=='54' and c['strict_minimum']=='1/5' and c['finite_coupon_increment']=='0'
    b=data[B12];assert b['old_cost']==64 and b['new_cost']==63 and b['removed_pair']==[0,9]
    assert len(b['literal_cases'])==6 and all(F(r['total'])>=1 for r in b['literal_cases'])
    b_units=[]
    for row in b['literal_cases']:
        for corner in row['corners']:
            if corner.get('kind')=='whole_affine_D11':
                assert corner['marked_distance']==8 and corner['payment']=='1'
                g=graph('D',12);ds=distances(g)
                centers=[i for i,r in enumerate(g) if sum(bool(x) for x in r)==3]
                leaves=[i for i,r in enumerate(g) if sum(bool(x) for x in r)==1]
                ps=[tuple(sorted((i,j))) for i in centers for j in leaves if ds[i][j]==8]
                assert len(ps)==4 and all(weights[12,'D'].get(p,F())==1 for p in ps)
                b_units.append(dict(kind='D',order=12,marked_role='fork_center/short_leaf',distance=8,all_equivalent_pairs=ps))
            if corner.get('reason')=='whole_D12_center_long_endpoint_unit':
                assert corner['payment']=='1'
                assert any(p['type']=='D12' and set(corner['marks'])<=set(p['nodes']) for p in corner['components'])
    d=data[D13];assert d['results']['new_affine_total']=='147/2' and d['new_weight']=='1/2'
    assert len(d['immediate_affine_unit_corners'])==4 and len(d['finite_coupon_corners'])==8
    for row in d['immediate_affine_unit_corners']:
        kind=row['whole_factor'][7];pair=tuple(row['canonical_pair'])
        assert kind in 'BD' and weights[11,kind][pair]==1
        assert row['retained_payment']=='1'
    assert all(r['coupon_payment']==r['original_ideal_payment']=='1/2' and r['total']=='1' for r in d['finite_coupon_corners'])
    assert all(r['opposite_black_path']==4 and r['payment']=='1' for r in d['mixed_energy_one_cases'])
    finite=data[FINITE]
    assert finite['baseline_ideal_rule']==CURRENT and finite['ideal_weights_unchanged'] and finite['current_coupons_retained']
    assert finite['new_finite_component_costs']==dict(A16=114,B16=104,D16=114,D17=125)
    for row in finite['rules']:
        ds=distances(row['target']);ws={tuple(p['pair']):F(p['weight']) for p in row['weighted_pairs']}
        assert str(sum(ws.values(),F()))==row['new_cost']
        assert all(ds[i][j]>=10 and ws[i,j]==0 for i,j in row['deleted_pairs'])
        assert all(ws[tuple(p)]==1 for p in row['retained_basic_pairs'])
    assert finite['literal_interfaces']['complete_configurations']==216
    assert finite['literal_interfaces']['minimum_total']=='1'
    assert len(finite['proper_completion_cases'])==6 and len(finite['radical_cases'])==10
    budget=data[BUDGET]
    assert (budget['finite_maximum'],budget['finite_curvature'],budget['finite_slack'])==('449','899/2','1/2')
    assert budget['actual_attaining_forest']==['D17','D17','D29'] and budget['actual_attainer_cost']=='449'
    assert budget['complete_large_vectors']==1329 and budget['complete_keys']==993
    expected_new=[(['A4','A12'],'4'),(['A4','A13'],'4'),(['D5','A12'],'2'),
                  (['A2','A13'],'2'),(['A2','A14'],'2'),(['E6','A11'],'4'),
                  (['A2','B9'],'1'),(['A2','D10'],'1'),(['A3','D10'],'1/2')]
    actual_new=[([kind+str(rank) for kind,rank in r['types']],r['added_component_pair_cost']) for r in budget['actual_cross_coupon_lists']]
    assert actual_new==expected_new
    assert budget['representative_allowance']['actual_total']=='2' and budget['representative_allowance']['upper_total']=='4'
    assert [r['pair_cost'] for r in finite['coupons']]==[r[1] for r in expected_new[:6]]
    assert b['new_finite_coupon']['component_pair_cost']==expected_new[6][1]
    assert [r['component_pair_increment'] for r in d['new_coupons']]==[r[1] for r in expected_new[7:]]
    tri=data[TRIANGLE]
    assert tri['finite_distance_at_most']==7 and tri['ideal_distance_at_most']==5
    assert tri['ideal_feature_pair_distance_at_most']==7 and tri['whole_ideal_factor_order_at_most']==9
    assert tri['triangle_minimum']==2 and data[ZERO]['changed_weight']=='0'
    ideal=data[IDEAL];assert ideal['complete_terminal_profiles']==2204 and ideal['all_profiles_positive']
    assert ideal['required_new_total_cost_envelopes']==envelopes
    assert ideal['minimum_normalized_margin']['normalized_margin']=='3505/201344'
    assert ideal['minimum_absolute_margin']['absolute_margin']=='139/46'
    assert ideal['minimum_absolute_margin']['ranks']==[62]
    old_env={r['rank']:F(r['cost']) for r in data[OLD_IDEAL]['envelopes']}
    new_env={r['rank']:F(r['cost']) for r in ideal['envelopes']}
    for rank,value in new_env.items():
        assert value==(F(envelopes[str(rank+1)])/(rank+1) if rank in (10,11,12) else old_env[rank])
    assert ideal['global_additive_constant']=='4960'
    n=63;beta=F(5*sum(min(j,n-1-j) for j in range(n))-n*(n-1),2)
    assert beta==F(899,2) and beta-F(budget['finite_maximum'])==F(1,2)
    extrema_gain=5*sum(2*j-(n-1) for j in range(n) if 2*j>n-1)
    assert extrema_gain==5*(n*n//4)==4960
    current=data[CURRENT]
    assert current['dimension_bound']==63 and current['unrestricted_finite_volume'] and not current['arithmetic_hypothesis']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{tag}.{ext}' for tag in STEMS
                           for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md',
              ROOT/'results/fresh_cusp_link_shape_curvature_20261017.md']
    return dict(status='PASS_SPARSE_COUPON_ORDER11_13_UNRESTRICTED_DIMENSION62',
                sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
                dimension_bound=62,excluded_dimension=63,inherited_higher_tail_bound=63,
                unrestricted_finite_volume=True,arithmetic_hypothesis=False,
                triangle_minimum=2,quadrilateral_minimum=1,new_factor_interfaces=factors,
                new_total_cost_envelopes=envelopes,B12_retained_affine_unit_interfaces=b_units,
                finite_component_costs=finite['new_finite_component_costs'],new_finite_coupon_increments=expected_new,
                finite_maximum='449',finite_curvature='899/2',finite_surplus='1/2',
                all_cusp_absolute_surplus='139/46',strict_global_constant='4960',
                impossible_inequality='0 >= (1/2)*v_f + (139/46)*v_i + 4960',
                quadrilateral_case_split=[
                    'Any changed C11 endpoint corner: independent distinct path1/2 unit.',
                    'Otherwise any changed B12 center/terminal corner: retained affine D11 center/leaf unit, whole finite D12 unit, or two typed finite halves.',
                    'Otherwise any changed D13 center pair: retained original half plus whole A2/A3-middle-D10-center half, or distinct retained B/D11 center or short4 unit.',
                    'No changed ideal corner: apply complete sparse finite theorem under old63 ideal weights; actual ideal weights agree pointwise, and B12/D13 coupons only increase finite payments.'],
                geometric_coupon_budget_hypotheses_discharged=True,
                unchanged_compact_arithmetic_bound=21,unchanged_arbitrary_compact_bound=29,
                scope='Complete written actual-corner local-to-global composition. Exact frozen coupon/factor interfaces, complete rank63 finite maximum449 and independent2204-profile ideal theorem are combined. Original geometric censuses, all-rank classification and concentration proofs remain inherited inputs. No new example, external review or formal verification claimed.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    output=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert output==json.loads(path.read_text())
    else:path.write_text(json.dumps(output,indent=2)+'\n')
    print(output['status'],'bound',output['dimension_bound'],'factor interfaces',len(output['new_factor_interfaces']),flush=True)
