#!/usr/bin/env python3
"""Compose complete finite9/ideal5 payments without the distance8 overlap domain."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_ideal5_small_factor_payments_20261017 as S

ROOT = S.ROOT
TAG = 'fresh_ideal5_combined_payments_20261017'
TRI = 'fresh_ideal5_finite8_triangle_payments_20261017'
MIXED = 'fresh_finite8_ideal5_local_replacement_20261017'
FINITE = 'fresh_finite8_a17_d11_coupon_budget_20261017'
LOW = 'fresh_distance9_low_core_20261017'
ENVELOPE = {'1_to_12':'r/2', '13':'87/14', '14':'20/3', '15':'111/16',
            '16':'127/17', '17_and_larger':'163/23'}


def digest(p):
    h = sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20), b''):
            h.update(chunk)
    return h.hexdigest()


def finite_increment():
    rows = []
    for rank in range(18, 81):
        for kind in ('A_or_B', 'D'):
            g = [[0]*rank for _ in range(rank)]
            edges = ([(i,i+1) for i in range(rank-1)] if kind == 'A_or_B' else
                     [(i,i+1) for i in range(rank-3)]+[(0,rank-2),(0,rank-1)])
            for a,b in edges:
                g[a][b] = g[b][a] = 3
            adj,ds,_ = S.P.distance_data(g)
            added = [list(p) for p in combinations(range(rank),2)
                     if ds[p[0]][p[1]] == 9 and len(adj[p[0]]) != 1 and len(adj[p[1]]) != 1]
            expected = rank-(11 if kind == 'A_or_B' else 12)
            assert len(added) == expected < rank
            rows.append(dict(type=kind, rank=rank, added_pairs=added, count=expected))
    return dict(rows=rows, finite8_budget='8n', additional_bound='at most n',
        result='9n', dimension_maximum=80,
        proof='Ranks<=17 already pay all pairs. At ranks>=18 only finite A/B/D occur; all distance9 pairs number r-9. Exactly2 contain a leaf in A/B and3 in D, leaving r-11 or r-12 new pairs, at most r. Sum over irreducible components.')


def run():
    stems = (S.TAG, TRI, MIXED, S.D13, S.D14, FINITE, S.OLD16, S.P.TAG,
             'fresh_cycle14_ideal5_transfer_20261017',
             'fresh_cycle16_half_coupon_transfer_20261017', S.P.OLD, LOW)
    inputs = {stem:S.read(stem) for stem in stems}
    for stem,data in inputs.items():
        assert data['status'].startswith('PASS')
        for path,h in data['sources'].items():
            assert digest(ROOT/path) == h, (stem,path)
    assert len(inputs[TRI]['rows']) == 1025
    assert len(inputs[S.TAG]['patterns']) == 571
    assert len(inputs[S.P.OLD]['isolation']['cases']) == 949
    assert max(r['selected_corner']['component_order_bound'] for r in inputs[LOW]['extensions']['cases']) <= 12
    assert inputs[FINITE]['dimension_maximum'] == 80
    assert all(F(r['excess']) <= 0 for r in inputs[FINITE]['dimension_endpoints'])
    records = inputs[S.TAG]['records']
    triangle_checks = []
    for record in records:
        target = record['target']
        adj,ds,features = S.P.distance_data(target)
        near = {i for i in range(len(target)) if any(ds[i][f] <= 5 for f in features)}
        required = {p for p in combinations(range(len(target)),2)
                    if ds[p[0]][p[1]] <= 5 or (ds[p[0]][p[1]] <= 8 and set(p)&near)}
        assert required <= set(map(tuple,record['paid_pairs']))
        triangle_checks.append(dict(order=record['order'], kind=record['kind'],
            elementary_triangle_pairs=[list(p) for p in sorted(required)], all_retained=True))
    # Immediate whole-affine alternatives are order12/13, not order14.
    d13small = inputs[S.D13]['formulas']['whole_affine12_payments']
    d14small = inputs[S.D14]['profile_formulas']['small_affine_rows']
    assert len(d13small) == len(d14small) == 2
    assert {r['order'] for r in d13small} == {12}
    assert {r['order'] for r in d14small} == {13}
    for stem,key,unit_key in ((S.D13,'extreme_equality','units'),(S.D14,'equality','finite_units')):
        assert len(inputs[stem][key]['rows']) == 72
        assert all(r[unit_key]+F(r['coupon_total']) >= 1 for r in inputs[stem][key]['rows'])
    old16 = inputs[S.OLD16]['exceptional_D15']['cases']
    old17 = inputs[S.P.TAG]['order17']['replacement_cases']
    assert len(old16) == len(old17) == 9
    for row in old16+old17:
        assert sum(F(c['weight']) for c in row['corners'] if c['finite']) >= 1
        assert all(F(c['weight']) == 0 for c in row['corners'] if not c['finite'])
    costs = {}
    for row in records:
        costs.setdefault(str(row['order']),{})[row['kind']] = row['cost']
    assert costs == {'14':dict(A=70,B=87,C=84,D=82), '15':dict(A=75,B=100,C=95,D=96),
                     '16':dict(A=80,B=107,C=103,D=111), '17':dict(A=85,B=111,C=108,D=127)}
    maxima = {m:str(F(max(v.values()),int(m))) for m,v in costs.items()}
    assert maxima == {'14':'87/14','15':'20/3','16':'111/16','17':'127/17'}
    assert inputs[S.TAG]['basic_affine_cost']['normalized_common_bound'] == '163/23'
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in stems
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths.append(ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md')
    return dict(status='PASS_COMPLETE_FINITE9_IDEAL5_LOCAL_PAYMENT_AND_BUDGET',
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        small_order_pair_sets=records, small_order_costs=costs,
        small_order_normalized_maxima=maxima, triangle_pair_retention=triangle_checks,
        ideal_rank_envelope=ENVELOPE, finite_vertex_budget='9n', finite_budget_dimension_maximum=80,
        finite9_increment=finite_increment(), triangle_minimum='2', quadrilateral_minimum='1',
        immediate_ideal_alternatives=dict(D13=d13small,D14=d14small,whole_factor_cutoff=13),
        old_order16_17_replacements=dict(order16_cases=[r['case'] for r in old16],
            order17_cases=[r['case'] for r in old17], all_terminal_payments_finite=True,
            nonorthogonal_support_alternative='marked path of length2'),
        low_core_maximum_factor_order=12, inherited_high_overlap_cases=949,
        large_order_linear_branches=inputs[S.TAG]['basic_affine_cost']['all_rank_branches'],
        scope='Complete finite9/ideal5 local payment composition, with actual cusp-factor multiplicity and piecewise factor costs. Finite8 whole-overlap coverage is not asserted. Ambient comparison is separate; published classification and inherited geometric proofs remain explicit inputs.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify',action='store_true')
    args = ap.parse_args()
    data = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['small_order_costs'],flush=True)
