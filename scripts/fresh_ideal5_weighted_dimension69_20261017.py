#!/usr/bin/env python3
"""Two affine-B entry obstructions and the unrestricted ceiling69."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_ideal5_weighted_dimension70_20261017 as B

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_ideal5_weighted_dimension69_20261017'
LOCAL = 'fresh_b14_b15_singleton_completion_transfer_20261017'
SMALL = 'fresh_distance8_ideal5_small_factor_20261017'
FINITE = 'fresh_finite8_d12_b10_coupon_budget_20261017'
ENVELOPE = dict(B.ENVELOPE, **{'13': '43/7', '14': '32/5'})


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def budget(rank):
    return F(rank, 2) if rank <= 11 else F(ENVELOPE[str(rank)]) if rank <= 16 else F(163, 23)


def profile(ranks):
    row = B.profile(ranks)
    cost = sum(map(budget, ranks), F())
    assert F(row['normalized_cost'])-cost == F(ranks.count(13), 14)+F(ranks.count(14), 15)
    row['normalized_cost'] = str(cost)
    row['margin'] = str(F(row['normalized_curvature'])-cost)
    return row


def scalar_comparison():
    rows = [profile(ranks) for ranks in B.B.terminal_shapes(69)]
    assert len(rows) == 1339 and all(F(r['margin']) > 0 for r in rows)
    minimum = min(rows, key=lambda r: F(r['margin']))
    assert minimum['margin'] == '1819/51840' and minimum['ranks'] == [11, 11, 11, 11, 11, 14]
    beta = F(70*64, 8)
    # Equality is sufficient because the global shape inequality adds345.
    assert beta-8*70 == 0
    gain = budget(11)+budget(13)-2*budget(12)
    assert gain == F(19, 182) > 0
    boundary = profile((11,)*5+(13,))
    assert boundary['link_dimension'] == 68 and boundary['margin'] == '-2029/13608'
    return dict(reduced_profiles=rows, dimension_minima=[dict(ambient_dimension=70,
        profiles=len(rows), minimum_margin=minimum['margin'], minimizing_ranks=minimum['ranks'],
        finite_slack='0', strict_global_slack=345)], boundary_at_n69=boundary,
        concentration=dict(repeated_cap_ranks=[11, 13, 14, 15, 16, 17],
            at_most_one_low_part=list(range(1, 11)), at_most_one_each_middle_part=[12],
            at_most_one_high_part_minimum=18, double12_budget_increase=str(gain),
            cost_direction='preserved or increased; normalized curvature decreases',
            termination='Each nontrivial move strictly increases the bounded integer sum of squared parts.',
            scope='The same complete terminal domain is retained; only the duplicate12 cost gain changes.'))


def local_interface(local, small):
    assert not local['added_finite_coupons']
    assert local['triangle_minimum'] == 2 and local['quadrilateral_minimum'] == 1
    rows = []
    for rec in local['records']:
        order = rec['order']
        assert order in (14, 15)
        old_cost, new_cost, drop, count = (87, 86, [0, 10], 46) if order == 14 else (97, 96, [1, 11], 18)
        assert rec['inherited_cost'] == old_cost and rec['new_cost'] == new_cost
        assert rec['removed_pair'] == drop
        paid = set(map(tuple, rec['new_paid_pairs']))
        assert len(paid) == new_cost and tuple(drop) not in paid
        for key in ('triangle_pairs_retained', 'short8_pairs_retained'):
            assert set(map(tuple, rec[key])) <= paid
        assert rec['basic_pairs_retained']
        old = next(r for r in small['records'] if r['order'] == order and r['kind'] == 'B')
        assert old['target'] == rec['target']
        support = [r['case'] for r in old['case_embeddings'] if drop in r['marked_images']]
        assert support == rec['affected_cases'] == [r['index'] for r in rec['cases']]
        maps = [mp for case in rec['cases'] for mp in case['maps']]
        assert len(maps) == rec['affected_map_count'] == count
        assert rec['map_outcomes'] == {'projection_exclusion': count}
        for case in rec['cases']:
            end_proof = case['end_projection_proof']
            assert end_proof['both_ends_shortcut_bound'] <= 8
            assert all(F(x) > 0 for block in end_proof['inverse_end_blocks'] for row in block for x in row)
            for mp in case['maps']:
                witness = mp['witness']
                assert mp['outcome'] == 'projection_exclusion'
                assert witness in mp['forced_entries']
                assert witness['new_target_node'] in mp['missing_target_nodes']
                assert witness['norm_two_contact_square'] == '2' and witness['coxeter_label'] == 4
                assert F(witness['forced_projection_lower_bound']) == 2*F(witness['full_core_inverse_diagonal']) > 2
        rows.append(dict(order=order, type=rec['affine_type'], old_cost=old_cost, new_cost=new_cost,
            removed_pair=drop, affected_cases=support, affected_maps=count,
            minimum_forced_projection=str(min(F(mp['witness']['forced_projection_lower_bound']) for mp in maps)),
            new_finite_cost='0'))
    assert [r['order'] for r in rows] == [14, 15]
    return rows


def run():
    stems = (B.TAG, B.TAG+'_review', LOCAL, SMALL, FINITE)
    inputs = {}
    for stem in stems:
        data = json.loads((ROOT/f'results/{stem}.json').read_text())
        assert data['status'].startswith('PASS'), stem
        for name, checksum in data['sources'].items():
            assert digest(ROOT/name) == checksum, (stem, name)
        inputs[stem] = data
    base = inputs[B.TAG]
    assert base['dimension_bound'] == inputs[B.TAG+'_review']['dimension_bound'] == 70
    assert base['unrestricted_finite_volume'] and not base['arithmetic_hypothesis'] and base['cusp_restriction'] is None
    finite = inputs[FINITE]
    assert finite['finite_vertex_budget'] == '8n' and finite['dimension_maximum'] == 75
    interface = local_interface(inputs[LOCAL], inputs[SMALL])
    costs = {order: dict(row) for order, row in base['small_order_costs'].items()}
    for row in interface:
        assert costs[str(row['order'])]['B'] == row['old_cost']
        costs[str(row['order'])]['B'] = row['new_cost']
    for order, row in costs.items():
        assert F(max(row.values()), int(order)) == budget(int(order)-1)
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in stems
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    out = scalar_comparison()
    out.update(status='PASS_UNRESTRICTED_WEIGHTED_IDEAL5_DIMENSION69', dimension_bound=69,
        inherited_dimension_bound=70, excluded_dimensions=[70], unrestricted_finite_volume=True,
        arithmetic_hypothesis=False, cusp_restriction=None, finite_vertex_budget='8n',
        finite_budget_dimension_maximum=75, ideal_rank_envelope=ENVELOPE,
        small_order_costs=costs, affine_B_interfaces=interface,
        unchanged_compact_arithmetic_bound=21, unchanged_arbitrary_compact_bound=29,
        effective_compact_facet_ambient_bound=69, older_auxiliary_compact_facet_bound=79,
        sources={str(path.relative_to(ROOT)): digest(path) for path in paths},
        scope='Complete two-singleton affine-B composition, unchanged finite8 budget, all1339 reduced cusp profiles at70 and the completed70 higher tail. Zero finite slack is sufficient because the inherited global additive constant is345. No new example, external refereeing or formalization claimed.')
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], len(result['reduced_profiles']), result['dimension_minima'], flush=True)
