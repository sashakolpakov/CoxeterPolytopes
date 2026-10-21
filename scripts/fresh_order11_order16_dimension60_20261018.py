#!/usr/bin/env python3
"""Compose proved local transfers, finite420 and complete positive cusp margins.

This binder reconstructs the terminal partition domain and exact saved-row
scalar comparisons. It does not rerun geometric censuses, convolution or DP.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import fresh_rank61_complete_finite420_budget_20261018 as C

ROOT = C.ROOT
TAG = 'fresh_order11_order16_dimension60_20261018'
BC = 'fresh_order11_bc_zero_transfer_20261018'
D105 = 'fresh_order16_d_distance9_half_transfer_20261018'
FRONT = 'fresh_ambient61_payment_frontier_20261018'
FRONT105 = 'fresh_order16_d105_ambient61_frontier_20261018'
FINITE = C.TAG
FINITE_LOCAL = C.D16
TAIL = 'fresh_rank15_rank16_dimension61_20261018'
IDEAL_GEOMETRY = 'fresh_order14_ambient62_ideal_margin_closure_20261018'
DEPTH = 'fresh_depth8_finite8_zero_transfer_20261018'
STEMS = (BC, D105, FRONT, FRONT105, FINITE, FINITE_LOCAL, TAIL, IDEAL_GEOMETRY, DEPTH)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def terminal_domain(total):
    rows = []
    def recurse(left, minimum, prefix):
        if left == 0:
            rows.append(prefix)
            return
        for r in range(minimum, left+1):
            recurse(left-r, r, prefix+[r])
    recurse(total, 9, [])
    for r in range(1, 9):
        recurse(total-r, 9, [r])
    return sorted(rows)


def key(row):
    return row['source'], json.dumps(row['case'], sort_keys=True)


def run():
    docs = {stem: json.loads((ROOT/f'results/{stem}.json').read_text()) for stem in STEMS}
    for stem, data in docs.items():
        assert data['status'].startswith('PASS'), stem
        for path, expected in data['sources'].items():
            assert digest(ROOT/path) == expected, (stem, path)
    bc, ideal16, finite = (docs[s] for s in (BC, D105, FINITE))
    assert bc['new_family_totals'] == dict(A='99/2', B='53', C='52', D='53')
    assert bc['new_maximum'] == '53' and bc['new_normalized_rank10_envelope'] == '53/11'
    assert bc['coupon_upgrade'] == dict(small='whole A1', large='whole B10',
        large_role='doubled terminal', old_marked_weight='1/2', new_marked_weight='1',
        added_component_pair_cost='1/2')
    assert ideal16['new_family_totals'] == dict(A='80', B='105', C='103', D='105')
    assert ideal16['new_order16_maximum'] == '105'
    assert ideal16['new_finite_coupons'] == [] and ideal16['all_finite_weights_unchanged']
    for stem in (BC, D105, FINITE_LOCAL, DEPTH):
        assert docs[stem]['triangle_minimum'] == 2 and docs[stem]['quadrilateral_minimum'] == 1
    assert docs[DEPTH]['feature_depth'] == 8
    assert docs[FINITE_LOCAL]['new_costs'] == dict(D16='111', D17='119')
    assert finite['geometric_inputs_final'] and finite['new_residual_tables']
    assert finite['finite_maximum'] == finite['finite_curvature'] == '420'
    assert finite['finite_slack'] == '0'
    assert finite['actual_component_costs'] == C.ACTUAL
    assert (finite['complete_large_vectors'], finite['interaction_keys']) == (1132, 854)
    assert finite['optional_A1_B10_half_to_unit_absorption'] and not finite['optional_A1_B10_rule_assumed']
    # The only extra actual finite cost beyond the complete420 theorem is
    # half on A1/B10. All other actual and representative costs are identical.
    assert C.actual_cross(('A', 1), ('B', 10)) == 1
    assert C.upper_cross(('A', 1), ('D', 10)) == 2
    assert C.B.representative(('B', 10)) == ('D', 10)
    assert C.B.representative(('A', 1)) == ('A', 1)
    assert C.actual_cross(('A', 1), ('B', 10))+1 <= C.upper_cross(('A', 1), ('D', 10))
    # The D105 replacement spends this unchanged whole-component coupon.
    assert C.actual_cross(('D', 5), ('D', 11)) == 15
    assert ideal16['existing_coupon']['unchanged_total_for_all_three_D11_fork_roles'] == '15/2'
    assert ideal16['existing_coupon']['marked_weight'] == '1/2'
    for forest in finite['actual_attaining_forests']:
        ts = [(x[0], int(x[1:])) for x in forest]
        new_extra = ts.count(('A', 1))*ts.count(('B', 10))
        assert new_extra == 0 and C.cost(ts)+new_extra == 840
    # The latest finite literal audit and the D105 no-debit audit refer to
    # the same216 actual records, despite their different historical baselines.
    finite_rows = {key(r): r for r in docs[FINITE_LOCAL]['literal_rows']}
    ideal_rows = {key(r): r for r in ideal16['literal_rows']}
    assert len(finite_rows) == len(ideal_rows) == 216 and finite_rows.keys() == ideal_rows.keys()
    for k, row in finite_rows.items():
        other = ideal_rows[k]
        assert [c['marks'] for c in row['corners']] == [c['marks'] for c in other['corners']]
        assert all(F(c['new_debit']) == 0 for c in other['corners'])
        assert all(F(c['retained_payment']) >= 0 for c in row['corners'])
        assert sum(F(c['retained_payment']) for c in row['corners']) == F(row['total']) >= 1

    front = docs[FRONT]
    assert front['ambient_dimension'] == 61 and front['link_dimension'] == 60
    old = front['ideal']
    assert old['low9_concentration_only'] and not old['high_factor_concentration_used']
    expected = terminal_domain(60)
    assert len(expected) == old['complete_profiles'] == 2365
    assert sorted(r['ranks'] for r in old['rows']) == expected
    assert len({tuple(r['ranks']) for r in old['rows']}) == 2365
    envelopes = {int(k): F(v) for k, v in old['small_rank_envelopes'].items()}
    envelopes.update({r['order']-1: F(r['envelope']) for r in old['exact_family_tail']})
    assert set(range(1, 61)) <= envelopes.keys()
    assert all(envelopes[r] == F(r, 2) for r in range(1, 10))
    assert envelopes[10] == F(54, 11) and envelopes[15] == F(106, 16)
    new_envelopes = dict(envelopes)
    new_envelopes[10] = F(bc['new_normalized_rank10_envelope'])
    new_envelopes[15] = F(ideal16['new_normalized_rank15_envelope'])
    after105 = {tuple(r['ranks']): r for r in docs[FRONT105]['rows']}
    assert set(after105) == {tuple(r) for r in expected}
    rows = []
    for row in old['rows']:
        ranks = row['ranks']
        h = row['h_polynomial']
        vertices = prod(r+1 for r in ranks)
        assert len(h) == 61 and h == h[::-1] and all(x >= 1 for x in h)
        assert sum(h) == row['link_vertices'] == vertices
        moment = sum(min(j, 60-j)*x for j, x in enumerate(h))
        assert moment == row['G']
        curvature = F(5*moment-60*vertices, 2)
        assert curvature == F(row['curvature'])
        old_cost = sum(envelopes[r] for r in ranks)
        assert old_cost == F(row['normalized_cost'])
        old_margin = curvature/vertices-old_cost
        assert old_margin == F(row['normalized_margin'])
        assert vertices*old_margin == F(row['absolute_margin'])
        saving = F(ranks.count(10), 11)+F(ranks.count(15), 16)
        new_cost = sum(new_envelopes[r] for r in ranks)
        margin = curvature/vertices-new_cost
        assert margin == old_margin+saving
        assert margin == F(after105[tuple(ranks)]['normalized_margin'])+F(ranks.count(10), 11)
        assert margin > 0
        rows.append(dict(ranks=ranks, link_vertices=vertices, curvature=str(curvature),
            normalized_cost=str(new_cost), normalized_margin=str(margin), absolute_margin=str(vertices*margin)))
    normalized = min(F(r['normalized_margin']) for r in rows)
    absolute = min(F(r['absolute_margin']) for r in rows)
    assert normalized == F(3, 250) and absolute == 64
    d = 60
    G = sum(min(j, d-j) for j in range(d+1))
    curvature = F(5*G-d*(d+1), 2)
    extremum_difference = sum(2*j-d for j in range(d+1) if 2*j > d)
    constant = 5*extremum_difference
    assert (G, curvature, extremum_difference, constant) == (900, F(420), 930, 4650)
    assert constant == F(front['strict_global_constant'])
    assert docs[TAIL]['dimension_bound'] == 61
    assert docs[IDEAL_GEOMETRY]['triangle_minimum'] == 2
    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{ext}' for stem in STEMS
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    return dict(status='PASS_COMPLETE_ORDER11_ORDER16_DIMENSION60_COMPOSITION',
        sources={str(p.relative_to(ROOT)): digest(p) for p in paths},
        excluded_dimension=61, inherited_dimension_bound=61, dimension_bound=60,
        hypotheses=['finite-sided real hyperbolic Coxeter polytope', 'finite volume'],
        arithmeticity_required=False, compactness_required=False,
        simple_cusps_required=False, facet_count_bound_required=False,
        actual_component_costs=C.ACTUAL, last_extra_finite_coupon_cost='1/2',
        last_extra_finite_coupon='whole A1 / whole B10 doubled terminal, half upgraded to unit',
        upper_finite_rule_unchanged=True, finite_maximum='420', finite_curvature='420', finite_surplus='0',
        current_literal_interfaces=216, complete_terminal_ideal_profiles=2365,
        ideal_order11_costs=bc['new_family_totals'], ideal_order16_costs=ideal16['new_family_totals'],
        rows=rows, normalized_ideal_minimum=str(normalized), actual_ideal_absolute_surplus=str(absolute),
        normalized_minimizers=[r['ranks'] for r in rows if F(r['normalized_margin']) == normalized],
        absolute_minimizers=[r['ranks'] for r in rows if F(r['absolute_margin']) == absolute],
        low9_concentration_unchanged=True, high_factor_concentration_used=False,
        strict_global_constant=str(constant), contradiction='0 >= 0*v_f + 64*v_i + 4650 > 0',
        local_composition='Independent changed C11/B11 actual-corner proofs first. The D105 actual-corner transfer uses the unchanged wholeD5/D11 half and independent short8/depth8 escape. Latest finite420 alternatives and all216 literal records remain compatible; all finite additions including A1/B10 are counted. Triangles retain2 and quadrilaterals1.',
        higher_dimension_tail='The completed bound61 excludes every dimension>=62; the present contradiction excludes61.',
        new_finite_volume_example=False, external_refereeing=False, formal_proof_assistant_verified=False,
        scope='Final mathematical composition, independently rebuilt terminal partition domain and exact saved-polynomial scalar checks. Local geometry, classified factor costs, original convolution, finite Bellman coverage, low9 concentration, factor fibers and global shape-curvature inequality remain explicit completed theorem inputs. No historical census, DP, convolution or optimization is rerun.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        assert not path.exists()
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['contradiction'], flush=True)
