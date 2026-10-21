#!/usr/bin/env python3
"""Bind completed local transfers, rank62 finite count and positive cusp margins.

This is the final composition/scalar check, not a historical geometric census
or a second implementation of the inherited finite-forest verifier.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_rank15_rank16_dimension61_20261018'
R15 = 'fresh_finite_rank15_existing_half_transfer_20261018'
R16 = 'fresh_rank16_existing_payment_exchange_20261018'
FINITE = 'fresh_finite_rank15_rank16_rank62_budget_20261018'
IDEAL = 'fresh_order14_ambient62_ideal_margin_closure_20261018'
D17 = 'fresh_finite_d17_five_unit_reduction_20261018'
TAIL = 'fresh_sparse_coupon_order11_13_dimension62_20261018'
STEMS = (R15, R16, FINITE, IDEAL, D17, TAIL)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run():
    docs = {s: json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem, data in docs.items():
        assert data['status'].startswith('PASS'), stem
        for name, expected in data['sources'].items():
            assert digest(ROOT/name) == expected, (stem, name)
    r15, r16, finite, ideal = (docs[s] for s in (R15, R16, FINITE, IDEAL))
    assert r15['new_costs'] == dict(A15='104', B15='104', D15='104')
    assert r16['new_costs'] == dict(A16='113', D16='225/2', A17='115')
    for local in (r15, r16):
        assert local['new_cross_component_coupons'] == []
        assert local['triangle_minimum'] == 2
        assert local['quadrilateral_minimum'] == 1
    assert r16['A16_half_pairs'] == [[0,11], [4,15]]
    assert r16['D16_half_pairs'] == [[2,13], [2,14], [2,15]]
    assert r16['A17_restored_half_pairs'] == [[1,12], [4,15]]
    assert r16['literal_configurations'] == 216 and r16['literal_changed_corners'] == {}
    assert r15['simultaneous_rank16_exchange'] == r16['new_costs']
    assert r15['preserved_D16_alternative_pairs'] == [[1,14], [1,15]]
    assert r15['preserved_D16_alternative_weight'] == '1/2'
    assert not (set(map(tuple,r15['preserved_D16_alternative_pairs'])) &
                set(map(tuple,r16['D16_half_pairs'])))
    assert r15['literal_configurations'] == 216
    assert len(r15['literal_rows']) == 216
    assert r15['literal_changed_corners'] == {'D15':16} and r15['rank16_literal_changes'] == {}
    for row in r15['literal_rows']:
        assert all(F(c['retained_payment']) >= 0 for c in row['corners'])
        assert sum(F(c['retained_payment']) for c in row['corners']) == F(row['total']) >= 1
    # The new selected replacements use disjoint contexts and retained payments.
    a = {r['index'] for r in r15['selected_proofs']}
    b = {r['index'] for r in r16['selected_proofs']}
    assert len(a) == 6 and len(b) == 14 and a.isdisjoint(b)
    for local in (r15, r16):
        for row in local['selected_proofs']:
            assert sum(F(c['payment']) for c in row['corners']) == F(row['total']) == 1
    assert docs[D17]['new_cost'] == '120' and docs[D17]['new_finite_coupons'] == []
    assert docs[TAIL]['dimension_bound'] == 62
    assert ideal['ambient_dimension'] == finite['ambient_dimension'] == 62
    assert ideal['inherited_unrestricted_ceiling'] == 62
    assert ideal['finite_weights_unchanged'] and ideal['new_finite_coupons'] == []
    assert ideal['triangle_minimum'] == 2 and ideal['quadrilateral_minimum'] == 1
    assert ideal['order14_family_costs'] == dict(A=70, B=82, C=83, D=82)
    assert finite['finite_maximum'] == finite['finite_curvature'] == '434'
    assert finite['finite_slack'] == '0'
    assert finite['geometric_inputs_final'] and finite['new_cross_component_coupons'] == []
    assert finite['complete_large_vectors'] == 1227 and finite['inherited_residual_tables'] == 921
    assert finite['actual_attaining_forest'] == ['D17', 'D45']
    # The successor must bind both actual local rules, including the A17 rise.
    for stem in (R15, R16):
        for folder, suffix in (('scripts','py'), ('results','json'), ('results','md')):
            name = f'{folder}/{stem}.{suffix}'
            assert finite['sources'][name] == digest(ROOT/name)

    rows = ideal['rows']
    assert len(rows) == ideal['complete_profiles'] == 2641
    assert len({tuple(r['ranks']) for r in rows}) == len(rows)
    for row in rows:
        ranks = row['ranks']
        assert ranks == sorted(ranks) and sum(ranks) == 61
        assert all(r >= 1 for r in ranks) and sum(r < 9 for r in ranks) <= 1
        assert row['link_vertices'] == prod(r+1 for r in ranks)
        assert F(row['absolute_margin']) == row['link_vertices']*F(row['normalized_margin'])
        assert F(row['normalized_margin']) > 0
    normalized = min(F(r['normalized_margin']) for r in rows)
    absolute = min(F(r['absolute_margin']) for r in rows)
    assert normalized == F(2741,371712) and absolute == 73
    assert ideal['all_actual_cusp_absolute_surplus'] == '73'
    # The ideal theorem supplies the low9 concentration and actual factor fibers.
    # The current finite rule has zero deficit, not a strictly positive surplus.
    d = 61
    finite_link_vertices = d+1
    finite_link_G = sum(min(j,d-j) for j in range(d+1))
    finite_curvature = F(5*finite_link_G-d*finite_link_vertices,2)
    extremum_difference = sum(2*j-d for j in range(d+1) if 2*j>d)
    assert finite_curvature == F(finite['finite_curvature']) == 434
    assert extremum_difference == 961
    finite_surplus = finite_curvature-F(finite['finite_maximum'])
    strict_constant = F(ideal['strict_global_constant'])
    assert finite_surplus == 0 and strict_constant == 5*extremum_difference == 4805
    assert absolute > 0 and strict_constant > 0

    paths = [Path(__file__)] + [ROOT/f'{folder}/{stem}.{suffix}'
        for stem in STEMS for folder, suffix in (('scripts','py'), ('results','json'), ('results','md'))]
    return dict(status='PASS_COMPLETE_RANK15_RANK16_DIMENSION61_COMPOSITION',
        sources={str(p.relative_to(ROOT)):digest(p) for p in paths},
        excluded_dimension=62, inherited_dimension_bound=62, dimension_bound=61,
        hypotheses=['finite-sided real hyperbolic Coxeter polytope', 'finite volume'],
        arithmeticity_required=False, compactness_required=False,
        simple_cusps_required=False, facet_count_bound_required=False,
        actual_component_costs=dict(A15='104', B15='104', D15='104',
            A16='113', B16='104', D16='225/2', A17='115', B17='105', D17='120'),
        finite_maximum='434', finite_curvature='434', finite_surplus='0',
        actual_finite_attaining_forest=['D17','D45'],
        complete_terminal_ideal_profiles=2641,
        normalized_ideal_minimum=str(normalized), actual_ideal_absolute_surplus=str(absolute),
        strict_global_constant=str(strict_constant),
        contradiction='0 >= 0*v_f + 73*v_i + 4805 > 0',
        local_composition='Current ideal priorities and independent short8 route first. Rank15 uses unchanged whole A2/A13 endpoint and D16(1,14)/(1,15) halves. Rank16 uses unchanged whole-small/D13 halves and restored A17 halves. D17 uses unchanged whole-small/D12 halves. Selected changed contexts are patched simultaneously; rank16 has no current literal debit, and rank15 checks all216 literal obligations after D17=120. All other nonselected, triangle and later outputs are retained by the local theorems.',
        higher_dimension_tail='The completed bound62 excludes all dimensions>=63; the present contradiction excludes62.',
        new_finite_volume_example=False, external_refereeing=False, formal_proof_assistant_verified=False,
        scope='Final source-bound mathematical composition and exact scalar inequalities. Local geometry, complete finite-type/Bellman coverage, cusp-domain completeness, low9 concentration, factor fibers and global shape-curvature inequality remain completed bound theorem inputs. No historical census, DP, scalar convolution or optimization is rerun.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        assert not path.exists(), 'Frozen output exists; use --verify.'
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], 'bound', out['dimension_bound'], out['contradiction'], flush=True)
