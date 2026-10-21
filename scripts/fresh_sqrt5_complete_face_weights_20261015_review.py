#!/usr/bin/env python3
"""Independent graph/DP verification of a rule including the full finite queue.

The frozen normalized-cross reviewer supplies a separate all-root graph
implementation and exact type-first forest DP. No weight producer or
optimizer is imported. Full mathematical coverage belongs to its binder.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_all_a1_cross_corrected_20261015_review as V

R, X = V.R, V.X
ROOT = Path(__file__).resolve().parents[1]
PRIMARY = 'fresh_sqrt5_complete_face_weights_20261015'
TAG = PRIMARY + '_review'
FINITE = 'fresh_sqrt5_nonsingular_rank29_fast_20261015'


def read(tag):
    return json.loads((ROOT / f'results/{tag}.json').read_text())


def run(dimension):
    primary = read(PRIMARY + '_n' + str(dimension))
    assert primary['dimension'] == dimension and primary['meets_target']
    rule = V.Rule(primary)
    assert primary['solution']['bonus_grid'] == 1000
    assert all((v * 1000).denominator == 1 for v in rule.x[:-1])
    old = read('fresh_sqrt5_parameter_weight_stress_20261015')
    queue = read(FINITE)
    assert queue['complete_queue'] and queue['state']['next_core'] == queue['marked_states'] == len(queue['records'])
    assert not queue['counts'].get('PENDING', 0)
    actual = [r['index'] for r in queue['records'] if r['result']['singular']]
    assert actual == queue['weight_obligations']
    faces = [('prior_catalogue', i, r['row']) for i, r in enumerate(old['faces'])
             if r['row']['nodes'] - 2 <= dimension]
    old_count = len(faces)
    faces += [('complete_finite_queue', r['index'], r['row']) for r in queue['records']
              if r['result']['singular'] and r['row']['nodes'] - 2 <= dimension]
    assert len(faces) == primary['faces']
    # The producer deliberately keeps literal repeated faces. No presumed
    # payment or deduplication is allowed to omit a singular stopping state.
    expected_connected = set()
    for rec in old['orbits']:
        if rec['rank'] <= dimension:
            g = R.graph(rec['nodes'], [(a, b, m) for (a, b), m in
                        zip(combinations(rec['nodes'], 2), rec['labels'])])
            expected_connected.add(R.tree_key(g, rec['marked']))
    expected_cross = set()
    records, minima = [], Counter()
    cases = 0
    minimum = None
    for family, index, row in faces:
        cs = V.corners(row)
        for c in cs:
            if c['connected'] is not None and c['baseline'] == 0:
                expected_connected.add(c['connected'])
            if c['cross'] is not None:
                expected_cross.add(c['cross'])
        payments = rule.payment(row)
        for c in payments:
            v = F(c['payment'])
            minima[str(v)] += 1
            minimum = v - 1 if minimum is None else min(minimum, v - 1)
        cases += len(payments)
        records.append(dict(family=family, index=index, nodes=row['nodes'],
                            cases=[dict(external_A1=c['external_A1'], payment=c['payment'],
                                        corner_weights=c['corner_weights']) for c in payments]))
    assert expected_connected == set(rule.connected) and expected_cross == set(rule.cross)
    assert cases == primary['face_cases']
    assert minimum == F(primary['exact_certificate']['minimum_face_slack'])
    # Original literal inputs and exact dotted parameters remain bound even
    # though the final payment calculation uses their reviewed quotient.
    raw = X.inputs(); assert len(raw) == len(old['literal_maps']) == 2397
    literal_count = 0
    for (family, index, row, ratio), mp in zip(raw, old['literal_maps']):
        assert (family, index) == (mp['family'], mp['index'])
        target = old['faces'][mp['representative']]
        assert ratio == target['dotted_ratio']
        R.check_bijection(R.label_matrix(row), R.label_matrix(target['row']), mp['to_representative'])
        assert (row['nodes'] <= dimension + 2) == (target['row']['nodes'] <= dimension + 2)
        literal_count += row['nodes'] <= dimension + 2
    models, totals, fronts, connected_mult, cross_mult = V.component_bounds(rule)
    assert len(fronts) == len(primary['exact_certificate']['fronts']) == dimension
    for ours, saved in zip(fronts, primary['exact_certificate']['fronts']):
        assert ours['context'] == saved['context'] and F(ours['maximum']) == F(saved['maximum'])
        a = ours['context']; forest = [R.alias(t) for t in saved['maximizing_forest']]
        assert sum(len(models[t]) for t in forest) == dimension and forest.count('A1') == a
        assert sum(totals[t][a] for t in forest) == F(ours['maximum'])
    maximum = max(F(r['maximum']) for r in fronts)
    target = F(dimension*(dimension-6), 8) if dimension % 2 == 0 else F((dimension-1)*(dimension-5), 8)
    assert maximum == rule.x[-1] == F(primary['solution']['feasible_maximum'])
    assert maximum == F(primary['exact_certificate']['maximum']) <= target == F(primary['target'])
    triangle = read('fresh_sqrt5_triangle_weights_20261015_review')
    assert len(triangle['markings']) == 76
    assert all(sum(c['h4_only_weight'] for c in r['corners']) >= 2 for r in triangle['markings'])
    paths = [f'scripts/{TAG}.py', f'scripts/{V.TAG}.py', f'results/{V.TAG}.json',
             f'scripts/{PRIMARY}.py', f'results/{PRIMARY}_n{dimension}.json',
             f'results/{FINITE}.json', 'results/fresh_sqrt5_parameter_weight_stress_20261015.json',
             'results/fresh_sqrt5_parameter_weight_stress_20261015_cross_review.json',
             'results/fresh_sqrt5_triangle_weights_20261015_review.json']
    for p in paths:
        if p.endswith('.json'):
            for q, digest in json.loads((ROOT / p).read_text()).get('sources', {}).items():
                assert sha256((ROOT / q).read_bytes()).hexdigest() == digest, (p, q)
    return dict(status='PASS_INDEPENDENT_EXACT_WEIGHTS_ON_PRIOR_CATALOGUE_AND_ALL_FINITE_SINGULAR_STOPS',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                dimension=dimension, original_literal_inputs=2397, dimension_feasible_literals=literal_count,
                old_parameter_marked_faces=old_count, finite_queue_states=queue['marked_states'],
                all_finite_singular_stops=len(actual), dimension_feasible_finite_stops=len(faces)-old_count,
                face_rows=len(faces), face_cases=cases, minimum_face_slack=str(minimum),
                variables=primary['variables'], connected_orbits=len(rule.connected), rooted_cross_orbits=len(rule.cross),
                finite_component_types=len(models), exact_A1_contexts=len(fronts),
                vertex_maximum=str(maximum), strict_face_threshold=str(target), threshold_margin=str(target-maximum),
                forest_maxima=fronts, face_payments=records,
                connected_pair_multiplicities=[dict(label=rule.connected[k]['label'], multiplicity=v)
                                               for k, v in sorted(connected_mult.items())],
                cross_root_multiplicities=[dict(label=rule.cross[k]['label'], multiplicity=v)
                                           for k, v in sorted(cross_mult.items())],
                method='No optimizer or producer imports. Frozen separate all-root labelled trees, intrinsic type recognition, exact rational face payments, explicit A1/root normalization, and reverse-type forest DP. Every singular stopping state is evaluated, including previously paid and expanded states.',
                scope='Weight/payment review only. The complete arithmetic attachment queues, minimal-union reductions, arbitrary-white sealing and the global coverage composition remain separate reviewed dependencies. No optimality or dimension theorem is asserted here.')


if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('--dimension', type=int, required=True, choices=(25, 26, 27));p.add_argument('--verify', action='store_true');a = p.parse_args()
    out = json.loads(json.dumps(run(a.dimension)))
    path = ROOT / f'results/{TAG}_n{a.dimension}.json'
    if a.verify:
        assert out == json.loads(path.read_text());print('PASS saved full-catalogue weight review')
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], a.dimension, out['vertex_maximum'], out['face_cases'])
