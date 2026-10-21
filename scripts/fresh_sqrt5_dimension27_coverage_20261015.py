#!/usr/bin/env python3
"""Compose the reviewed Qsqrt5 face reductions and current cross weights.

This binder does not enumerate a new arithmetic domain. It recomputes every
required zero-touch payment from the literal graphs, checks the finite
reduction links, and requires complete independently reviewed finite/dotted
queues. --stage checks only the already frozen branches and writes no result.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_parameter_weight_stress_20261015_cross_review as X

R = X.R
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_dimension27_coverage_20261015'
FINITE = 'fresh_sqrt5_nonsingular_rank29_fast_20261015'
FINITE_PREDECESSOR = 'fresh_sqrt5_nonsingular_rank29_exhaustion_20261015'
DOTTED = 'fresh_sqrt5_dotted_rank29_exhaustion_20261015'
WEIGHTS = 'fresh_sqrt5_parameter_weight_stress_20261015_cross'


class Binder:
    def __init__(self):
        self.paths = {f'scripts/{TAG}.py'}
        self.manifests = set()
        self.references = 0
        self.payments = []
        self.maps = 0
        self.primary = self.load(WEIGHTS)
        self.review(WEIGHTS + '_review')
        self.paths.update([f'scripts/{X.TAG}.py', f'scripts/{R.TAG}.py'])
        self.values = list(map(F, self.primary['solution']['exact_primal']))
        assert self.values[-1] == F(230, 3) and min(self.values) >= 0
        labels = {}
        for rec in self.primary['connected_orbits']:
            g = R.graph(rec['nodes'], [(i, j, m) for (i, j), m in
                        zip(combinations(rec['nodes'], 2), rec['labels'])])
            labels[rec['key']] = ('connected', R.tree_key(g, rec['marked']))
        for rec in self.primary['cross_orbits']:
            g = R.graph(rec['nodes'], [(i, j, m) for (i, j), m in
                        zip(combinations(rec['nodes'], 2), rec['labels'])])
            labels[rec['key']] = ('cross', R.tree_key(g, [rec['marked_root']]))
        self.weights = {}
        for i, rec in enumerate(self.primary['variable_order']):
            kind, key = labels[rec['orbit']]
            item = kind, key, rec['A1_count_capped2']
            assert item not in self.weights
            self.weights[item] = self.values[i]
        assert len(self.weights) == 922

    def load(self, tag):
        path = f'results/{tag}.json'
        self.paths.add(path)
        data = json.loads((ROOT / path).read_text())
        assert 'sources' in data, path
        return data

    def review(self, tag):
        data = self.load(tag)
        assert data['status'].startswith('PASS'), (tag, data['status'])
        return data

    def bind_manifests(self):
        def check(path):
            if path in self.manifests:
                return
            self.manifests.add(path)
            data = json.loads((ROOT / path).read_text())
            for p, digest in data.get('sources', {}).items():
                self.references += 1
                assert sha256((ROOT / p).read_bytes()).hexdigest() == digest, (path, p)
                if p.endswith('.json'):
                    check(p)
        for path in sorted(self.paths):
            if path.endswith('.json'):
                check(path)
        return {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in sorted(self.paths)}

    def payment(self, row, origin):
        cs = X.corners(row)
        rank = row['nodes'] - 2
        assert rank <= 29
        records = []
        for dimension in (28, 29):
            residual = dimension - rank
            if residual < 0:
                continue
            feasible = sorted({min(2, a) for a in range(residual + 1)
                               if residual - a == 0 or residual - a >= 2})
            for external in feasible:
                weights = []
                for c in cs:
                    context = min(2, c['A1'] + external)
                    value = c['baseline']
                    if c['key'] is not None:
                        value += self.weights.get(('connected', c['key'], context), F(0))
                    if c['cross'] is not None and c['A1'] + external == 1:
                        value += self.weights.get(('cross', c['cross'], 1), F(0))
                    weights.append(value)
                total = sum(weights)
                assert total >= 1, (origin, dimension, external, total)
                records.append(dict(dimension=dimension, external_A1_capped2=external,
                                    corners=list(map(str, weights)), total=str(total)))
        assert records
        self.payments.append(dict(origin=origin, nodes=row['nodes'], cases=records))

    def check_map(self, row, target, mapping):
        R.check_bijection(R.label_matrix(row), R.label_matrix(target), mapping)
        self.maps += 1


def extension(row, contacts):
    n = row['nodes']
    assert len(contacts) == n
    old = R.label_matrix(row)
    return dict(nodes=n + 1, white=list(range(4, n + 1)),
                labels=[contacts[i] if j == n else old[i][j]
                        for i, j in combinations(range(n + 1), 2)])


def induced(row, ids):
    assert ids[:4] == list(range(4))
    old = R.label_matrix(row)
    return dict(nodes=len(ids), white=list(range(4, len(ids))),
                labels=[old[i][j] for i, j in combinations(ids, 2)])


def triangles_and_vertices(B):
    triangle = B.review('fresh_sqrt5_triangle_weights_20261015_review')
    B.review('fresh_compact_global_bound_20261015_review')
    B.paths.add('results/arithmetic_dimension_bounds_and_effective_high_dimension_search.md')
    assert len(triangle['markings']) == 76
    totals = []
    for rec in triangle['markings']:
        n = rec['order']
        edges = [(i, j, m) for (i, j), m in
                 zip(combinations(range(n), 2), rec['labels'])]
        face_total = F(0)
        for omit in rec['black']:
            marked = [i for i in rec['black'] if i != omit]
            parts = R.components(R.graph([i for i in range(n) if i != omit], edges))
            g = next(h for h in parts if marked[0] in h)
            face_total += R.baseline(g, *marked) if marked[1] in g else 0
        assert face_total >= 2
        totals.append(str(face_total))
    component_totals = {}
    for name, g in R.MODELS.items():
        values = [F(0)] * 3
        for a, b in combinations(g, 2):
            key = R.tree_key(g, (a, b))
            for context in range(3):
                values[context] += R.baseline(g, a, b) + B.weights.get(
                    ('connected', key, context), F(0))
        if name == 'E7':
            for node in g:
                values[1] += B.weights.get(('cross', R.tree_key(g, [node]), 1), F(0))
        component_totals[name] = values
    dp = R.forest_dp(component_totals)
    maxima = {str(d): [str(dp[c][d, c]) for c in range(3)] for d in (28, 29)}
    assert maxima == {'28': ['230/3'] * 3, '29': ['239/3', '84', '236/3']}
    hvector = []
    for n, maximum in ((28, F(230, 3)), (29, F(84))):
        paired = [F((n-i)*(n-i-1) + i*(i-1), 4) for i in range(n + 1)]
        minimum = min(paired)
        assert paired[0] > minimum  # h0=hn=1 makes the ratio strictly larger.
        threshold = 5 * minimum - F(n*(n-1), 2)
        assert threshold == {28: 77, 29: 84}[n] and maximum <= threshold
        hvector.append(dict(dimension=n, minimum_paired_f2_coefficient=str(minimum),
                            strict_face_weight_threshold=str(threshold),
                            vertex_upper_bound=str(maximum),
                            strictness='h_0=h_n=1 and the extreme paired coefficient is larger'))
    return dict(markings=76, baseline_totals=totals, component_types=92,
                forest_maxima=maxima, hvector_contradictions=hvector)


def old_branches(B):
    union = B.load('fresh_sqrt5_quad_exhaustion_20261015')
    B.review('fresh_sqrt5_quad_exhaustion_20261015_incidence_review')
    base = union['labelled_union_supports']
    assert len(base) == 5440
    filters = B.load('fresh_sqrt5_quad_filters_20261015')
    B.review('fresh_sqrt5_quad_filters_20261015_review')
    assert len(filters['records']) == 1180
    assert filters['counts'] == {'DOTTED_PARAMETER_UNFILTERED': 34,
        'REJECT_TWO_NEGATIVE_DIRECTIONS': 617, 'PASS_PRINCIPAL_UNION_NECESSARY_TESTS': 362,
        'REJECT_CONJUGATE_NOT_PSD': 78, 'REJECT_CYCLE_FIELD': 89}
    contacts = B.load('fresh_sqrt5_quad_contacts_20261015')
    B.review('fresh_sqrt5_quad_contacts_20261015_review')
    assert len(contacts['records']) == 396
    survivor_ids = {r['class_index'] for r in filters['records'] if r['status'] in
                     ('DOTTED_PARAMETER_UNFILTERED', 'PASS_PRINCIPAL_UNION_NECESSARY_TESTS')}
    assert survivor_ids == {r['class_index'] for r in contacts['records']}
    assert contacts['counts'] == {'ALREADY_PAID_BY_ENTIRE_H4': 84,
        'UNPAID_DOTTED_CONTACT_DOMAIN': 24, 'REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM': 2,
        'UNPAID_FINITE_CONTACT_DOMAIN': 286}
    sealed = {}
    accepted = {}
    finite_rank_rejects = []
    dotted_ids = []
    for rec in contacts['records']:
        row = base[rec['union_index']]
        if rec['status'] == 'ALREADY_PAID_BY_ENTIRE_H4':
            assert sum(c['baseline'] for c in X.corners(row)) >= 1
        elif rec['status'] == 'UNPAID_DOTTED_CONTACT_DOMAIN':
            dotted_ids.append(rec['class_index'])
        elif rec['status'] == 'UNPAID_FINITE_CONTACT_DOMAIN':
            if not rec['accepted']:
                if rec['singular']:
                    sealed[rec['class_index']] = row
                    B.payment(row, ['old_sealed', rec['class_index']])
                else:
                    finite_rank_rejects.append(rec['class_index'])
            else:
                assert not rec['singular']
                for ti, trial in enumerate(rec['trials']):
                    if trial['status'] == 'PASS_EXTRA_WHITE_NECESSARY_TESTS':
                        key = rec['class_index'], ti
                        assert key not in accepted
                        accepted[key] = extension(row, trial['contacts'])
    assert len(sealed) == 112 and len(finite_rank_rejects) == 44 and len(accepted) == 458
    sub = B.load('fresh_sqrt5_sealed_subcore_closure_20261015')
    B.review('fresh_sqrt5_sealed_subcore_closure_20261015_review')
    assert sub['counts'] == {'OPEN_NONSINGULAR': 311, 'OPEN_SINGULAR': 61,
                            'PAID_BY_REVIEWED_SEALED_SUBCORE': 86}
    assert {(r['source_class_index'], r['trial_index']) for r in sub['records']} == set(accepted)
    open_singular, open_nonsingular = [], []
    for rec in sub['records']:
        row = accepted[rec['source_class_index'], rec['trial_index']]
        if rec['singular']:
            B.payment(row, ['singular_first_extension', rec['source_class_index'], rec['trial_index']])
        if rec['status'] == 'PAID_BY_REVIEWED_SEALED_SUBCORE':
            for match in rec['sealed_subcores']:
                B.check_map(induced(row, match['nodes_in_extension']),
                            sealed[match['class_index']], match['to_sealed_representative'])
            assert rec['sealed_subcores']
        elif rec['status'] == 'OPEN_SINGULAR':
            open_singular.append(row)
        else:
            assert rec['status'] == 'OPEN_NONSINGULAR'
            open_nonsingular.append(row)
    first = B.load('fresh_sqrt5_singular_extension_frontier_20261015')
    B.review('fresh_sqrt5_extended_forest_weights_20261015_review')
    assert len(first['records']) == 51 and first['input_literal_extensions'] == 61
    reached = set()
    for rec in first['records']:
        B.payment(rec['row'], ['singular51', rec['index']])
        for mp in rec['origins']:
            i = mp['open_index']; assert i not in reached; reached.add(i)
            B.check_map(open_singular[i], rec['row'], mp['to_representative'])
    assert reached == set(range(61))
    closure = B.load('fresh_sqrt5_singular_attachment_closure_20261015')
    B.review('fresh_sqrt5_singular_attachment_closure_20261015_review')
    assert closure['complete'] and closure['processed'] == closure['marked_states'] == 17
    assert closure['counts'] == {'REDUCES_TO_KNOWN_158_SEALED_CATALOGUE': 6, 'NEW_SEALED_PAID': 11}
    for rec in closure['records']:
        B.payment(rec['row'], ['singular_closure17', rec['index']])
        for origin in rec['origins']:
            if origin['kind'] == 'seven_continuations':
                source = first['records'][origin['source_index']]
                assert origin['contacts'] in source['accepted']
                B.check_map(extension(source['row'], origin['contacts']), rec['row'],
                            origin['to_representative'])
        if rec['status'] == 'REDUCES_TO_KNOWN_158_SEALED_CATALOGUE':
            for mp in rec['known_subcores']:
                target = sealed[mp['index']] if mp['family'] == 'old' else first['records'][mp['index']]['row']
                B.check_map(induced(rec['row'], mp['nodes']), target, mp['to_target'])
    return dict(no_short_literal_unions=5440, minimal_union_classes=1180,
                principal_survivors=396, H4_paid=84, affine_excluded=2,
                original_sealed=112, nonsingular_sealed_rank_excluded=44,
                first_extensions=458, singular_first_extensions=147,
                singular_subcore_reductions=86, singular_literal_frontier=61,
                singular_marked_frontier=51, singular_closure_states=17,
                nonsingular_literal_roots=311, dotted_seed_classes=24), open_nonsingular, set(dotted_ids)


def dotted_branch(B, expected):
    parameters = B.load('fresh_sqrt5_dotted_parameters_20261015')
    B.review('fresh_sqrt5_dotted_parameters_20261015_review')
    assert {r['source_class_index'] for r in parameters['records']} == expected
    assert sum(parameters['counts'].values()) == 532
    assert parameters['counts']['ACCEPT_SINGULAR'] == 18 and parameters['counts']['ACCEPT_NONSINGULAR'] == 14
    queue = B.load(DOTTED)
    B.review(DOTTED + '_review')
    assert queue['status'] == 'COMPLETE_DOTTED_ARITHMETIC_ATTACHMENT_QUEUE_THROUGH_DIMENSION29'
    assert queue['processed'] == queue['marked_states'] == len(queue['records']) == 438
    actual = []
    for rec in queue['records']:
        if rec['result']['singular']:
            actual.append(rec['index'])
            B.payment(rec['row'], ['dotted', rec['index'], rec['ratio']])
    assert actual == queue['weight_obligations'] and len(actual) == 264
    return dict(parameter_classes=24, arithmetic_parameter_trials=532, exact_seeds=32,
                queue_states=438, singular_zero_touch_obligations=264)


def finite_branch(B, roots):
    queue = B.load(FINITE)
    assert queue['status'] == 'COMPLETE_ACCELERATED_FINITE_ATTACHMENT_QUEUE_DIMENSION29'
    assert queue['complete_queue'] and not queue['counts'].get('PENDING', 0)
    review = B.review(FINITE + '_review')
    assert review['status'] == 'PASS_INDEPENDENT_COMPLETE_FINITE_ATTACHMENT_QUEUE_DIMENSION29'
    assert review['complete_queue'] and review['pending'] == 0
    assert review['processed'] == review['marked_states'] == queue['marked_states']
    assert review['counts'] == queue['counts'] and review['trial_counts'] == queue['trial_counts']
    assert review['weight_obligations'] == len(queue['weight_obligations'])
    predecessor = B.review(FINITE_PREDECESSOR + '_review')
    assert not predecessor['complete_queue'] and predecessor['processed'] == 1433
    assert queue['state']['inherited_processed'] == 1433
    assert queue['state']['inherited_marked_states'] == 1649
    bound = review.get('sources', {}).get(f'results/{FINITE}.json')
    assert bound == sha256((ROOT / f'results/{FINITE}.json').read_bytes()).hexdigest()
    assert queue['state']['next_core'] == queue['marked_states'] == len(queue['records'])
    assert len(queue['root_maps']) == len(roots) == 311
    for row, mp in zip(roots, queue['root_maps']):
        # Schema is checked below against the producer's literal root map.
        target = mp.get('target', mp.get('index'))
        mapping = mp.get('map', mp.get('to_representative'))
        assert target is not None and mapping is not None, mp
        B.check_map(row, queue['records'][target]['row'], mapping)
    actual = []
    for rec in queue['records']:
        if rec['result']['singular']:
            actual.append(rec['index'])
            B.payment(rec['row'], ['finite_queue', rec['index']])
    assert actual == queue['weight_obligations']
    return dict(literal_roots=311, queue_states=len(queue['records']),
                singular_zero_touch_obligations=len(actual), counts=queue['counts'])


def run(stage=False):
    if not stage:
        # Refuse to create a complete certificate while either required file
        # is absent, active, or not bound by its separate review.
        q = json.loads((ROOT / f'results/{FINITE}.json').read_text())
        assert q['complete_queue'], 'Finite queue is not complete; use --stage only.'
        assert (ROOT / f'results/{FINITE}_review.json').exists(), 'Finite queue review is missing.'
    B = Binder()
    totals = triangles_and_vertices(B)
    old, roots, dotted_ids = old_branches(B)
    dotted = dotted_branch(B, dotted_ids)
    finite = None if stage else finite_branch(B, roots)
    sources = B.bind_manifests()
    return dict(status='STAGED_FROZEN_BRANCHES_ONLY_NO_DIMENSION_CLAIM' if stage else
                'PASS_COMPLETE_QSQRT5_FACE_COVERAGE_EXCLUDES_COMPACT_DIMENSIONS28_AND29',
                arithmetic_ground_field='Q(sqrt5)',
                compact_dimension_upper_bound=None if stage else 27,
                excluded_dimensions=[] if stage else [28, 29],
                complete_coverage=not stage,
                sources=sources, recursive_manifests=len(B.manifests), source_references=B.references,
                baseline_and_vertex_bounds=totals, old_reductions=old, dotted=dotted, finite=finite,
                literal_map_checks=B.maps, zero_touch_rows=len(B.payments),
                zero_touch_cases=sum(len(r['cases']) for r in B.payments), payments=B.payments,
                scope='Source-bound proof composition with direct exact graph payments, not a new arithmetic enumeration. Complete mode requires frozen independently reviewed finite and dotted queues. The conclusion concerns compact arithmetic hyperbolic Coxeter polytopes with ground field Q(sqrt5); no other field, noncompact class, or construction claim. The known compact upper bound29 reduces the possible improvement to dimensions28 and29.')


if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('--verify', action='store_true');p.add_argument('--stage', action='store_true');a=p.parse_args()
    assert not (a.verify and a.stage)
    result = json.loads(json.dumps(run(a.stage)))
    path = ROOT / f'results/{TAG}.json'
    if a.stage:
        print(result['status'], result['zero_touch_rows'], result['zero_touch_cases'])
    elif a.verify:
        assert result == json.loads(path.read_text());print('PASS saved coverage replay')
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'])
