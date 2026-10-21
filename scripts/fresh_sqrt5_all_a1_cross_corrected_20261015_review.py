#!/usr/bin/env python3
"""Separate literal-graph/Fraction verification of normalized A1 cross weights.

No optimizer or producer is imported. The frozen independent 224-graph
engine supplies all-root canonicalization and independently built Dynkin
models. This checks a feasible rule, not an optimality certificate.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_enlarged_a1_weights_20261015_review as R
import fresh_sqrt5_parameter_weight_stress_20261015_cross_review as X

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = 'fresh_sqrt5_all_a1_cross_corrected_20261015'
TAG = PRIMARY + '_review'
BOUND = 27


def read(tag):
    return json.loads((ROOT / f'results/{tag}.json').read_text())


def corners(row):
    matrix = R.label_matrix(row)
    edges = [(i, j, matrix[i][j]) for i, j in combinations(range(row['nodes']), 2)]
    result = []
    for a, b in product((0, 2), (1, 3)):
        parts = R.components(R.graph(row['white'] + [a, b], edges))
        left = next(g for g in parts if a in g)
        right = next(g for g in parts if b in g)
        key = cross = None
        baseline = F(0)
        if b in left:
            key = R.tree_key(left, [a, b])
            baseline = R.baseline(left, a, b)
        elif len(left) == 1 and len(right) > 1:
            cross = R.tree_key(right, [b])
        elif len(right) == 1 and len(left) > 1:
            cross = R.tree_key(left, [a])
        result.append(dict(marked=[a, b], forest=sorted(R.type_of(g) for g in parts),
                           A1=sum(len(g) == 1 for g in parts), baseline=baseline,
                           connected=key, cross=cross))
    return result


class Rule:
    def __init__(self, primary):
        self.primary = primary
        self.bound = primary['dimension']
        assert 1 <= self.bound <= 29
        self.x = list(map(F, primary['solution']['exact_primal']))
        assert len(self.x) == primary['variables'] + 1 and min(self.x) >= 0
        self.connected = {}
        self.cross = {}
        labels = {}
        for kind, records, lookup in [('connected', primary['connected_orbits'], self.connected),
                                      ('cross', primary['cross_orbits'], self.cross)]:
            for rec in records:
                g = R.graph(rec['nodes'], [(i, j, m) for (i, j), m in
                            zip(combinations(rec['nodes'], 2), rec['labels'])])
                marks = rec['marked'] if kind == 'connected' else [rec['marked_root']]
                key = R.tree_key(g, marks)
                assert key not in lookup and R.type_of(g) == R.alias(rec['type'])
                assert len(g) == rec['rank'] <= self.bound
                if kind == 'connected':
                    assert len(marks) == 2 and R.baseline(g, *marks) == 0
                else:
                    assert len(g) > 1
                lookup[key] = dict(label=rec['key'], rank=len(g), graph=g)
                assert rec['key'] not in labels
                labels[rec['key']] = kind, key
        self.index = {}
        for i, rec in enumerate(primary['variable_order']):
            kind, key = labels[rec['orbit']]
            count = rec['exact_A1_count']
            item = kind, key, count
            rank = (self.connected if kind == 'connected' else self.cross)[key]['rank']
            assert int(kind == 'cross') <= count <= self.bound - rank and item not in self.index
            self.index[item] = i
        wanted = {('connected', k, a) for k, v in self.connected.items()
                  for a in range(self.bound - v['rank'] + 1)}
        wanted |= {('cross', k, a) for k, v in self.cross.items()
                   for a in range(1, self.bound - v['rank'] + 1)}
        assert set(self.index) == wanted and len(self.index) == primary['variables']

    def bonus(self, kind, key, count):
        i = self.index.get((kind, key, count))
        return self.x[i] if i is not None else F(0)

    def payment(self, row):
        residual = self.bound - (row['nodes'] - 2)
        assert residual >= 0
        data = corners(row)
        answer = []
        for external in range(residual + 1):
            if residual - external == 1:
                continue
            terms = []
            coefficients = Counter()
            baseline = sum(c['baseline'] for c in data)
            for c in data:
                count = c['A1'] + external
                value = c['baseline']
                if c['connected'] in self.connected:
                    assert c['baseline'] == 0
                    i = self.index['connected', c['connected'], count]
                    coefficients[i] += 1
                    value += self.x[i]
                if c['cross'] in self.cross:
                    assert count > 0
                    i = self.index['cross', c['cross'], count]
                    coefficients[i] += F(1, count)
                    value += self.x[i] / count
                terms.append(value)
            total = sum(terms)
            assert total >= 1, (row['nodes'], external, total)
            assert total == baseline + sum(v * self.x[i] for i, v in coefficients.items())
            answer.append(dict(external_A1=external, baseline=str(baseline),
                               coefficients=[[i, str(v)] for i, v in sorted(coefficients.items())],
                               corner_weights=list(map(str, terms)), payment=str(total)))
        return answer


def component_bounds(rule):
    bound = rule.bound
    models = {name: g for name, g in R.MODELS.items() if len(g) <= bound}
    totals = {}
    cross_multiplicity = Counter()
    connected_multiplicity = Counter()
    for name, g in models.items():
        base = sum(R.baseline(g, a, b) for a, b in combinations(g, 2))
        values = [base for _ in range(bound + 1)]
        for a, b in combinations(g, 2):
            key = R.tree_key(g, [a, b])
            if key not in rule.connected:
                continue
            connected_multiplicity[key] += 1
            for count in range(bound - len(g) + 1):
                values[count] += rule.bonus('connected', key, count)
        if name != 'A1':
            for node in g:
                key = R.tree_key(g, [node])
                if key not in rule.cross:
                    continue
                cross_multiplicity[key] += 1
                for count in range(1, bound - len(g) + 1):
                    # There are count distinct A1/root pairs. Verify the
                    # normalization explicitly before using additivity.
                    u = rule.bonus('cross', key, count)
                    assert sum((u / count for _ in range(count)), F(0)) == u
                    values[count] += u
        totals[name] = values
    assert set(cross_multiplicity) == set(rule.cross)
    assert set(connected_multiplicity) == set(rule.connected)
    assert all(v == 0 for v in totals['A1'])  # No A1--A1 bonus.
    # Type-first unbounded coin recurrence, distinct from rank-first producer.
    fronts = []
    for count in range(bound + 1):
        remaining = bound - count
        if remaining == 1:
            continue
        dp = {0: F(0)}
        for name in sorted(models, reverse=True):
            if name == 'A1':
                continue
            rank = len(models[name])
            for n in range(rank, remaining + 1):
                if n - rank in dp:
                    v = dp[n-rank] + totals[name][count]
                    if n not in dp or v > dp[n]:
                        dp[n] = v
        assert remaining in dp
        fronts.append(dict(context=count, maximum=str(dp[remaining]),
                           all_residual_rank_maxima=[str(dp[n]) if n in dp else None
                                                    for n in range(remaining + 1)]))
    return models, totals, fronts, connected_multiplicity, cross_multiplicity


def run():
    primary = read(PRIMARY)
    assert primary['dimension'] == BOUND and primary['meets_target']
    rule = Rule(primary)
    assert primary['solution']['bonus_grid'] == 1000
    assert all((v * 1000).denominator == 1 for v in rule.x[:-1])
    frozen = read('fresh_sqrt5_parameter_weight_stress_20261015')
    frozen_review = read('fresh_sqrt5_parameter_weight_stress_20261015_cross_review')
    assert frozen_review['status'] == 'PASS_SEPARATE_UNIQUE_A1_E7_GRAPH_DUAL_DP_REPLAY'
    faces = [(i, rec) for i, rec in enumerate(frozen['faces'])
             if rec['row']['nodes'] - 2 <= BOUND]
    assert len(faces) == primary['faces']
    connected_expected = {r['key'] for r in frozen['orbits'] if r['rank'] <= BOUND}
    assert connected_expected == {r['key'] for r in primary['connected_orbits']}
    seen_cross = set()
    payments = []
    minimum = None
    case_count = 0
    for old_index, rec in faces:
        for c in corners(rec['row']):
            if c['cross'] is not None:
                seen_cross.add(c['cross'])
        cases = rule.payment(rec['row'])
        for case in cases:
            slack = F(case['payment']) - 1
            minimum = slack if minimum is None else min(minimum, slack)
        case_count += len(cases)
        payments.append(dict(frozen_face_index=old_index, cases=cases))
    assert seen_cross == set(rule.cross) and case_count == primary['face_cases']
    assert minimum == F(primary['exact_certificate']['minimum_face_slack'])
    # Check all literal provenance, including the dimension filter, rather
    # than treating 1343 distinct representatives as unbound input data.
    raw = X.inputs()
    assert len(raw) == len(frozen['literal_maps']) == 2397
    literal_count = 0
    used = Counter()
    for (family, index, row, ratio), mp in zip(raw, frozen['literal_maps']):
        assert (family, index) == (mp['family'], mp['index'])
        target = frozen['faces'][mp['representative']]
        assert ratio == target['dotted_ratio']
        R.check_bijection(R.label_matrix(row), R.label_matrix(target['row']), mp['to_representative'])
        assert (row['nodes'] <= BOUND + 2) == (target['row']['nodes'] <= BOUND + 2)
        if row['nodes'] <= BOUND + 2:
            literal_count += 1
            used[mp['representative']] += 1
    assert set(used) == {i for i, _ in faces}
    models, totals, fronts, connected_mult, cross_mult = component_bounds(rule)
    saved = primary['exact_certificate']['fronts']
    assert len(saved) == len(fronts)
    for ours, rec in zip(fronts, saved):
        assert ours['context'] == rec['context'] and F(ours['maximum']) == F(rec['maximum'])
        forest = [R.alias(t) for t in rec['maximizing_forest']]
        count = ours['context']
        assert sum(len(models[t]) for t in forest) == BOUND and forest.count('A1') == count
        assert sum(totals[t][count] for t in forest) == F(ours['maximum'])
    maximum = max(F(r['maximum']) for r in fronts)
    assert maximum == rule.x[-1] == F(primary['solution']['feasible_maximum'])
    assert maximum == F(primary['exact_certificate']['maximum']) <= F(143, 2)
    predecessor = read('fresh_sqrt5_all_a1_cross_20261015')
    assert predecessor['variable_order'] == primary['variable_order']
    assert predecessor['solution']['exact_primal'] == primary['solution']['exact_primal']
    corrections = []
    for ours, old in zip(fronts, predecessor['exact_certificate']['fronts']):
        assert ours['context'] == old['context']
        if F(ours['maximum']) != F(old['maximum']):
            corrections.append(dict(exact_A1_count=ours['context'], previous=str(old['maximum']),
                                    corrected=ours['maximum']))
    assert corrections == [dict(exact_A1_count=20, previous='6687/200', corrected='12359/250')]
    triangle = read('fresh_sqrt5_triangle_weights_20261015_review')
    assert len(triangle['markings']) == 76
    for rec in triangle['markings']:
        assert sum(c['h4_only_weight'] for c in rec['corners']) >= 2
    paths = [f'scripts/{TAG}.py', f'scripts/{R.TAG}.py', f'scripts/{X.TAG}.py']
    for tag in (PRIMARY, 'fresh_sqrt5_parameter_weight_stress_20261015',
                'fresh_sqrt5_parameter_weight_stress_20261015_cross_review',
                'fresh_sqrt5_triangle_weights_20261015_review',
                'fresh_sqrt5_all_a1_cross_20261015'):
        paths.append(f'results/{tag}.json')
    paths.append(f'scripts/{PRIMARY}.py')
    for p in paths:
        if p.endswith('.json'):
            for q, digest in json.loads((ROOT / p).read_text()).get('sources', {}).items():
                assert sha256((ROOT / q).read_bytes()).hexdigest() == digest, (p, q)
    return dict(status='PASS_SEPARATE_EXACT_A1_NORMALIZED_CROSS_FEASIBILITY_REVIEW',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                dimension=BOUND, original_literal_inputs=2397, dimension_feasible_literals=literal_count,
                distinct_faces=len(faces), face_cases=case_count, variables=primary['variables'],
                connected_orbits=len(rule.connected), rooted_cross_orbits=len(rule.cross),
                minimum_face_slack=str(minimum), feasible_vertex_maximum=str(maximum),
                target='143/2', target_margin=str(F(143, 2) - maximum),
                finite_component_types=len(models), forest_maxima=fronts, face_payments=payments,
                predecessor_frontier_corrections=corrections,
                cross_root_multiplicities=[dict(label=rule.cross[k]['label'], multiplicity=v)
                                           for k, v in sorted(cross_mult.items())],
                connected_pair_multiplicities=[dict(label=rule.connected[k]['label'], multiplicity=v)
                                               for k, v in sorted(connected_mult.items())],
                proof='At exact A1 count a>0, each non-A1 root makes exactly a distinct A1/root pairs, each carrying u_a/a, so its total is u_a. A1/A1 pairs carry no bonus. Fixing a leaves an additive non-A1 forest partition of rank27-a. Reverse-type exact DP covers every such partition. The common external count in each face is evaluated without capping.',
                scope='Separate graph and rational feasibility implementation, no producer/optimizer arithmetic imports and no optimality claim. Frozen literal catalogue and prior sealing/classification remain dependencies. Every new finite-queue zero-touch obligation must be checked before any dimension26 conclusion; this partial certificate establishes no dimension bound.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser();ap.add_argument('--verify', action='store_true');args = ap.parse_args()
    out = json.loads(json.dumps(run()));path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text());print('PASS saved normalized-cross review')
    else:
        path.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], out['feasible_vertex_maximum'], out['face_cases'])
