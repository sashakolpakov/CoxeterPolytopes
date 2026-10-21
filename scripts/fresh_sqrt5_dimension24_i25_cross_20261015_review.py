#!/usr/bin/env python3
"""Independent literal-graph and rational review of the rank24 hub rule.

No optimizer or producer import. The frozen independent graph engine provides
all-root canonical keys and separately constructed finite Dynkin diagrams.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import fresh_sqrt5_all_a1_cross_corrected_20261015_review as V

R = V.R
ROOT = Path(__file__).resolve().parents[1]
PRIMARY = 'fresh_sqrt5_dimension24_i25_cross_20261015'
TAG = PRIMARY + '_review'
BOUND = 24


def read(tag):
    return json.loads((ROOT / f'results/{tag}.json').read_text())


def corners(row):
    matrix = R.label_matrix(row)
    edges = [(i, j, matrix[i][j]) for i, j in combinations(range(row['nodes']), 2)]
    out = []
    for left, right in product((0, 2), (1, 3)):
        parts = R.components(R.graph(row['white'] + [left, right], edges))
        lc = next(g for g in parts if left in g)
        rc = next(g for g in parts if right in g)
        connected = a1cross = hubcross = None
        baseline = F(0)
        if right in lc:
            connected = R.tree_key(lc, [left, right])
            baseline = R.baseline(lc, left, right)
        else:
            if len(lc) == 1 and len(rc) > 1:
                a1cross = R.tree_key(rc, [right])
            elif len(rc) == 1 and len(lc) > 1:
                a1cross = R.tree_key(lc, [left])
            if R.type_of(lc) == 'I2(5)' and R.type_of(rc) not in ('A1', 'I2(5)'):
                hubcross = R.tree_key(rc, [right])
            elif R.type_of(rc) == 'I2(5)' and R.type_of(lc) not in ('A1', 'I2(5)'):
                hubcross = R.tree_key(lc, [left])
        forest = [R.type_of(g) for g in parts]
        out.append(dict(baseline=baseline, connected=connected, a1cross=a1cross,
                        hubcross=hubcross, a=forest.count('A1'), b=forest.count('I2(5)')))
    return out


class Rule:
    def __init__(self, primary):
        assert primary['dimension'] == BOUND and primary['split_existing_bonus_by_hub_presence']
        assert primary['certificate']['kind'] == 'exact_feasible_rule'
        self.x = list(map(F, primary['certificate']['exact_primal']))
        assert len(self.x) == primary['variables'] + 1 and min(self.x) >= 0
        assert all((v * 10000).denominator == 1 for v in self.x[:-1])
        catalogue = read('fresh_sqrt5_complete_face_weights_20261015_n25')
        labels = {}
        for kind, records in [('connected', catalogue['connected_orbits']),
                              ('a1cross', catalogue['cross_orbits']),
                              ('hubcross', primary['new_cross_orbits'])]:
            for rec in records:
                graph = R.graph(rec['nodes'], [(i, j, m) for (i, j), m in
                                zip(combinations(rec['nodes'], 2), rec['labels'])])
                marks = rec['marked'] if kind == 'connected' else [rec['marked_root']]
                assert R.type_of(graph) == R.alias(rec['type']) and len(graph) == rec['rank']
                if kind == 'connected':
                    assert R.baseline(graph, *marks) == 0
                elif kind == 'a1cross':
                    assert len(graph) > 1
                else:
                    assert R.type_of(graph) not in ('A1', 'I2(5)')
                assert rec['key'] not in labels
                labels[rec['key']] = kind, R.tree_key(graph, marks), len(graph)
        self.index = {}
        self.orbits = {k: {} for k in ('connected', 'a1cross', 'hubcross')}
        for i, variable in enumerate(primary['variable_order']):
            label = variable['orbit']
            present = None
            for suffix, flag in [('@I25absent', False), ('@I25present', True)]:
                if label.endswith(suffix):
                    label, present = label[:-len(suffix)], flag
                    break
            kind, key, rank = labels[label]
            assert (present is None) == (kind == 'hubcross')
            a = variable['A1_count']
            assert 0 <= a <= BOUND - rank
            item = kind, key, a, present
            assert item not in self.index
            self.index[item] = i
            self.orbits[kind][key] = dict(rank=rank, label=label)
        wanted = set()
        for kind, orbits in self.orbits.items():
            for key, rec in orbits.items():
                for a in range(int(kind == 'a1cross'), BOUND - rec['rank'] + 1 - 2 * (kind == 'hubcross')):
                    for present in ([None] if kind == 'hubcross' else [False, True]):
                        wanted.add((kind, key, a, present))
        assert wanted == set(self.index) and len(wanted) == primary['variables']

    def bonus(self, kind, key, a, present):
        index = self.index.get((kind, key, a, None if kind == 'hubcross' else present))
        return F(0) if index is None else self.x[index]

    def payment(self, row):
        remaining = BOUND - row['nodes'] + 2
        assert remaining >= 0
        cs = corners(row)
        out = []
        for external_a in range(remaining + 1):
            for external_b in range((remaining - external_a) // 2 + 1):
                if remaining - external_a - 2 * external_b == 1:
                    continue
                values = []
                for c in cs:
                    a, b = c['a'] + external_a, c['b'] + external_b
                    value = c['baseline'] + self.bonus('connected', c['connected'], a, b > 0)
                    if c['a1cross'] in self.orbits['a1cross']:
                        assert a > 0
                        value += self.bonus('a1cross', c['a1cross'], a, b > 0) / a
                    if c['hubcross'] in self.orbits['hubcross']:
                        assert b > 0
                        value += self.bonus('hubcross', c['hubcross'], a, b > 0) / (2 * b)
                    values.append(value)
                assert min(values) >= 0 and sum(values) >= 1, (row['nodes'], external_a, external_b, values)
                out.append(dict(external_A1=external_a, external_I25=external_b,
                                corner_weights=list(map(str, values)), payment=str(sum(values))))
        return out


def component_bounds(rule):
    models = {name: graph for name, graph in R.MODELS.items() if len(graph) <= BOUND}
    profiles = {}
    observed = {kind: Counter() for kind in rule.orbits}
    for name, graph in models.items():
        base = sum(R.baseline(graph, u, v) for u, v in combinations(graph, 2))
        connected = Counter(R.tree_key(graph, [u, v]) for u, v in combinations(graph, 2))
        rooted = Counter(R.tree_key(graph, [u]) for u in graph)
        profiles[name] = base, connected, rooted
        for key, count in connected.items():
            if key in rule.orbits['connected']:
                observed['connected'][key] += count
        for kind in ('a1cross', 'hubcross'):
            for key, count in rooted.items():
                if key in rule.orbits[kind]:
                    observed[kind][key] += count
    assert all(set(observed[kind]) == set(rule.orbits[kind]) for kind in observed)
    fronts, totals = [], {}
    for a in range(BOUND + 1):
        for b in range((BOUND - a) // 2 + 1):
            remaining = BOUND - a - 2 * b
            if remaining == 1:
                continue
            values = {}
            for name, (base, connected, rooted) in profiles.items():
                value = base + sum(count * rule.bonus('connected', key, a, b > 0)
                                   for key, count in connected.items())
                if name != 'A1' and a:
                    for key, count in rooted.items():
                        bonus = rule.bonus('a1cross', key, a, b > 0)
                        assert sum((bonus / a for _ in range(a)), F(0)) == bonus
                        value += count * bonus
                if name not in ('A1', 'I2(5)') and b:
                    for key, count in rooted.items():
                        bonus = rule.bonus('hubcross', key, a, True)
                        assert sum((bonus / (2 * b) for _ in range(2 * b)), F(0)) == bonus
                        value += count * bonus
                values[name] = value
            assert values['A1'] == 0
            dp = {0: F(0)}
            # Reverse-type unbounded coin DP, independent of the producer's
            # rank-first recurrence. Exclude exactly the two fixed hub types.
            for name in sorted(models, reverse=True):
                if name in ('A1', 'I2(5)'):
                    continue
                rank = len(models[name])
                for n in range(rank, remaining + 1):
                    if n - rank in dp:
                        v = dp[n - rank] + values[name]
                        dp[n] = max(dp.get(n, F(-1)), v)
            assert remaining in dp
            maximum = dp[remaining] + b * values['I2(5)']
            totals[a, b] = values
            fronts.append(dict(context=[a, b], maximum=str(maximum),
                               residual_rank_maxima=[str(dp[n]) if n in dp else None
                                                     for n in range(remaining + 1)]))
    return models, totals, fronts, observed


def run():
    primary = read(PRIMARY + '_hub_context')
    rule = Rule(primary)
    old = read('fresh_sqrt5_parameter_weight_stress_20261015')
    queue = read('fresh_sqrt5_nonsingular_rank29_fast_20261015')
    assert queue['complete_queue'] and queue['state']['next_core'] == queue['marked_states'] == len(queue['records'])
    assert not queue['counts'].get('PENDING', 0)
    singular = [r for r in queue['records'] if r['result']['singular']]
    assert [r['index'] for r in singular] == queue['weight_obligations']
    faces = [('prior_catalogue', i, r['row']) for i, r in enumerate(old['faces']) if r['row']['nodes'] <= 26]
    old_count = len(faces)
    faces += [('finite_queue', r['index'], r['row']) for r in singular if r['row']['nodes'] <= 26]
    assert len(faces) == primary['faces']
    expected = {kind: set() for kind in rule.orbits}
    for rec in old['orbits']:
        if rec['rank'] <= BOUND:
            graph = R.graph(rec['nodes'], [(i, j, m) for (i, j), m in
                            zip(combinations(rec['nodes'], 2), rec['labels'])])
            expected['connected'].add(R.tree_key(graph, rec['marked']))
    payments, case_count, minimum = [], 0, None
    for family, index, row in faces:
        for c in corners(row):
            for kind in expected:
                if c[kind] is not None and (kind != 'connected' or c['baseline'] == 0):
                    expected[kind].add(c[kind])
        cases = rule.payment(row)
        for case in cases:
            slack = F(case['payment']) - 1
            minimum = slack if minimum is None else min(minimum, slack)
        case_count += len(cases)
        payments.append(dict(family=family, index=index, nodes=row['nodes'], cases=cases))
    assert all(expected[kind] == set(rule.orbits[kind]) for kind in expected)
    assert case_count == primary['face_cases']
    raw = V.X.inputs()
    assert len(raw) == len(old['literal_maps']) == 2397
    literals = 0
    for (family, index, row, ratio), mp in zip(raw, old['literal_maps']):
        assert (family, index) == (mp['family'], mp['index'])
        target = old['faces'][mp['representative']]
        assert ratio == target['dotted_ratio']
        R.check_bijection(R.label_matrix(row), R.label_matrix(target['row']), mp['to_representative'])
        assert (row['nodes'] <= 26) == (target['row']['nodes'] <= 26)
        literals += row['nodes'] <= 26
    models, totals, fronts, multiplicities = component_bounds(rule)
    saved = primary['verification']['fronts']
    assert len(saved) == len(fronts) == primary['verification']['contexts']
    for ours, theirs in zip(fronts, saved):
        assert ours['context'] == theirs['context'] and F(ours['maximum']) == F(theirs['maximum'])
        a, b = ours['context']
        forest = [R.alias(name) for name in theirs['types']]
        assert sum(len(models[name]) for name in forest) == BOUND
        assert forest.count('A1') == a and forest.count('I2(5)') == b
        assert sum(totals[a, b][name] for name in forest) == F(ours['maximum'])
    maximum = max(F(row['maximum']) for row in fronts)
    assert maximum == rule.x[-1] == F(primary['certificate']['maximum']) == F(primary['verification']['maximum']) <= 54
    triangle = read('fresh_sqrt5_triangle_weights_20261015_review')
    assert len(triangle['markings']) == 76
    assert all(sum(c['h4_only_weight'] for c in row['corners']) >= 2 for row in triangle['markings'])
    paths = [f'scripts/{TAG}.py', f'scripts/{V.TAG}.py', f'results/{V.TAG}.json',
             f'scripts/{PRIMARY}.py', f'results/{PRIMARY}_hub_context.json',
             'results/fresh_sqrt5_complete_face_weights_20261015_n25.json',
             'results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json',
             'results/fresh_sqrt5_parameter_weight_stress_20261015.json',
             'results/fresh_sqrt5_parameter_weight_stress_20261015_cross_review.json',
             'results/fresh_sqrt5_triangle_weights_20261015_review.json']
    for path in paths:
        if path.endswith('.json'):
            for dep, digest in json.loads((ROOT / path).read_text()).get('sources', {}).items():
                assert sha256((ROOT / dep).read_bytes()).hexdigest() == digest, (path, dep)
    return dict(status='PASS_INDEPENDENT_EXACT_DIMENSION24_HUB_CONTEXT_FACE_AND_FOREST_REVIEW',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                dimension=BOUND, original_literal_inputs=len(raw), dimension_feasible_literals=literals,
                old_face_rows=old_count, all_finite_singular_stops=len(singular),
                dimension_feasible_finite_stops=len(faces)-old_count, face_rows=len(faces),
                face_cases=case_count, minimum_face_slack=str(minimum), variables=primary['variables'],
                orbit_counts={kind: len(orbits) for kind, orbits in rule.orbits.items()},
                finite_component_types=len(models), exact_A1_I25_contexts=len(fronts),
                vertex_maximum=str(maximum), strict_face_threshold='54', threshold_margin=str(54-maximum),
                forest_maxima=fronts, face_payments=payments,
                multiplicities={kind: [dict(label=rule.orbits[kind][k]['label'], multiplicity=v)
                                      for k, v in sorted(counts.items())]
                                for kind, counts in multiplicities.items()},
                proof='Each destination root has exactly a A1/root pairs and exactly 2b I2(5)/root pairs. The respective denominators cancel in vertex sums. Hub destinations exclude A1 and I2(5), so the two cross families are disjoint. Fixing exact(a,b) and b>0 leaves additive components; reverse-type rational DP covers every residual forest. Faces use the same external(a,b) at all four corners and retain the fixed short-path/H4 baseline.',
                scope='Independent weight verification only, without producer or optimizer imports. Complete geometric face coverage and field/dimension composition require separate certificates. No optimality or priority claim.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser();ap.add_argument('--verify', action='store_true');args = ap.parse_args()
    result = run();path = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['vertex_maximum'], result['face_cases'])
