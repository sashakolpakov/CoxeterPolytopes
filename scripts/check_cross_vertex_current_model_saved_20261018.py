#!/usr/bin/env python3
"""Independent exact check of new exports and saved model certificate fields.

No producer import, matrix reconstruction or optimization. The full historical
matrix inequalities remain the scope of the separate full-model replay.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEM = 'fresh_cross_vertex_current_full_model_20261018'


def read(stem):
    return json.loads((ROOT/'results'/f'{stem}.json').read_text())


def rows(items):
    result = Counter()
    for key, value in items:
        result[key] += Fraction(value)
    return {key:value for key, value in result.items() if value}


def main():
    saved = read(STEM)
    assert saved['status'].startswith('PASS_EXACT') and saved['solver_invocations'] == 1
    p = saved['prepared']; counts = p['counts']
    assert counts['exported_columns'] == 213 and counts['distinct_added_columns'] == 209
    assert counts['several_vertex_columns'] == 14 and counts['variables'] == 439619
    assert counts['inherited_variables'] == 439410 and counts['duplicate_exports'] == 4
    for path, expected in p['sources'].items():
        assert sha256((ROOT/path).read_bytes()).hexdigest() == expected, path
    parent = read('fresh_weighted_raw110_full_model_20261017')['prepared']
    assert parent['model_sha256'] == p['inherited_model_sha256']
    assert parent['provenance_to_face'] == p['provenance_to_face']
    inputs = {stem:read(stem) for stem, _ in p['input_column_counts']}
    expected_origins = {(stem,index) for stem, n in p['input_column_counts'] for index in range(n)}
    assert len(expected_origins) == 213 and len(inputs) == 12
    for stem, n in p['input_column_counts']:
        assert len(inputs[stem]['columns']) == n
    seen, keys, columns = set(), set(), []
    for offset, column in enumerate(p['new_columns']):
        j = column['variable']; assert j == 439410+offset
        face = rows(column['integer_face_coefficients'])
        cost = rows(column['integer_vertex_capacities'])
        assert len(face) == len(column['integer_face_coefficients'])
        assert len(cost) == len(column['integer_vertex_capacities'])
        assert all(0 <= f < 82806 and a > 0 and a.denominator == 1 for f, a in face.items())
        assert all(0 <= v < 6558 and a > 0 and a.denominator == 1 for v, a in cost.items())
        assert gcd(*(int(x) for x in list(face.values())+list(cost.values()))) == 1
        key = tuple(sorted(face.items())), tuple(sorted(cost.items()))
        assert key not in keys; keys.add(key)
        for origin in column['origins']:
            pair = origin['package'], origin['index']; assert pair not in seen; seen.add(pair)
            raw = inputs[pair[0]]['columns'][pair[1]]
            old_face = rows(raw['positive_face_coefficients'])
            if 'vertex_capacities' in raw:
                old_cost = rows(raw['vertex_capacities'])
            elif 'source_vertex' in raw:
                old_cost = rows([[raw['source_vertex'],1], [raw['target_vertex'],raw['target_vertex_cost']]])
            else:
                old_cost = rows([[raw['vertex'],raw['vertex_cost']]])
            assert set(face) == set(old_face) and set(cost) == set(old_cost)
            anchor = min(face); ratio = face[anchor]/old_face[anchor]
            assert ratio > 0
            assert all(face[f] == ratio*a for f, a in old_face.items())
            assert all(cost[v] == ratio*a for v, a in old_cost.items())
            scale = origin['scaling']
            assert ratio == Fraction(scale['denominator'],scale['integer_gcd'])
        columns.append((j,face,cost))
    assert seen == expected_origins and len(columns) == 209
    assert sum(len(cost)>1 for _,_,cost in columns) == 14
    primal = next(c for c in saved['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_PRIMAL')
    dual = next(c for c in saved['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_DUAL')
    alpha, mu = rows(dual['face_multipliers']), rows(dual['vertex_multipliers'])
    assert len(alpha) == len(dual['face_multipliers']) and len(mu) == len(dual['vertex_multipliers'])
    assert all(0 <= f < 82806 and a > 0 for f,a in alpha.items())
    assert all(0 <= v < 6558 and a > 0 for v,a in mu.items())
    for _, face, cost in columns:
        lhs = sum((a*alpha.get(f,0) for f,a in face.items()),Fraction())
        rhs = sum((a*mu.get(v,0) for v,a in cost.items()),Fraction())
        assert lhs <= rhs
    x = rows(primal['rounded_variable_numerators']); den = primal['integer_denominator']
    assert len(x) == len(primal['rounded_variable_numerators']) and den > 0
    assert all(0 <= j < 439619 and a >= 0 and a.denominator == 1 for j,a in x.items())
    assert Fraction(primal['repair_delta']) >= 0
    assert len(set(primal['repair_variable_indices'])) == len(primal['repair_variable_indices'])
    assert all(0 <= j < 439410 for j in primal['repair_variable_indices'])
    budgets = rows(primal['vertex_budgets'])
    assert len(budgets) == len(primal['vertex_budgets'])
    assert all(0 <= v < 6558 and a > 0 for v,a in budgets.items())
    assert max(budgets.values()) == Fraction(primal['maximum'])
    added = Counter()
    for j, _, cost in columns:
        for v,a in cost.items():
            added[v] += a*x.get(j,0)/den
    assert all(a <= budgets.get(v,0) for v,a in added.items())
    normalizer = sum(mu.values(),Fraction())
    assert normalizer == Fraction(dual['normalizer']) and normalizer > 0
    lower = Fraction(dual['objective_numerator'])/normalizer
    upper = Fraction(primal['maximum'])
    assert lower == Fraction(dual['bound']) and lower <= upper
    assert dual['excludes_target12'] == (lower > 12)
    assert primal['meets_target12'] == (upper <= 12)
    print('PASS_SAVED_CROSS_VERTEX_CURRENT_MODEL_EXPORTS',len(p['sources']),'bindings',
          len(seen),'exports',len(columns),'columns','14 several-vertex',
          'bounds',float(lower),float(upper))


if __name__ == '__main__':
    main()
