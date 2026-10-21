#!/usr/bin/env python3
"""Exact, partial dimension-27 stress test of the frozen unique-A1/E7 rule.

Only dimension-feasible spectator cases are imposed. This is an optimization
diagnostic, not complete face coverage or a new dimension theorem.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_parameter_weight_stress_20261015_cross as X

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_rank27_cross_diagnostic_20261015'
BOUND = 27


def build():
    data = X.build()
    cases, unique = [], {}
    for fc in data['faceconstraints']:
        face = data['faces'][fc['face']]
        residual = BOUND - (face['row']['nodes'] - 2)
        feasible = {min(2, a) for a in range(residual + 1)
                    if residual - a == 0 or residual - a >= 2}
        if fc['external_A1_count_capped2'] not in feasible:
            continue
        co = tuple(tuple(v) for v in fc['coefficients'])
        rhs = F(fc['baseline']) - 1
        unique.setdefault((co, rhs), []).append(dict(
            face=fc['face'], origins=face['origins'],
            external_A1_count_capped2=fc['external_A1_count_capped2']))
        cases.append(fc)
    rows, rhs, meta = [], [], []
    for (co, b), witnesses in sorted(unique.items()):
        rows.append({j: -v for j, v in co})
        rhs.append(b)
        meta.append(dict(kind='face', witnesses=witnesses))
    data.update(faceconstraints=cases, rows=rows, rhs=rhs, meta=meta,
                components=[c for c in data['components'] if c['rank'] <= BOUND])
    return data


def run(saved=None):
    data = build()
    # Parameterize the imported exact solver only inside this process. Its
    # original file and every frozen certificate remain unchanged.
    engine = X.C
    old_bound, old_maxforests = engine.BOUND, engine.maxforests
    engine.BOUND = BOUND
    engine.maxforests = lambda d, x, bound=BOUND: old_maxforests(d, x, bound)
    try:
        solution = engine.solve(data) if saved is None else saved['solution']
        certificate = engine.verify(data, solution)
    finally:
        engine.BOUND, engine.maxforests = old_bound, old_maxforests
    target = F((BOUND - 1) * (BOUND - 5), 8)
    paths = [f'scripts/{TAG}.py',
             'scripts/fresh_sqrt5_a1_count_weights_20261015.py',
             'scripts/fresh_sqrt5_parameter_weight_stress_20261015.py',
             'scripts/fresh_sqrt5_parameter_weight_stress_20261015_cross.py',
             'results/fresh_sqrt5_parameter_weight_stress_20261015_cross.json']
    return dict(status='PASS_EXACT_PARTIAL_RANK27_WEIGHT_DIAGNOSTIC',
                sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths},
                dimension=BOUND, face_cases=len(data['faceconstraints']),
                distinct_constraints=len(data['rows']), variables=data['zi'],
                components=len(data['components']), solution=solution,
                exact_certificate=certificate, target=str(target),
                meets_target=F(solution['optimum']) <= target,
                scope='Only the frozen partial parameter-family input catalogue. '
                      'Exact optimization in the specified connected and unique-A1/E7 '
                      'weight class; no complete finite-attachment coverage and no '
                      'global nonexistence or method-independent obstruction claim.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    path = ROOT / f'results/{TAG}.json'
    saved = json.loads(path.read_text()) if args.verify else None
    result = json.loads(json.dumps(run(saved)))
    if args.verify:
        assert result == saved
    else:
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['solution']['optimum'],
          'target', result['target'], 'meets', result['meets_target'])
