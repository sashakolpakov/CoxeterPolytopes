#!/usr/bin/env python3
"""Numerical polar-truncation bank with optional old-wall replacement.

All numerical predicates are proposal gates. Neither a failed SAT search
nor LP scouts are exact geometric or nonexistence certificates.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.optimize import linprog
import z3

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_makarov_polar_closure_20261013'
SOURCE = ROOT/'results/fresh_h9_makarov_joint_20261013_exact.json'
SCOUT = ROOT/'results/fresh_h9_makarov_joint_20261013_new_frontier_scout.json'


def f(x):
    return float(Fraction(x[0]))+np.sqrt(5)*float(Fraction(x[1]))


def setup():
    source = json.loads(SOURCE.read_text())
    scout = json.loads(SCOUT.read_text())
    Q = np.array([[f(x) for x in row] for row in source['form']])
    R = [[f(x) for x in row] for row in source['seed_normals']]
    R += [[f(a)+np.sqrt(f(row['radicand']))*f(b) for a, b in row['normal']]
          for row in source['caps']]
    rows = [row for row in scout['vertices'] if row['type'] == 'spacelike']
    R += [row['point'] for row in rows]
    R = np.array(R)
    norm = np.einsum('ij,jk,ik->i', R, Q, R)
    assert min(norm) > 0
    R /= np.sqrt(norm)[:, None]
    c = np.array([f(x) for x in source['strict_center']])
    c /= np.sqrt(-c@Q@c)
    assert max(R@Q@c) < 0
    refs = [{'old_wall': i} for i in range(27)]+[
        {'polar_vertex': row['id'], 'support': row['active']} for row in rows]
    return source, Q, R, c, refs, np.array([row['point'] for row in rows])


def volume(Q, R, c, objectives, rng):
    A = R@Q
    A /= np.linalg.norm(A, axis=1)[:, None]
    best = (-1e100, None)
    for i in range(objectives):
        direction = Q@best[1] if i % 3 == 0 and best[1] is not None else rng.normal(size=10)
        result = linprog(-direction, A_ub=A, b_ub=np.zeros(len(R)),
                         A_eq=[c@Q], b_eq=[-1], bounds=[(None, None)]*10, method='highs')
        if result.status == 3:
            # A finite Euclidean box supplies a concrete positive witness
            # rather than treating unboundedness as a missing point.
            result = linprog(-direction, A_ub=A, b_ub=np.zeros(len(R)),
                             A_eq=[c@Q], b_eq=[-1], bounds=[(-1e5, 1e5)]*10, method='highs')
        if not result.success:
            continue
        square = float(result.x@Q@result.x)
        if square > best[0]:
            best = square, result.x
    return best


def run(rounds=500, cap_budget=100):
    start = time.time()
    source, Q, R, c, refs, points = setup()
    p = R@Q@R.T
    finite = np.array([-np.cos(np.pi/m) for m in range(2, 41)])
    ok = p < -1+2e-6
    for a in finite:
        ok |= abs(p-a) < 2e-6
    np.fill_diagonal(ok, True)
    base_compatible = [i for i in range(27, len(R)) if all(ok[i, :27])]
    frozen_escape = np.array(source['escape']['point'], dtype=float)
    frozen_escape /= max(abs(frozen_escape))
    rng = np.random.default_rng(20431013)
    whole = volume(Q, R, c, 100, rng)
    s = z3.Solver()
    s.set(timeout=10000)
    variables = [z3.Bool('polar_'+str(i)) for i in range(len(R))]
    if cap_budget:
        s.add(z3.PbLe([(v, 1) for v in variables], cap_budget))
    coverage = R@Q@points.T > 1e-7
    clause_keys = set()
    for mask in coverage.T:
        ids = tuple(map(int, np.flatnonzero(mask)))
        if ids not in clause_keys:
            clause_keys.add(ids)
            s.add(z3.Or(*[variables[i] for i in ids]))
    trace = []
    conflicts = set()
    witnesses = []
    answer = 'ROUND_LIMIT'
    selected = []
    for step in range(rounds):
        status = s.check()
        if status != z3.sat:
            answer = 'NUMERICAL_'+str(status).upper()
            break
        model = s.model()
        selected = [i for i, v in enumerate(variables) if z3.is_true(model.eval(v, model_completion=True))]
        bad = [(a, b) for k, a in enumerate(selected) for b in selected[:k] if not ok[a, b]]
        trace.append({'round': step, 'selected': len(selected), 'new_conflicts': len(bad)})
        if bad:
            for a, b in bad:
                key = tuple(sorted([a, b]))
                assert key not in conflicts
                conflicts.add(key)
                s.add(z3.Or(z3.Not(variables[a]), z3.Not(variables[b])))
        else:
            escape_square, point = volume(Q, R[selected], c, 60, rng)
            trace[-1]['maximum_scout_square'] = escape_square
            if point is None:
                answer = 'LP_UNRESOLVED'
                break
            if escape_square <= 1e-7:
                answer = 'NUMERICAL_COMPATIBLE_NO_ESCAPE_IN_SCOUT'
                break
            ids = list(map(int, np.flatnonzero(R@Q@point > 1e-7)))
            witnesses.append({'square': escape_square, 'point': point.tolist(), 'cutters': ids})
            if not ids:
                answer = 'NUMERICAL_WHOLE_BANK_ESCAPE_FOUND'
                break
            s.add(z3.Or(*[variables[i] for i in ids]))
        if step % 20 == 0:
            print(step, len(selected), 'conflicts', len(conflicts), 'new_gaps', len(witnesses), flush=True)
    out = {'status': answer, 'numerical_only': True,
           'sources': {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest() for p in [SOURCE, SCOUT]},
           'bank_size': len(R), 'old_walls_optional': True, 'maximum_selected_facets': cap_budget,
           'initial_spacelike_witnesses': len(points), 'distinct_initial_coverage_clauses': len(clause_keys),
           'polar_caps_compatible_with_all_27': [refs[i] for i in base_compatible],
           'number_of_these_cutting_frozen_escape': sum(float(R[i]@Q@frozen_escape) > 1e-7 for i in base_compatible),
           'whole_bank_scout_maximum_square': whole[0],
           'whole_bank_scout_point': None if whole[1] is None else whole[1].tolist(),
           'learned_conflicts': len(conflicts), 'trace': trace, 'new_global_witnesses': witnesses,
           'selected_references': [refs[i] for i in selected], 'selected_normals': R[selected].tolist(),
           'seconds': time.time()-start,
           'scope': 'Numerical polar bank, approximate compatibility and bounded SAT/LP scout. No exact infeasibility, exhaustive geometric census, compactness, or finite-volume claim.'}
    (ROOT/'results'/f'{TAG}.json').write_text(json.dumps(out, indent=2)+'\n')
    print(answer, 'conflicts', len(conflicts), 'seconds', time.time()-start, flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--rounds', type=int, default=500)
    ap.add_argument('--cap-budget', type=int, default=100)
    args = ap.parse_args()
    run(args.rounds, args.cap_budget)
