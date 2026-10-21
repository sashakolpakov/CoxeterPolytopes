#!/usr/bin/env python3
"""Exact dual recovery on numerically suggested active triangle support, not a new full LP."""
import argparse, hashlib, json
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017 as B
import fresh_sqrt5_rank13_triangle_dual_core_20261017 as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_rank13_active_triangle_dual_20261017'
INPUT = C.INPUT
EPSILON = F(1, 10**6)

def domain():
    original = B.read(INPUT)
    primal = next(c for c in original['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_PRIMAL')
    B.verify_certificate(original, primal)
    weights = {i: F(a) for i, a in primal['weights']}
    maximum = F(primal['maximum'])
    active = {i for i, a in enumerate(primal['vertex_budgets']) if maximum-F(a) <= EPSILON}
    indices = []
    for face in original['faces']:
        if face['required'] != 2:
            continue
        if any(original['variables'][j]['vertex'] not in active for j, m in face['coefficients']):
            continue
        total = sum((m*weights.get(j, F(0)) for j, m in face['coefficients']), F(0))
        if total-2 <= EPSILON:
            indices.append(face['index'])
    vertex_ids = sorted({original['variables'][j]['vertex'] for i in indices for j, m in original['faces'][i]['coefficients']})
    variables, vertices, lookup = [], [], {}
    for old_vertex in vertex_ids:
        v = original['vertices'][old_vertex]
        vid, vv = len(vertices), []
        for old, m in v['variables']:
            j = len(variables)
            lookup[old] = j
            var = dict(original['variables'][old])
            var.update(index=j, vertex=vid, original_variable=old)
            variables.append(var)
            vv.append([j, m])
        assert sum(m for j, m in vv) == 78
        vertices.append(dict(index=vid, forest=v['forest'], variables=vv, original_vertex=old_vertex))
    faces = []
    for old in indices:
        face = original['faces'][old]
        provenance = face['provenance'][0]
        assert 'triangle' in provenance['kind']
        faces.append(dict(index=len(faces), required=2,
            coefficients=[[lookup[j], m] for j, m in face['coefficients']],
            provenance=[provenance], original_face=old))
    paths = [f'scripts/{TAG}.py', f'scripts/{B.TAG}.py', f'scripts/{C.TAG}.py', f'results/{INPUT}.json']
    return dict(sources={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        selection=dict(epsilon=str(EPSILON), saved_primal_maximum=str(maximum),
            eligible_active_vertex_count=len(active), triangle_only=True,
            scope='The active-set heuristic only selects a subset. Every final dual is verified exactly without assuming complementary slackness or numerical accuracy.'),
        variables=variables, vertices=vertices, faces=faces,
        counts=dict(vertices=len(vertices), variables=len(variables), distinct_face_constraints=len(faces)))

def run(verify):
    B.dual_bound = C.sparse_dual_bound
    d = domain()
    out = ROOT/f'results/{TAG}.json'
    if verify:
        saved = json.loads(out.read_text())
        for k, v in d.items():
            assert saved[k] == v, k
        for c in saved['certificates']:
            B.verify_certificate(d, c)
        return saved
    B.recover_dual = C.recover_unminimized
    solver, certificates = B.optimize(d)
    assert solver['success'] and certificates
    d.update(status='PASS_EXACT_ACTIVE_TRIANGLE_DUAL', solver=solver, certificates=certificates,
        scope='Only retained local triangle inequalities and complete78-pair budgets of their whole forests are imposed. This exact obstruction concerns nonnegative whole-forest corner-orbit payments, not a dimension theorem or simultaneous global realization.')
    out.write_text(json.dumps(d, indent=2)+'\n')
    return d

if __name__ == '__main__':
    p=argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); a=p.parse_args()
    d=run(a.verify)
    print(d['status'], d['solver'], [(c['kind'], c.get('maximum', c.get('bound'))) for c in d['certificates']], flush=True)
