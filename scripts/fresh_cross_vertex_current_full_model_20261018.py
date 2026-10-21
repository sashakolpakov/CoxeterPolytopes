#!/usr/bin/env python3
"""Append all selected post-raw110 exports, including several-vertex costs.

One bounded solve of a changed finite matrix. Exact recovery uses the actual
budget matrix: a cross-vertex column is never treated as a zero-cost transfer.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from math import gcd, lcm
from pathlib import Path
import shutil
import time

import fresh_weighted_raw110_full_model_20261017 as BASE

ROOT = BASE.ROOT
ENGINE = BASE.ENGINE
TAG = 'fresh_cross_vertex_current_full_model_20261018'
OUT = ROOT/f'results/{TAG}.json'
INPUTS = (
    ('fresh_current_historical_weighted_corners_20261017', 78),
    ('fresh_unowned_cover_pair_geometry_20261017', 5),
    ('fresh_historical22_corner_cones_20261017', 17),
    ('fresh_historical22_eight_pair_geometry_20261017', 1),
    ('fresh_raw125_raw205_actual_corner_frontier_20261017', 82),
    ('fresh_290_class_shortlist_pair_geometry_20261017', 8),
    ('fresh_2599_raw110_raw24_pair_geometry_20261017', 1),
    ('fresh_forest1325_cross_vertex_incidence_20261017', 1),
    ('fresh_raw23_raw25_cross_vertex_family_20261017', 14),
    ('fresh_shared_prism_tetrahedron_incidence_20261017', 2),
    ('fresh_raw39_shared_edge_target_incidence_20261018', 2),
    ('fresh_raw39_additional_target_frontier_20261018', 2),
)


def read(stem):
    return json.loads((ROOT/f'results/{stem}.json').read_text())


def canonical(column):
    faces = Counter()
    for f, value in column['positive_face_coefficients']:
        assert F(value) > 0
        faces[f] += F(value)
    if 'vertex_capacities' in column:
        raw = column['vertex_capacities']
    elif 'source_vertex' in column:
        raw = [[column['source_vertex'], 1],
               [column['target_vertex'], column['target_vertex_cost']]]
    else:
        raw = [[column['vertex'], column['vertex_cost']]]
    costs = Counter()
    for v, value in raw:
        assert F(value) >= 0
        costs[v] += F(value)
    costs = {v: a for v, a in costs.items() if a}
    assert faces and costs
    den = lcm(*(a.denominator for a in list(faces.values())+list(costs.values())))
    div = gcd(*(int(a*den) for a in list(faces.values())+list(costs.values())))
    assert den > 0 and div > 0
    terms = tuple(sorted((f, int(a*den)//div) for f, a in faces.items()))
    caps = tuple(sorted((v, int(a*den)//div) for v, a in costs.items()))
    assert all(F(m*div, den) == faces[f] for f, m in terms)
    assert all(F(m*div, den) == costs[v] for v, m in caps)
    return terms, caps, dict(denominator=den, integer_gcd=div)


def build():
    assert shutil.disk_usage(ROOT).free > 200*1024**2
    d, parent = BASE.build()
    assert json.loads(BASE.dump(parent)) == read(BASE.TAG)['prepared']
    nv, nf, nb, nz = map(len, (d['variables'], d['faces'], d['budgets'], d['zeros']))
    assert (nv, nf, nb, nz) == (439410, 82806, 6558, 96176)
    old_variables = BASE.digest(d['variables'])
    old_zeros = BASE.digest(d['zeros'])
    old_faces = [[list(x) for x in f['coefficients']] for f in d['faces']]
    old_budgets = [[list(x) for x in b['coefficients']] for b in d['budgets']]
    old_meta = [(f['index'], f['kind'], f['required']) for f in d['faces']]
    old_active = set(d['active'])
    final = read(INPUTS[-1][0])
    alpha = {f: F(a) for f, a in final['face_multipliers']}
    mu = {v: F(a) for v, a in final['repairs'][0]['vertex_multipliers']}
    parent_dual = next(c for c in read(BASE.TAG)['certificates']
                       if c['kind'] == 'EXACT_PARTIAL_MODEL_DUAL')
    assert final['face_multipliers'] == parent_dual['face_multipliers']
    assert all(mu.get(v, F()) >= F(a) for v, a in parent_dual['vertex_multipliers'])
    assert F(final['repairs'][0]['bound']) < 12
    columns, bindings, counts, records = {}, {}, [], []
    for stem, expected in INPUTS:
        doc = read(stem)
        assert doc['status'].startswith('PASS') and len(doc['columns']) == expected
        assert doc['face_multipliers'] == final['face_multipliers']
        for path, expected_hash in doc['sources'].items():
            if path not in bindings:
                bindings[path] = sha256((ROOT/path).read_bytes()).hexdigest()
            assert bindings[path] == expected_hash, (stem, path)
        counts.append([stem, expected])
        for index, column in enumerate(doc['columns']):
            terms, caps, scale = canonical(column)
            assert all(0 <= f < nf and value > 0 for f, value in terms)
            assert all(0 <= v < nb and value > 0 for v, value in caps)
            load = sum((value*alpha.get(f, F()) for f, value in terms), F())
            cost = sum((value*mu.get(v, F()) for v, value in caps), F())
            assert load <= cost, (stem, index, load, cost)
            columns.setdefault((terms, caps), []).append(dict(package=stem, index=index, scaling=scale))
    assert sum(n for _, n in counts) == 213
    for index, ((terms, caps), origins) in enumerate(sorted(columns.items())):
        j = len(d['variables'])
        # Metadata retains all owners. Optimization and recovery use budgets,
        # not a fictitious single owner or the metadata's null vertex field.
        d['variables'].append(dict(index=j, kind='current_cross_vertex_capacity',
            vertex=caps[0][0] if len(caps) == 1 else None,
            cost=caps[0][1] if len(caps) == 1 else 0,
            vertex_capacities=[list(x) for x in caps], source_index=index))
        for f, value in terms:
            d['faces'][f]['coefficients'].append([j, value])
        for v, value in caps:
            d['budgets'][v]['coefficients'].append([j, value])
        records.append(dict(variable=j, integer_face_coefficients=[list(x) for x in terms],
            integer_vertex_capacities=[list(x) for x in caps], origins=origins))
    assert BASE.digest(d['variables'][:nv]) == old_variables
    assert BASE.digest(d['zeros']) == old_zeros
    assert (len(d['faces']), len(d['budgets']), len(d['zeros'])) == (nf, nb, nz)
    for i, f in enumerate(d['faces']):
        assert (f['index'], f['kind'], f['required']) == old_meta[i]
        assert [[j, a] for j, a in f['coefficients'] if j < nv] == old_faces[i]
        assert all(a > 0 for j, a in f['coefficients'] if j >= nv)
    for i, b in enumerate(d['budgets']):
        assert [[j, a] for j, a in b['coefficients'] if j < nv] == old_budgets[i]
        assert all(a > 0 for j, a in b['coefficients'] if j >= nv)
    d['active'] = sorted({j for f in d['faces'] for j, a in f['coefficients'] if a > 0})
    assert old_active <= set(d['active'])
    eligible = set(d['repair_indices'])
    assert all(j < nv for j in eligible)
    assert all(sum(a for j, a in f['coefficients'] if j in eligible) > 0 for f in d['faces'])
    assert all(a > 0 for f in d['faces'] for j, a in f['coefficients'] if j in eligible)
    model_hash = sha256()
    for group in (d['variables'], d['faces'], d['budgets'], d['zeros']):
        model_hash.update(BASE.digest(group).encode()); model_hash.update(b'\n')
    paths = [Path(__file__), Path(BASE.__file__), Path(ENGINE.__file__),
             ROOT/f'results/{BASE.TAG}.json', ROOT/f'results/{BASE.TAG}.md']
    for stem, _ in INPUTS:
        paths += [ROOT/f'{folder}/{stem}.{ext}' for folder, ext in
                  (('scripts', 'py'), ('results', 'json'), ('results', 'md'))]
    prepared = dict(sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        dimension=13, ground_field='Q(sqrt5)', target=12,
        inherited_model_sha256=parent['model_sha256'], model_sha256=model_hash.hexdigest(),
        provenance_to_face=parent['provenance_to_face'], new_columns=records, input_column_counts=counts,
        counts=dict(inherited_variables=nv, variables=len(d['variables']), faces=nf, vertex_forests=nb,
            active_variables=len(d['active']), zero_budget_signed_variables=nz,
            inherited_vertex_cost_signed_variables=52, exported_columns=213,
            distinct_added_columns=len(records), duplicate_exports=213-len(records),
            several_vertex_columns=sum(len(r['integer_vertex_capacities']) > 1 for r in records),
            face_matrix_nonzeros=sum(len(f['coefficients']) for f in d['faces'])),
        latest_fixed_mass_ratio=final['repairs'][0]['bound'],
        preservation='Every439410 inherited variable,82806 face row,6558 budget and96176 zero-budget record keeps its exact coefficients. All52 positive-cost signed columns stay budgeted. All213 selected actual exports are included, modulo exact positive scaling duplicates; hypothetical zero-event exports are excluded.',
        recovery='Exact primal repair uses inherited positive old columns. Exact dual recovery uses the full budget matrix and nonnegative multipliers; several-vertex columns remain charged to every stated vertex. Old zero-budget transfers receive their inherited simplex-mass repair.',
        scope='Enlarged finite exported compact-arithmetic H13 model over fixed Q(sqrt5). Complete implicit packing cones and missing face/forest/field coverage are not supplied. Neither a target12 primal nor a dual obstruction alone is a dimension theorem or geometric realization.')
    return d, json.loads(BASE.dump(prepared))


def cost_rows(d):
    costs = [[] for _ in d['variables']]
    for v, budget in enumerate(d['budgets']):
        assert budget['index'] == v
        for j, value in budget['coefficients']:
            assert 0 <= j < len(costs) and value > 0
            costs[j].append((v, value))
    assert {j for j, row in enumerate(costs) if not row} == {z['variable'] for z in d['zeros']}
    return costs


def loads(d, alpha):
    values = [F() for _ in d['variables']]
    for f, value in alpha.items():
        for j, coefficient in d['faces'][f]['coefficients']:
            values[j] += coefficient*value
    return values


def verify_dual(d, certificate):
    alpha = {f:F(a) for f, a in certificate['face_multipliers']}
    mu = {v:F(a) for v, a in certificate['vertex_multipliers']}
    assert len(alpha) == len(certificate['face_multipliers'])
    assert len(mu) == len(certificate['vertex_multipliers'])
    assert all(0 <= f < len(d['faces']) and a > 0 for f, a in alpha.items())
    assert all(0 <= v < len(d['budgets']) and a > 0 for v, a in mu.items())
    values = loads(d, alpha)
    for j, row in enumerate(cost_rows(d)):
        assert values[j] <= sum((a*mu.get(v, F()) for v, a in row), F()), j
    numerator = sum((d['faces'][f]['required']*a for f, a in alpha.items()), F())
    denominator = sum(mu.values(), F())
    assert denominator > 0
    assert str(numerator) == certificate['objective_numerator']
    assert str(denominator) == certificate['normalizer']
    assert str(numerator/denominator) == certificate['bound']
    assert certificate['excludes_target12'] == (numerator > 12*denominator)


def recover_dual(d, marginals):
    denominator = 10**10
    alpha = {f:F(max(0, round(-float(marginals[f])*denominator)), denominator)
             for f in range(len(d['faces']))}
    # The inherited signed transfer order provides the same simplex repair
    # as the old engine; all resulting zero-cost inequalities are checked below.
    for zero in d['zeros']:
        need = sum((m*alpha.get(f, F()) for f, m in zero['positive_coefficients']), F())/zero['capacity']
        f = zero['simplex_face']; alpha[f] = max(alpha.get(f, F()), need)
    alpha = {f:a for f, a in alpha.items() if a}
    mu = {v:F(max(0, round(-float(marginals[len(d['faces'])+v])*denominator)), denominator)
          for v in range(len(d['budgets']))}
    values = loads(d, alpha)
    repairs = []
    for j, row in enumerate(cost_rows(d)):
        if not row:
            assert values[j] <= 0, j
            continue
        missing = values[j]-sum((a*mu[v] for v, a in row), F())
        if missing > 0:
            # Increasing one nonnegative vertex allowance cannot invalidate
            # any other column. Largest coefficient minimizes this increment.
            v, coefficient = min(row, key=lambda x:(-x[1], x[0]))
            increase = missing/coefficient; mu[v] += increase
            repairs.append([j, v, str(increase)])
    numerator = sum((d['faces'][f]['required']*a for f, a in alpha.items()), F())
    normalizer = sum(mu.values(), F())
    assert normalizer > 0
    certificate = dict(kind='EXACT_PARTIAL_MODEL_DUAL',
        face_multipliers=[[f,str(a)] for f, a in sorted(alpha.items())],
        vertex_multipliers=[[v,str(a)] for v, a in sorted(mu.items()) if a],
        objective_numerator=str(numerator), normalizer=str(normalizer), bound=str(numerator/normalizer),
        excludes_target12=numerator > 12*normalizer, allowance_repairs=repairs,
        recovery='Rounded solver face/budget multipliers, inherited simplex repair, monotone exact full-budget repair.')
    verify_dual(d, certificate)
    return certificate


def verify(d, certificate):
    if certificate['kind'] == 'EXACT_PARTIAL_MODEL_PRIMAL':
        maximum = ENGINE.verify_primal(d, certificate)
        assert certificate['meets_target12'] == (maximum <= 12)
    else:
        assert certificate['kind'] == 'EXACT_PARTIAL_MODEL_DUAL'
        verify_dual(d, certificate)


def optimize(d):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    active = d['active']; indices = {j:i for i, j in enumerate(active)}; t = len(active)
    rr, cc, vv, rhs = [], [], [], []
    for face in d['faces']:
        row = len(rhs); rhs.append(-face['required'])
        for j, value in face['coefficients']:
            if j in indices:
                rr.append(row); cc.append(indices[j]); vv.append(-value)
    for budget in d['budgets']:
        row = len(rhs); rhs.append(0)
        for j, value in budget['coefficients']:
            if j in indices:
                rr.append(row); cc.append(indices[j]); vv.append(value)
        rr.append(row); cc.append(t); vv.append(-1)
    matrix = coo_matrix((np.asarray(vv,dtype=float),(rr,cc)),shape=(len(rhs),t+1)).tocsr()
    objective = np.zeros(t+1); objective[t] = 1
    print('ONE_CHANGED_MATRIX_HIGHS_IPM_BEGIN', matrix.shape, matrix.nnz, flush=True)
    start = time.monotonic()
    result = linprog(objective,A_ub=matrix,b_ub=np.asarray(rhs,dtype=float),bounds=(0,None),
                     method='highs-ipm',options={'time_limit':300})
    solver = dict(status=int(result.status), message=str(result.message), success=bool(result.success),
        objective=None if result.fun is None else float(result.fun), elapsed_seconds=time.monotonic()-start,
        method='highs-ipm', time_limit_seconds=300, iterations=int(getattr(result,'nit',0)))
    print('ONE_CHANGED_MATRIX_HIGHS_IPM_END', solver, flush=True)
    certificates = []
    if result.x is not None:
        certificates.append(ENGINE.recover_primal(d,result.x))
    marginals = getattr(getattr(result,'ineqlin',None),'marginals',None)
    if marginals is not None:
        certificates.append(recover_dual(d,marginals))
    for certificate in certificates:
        verify(d,certificate)
    return solver, certificates


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare',action='store_true')
    mode.add_argument('--solve',action='store_true')
    mode.add_argument('--verify',action='store_true')
    args = ap.parse_args()
    d, prepared = build()
    print('PREPARED_CROSS_VERTEX_CURRENT_MODEL',prepared['counts'],prepared['model_sha256'],flush=True)
    if args.prepare:
        assert not OUT.exists(), 'Existing package; no overwrite.'
        OUT.write_text(BASE.dump(dict(status='PREPARED_CROSS_VERTEX_CURRENT_MODEL',prepared=prepared,
            solver_invocations=0,solver=None,certificates=[]))+'\n')
    elif args.verify:
        saved = json.loads(OUT.read_text()); assert saved['prepared'] == prepared
        for certificate in saved['certificates']:
            verify(d,certificate)
        print('PASS_SAVED_CROSS_VERTEX_CURRENT_MODEL',saved['status'],flush=True)
    else:
        saved = json.loads(OUT.read_text())
        assert saved['prepared'] == prepared and saved['solver_invocations'] == 0
        saved.update(status='ONE_BOUNDED_SOLVE_STARTED',solver_invocations=1)
        OUT.write_text(BASE.dump(saved)+'\n')
        solver, certificates = optimize(d)
        status = 'PASS_EXACT_CROSS_VERTEX_MODEL_CERTIFICATES' if certificates else 'UNRESOLVED_SINGLE_BOUNDED_SOLVE'
        if any(c.get('meets_target12') for c in certificates):
            status = 'PASS_EXACT_PARTIAL_PRIMAL_TARGET12'
        saved.update(status=status,solver=solver,certificates=certificates)
        OUT.write_text(BASE.dump(saved)+'\n')
        print(status,[(c['kind'],c.get('maximum',c.get('bound'))) for c in certificates],flush=True)
