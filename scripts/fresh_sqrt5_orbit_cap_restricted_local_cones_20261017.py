#!/usr/bin/env python3
"""Fresh complete local covers after additional cross-coface orbit caps.

Only bounded local packing LPs are used. Saved replay invokes none.
The cap package is an explicit frozen input selected at generation.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import time

import fresh_sqrt5_cross_coface_cap_full_model_20261017 as N
import fresh_sqrt5_restricted_flag_local_cones_20261017 as R

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_orbit_cap_restricted_local_cones_20261017'
SUPPORT = 'fresh_sqrt5_cross_coface_cap_dual_support_20261017'
PROGRESS = Path('/private/tmp')/(TAG+'_progress.jsonl')


def dump(x): return json.dumps(x, sort_keys=True, separators=(',', ':'))


def verify_columns(d, alpha, mu, caps):
    loads = defaultdict(F)
    for fi, a in alpha.items():
        assert a > 0
        for j, m in d['faces'][fi]['coefficients']: loads[j] += a*m
    for col in d['variables']:
        if col['vertex'] is None: assert loads[col['index']] <= 0
        else: assert loads[col['index']] <= col['cost']*mu.get(col['vertex'], F())
    for c in caps:
        load = sum((m*alpha.get(fi, F()) for fi, m in c['positive_face_coefficients']), F())
        assert load <= c['capacity']*alpha.get(c['simplex_face'], F())


def run(cap_tag, saved=None):
    support = N.read(SUPPORT); capdata = N.read(cap_tag)
    d, prepared = N.build()
    assert json.loads(json.dumps(prepared)) == N.read(N.TAG)['prepared']
    test = next(t for t in support['tests'] if t['name'] == 'all_S_and_selected_actual_Q')
    assert not test['restored_S_outside_seed_domain']
    oldcert = test['repaired_dual']; N.N.verify(d, oldcert)
    olda = {i: F(a) for i, a in oldcert['face_multipliers']}; alpha = dict(olda)
    pmap = dict(prepared['provenance_to_face'])
    caploads = []
    for c in capdata['columns']:
        positive = [[pmap[j], m] for j, m in c['provenance_coefficients']]
        assert positive == c['positive_face_coefficients']
        assert c['capacity'] > 0
        load = sum((m*olda.get(i, F()) for i, m in positive), F())
        fi = c['simplex_face']; need = load/c['capacity']
        alpha[fi] = max(alpha.get(fi, F()), need)
        caploads.append([c['index'], str(load), str(c['capacity']*olda.get(fi, F()))])
    alpha = {i: a for i, a in alpha.items() if a}
    initial = N.N.M.dual_certificate(d, alpha); N.N.verify(d, initial)
    alpha = {i: F(a) for i, a in initial['face_multipliers']}
    mu = {v: F(a) for v, a in initial['vertex_multipliers']}
    assert all(alpha.get(i, F()) >= a for i, a in olda.items())
    qrows = {f['index'] for f in d['faces'] if f['kind'] == 'canonical_dotted_quad_aggregate'}
    assert all(alpha.get(i, F()) == olda.get(i, F()) for i in qrows)
    assert sum(i in qrows for i in alpha) == 336
    verify_columns(d, alpha, mu, capdata['columns'])
    print('ORBIT_CAP_INITIAL', len(alpha), len(capdata['columns']), float(F(initial['bound'])), flush=True)
    columns, features, dimensions, corners = R.A.corner_features(d, alpha)
    columns, features = json.loads(dump(columns)), json.loads(dump(features))
    literals = R.literal_orbits(d, columns)
    weights = defaultdict(F); features_by_vertex = defaultdict(list)
    vmap = {c['index']: c['vertex'] for c in columns}
    for f in features:
        touched = set()
        for j, m in f['block_orbits']:
            weights[j] += m*alpha[f['face']]; touched.add(vmap[j])
        for v in touched: features_by_vertex[v].append(f)
    grouped = defaultdict(list)
    for c in literals: grouped[c['vertex']].append(c)
    for cs in grouped.values(): cs.sort(key=lambda c: c['index'])
    order = sorted(grouped, key=lambda v: (-mu.get(v, F()), -len(grouped[v]), v))
    key = sha256(dump([cap_tag, capdata['sources'], initial['face_multipliers']]).encode()).hexdigest()
    checkpoint = {}
    if saved is None and PROGRESS.exists():
        for line in PROGRESS.read_text().splitlines():
            item = json.loads(line); assert item['key'] == key
            assert item['vertex'] not in checkpoint
            checkpoint[item['vertex']] = item['record']
    oldrows = {} if saved is None else {r['vertex']: r for r in saved['local_cones']}
    repaired = dict(mu); records = []; cuts = []; open_cases = []; start = time.monotonic()
    for serial, v in enumerate(order):
        ws = [weights[c['index']] for c in grouped[v]]
        try: pack = R.L.Packings(grouped[v])
        except R.L.ResourceStop as exc:
            rec = dict(vertex=v, forest=d['budgets'][v]['forest'], status='OPEN_ENUMERATION', reason=str(exc))
            open_cases.append(v); records.append(rec)
            if saved is not None: assert rec == oldrows[v]
            elif v not in checkpoint:
                with PROGRESS.open('a') as stream: stream.write(dump(dict(key=key, vertex=v, record=rec))+'\n')
            continue
        if saved is not None: certificate = oldrows[v]['certificate']
        elif v in checkpoint:
            assert 'certificate' in checkpoint[v], 'An old open resource stop requires a fresh successor.'
            certificate = checkpoint[v]['certificate']
        else: certificate = R.L.solve_and_recover(pack, ws, max(mu.get(v, F()), max(ws)))
        lower, upper, slack = R.L.verify_certificate(pack, ws, certificate)
        repaired[v] = max(mu.get(v, F()), upper)
        if certificate['status'].startswith('OPEN'): open_cases.append(v)
        if lower > mu.get(v, F()):
            q = {j: F(a) for j, a in certificate['charges']}; terms = []
            for f in features_by_vertex[v]:
                a = sum((m*q.get(j, F()) for j, m in f['block_orbits']), F())
                if a: terms.append([f['face'], str(a)])
            assert sum((alpha[i]*F(a) for i, a in terms), F()) == lower
            cuts.append(dict(index=len(cuts), vertex=v, vertex_cost=1,
                positive_simplex_coefficients=terms, exact_load=str(lower), old_vertex_mass=str(mu.get(v, F())),
                excess=str(lower-mu.get(v, F())), orbit_charges=certificate['charges']))
        rec = dict(vertex=v, forest=d['budgets'][v]['forest'], status=certificate['status'],
            counts=dict(block_orbits=len(pack.ids), literal_blocks=len(pack.masks), packing_count_vectors=len(pack.vectors),
                cached_vector_entries=pack.entries, positive_mixture_packings=len(certificate['packing_mixture'])),
            block_orbit_marginals=[[j, str(a)] for j, a in zip(pack.ids, ws)], literal_block_orbits=grouped[v],
            count_vector_digest=sha256(dump(pack.vectors).encode()).hexdigest(), certificate=certificate,
            primal_cover_slacks=slack, old_vertex_mass=str(mu.get(v, F())), repaired_vertex_mass=str(repaired[v]))
        records.append(rec)
        if saved is None and v not in checkpoint:
            with PROGRESS.open('a') as stream: stream.write(dump(dict(key=key, vertex=v, record=rec))+'\n')
        if serial % 40 == 0 or lower > mu.get(v, F()):
            print('ORBIT_CAP_LOCAL', serial+1, '/', len(order), 'cuts', len(cuts), 'open', len(open_cases),
                  'seconds', round(time.monotonic()-start, 2), flush=True)
        pack.dp.cache_clear()
    verify_columns(d, alpha, repaired, capdata['columns'])
    num = sum((a*d['faces'][i]['required'] for i, a in alpha.items()), F()); den = sum(repaired.values(), F())
    assert num == F(initial['objective_numerator'])
    paths = [f'scripts/{TAG}.py', f'scripts/{N.TAG}.py', f'results/{N.TAG}.json',
        f'scripts/{R.TAG}.py', f'scripts/{R.A.TAG}.py', f'scripts/{R.L.TAG}.py',
        f'results/{SUPPORT}.json', f'scripts/{cap_tag}.py', f'results/{cap_tag}.json']
    return dict(status='PASS_COMPLETE_NEW_ORBIT_CAP_LOCAL_CONES' if not open_cases else 'PASS_NEW_LOCAL_CONES_WITH_OPEN_CASES',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        cap_package=cap_tag, input_test='all_S_and_selected_actual_Q', inherited_model_sha256=prepared['model_sha256'],
        global_optimizer_invocations=0, saved_replay_optimizer_invocations=0,
        reused_old_fixed_mass_certificates=0, resource_guards=dict(count_vector_entries=R.L.ENTRY_LIMIT,
            enumeration_seconds_per_forest=R.L.ENUM_SECONDS, solver_seconds_per_forest=R.L.SOLVE_SECONDS),
        counts=dict(positive_faces=len(alpha), simplex_faces=len(features), actual_quad_faces=336,
            positive_simplex_forests=len(order), block_orbits=len(columns), simplex_corners=corners,
            raised_simplex_masses=sum(a > olda.get(i, F()) for i, a in alpha.items()),
            exact_local_optima=sum(r['status'] == 'EXACT_LOCAL_OPTIMUM' for r in records),
            exact_cover_brackets=sum(r['status'] == 'EXACT_COVER_BRACKET' for r in records),
            open_forests=len(open_cases), new_exact_cut_forests=len(cuts),
            positive_mixture_packings=sum(r.get('counts', {}).get('positive_mixture_packings', 0) for r in records),
            inherited_signed_columns=len(d['zeros']), additional_cap_columns=len(capdata['columns'])),
        cap_loads_before_repair=caploads, initial_after_cap_repair=initial,
        initial_restricted_bound=oldcert['bound'], bound=str(num/den), objective_numerator=str(num), normalizer=str(den),
        complete_positive_simplex_forest_coverage=not open_cases, open_forests=open_cases,
        fresh_simplex_dimension_counts=dimensions, priority_order=order,
        block_orbit_columns=columns, positive_simplex_features=features, literal_block_orbits=literals,
        local_cones=records, new_exact_charge_columns=cuts, face_multipliers=initial['face_multipliers'],
        vertex_multipliers=[[v, str(a)] for v, a in sorted(repaired.items()) if a],
        scope='Fresh changed-alpha complete local simplex packing test after the supplied exact orbit-cap family. Common336actual Q masses fixed; inherited94722 signed columns plus every supplied cap rechecked. Local solves are bounded and exact certificates replay without optimization. New whole joint3/4/5 families are not claimed beyond inherited finite model cuts. No global optimization, complete ambient realization, or dimension bound.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--caps'); p.add_argument('--verify', action='store_true'); a = p.parse_args()
    out = ROOT/f'results/{TAG}.json'; saved = json.loads(out.read_text()) if a.verify else None
    cap_tag = saved['cap_package'] if saved is not None else a.caps
    assert cap_tag and (not a.caps or a.caps == cap_tag)
    result = json.loads(dump(run(cap_tag, saved)))
    if a.verify: assert result == saved
    else: out.write_text(dump(result)+'\n')
    print(result['status'], result['counts'], result['bound'], float(F(result['bound'])), flush=True)
