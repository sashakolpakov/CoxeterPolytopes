#!/usr/bin/env python3
"""Complete local simplex packing for new, exactly restricted flag duals.

The entire frozen flag model is reconstructed and checked. Only local packing
LPs are allowed; fresh saved replay never optimizes. Old fixed-alpha covers
are not reused. A second restriction is tested only if the first cover ratio
is at most twelve. Resource stops remain explicit, not infeasibility claims.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import time

import fresh_sqrt5_flag_local_packing_full_model_20261017 as N
import fresh_sqrt5_all_subset_port_packing_20261017 as A
import fresh_sqrt5_all_local_simplex_packing_cones_20261017 as L

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_restricted_flag_local_cones_20261017'
SUPPORT = 'fresh_sqrt5_flag_local_packing_dual_support_20261017'
PROGRESS = Path('/private/tmp') / (TAG + '_progress.jsonl')


def dump(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'))


def sha(p):
    return sha256((ROOT / p).read_bytes()).hexdigest()


def read(tag):
    obj = json.loads((ROOT / f'results/{tag}.json').read_text())
    for doc in (obj, obj.get('prepared', {})):
        for p, h in doc.get('sources', {}).items():
            assert sha(p) == h, p
    return obj


def literal_orbits(d, columns):
    """Enumerate every root subset in the exact weighted-forest orbit."""
    models = dict(A.C.models())
    byvertex = defaultdict(list)
    for c in columns:
        byvertex[c['vertex']].append(c)
    result = []
    for serial, (v, cs) in enumerate(sorted(byvertex.items())):
        g = {}
        for name in d['budgets'][v]['forest']:
            offset = len(g)
            for i, row in models[name].items():
                g[offset + i] = {offset + j: m for j, m in row.items()}
        assert len(g) == 13
        pieces = A.C.pieces(g)
        wanted = {(c['block_size'], A.freeze(c['marked_forest_code'])): c['index'] for c in cs}
        assert len(wanted) == len(cs)
        found = defaultdict(list)
        for size in sorted({c['block_size'] for c in cs}):
            for block in combinations(range(13), size):
                j = wanted.get((size, A.code(pieces, block)))
                if j is not None:
                    found[j].append(sum(1 << p for p in block))
        assert set(found) == {c['index'] for c in cs}
        for c in cs:
            masks = found[c['index']]
            assert masks and len(masks) == len(set(masks))
            result.append(dict(index=c['index'], vertex=v, block_size=c['block_size'], literal_masks=masks))
        if serial % 100 == 0:
            print('LITERAL_ORBITS', serial + 1, '/', len(byvertex), flush=True)
    return sorted(result, key=lambda c: c['index'])


def verify_inherited(d, alpha, mu):
    """Directly check every old/new positive and signed column after repair."""
    loads = defaultdict(F)
    for fi, a in alpha.items():
        assert a > 0
        for j, m in d['faces'][fi]['coefficients']:
            loads[j] += a * m
    positive = signed = 0
    for c in d['variables']:
        if c['vertex'] is None:
            assert loads[c['index']] <= 0
            signed += 1
        else:
            assert loads[c['index']] <= c['cost'] * mu.get(c['vertex'], F())
            positive += 1
    assert signed == len(d['zeros']) == 94481
    return dict(positive_cost_columns=positive, signed_columns=signed,
                inherited_face_obligations=len(d['faces']), inherited_vertex_budgets=len(d['budgets']))


def run_case(d, test, input_digest, previous_features, saved=None):
    name = test['name']
    initial = test['repaired_dual']
    N.verify(d, initial)
    alpha = {i: F(a) for i, a in initial['face_multipliers']}
    mu = {v: F(a) for v, a in initial['vertex_multipliers']}
    assert not test['restored_S_outside_seed_domain']
    columns, features, dimensions, corners = A.corner_features(d, alpha)
    features = json.loads(dump(features))
    columns = json.loads(dump(columns))
    print('RECONSTRUCTED_BLOCKS', name, len(features), len(columns), corners, dimensions, flush=True)
    literals = literal_orbits(d, columns)
    weights = defaultdict(F)
    features_by_vertex = defaultdict(list)
    column_vertex = {c['index']: c['vertex'] for c in columns}
    for f in features:
        assert alpha[f['face']] > 0
        touched = set()
        for j, m in f['block_orbits']:
            assert isinstance(m, int) and m > 0
            weights[j] += alpha[f['face']] * m
            touched.add(column_vertex[j])
        for v in touched:
            features_by_vertex[v].append(f)
    assert set(weights) == set(column_vertex) and all(a > 0 for a in weights.values())
    grouped = defaultdict(list)
    for c in literals:
        grouped[c['vertex']].append(c)
    for cs in grouped.values():
        cs.sort(key=lambda c: c['index'])
    order = sorted(grouped, key=lambda v: (-mu.get(v, F()), -len(grouped[v]), v))
    key = sha256(dump([input_digest, name, initial['face_multipliers']]).encode()).hexdigest()
    existing = {}
    if saved is None and PROGRESS.exists():
        for line in PROGRESS.read_text().splitlines():
            rec = json.loads(line)
            if rec['name'] == name:
                assert rec['key'] == key
                assert rec['vertex'] not in existing
                existing[rec['vertex']] = rec['record']
    saved_rows = {} if saved is None else {r['vertex']: r for r in saved['local_cones']}
    repaired = dict(mu)
    records, cuts, open_cases = [], [], []
    start = time.monotonic()
    for serial, v in enumerate(order):
        ws = [weights[c['index']] for c in grouped[v]]
        try:
            pack = L.Packings(grouped[v])
        except L.ResourceStop as exc:
            rec = dict(vertex=v, forest=d['budgets'][v]['forest'], status='OPEN_ENUMERATION', reason=str(exc))
            open_cases.append(v)
            records.append(rec)
            if saved is not None:
                assert rec == saved_rows[v]
            elif v not in existing:
                with PROGRESS.open('a') as stream:
                    stream.write(dump(dict(name=name, key=key, vertex=v, record=rec)) + '\n')
            continue
        if saved is not None:
            certificate = saved_rows[v]['certificate']
        elif v in existing:
            assert 'certificate' in existing[v], 'Previously open forest requires a fresh successor.'
            certificate = existing[v]['certificate']
        else:
            scale = max(mu.get(v, F()), max(ws))
            assert scale > 0
            certificate = L.solve_and_recover(pack, ws, scale)
        lower, upper, slack = L.verify_certificate(pack, ws, certificate)
        repaired[v] = max(mu.get(v, F()), upper)
        if certificate['status'].startswith('OPEN'):
            open_cases.append(v)
        if lower > mu.get(v, F()):
            q = {j: F(a) for j, a in certificate['charges']}
            terms = []
            for f in features_by_vertex[v]:
                value = sum((m * q.get(j, F()) for j, m in f['block_orbits']), F())
                if value:
                    terms.append([f['face'], str(value)])
            assert sum((alpha[i] * F(a) for i, a in terms), F()) == lower
            cuts.append(dict(index=len(cuts), vertex=v, vertex_cost=1,
                positive_simplex_coefficients=terms, exact_load=str(lower), old_vertex_mass=str(mu.get(v, F())),
                excess=str(lower - mu.get(v, F())), orbit_charges=certificate['charges']))
        rec = dict(vertex=v, forest=d['budgets'][v]['forest'], status=certificate['status'],
            counts=dict(block_orbits=len(pack.ids), literal_blocks=len(pack.masks),
                packing_count_vectors=len(pack.vectors), cached_vector_entries=pack.entries,
                cached_available_masks=pack.dp.cache_info().currsize,
                positive_mixture_packings=len(certificate['packing_mixture'])),
            block_orbit_marginals=[[j, str(a)] for j, a in zip(pack.ids, ws)],
            literal_block_orbits=grouped[v], count_vector_digest=sha256(dump(pack.vectors).encode()).hexdigest(),
            certificate=certificate, primal_cover_slacks=slack, old_vertex_mass=str(mu.get(v, F())),
            repaired_vertex_mass=str(repaired[v]), lower_excess=str(max(F(), lower - mu.get(v, F()))),
            cover_mass_increase=str(repaired[v] - mu.get(v, F())))
        records.append(rec)
        if saved is None and v not in existing:
            with PROGRESS.open('a') as stream:
                stream.write(dump(dict(name=name, key=key, vertex=v, record=rec)) + '\n')
        if serial % 40 == 0 or lower > mu.get(v, F()):
            print('RESTRICTED_LOCAL_CONES', name, serial + 1, '/', len(order),
                  'cuts', len(cuts), 'open', len(open_cases), 'seconds', round(time.monotonic() - start, 2), flush=True)
        pack.dp.cache_clear()
    inherited = verify_inherited(d, alpha, repaired)
    assert all(repaired.get(v, F()) >= a for v, a in mu.items())
    num = sum((a * d['faces'][i]['required'] for i, a in alpha.items()), F())
    assert num == F(initial['objective_numerator'])
    den = sum(repaired.values(), F())
    complete = not open_cases
    result = dict(name=name, status='PASS_COMPLETE_RESTRICTED_LOCAL_CONE_COVER' if complete else 'PASS_LOCAL_CONES_WITH_OPEN_CASES',
        dimension=13, target=12, input_bound=initial['bound'], complete_positive_simplex_forest_coverage=complete,
        input_face_mass_digest=sha256(dump(initial['face_multipliers']).encode()).hexdigest(),
        inherited_columns_checked=inherited, simplex_dimension_counts=dimensions,
        counts=dict(positive_faces=len(alpha), positive_simplex_faces=len(features),
            positive_quad_faces=len(alpha)-len(features), positive_vertex_forests=len(mu),
            positive_simplex_forests=len(order), occupied_block_orbits=len(columns),
            new_positive_simplex_faces_outside_old_full_block_export=len({f['face'] for f in features}-previous_features),
            corner_occurrences=corners, open_forests=len(open_cases), processed_forests=len(records),
            exact_local_optima=sum(r['status']=='EXACT_LOCAL_OPTIMUM' for r in records),
            exact_cover_brackets=sum(r['status']=='EXACT_COVER_BRACKET' for r in records),
            packing_count_vectors=sum(r.get('counts',{}).get('packing_count_vectors',0) for r in records),
            positive_mixture_packings=sum(r.get('counts',{}).get('positive_mixture_packings',0) for r in records),
            new_exact_cut_forests=len(cuts), increased_vertex_masses=sum(repaired[v]>mu.get(v,F()) for v in repaired)),
        generation_local_solver_invocations=sum('certificate' in r for r in records),
        reused_old_fixed_alpha_certificates=0, saved_replay_solver_invocations=0,
        priority_order=order, open_forests=open_cases, block_orbit_columns=columns,
        positive_simplex_features=features, literal_block_orbits=literals,
        local_cones=records, new_exact_charge_columns=cuts,
        face_multipliers=initial['face_multipliers'], initial_vertex_multipliers=initial['vertex_multipliers'],
        vertex_multipliers=[[v,str(a)] for v,a in sorted(repaired.items()) if a],
        objective_numerator=str(num), normalizer=str(den), bound=str(num/den),
        excludes_target12=complete and num/den>12)
    print('RESTRICTED_CONE_RESULT', name, result['counts'], result['bound'], float(num/den), flush=True)
    return result


def run(saved=None):
    support, model, prior = read(SUPPORT), read(N.TAG), read(A.TAG)
    d, prepared = N.build()
    assert json.loads(dump(prepared)) == model['prepared']
    assert len(d['faces']) == 82806 and len(d['zeros']) == 94481
    paths = [f'scripts/{TAG}.py', f'scripts/{SUPPORT}.py', f'results/{SUPPORT}.json',
        f'scripts/{N.TAG}.py', f'results/{N.TAG}.json', f'scripts/{A.TAG}.py', f'results/{A.TAG}.json',
        f'scripts/{L.TAG}.py', f'scripts/{L.H.TAG}.py',
        'results/fresh_simplex_disjoint_ports_20261017.md']
    sources = {p: sha(p) for p in paths}
    input_digest = sha256(dump({p:h for p,h in sources.items() if p!=f'scripts/{TAG}.py'}).encode()).hexdigest()
    previous = {f['face'] for f in prior['positive_simplex_features']}
    tests = {t['name']:t for t in support['tests']}
    outputs = []
    saved_tests = {} if saved is None else {t['name']:t for t in saved['tests']}
    name = 'legacy_export_seed'
    outputs.append(run_case(d, tests[name], input_digest, previous, saved_tests.get(name)))
    if outputs[0]['complete_positive_simplex_forest_coverage'] and F(outputs[0]['bound']) <= 12:
        name = 'remove_unselected_Q_keep_all_S'
        outputs.append(run_case(d, tests[name], input_digest, previous, saved_tests.get(name)))
    return dict(status='PASS_RESTRICTED_FLAG_DUAL_COMPLETE_LOCAL_CONES' if all(t['complete_positive_simplex_forest_coverage'] for t in outputs)
                else 'PASS_RESTRICTED_FLAG_DUAL_LOCAL_CONES_OPEN_CASES',
        sources=sources, input_model_sha256=prepared['model_sha256'], dimension=13, target=12,
        global_optimizer_invocations=0, saved_replay_optimizer_invocations=0,
        resource_limits=dict(cached_vector_entries_per_forest=L.ENTRY_LIMIT,
            enumeration_seconds_per_forest=L.ENUM_SECONDS, local_solver_seconds=L.SOLVE_SECONDS),
        tests=outputs,
        theorem='For each finite forest, explicit rational mixtures of disjoint literal simplex-port blocks cover every orbit marginal at the NEW face masses. Averaging each whole packing over forest automorphisms gives uniform literal marginals. For any nonnegative orbit-local charge q, marginal load is bounded by mixture mass times the maximum disjoint-packing cost of q. Thus the repaired vertex masses simultaneously satisfy the whole nonnegative simplex-block charge cone, not just the saved separating charges. Full inherited positive and signed columns are checked directly; face and actual-provenance masses remain fixed.',
        scope='All82806 inherited face obligations and94481 signed columns are preserved in the verified model. No global optimization, new face enumeration, old-alpha cover reuse, or primal feasibility claim. This is an exact dual obstruction only when the complete repaired ratio exceeds12. It does not cover joint geometric restrictions on simplex/prism configurations, omitted face domains, or an ambient dimension bound.')


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    a=p.parse_args()
    out=ROOT/f'results/{TAG}.json'
    saved=json.loads(out.read_text()) if a.verify else None
    result=json.loads(dump(run(saved)))
    if a.verify:
        assert result==saved
    else:
        out.write_text(dump(result)+'\n')
    print(result['status'], [(t['name'],t['bound']) for t in result['tests']], flush=True)
