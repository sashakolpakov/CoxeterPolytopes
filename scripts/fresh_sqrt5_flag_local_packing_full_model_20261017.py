#!/usr/bin/env python3
"""One new bounded full-model solve: exact flags and 175 finite packing cuts.

All inherited rows, variables, budgets and transfers are retained. This is
an explicitly finite selected-charge model, not the entire local cone.
"""
import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction as F
from math import lcm
from pathlib import Path

import fresh_sqrt5_labelled_packing_full_model_20261017 as M

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_flag_local_packing_full_model_20261017'
OUT = ROOT / f'results/{TAG}.json'
FLAGS = 'fresh_sqrt5_flag_refined_base_capacity_20261017'
SUBSETS = 'fresh_sqrt5_all_subset_port_packing_20261017'
CLIQUES = 'fresh_sqrt5_intersecting_block_cliques_20261017'
LOCAL = 'fresh_sqrt5_all_local_simplex_packing_cones_20261017'


def digest(records):
    h = hashlib.sha256()
    for r in records:
        h.update(M.R.B.dump(r).encode())
        h.update(b'\n')
    return h.hexdigest()


def build():
    d, inherited = M.build()
    saved = json.loads((ROOT / f'results/{M.TAG}.json').read_text())
    assert json.loads(json.dumps(inherited)) == saved['prepared']
    flags = M.R.P.read(FLAGS)
    subset = M.R.P.read(SUBSETS)
    clique = M.R.P.read(CLIQUES)
    local = M.R.P.read(LOCAL)
    assert flags['inherited_model_sha256'] == inherited['model_sha256']
    assert len(flags['flag_transfer_columns']) == 77391
    oldvars, oldfaces, oldzeros = len(d['variables']), len(d['faces']), len(d['zeros'])
    assert (oldvars, oldfaces, oldzeros) == (350744, 82806, 17090)
    oldvarhash = digest(d['variables'])
    oldzerohash = digest(d['zeros'])
    original_faces = [[list(x) for x in f['coefficients']] for f in d['faces']]
    original_budgets = [[list(x) for x in b['coefficients']] for b in d['budgets']]
    old_meta = [(f['index'], f['kind'], f['required']) for f in d['faces']]
    oldactive = set(d['active'])
    pmap = dict(inherited['provenance_to_face'])
    blocks, scalings = [], []
    families = [
        ('all_subset_critical', SUBSETS, subset['critical_columns'], 148),
        ('intersecting_orbit_clique', CLIQUES, clique['critical_columns'], 4),
        ('complete_local_cone_cut', LOCAL, local['new_exact_charge_columns'], 23),
    ]
    for kind, source, columns, expected in families:
        assert len(columns) == expected
        start = len(d['variables'])
        for i, c in enumerate(columns):
            assert c['index'] == i
            values = [(fi, F(a)) for fi, a in c['positive_simplex_coefficients']]
            cost = F(c['vertex_cost'])
            assert cost > 0 and values and all(a > 0 for fi, a in values)
            assert len({fi for fi, a in values}) == len(values)
            den = lcm(cost.denominator, *(a.denominator for fi, a in values))
            j, v = len(d['variables']), c['vertex']
            scaled_cost = int(cost * den)
            assert F(scaled_cost, den) == cost and scaled_cost > 0
            d['variables'].append(dict(index=j, kind=kind, vertex=v, cost=scaled_cost, source_index=i))
            d['budgets'][v]['coefficients'].append([j, scaled_cost])
            for fi, a in values:
                assert d['faces'][fi]['kind'] in ('maximal_simplex_aggregate', 'maximal_edge_incidence')
                m = int(a * den)
                assert F(m, den) == a and m > 0
                d['faces'][fi]['coefficients'].append([j, m])
            scalings.append(dict(variable=j, kind=kind, source_index=i, vertex=v,
                positive_denominator=den, scaled_vertex_cost=scaled_cost,
                original_vertex_cost=str(cost), source_feature_rows=len(values)))
        blocks.append(dict(kind=kind, source=source, start=start, count=expected))
    start = len(d['variables'])
    positive_flag_terms = 0
    for i, c in enumerate(flags['flag_transfer_columns']):
        assert c['index'] == i
        j, fi, cap = len(d['variables']), c['simplex_face'], c['capacity']
        assert isinstance(cap, int) and cap > 0
        assert d['faces'][fi]['kind'] in ('maximal_simplex_aggregate', 'maximal_edge_incidence')
        positive = [[pmap[p], m] for p, m in c['provenance_coefficients']]
        assert positive == c['positive_face_coefficients']
        assert len({q for q, m in positive}) == len(positive)
        d['variables'].append(dict(index=j, kind='simplex_vertex_coface_flag_transfer',
                                   vertex=None, cost=0, source_index=i))
        d['faces'][fi]['coefficients'].append([j, -cap])
        for q, m in positive:
            assert isinstance(m, int) and m > 0
            assert d['faces'][q]['kind'] == 'canonical_dotted_quad_aggregate'
            d['faces'][q]['coefficients'].append([j, m])
            positive_flag_terms += 1
        d['zeros'].append(dict(variable=j, simplex_face=fi, capacity=cap,
                               positive_coefficients=positive))
    blocks.append(dict(kind='simplex_vertex_coface_flag_transfer', source=FLAGS,
                       start=start, count=len(flags['flag_transfer_columns'])))

    # Preserve the entire inherited matrix, not merely its positive support.
    assert len(d['faces']) == oldfaces and len(d['budgets']) == 6558
    assert digest(d['variables'][:oldvars]) == oldvarhash
    assert digest(d['zeros'][:oldzeros]) == oldzerohash
    for fi, f in enumerate(d['faces']):
        assert (f['index'], f['kind'], f['required']) == old_meta[fi]
        cc = Counter()
        for j, m in f['coefficients']:
            cc[j] += m
        f['coefficients'] = [[j, m] for j, m in sorted(cc.items()) if m]
        assert [[j, m] for j, m in f['coefficients'] if j < oldvars] == original_faces[fi]
        assert all(m > 0 for j, m in f['coefficients'] if d['variables'][j]['vertex'] is not None)
    for v, b in enumerate(d['budgets']):
        assert [[j, m] for j, m in b['coefficients'] if j < oldvars] == original_budgets[v]
    assert len(d['variables']) == 428310 and len(d['zeros']) == 94481
    d['active'] = sorted({j for f in d['faces'] for j, m in f['coefficients'] if m > 0})
    assert oldactive <= set(d['active'])
    eligible = set(d['repair_indices'])
    assert all(sum(m for j, m in f['coefficients'] if j in eligible) > 0 for f in d['faces'])
    h = hashlib.sha256()
    for group in (d['variables'], d['faces'], d['budgets'], d['zeros']):
        h.update(digest(group).encode())
        h.update(b'\n')
    paths = [f'scripts/{TAG}.py', f'scripts/{M.TAG}.py', f'results/{M.TAG}.json']
    for source in (FLAGS, SUBSETS, CLIQUES, LOCAL):
        paths += [f'scripts/{source}.py', f'results/{source}.json', f'results/{source}.md']
    prepared = dict(sources={p: M.R.P.sha(p) for p in paths}, dimension=13, target=12,
        inherited_model_sha256=inherited['model_sha256'], model_sha256=h.hexdigest(),
        counts=dict(faces=len(d['faces']), inherited_faces=oldfaces,
            vertex_forests=len(d['budgets']), variables=len(d['variables']), inherited_variables=oldvars,
            active_variables=len(d['active']), zero_only_variables=len(d['variables'])-len(d['active']),
            signed_variables=len(d['zeros']), new_flag_signed_variables=77391,
            new_finite_positive_cuts=175, all_subset_cuts=148, clique_cuts=4, local_cone_cuts=23,
            new_flag_positive_entries=positive_flag_terms,
            face_matrix_nonzeros=sum(len(f['coefficients']) for f in d['faces'])),
        added_column_blocks=blocks, exact_positive_column_scalings=scalings,
        provenance_to_face=inherited['provenance_to_face'],
        preservation='Every inherited variable, all82806 face obligations,6558 complete vertex budgets and17090 signed columns retain their exact coefficients and indices. New rational positive columns are scaled by a positive common denominator in both face entries and their vertex cost. Every flag positive entry is checked against the inherited1858 actual-provenance row map.',
        charge_scope='Finite successor:77391 flag signed columns and exactly148 previously maximizing all-subset columns,4 intersecting-orbit clique columns,23 new local-packing cut columns. It does not include every subset column or the entire local packing cone. New finite positive features are restricted to their exported simplex classes, with zero quad coefficients. New flag positives are restricted to the1858 selected actual quad provenances; their negative simplex terms are retained.',
        repair='All signed columns have positive quad coefficients and exactly one negative simplex coefficient. Raising simplex dual masses repairs their inequalities monotonically. Exact primal repair raises only inherited eligible positive pair/port variables; all face rows have a checked such entry. Inactive columns have no positive face entries and can be set to zero.',
        scope='Preserved partial face/forest domain with a strictly enlarged incidence family. Neither a numerical nor exact budget12 primal is a dimension theorem without completing the outstanding geometric domains. No certificate for changed alpha is borrowed from the old fixed-alpha local covers.')
    return d, prepared


def verify(d, certificate):
    if certificate['kind'] == 'EXACT_PARTIAL_MODEL_PRIMAL':
        entries = certificate['rounded_variable_numerators']
        assert len(dict(entries)) == len(entries)
        active = set(d['active'])
        assert all(j in active and isinstance(a, int) and a >= 0 for j, a in entries)
        maximum = M.verify_primal(d, certificate)
        assert certificate['meets_target12'] == (maximum <= 12)
    else:
        assert certificate['kind'] == 'EXACT_PARTIAL_MODEL_DUAL'
        for key, size in [('face_multipliers', len(d['faces'])), ('vertex_multipliers', len(d['budgets']))]:
            values = certificate[key]
            assert len(dict(values)) == len(values)
            assert all(0 <= i < size and F(a) > 0 for i, a in values)
        M.verify_dual(d, certificate)
        assert certificate['excludes_target12'] == (F(certificate['bound']) > 12)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true')
    mode.add_argument('--solve', action='store_true')
    mode.add_argument('--verify', action='store_true')
    a = p.parse_args()
    d, prepared = build()
    prepared = json.loads(json.dumps(prepared))
    print('PREPARED_FLAG_LOCAL_PACKING_MODEL', prepared['counts'], prepared['model_sha256'], flush=True)
    if a.verify:
        saved = json.loads(OUT.read_text())
        assert saved['prepared'] == prepared
        for c in saved.get('certificates', []):
            verify(d, c)
        print('PASS_SAVED_FLAG_LOCAL_PACKING_MODEL', saved['status'], flush=True)
    elif a.prepare:
        if OUT.exists():
            assert json.loads(OUT.read_text())['solver_invocations'] == 0
        OUT.write_text(M.R.B.dump(dict(status='PREPARED_FLAG_LOCAL_PACKING_MODEL', prepared=prepared,
            solver_invocations=0, solver=None, certificates=[])) + '\n')
    else:
        saved = json.loads(OUT.read_text())
        assert saved['prepared'] == prepared and saved['solver_invocations'] == 0
        saved.update(status='ONE_BOUNDED_SOLVE_STARTED', solver_invocations=1)
        OUT.write_text(M.R.B.dump(saved) + '\n')
        solver, certificates = M.optimize(d)
        for c in certificates:
            verify(d, c)
        status = 'PASS_EXACT_PARTIAL_PRIMAL_TARGET12' if any(c.get('meets_target12') for c in certificates) else \
            'PASS_EXACT_FLAG_LOCAL_PACKING_CERTIFICATES' if certificates else 'UNRESOLVED_SINGLE_BOUNDED_SOLVE'
        saved.update(status=status, solver=solver, certificates=certificates)
        OUT.write_text(M.R.B.dump(saved) + '\n')
        print(status, [(c['kind'], c.get('maximum', c.get('bound'))) for c in certificates], flush=True)
