#!/usr/bin/env python3
"""One bounded successor: cross-coface capacities and saved joint cuts."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import lcm
import json
from pathlib import Path

import fresh_sqrt5_flag_local_packing_full_model_20261017 as N

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_cross_coface_cap_full_model_20261017'
CAP = 'fresh_sqrt5_cross_coface_cap_capacity_20261017'
LOCAL = 'fresh_sqrt5_restricted_flag_local_cones_20261017'
THREE = 'fresh_sqrt5_joint_triangle_square_ports_20261017'
FOUR = 'fresh_sqrt5_joint_four_port_diagnostic_20261017'
FIVE = 'fresh_sqrt5_joint_five_port_diagnostic_20261017'
OUT = ROOT/f'results/{TAG}.json'


def read(tag):
    obj = json.loads((ROOT/f'results/{tag}.json').read_text())
    for doc in (obj, obj.get('prepared', {})):
        for p, h in doc.get('sources', {}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    return obj


def build():
    d, inherited = N.build()
    assert json.loads(json.dumps(inherited)) == read(N.TAG)['prepared']
    caps, local, three, four, five = map(read, [CAP, LOCAL, THREE, FOUR, FIVE])
    t = next(x for x in local['tests'] if x['name'] == 'legacy_export_seed')
    assert t['face_multipliers'] == three['face_multipliers'] == four['face_multipliers'] == five['face_multipliers']
    oldvars, oldfaces, oldzeros = len(d['variables']), len(d['faces']), len(d['zeros'])
    assert (oldvars, oldfaces, oldzeros) == (428310, 82806, 94481)
    hashes = (N.digest(d['variables']), N.digest(d['zeros']))
    original_faces = [[list(x) for x in f['coefficients']] for f in d['faces']]
    original_budgets = [[list(x) for x in b['coefficients']] for b in d['budgets']]
    old_meta = [(f['index'], f['kind'], f['required']) for f in d['faces']]
    oldactive = set(d['active']); pmap = dict(inherited['provenance_to_face'])
    families = [('restricted_local_packing_cut', LOCAL, t['new_exact_charge_columns'], 31),
                ('joint_three_port_cut', THREE, three['critical_columns'], 24),
                ('joint_four_port_cut', FOUR, four['critical_columns'], 47),
                ('joint_five_port_cut', FIVE, five['critical_columns'], 115)]
    blocks = []; scalings = []
    for kind, source, columns, expected in families:
        assert len(columns) == expected
        start = len(d['variables'])
        for i, c in enumerate(columns):
            assert c['index'] == i
            entries = c.get('positive_face_coefficients', c.get('positive_simplex_coefficients'))
            assert entries
            values = [(fi, F(a)) for fi, a in entries]
            cost = F(c['vertex_cost']); assert cost > 0 and all(a > 0 for fi, a in values)
            assert len({fi for fi, a in values}) == len(values)
            den = lcm(cost.denominator, *(a.denominator for fi, a in values))
            v, j = c['vertex'], len(d['variables'])
            scaled_cost = int(cost*den); assert F(scaled_cost, den) == cost
            d['variables'].append(dict(index=j, kind=kind, vertex=v, cost=scaled_cost, source_index=i))
            d['budgets'][v]['coefficients'].append([j, scaled_cost])
            for fi, a in values:
                assert 0 <= fi < oldfaces
                if kind == 'restricted_local_packing_cut':
                    assert d['faces'][fi]['kind'] in ('maximal_simplex_aggregate', 'maximal_edge_incidence')
                else:
                    assert d['faces'][fi]['kind'] in ('maximal_simplex_aggregate', 'maximal_edge_incidence', 'canonical_dotted_quad_aggregate')
                m = int(a*den); assert F(m, den) == a
                d['faces'][fi]['coefficients'].append([j, m])
            scalings.append(dict(variable=j, kind=kind, source_index=i, vertex=v,
                positive_denominator=den, original_cost=str(cost), scaled_cost=scaled_cost,
                exported_face_count=len(values)))
        blocks.append(dict(kind=kind, source=source, start=start, count=expected))
    start = len(d['variables'])
    assert len(caps['columns']) == 241
    for i, c in enumerate(caps['columns']):
        assert c['index'] == i and c['capacity'] == 1
        fi = c['simplex_face']; assert d['faces'][fi]['kind'] in ('maximal_simplex_aggregate', 'maximal_edge_incidence')
        positive = [[pmap[p], m] for p, m in c['provenance_coefficients']]
        assert positive == c['positive_face_coefficients']
        assert len({q for q, m in positive}) == len(positive)
        j = len(d['variables'])
        d['variables'].append(dict(index=j, kind='cross_coface_A1_cap_transfer', vertex=None, cost=0, source_index=i))
        d['faces'][fi]['coefficients'].append([j, -1])
        for q, m in positive:
            assert isinstance(m, int) and m > 0 and d['faces'][q]['kind'] == 'canonical_dotted_quad_aggregate'
            d['faces'][q]['coefficients'].append([j, m])
        d['zeros'].append(dict(variable=j, simplex_face=fi, capacity=1, positive_coefficients=positive))
    blocks.append(dict(kind='cross_coface_A1_cap_transfer', source=CAP, start=start, count=241))
    assert len(d['faces']) == oldfaces and len(d['budgets']) == 6558
    assert N.digest(d['variables'][:oldvars]) == hashes[0]
    assert N.digest(d['zeros'][:oldzeros]) == hashes[1]
    for i, face in enumerate(d['faces']):
        assert (face['index'], face['kind'], face['required']) == old_meta[i]
        co = Counter()
        for j, m in face['coefficients']: co[j] += m
        face['coefficients'] = [[j, m] for j, m in sorted(co.items()) if m]
        assert [[j, m] for j, m in face['coefficients'] if j < oldvars] == original_faces[i]
        assert all(m > 0 for j, m in face['coefficients'] if d['variables'][j]['vertex'] is not None)
    for v, b in enumerate(d['budgets']):
        assert [[j, m] for j, m in b['coefficients'] if j < oldvars] == original_budgets[v]
    assert len(d['variables']) == 428768 and len(d['zeros']) == 94722
    d['active'] = sorted({j for f in d['faces'] for j, m in f['coefficients'] if m > 0})
    assert oldactive <= set(d['active'])
    eligible = set(d['repair_indices'])
    assert all(sum(m for j, m in f['coefficients'] if j in eligible) > 0 for f in d['faces'])
    h = sha256()
    for group in (d['variables'], d['faces'], d['budgets'], d['zeros']):
        h.update(N.digest(group).encode()); h.update(b'\n')
    paths = [f'scripts/{TAG}.py', f'scripts/{N.TAG}.py', f'results/{N.TAG}.json']
    for source in (CAP, LOCAL, THREE, FOUR, FIVE):
        paths += [f'scripts/{source}.py', f'results/{source}.json', f'results/{source}.md']
    prepared = dict(sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        dimension=13, target=12, inherited_model_sha256=inherited['model_sha256'],
        model_sha256=h.hexdigest(), provenance_to_face=inherited['provenance_to_face'],
        counts=dict(faces=len(d['faces']), vertex_forests=len(d['budgets']), variables=len(d['variables']),
            inherited_variables=oldvars, active_variables=len(d['active']), signed_variables=len(d['zeros']),
            new_signed_variables=241, new_finite_positive_cuts=217, local_cuts=31,
            triple_cuts=24, quartet_cuts=47, quintet_cuts=115,
            face_matrix_nonzeros=sum(len(f['coefficients']) for f in d['faces'])),
        added_column_blocks=blocks, exact_positive_column_scalings=scalings,
        preservation='All82806 inherited face obligations,6558 complete vertex budgets,428310 old variables and94481 old signed transfers are retained coefficient by coefficient. Every new cap positive is checked against the same actual-provenance map.',
        charge_scope='Adds241 cross-coface cap capacities and exactly217 saved finite positive cuts:31 restricted-local,24 strengthened triple,47 four-port and115 five-port. New positive coefficients are zero outside their exported full-star classes. No entire local cone or all joint orbit family is asserted at a newly optimized alpha.',
        recovery='Uses inherited one-call highs-ipm with300-second limit and exact rational primal/dual recovery. All signed columns have positive Q terms and a negative S term, so raising simplex masses repairs signed dual loads monotonically. Every face retains positive inherited primal-repair coefficients. Inactive columns have no positive face entries and may be fixed to zero.',
        scope='A genuinely enlarged but partial face/forest and finite-charge model. A successful target12 primal would still require the missing full geometric face domains and complete charge separation before any dimension theorem.')
    return d, prepared


if __name__ == '__main__':
    p = argparse.ArgumentParser(); mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true'); mode.add_argument('--solve', action='store_true'); mode.add_argument('--verify', action='store_true')
    a = p.parse_args(); d, prepared = build(); prepared = json.loads(json.dumps(prepared))
    print('PREPARED_CROSS_COFACE_CAP_MODEL', prepared['counts'], prepared['model_sha256'], flush=True)
    if a.verify:
        saved = json.loads(OUT.read_text()); assert saved['prepared'] == prepared
        for cert in saved.get('certificates', []): N.verify(d, cert)
        print('PASS_SAVED_CROSS_COFACE_CAP_MODEL', saved['status'], flush=True)
    elif a.prepare:
        if OUT.exists(): assert json.loads(OUT.read_text())['solver_invocations'] == 0
        OUT.write_text(json.dumps(dict(status='PREPARED_CROSS_COFACE_CAP_MODEL', prepared=prepared,
            solver_invocations=0, solver=None, certificates=[]), sort_keys=True, separators=(',', ':'))+'\n')
    else:
        saved = json.loads(OUT.read_text()); assert saved['prepared'] == prepared and saved['solver_invocations'] == 0
        saved.update(status='ONE_BOUNDED_SOLVE_STARTED', solver_invocations=1)
        OUT.write_text(json.dumps(saved, sort_keys=True, separators=(',', ':'))+'\n')
        solver, certificates = N.M.optimize(d)
        for cert in certificates: N.verify(d, cert)
        status = 'PASS_EXACT_PARTIAL_PRIMAL_TARGET12' if any(c.get('meets_target12') for c in certificates) else \
                 'PASS_EXACT_CROSS_COFACE_CAP_MODEL_CERTIFICATES' if certificates else 'UNRESOLVED_SINGLE_BOUNDED_SOLVE'
        saved.update(status=status, solver=solver, certificates=certificates)
        OUT.write_text(json.dumps(saved, sort_keys=True, separators=(',', ':'))+'\n')
        print(status, [(c['kind'], c.get('maximum', c.get('bound'))) for c in certificates], flush=True)
