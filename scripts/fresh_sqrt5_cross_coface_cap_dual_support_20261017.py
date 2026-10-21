#!/usr/bin/env python3
"""Exact support restrictions of the new cross-coface dual; no optimizer."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import fresh_sqrt5_cross_coface_cap_full_model_20261017 as N
import fresh_sqrt5_all_subset_port_packing_20261017 as A

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_cross_coface_cap_dual_support_20261017'
LOCAL = 'fresh_sqrt5_restricted_flag_local_cones_20261017'
JOINT = 'fresh_sqrt5_joint_five_port_diagnostic_20261017'


def run():
    saved = N.read(N.TAG)
    d, prepared = N.build()
    assert json.loads(json.dumps(prepared)) == saved['prepared']
    original = next(c for c in saved['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_DUAL')
    N.N.verify(d, original)
    alpha = {i: F(a) for i, a in original['face_multipliers']}
    local, joint, caps = map(N.read, [LOCAL, JOINT, N.CAP])
    old = next(t for t in local['tests'] if t['name'] == 'legacy_export_seed')
    old_s = {r['face'] for r in old['positive_simplex_features']}
    old_q = {r['face'] for r in joint['actual_quad_features']}
    pmap = dict(prepared['provenance_to_face']); inverse = {fi: p for p, fi in pmap.items()}
    assert len(pmap) == len(inverse) == 1858 and len(old_s) == 1864 and len(old_q) == 355
    all_q = {f['index'] for f in d['faces'] if f['kind'] == 'canonical_dotted_quad_aggregate'}
    all_s = set(range(len(d['faces'])))-all_q
    assert old_s <= all_s and old_q <= set(inverse) <= all_q
    groups = Counter(); demand = defaultdict(F); rows = []
    for fi, a in sorted(alpha.items()):
        if fi in all_q:
            group = 'Q_previously_joint_scanned_actual' if fi in old_q else \
                    'Q_other_selected_actual' if fi in inverse else 'Q_unselected_merged'
        else:
            group = 'S_previous_complete_local_feature' if fi in old_s else 'S_needs_new_complete_local_feature'
        groups[group] += 1; demand[group] += a*d['faces'][fi]['required']
        rows.append(dict(face=fi, kind=d['faces'][fi]['kind'], mass=str(a), group=group,
            actual_provenance=inverse.get(fi), objective_contribution=str(a*d['faces'][fi]['required'])))
    assert sum(demand.values(), F()) == F(original['objective_numerator'])
    cap_rows = []
    for c in caps['columns']:
        load = sum((m*alpha.get(i, F()) for i, m in c['positive_face_coefficients']), F())
        capacity = c['capacity']*alpha.get(c['simplex_face'], F())
        assert load <= capacity
        cap_rows.append(dict(index=c['index'], base_type=c['base_type'], simplex_face=c['simplex_face'],
            positive_face_coefficients=c['positive_face_coefficients'], load=str(load), capacity=str(capacity),
            slack=str(capacity-load), positive=load > 0, tight=load == capacity))
    print('NEW_DUAL_SUPPORT', dict(groups), 'CAPS', sum(r['positive'] for r in cap_rows),
          sum(r['positive'] and r['tight'] for r in cap_rows), flush=True)
    cases = [('all_S_and_selected_actual_Q', all_s | set(inverse)),
             ('old_complete_S_and_selected_actual_Q', old_s | set(inverse)),
             ('old_complete_joint_domain', old_s | old_q),
             ('simplex_only', all_s)]
    tests = []
    for name, domain in cases:
        initial = {i: a for i, a in alpha.items() if i in domain}
        cert = N.N.M.dual_certificate(d, dict(initial)); N.N.verify(d, cert)
        repaired = {i: F(a) for i, a in cert['face_multipliers']}
        assert all(repaired.get(i, F()) >= a for i, a in initial.items())
        assert all(repaired.get(i, F()) == initial.get(i, F()) for i in all_q)
        restored = [i for i, a in repaired.items() if a > initial.get(i, F())]
        assert set(restored) <= all_s
        record = dict(name=name, seed_faces=len(initial), seed_Q=sum(i in all_q for i in initial),
            seed_S=sum(i in all_s for i in initial), increased_or_restored_S=len(restored),
            restored_S_outside_seed_domain=[i for i in repaired if i not in domain],
            S_missing_old_complete_features=[i for i in repaired if i in all_s and i not in old_s],
            Q_missing_old_joint_features=[i for i in repaired if i in all_q and i not in old_q],
            simplex_mass_changes=[[i, str(initial.get(i, F())), str(repaired[i])] for i in restored],
            repaired_dual=cert)
        tests.append(record)
        print('RESTRICTED_NEW_DUAL', name, len(repaired), 'RESTORED', len(restored),
              'OUTSIDE', len(record['restored_S_outside_seed_domain']), float(F(cert['bound'])), flush=True)
    # Full fresh feature reconstruction is independent of whether a face was
    # positive in any earlier dual. Literal orbits/covers are a successor step.
    columns, features, dimensions, corners = A.corner_features(d, alpha)
    assert {f['face'] for f in features} == set(alpha)&all_s
    weights = defaultdict(F)
    for f in features:
        assert sum(m for j, m in f['block_orbits']) == f['simplex_dimension']+1
        for j, m in f['block_orbits']: weights[j] += m*alpha[f['face']]
    assert set(weights) == {c['index'] for c in columns}
    paths = [f'scripts/{TAG}.py', f'scripts/{N.TAG}.py', f'results/{N.TAG}.json',
             f'scripts/{A.TAG}.py', f'results/{LOCAL}.json', f'results/{JOINT}.json',
             f'results/{N.CAP}.json']
    return dict(status='PASS_EXACT_NEW_CROSS_COFACE_DUAL_SUPPORT_AND_RESTRICTIONS',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        inherited_model_sha256=prepared['model_sha256'], optimizer_invocations=0,
        full_bound=original['bound'], positive_support_counts=dict(groups),
        objective_by_group={k: str(a) for k, a in demand.items()}, positive_face_support=rows,
        counts=dict(positive_faces=len(alpha), positive_forests=len(original['vertex_multipliers']),
            selected_actual_Q=sum(i in inverse for i in alpha), old_joint_actual_Q=sum(i in old_q for i in alpha),
            cap_columns=len(cap_rows), positive_cap_columns=sum(r['positive'] for r in cap_rows),
            positive_tight_caps=sum(r['positive'] and r['tight'] for r in cap_rows),
            all_face_obligations=len(d['faces']), all_signed_columns=len(d['zeros']),
            full_simplex_faces=len(features), fresh_block_orbits=len(columns), fresh_simplex_corners=corners),
        cap_loads=cap_rows, tests=tests, fresh_simplex_dimension_counts=dimensions,
        full_positive_simplex_features=features, full_block_orbits=columns,
        full_block_orbit_marginals=[[j, str(a)] for j, a in sorted(weights.items())],
        scope='Exact restrictions retain all94722 signed constraints through the new full-model repair. Restored simplex rows outside seed domains are explicitly listed. Full new simplex blocks are reconstructed but no changed-alpha packing cover or complete joint scan is claimed here. Unselected merged Q rows are never assigned a convenient actual provenance. No optimizer or dimension assertion.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true'); a = p.parse_args()
    result = json.loads(json.dumps(run())); out = ROOT/f'results/{TAG}.json'
    if a.verify: assert result == json.loads(out.read_text())
    else: out.write_text(json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n')
    print(result['status'], result['counts'], flush=True)
