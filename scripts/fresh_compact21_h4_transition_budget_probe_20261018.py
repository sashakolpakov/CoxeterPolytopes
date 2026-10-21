#!/usr/bin/env python3
"""Stream the saved primal and evaluate five explicit forests; no optimizer."""
from collections import Counter
import argparse
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
from itertools import combinations
import json
import os
from pathlib import Path
import resource
import signal
import sys
import fresh_compact21_priority_support_scout_20261018 as S
import fresh_sqrt5_forest_sensitive_weights_20261015_review as G

TAG = 'fresh_compact21_h4_transition_budget_probe_20261018'


def run():
    model = dict(G.models())
    families = {
        'H4': ['A1', 'H4', 'E8', 'E8'],
        'A4': ['A1', 'A4', 'E8', 'E8'],
        'A5': ['A5', 'E8', 'E8'],
        'B4': ['A1', 'B4', 'E8', 'E8'],
        'B5': ['B5', 'E8', 'E8'],
    }
    source_roles = {'A4': [(0, 3)], 'A5': [(0, 3), (1, 4)],
                    'B4': [(0, 3)], 'B5': [(1, 4)]}
    encode = lambda obj: json.dumps(obj, separators=(',', ':'))
    connected = {r['key'] for _, r in S.M.objects('connected_orbits')}
    cross = {r['key'] for _, r in S.M.objects('cross_orbits')}
    a4cross = {r['key'] for _, r in S.M.objects('a4_cross_orbits')}
    source_codes = sorted({encode(G.code(model['A4'], (r,))) for r in model['A4']})
    assert len(source_codes) == 2
    coefficients, baselines, context, roles = {}, {}, {}, {}
    for name, forest in families.items():
        a, a4 = forest.count('A1'), forest.count('A4')
        context[name] = (a, 0, 0, a4)
        co, baseline = Counter(), F()
        for family in forest:
            gg = model[family]
            baseline += sum((G.base(gg, u, v) for u, v in combinations(gg, 2)), F())
            for u, v in combinations(gg, 2):
                key = encode(G.code(gg, (u, v)))
                if key in connected:
                    co[key] += 1
            if a and family != 'A1':
                for r in gg:
                    key = 'A1x'+family+':'+encode(G.code(gg, (r,)))
                    if key in cross:
                        co[key] += 1
            if a4 and family not in ('A1', 'I2(5)', 'H3', 'A4'):
                for r in gg:
                    dc = encode(G.code(gg, (r,)))
                    for sc in source_codes:
                        key = 'A4sx'+family+':'+encode([sc, dc])
                        if key in a4cross:
                            co[key] += F(1, 2)
        coefficients[name], baselines[name] = co, baseline
        if name in source_roles:
            roles[name] = [encode(G.code(model[name], pair)) for pair in source_roles[name]]
    primal = json.loads((S.ROOT/'results/fresh_compact21_priority_support_scout_20261018.json').read_text())
    values = {j: F(v) for j, v in primal['positive_bonus_values']}
    upper = F(primal['upper_coordinate'])
    rows = {name: {} for name in families}
    role_values = {name: {} for name in source_roles}
    for j, record in S.M.objects('variable_order'):
        ctx = tuple(record[k] for k in ('A1_count', 'I25_context', 'H3_context', 'A4_context'))
        key = record['orbit']
        for name in families:
            if ctx != context[name]:
                continue
            if key in coefficients[name]:
                rows[name][j] = coefficients[name][key]
            if name in roles and key in roles[name]:
                role_values[name][key] = values.get(j, F())
    checked = []
    for _, term in S.M.objects('terms', 4):
        if term['metadata']['kind'] != 'forest':
            continue
        for name, forest in families.items():
            if sorted(forest) == sorted(term['metadata']['types']):
                saved = {j: F(v) for j, v in term['coefficients'] if j != 93045}
                assert rows[name] == saved and baselines[name] == -F(term['rhs'])
                checked.append(name)
    assert {'H4', 'A4'} <= set(checked)
    result = []
    for name in families:
        total = baselines[name]+sum((a*values.get(j, F()) for j, a in rows[name].items()), F())
        result.append({'forest': families[name], 'old_total': str(total),
                       'old_total_decimal': float(total), 'room_below_saved_upper': str(upper-total),
                       'source_angle_bonus_values': [str(role_values.get(name, {}).get(k, F()))
                                                     for k in roles.get(name, [])]})
    dual_path = 'results/fresh_compact21_priority_support_dual_repaired_20261018.json'
    dual = json.loads((S.ROOT/dual_path).read_text())
    masses = Counter()
    for r in dual['all_nonzero_forest_masses']:
        masses[tuple(sorted(r['types']))] += F(r['positive_mass'])
    source_capacities = {'A4': 1, 'A5': 2, 'B4': 1, 'B5': 1}
    h4_mass = masses[tuple(sorted(families['H4']))]
    source_mass = sum((cap*masses[tuple(sorted(families[name]))]
                       for name, cap in source_capacities.items()), F())
    tload = h4_mass-source_mass
    assert tload > 0
    return {'status': 'PASS_FIVE_FOREST_PRIMAL_BUDGET_PROBE',
            'saved_upper': str(upper), 'support_forest_rows_reconstructed': checked,
            'forest_rows': {name: {'coefficients': [[j, str(a)] for j, a in sorted(rows[name].items())],
                                   'baseline': str(baselines[name]), 'types': families[name]}
                            for name in families},
            'forests': result,
            'source_angle_capacities': source_capacities,
            'old_dual_H4_mass': str(h4_mass),
            'old_dual_weighted_source_mass': str(source_mass),
            'new_transfer_column_load': str(tload),
            'old_dual_violates_new_transfer_column': True,
            'old_dual_repaired_via_transfer_cap_row': str(F(dual['exact_lower_bound'])-tload),
            'new_dimension_bound': None}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, S.G.interrupt)
    signal.alarm(59)
    with (S.ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        S.G.resources()
        data = run()
        data['sources'] = {}
        for path in (f'scripts/{TAG}.py', S.M.MODEL,
                     'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
                     'results/fresh_compact21_priority_support_scout_20261018.json',
                     'results/fresh_compact21_priority_support_dual_repaired_20261018.json'):
            with (S.ROOT/path).open('rb') as stream:
                data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
        target = S.ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == data
        else:
            assert not target.exists()
            temporary = target.with_suffix('.json.tmp')
            try:
                with temporary.open('x') as stream:
                    json.dump(data, stream, indent=2)
                    stream.write('\n')
                    stream.flush()
                    os.fsync(stream.fileno())
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
        S.G.resources(start=False)
    signal.alarm(0)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    assert peak < 1024**3
    print(json.dumps({'status': data['status'], 'forest_totals': [r['old_total'] for r in data['forests']],
                      'peak_rss_bytes': peak, 'new_dimension_bound': None}))
