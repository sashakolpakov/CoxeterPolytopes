#!/usr/bin/env python3
"""One capped full-face/full-forest Bellman-epigraph research solve.

All finite forests are covered by an additive rank DP, with exact exclusions
for the two decreased physical forest types. Rational verification recomputes
minimal DP states from the bonus vector, never trusts numerical auxiliaries.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
import fcntl
from hashlib import file_digest
from itertools import product
import json
import os
from pathlib import Path
import resource
import selectors
import signal
import subprocess
import sys
import threading
import time
import warnings
import fresh_compact21_full_candidate_forest_separation_20261018 as D
import run_lean_batches_signal_safe_20261018 as SAFE

C, P, ROOT, Z, T = D.C, D.P, D.ROOT, D.Z, D.T
TAG = 'fresh_compact21_full_forest_epigraph_20261018'
ARTIFACT_LIMIT, OUTPUT_LIMIT = 8*1024**2, 512*1024
HUBS = D.HUBS
SPECIAL = {
    'H4': tuple(sorted(('A1', 'H4', 'E8', 'E8'))),
    'D5': tuple(sorted(('D5', 'E8', 'E8'))),
    'A4': tuple(sorted(('A1', 'A4', 'E8', 'E8'))),
    'A5': tuple(sorted(('A5', 'E8', 'E8'))),
    'B4': tuple(sorted(('A1', 'B4', 'E8', 'E8'))),
    'B5': tuple(sorted(('B5', 'E8', 'E8'))),
}


def groups():
    result = defaultdict(list)
    for a in range(22):
        for b in range((21-a)//2+1):
            for h in range((21-a-2*b)//3+1):
                for d in range((21-a-2*b-3*h)//4+1):
                    rank = 21-a-2*b-3*h-4*d
                    if rank != 1:
                        result[a, min(b, 2), min(h, 1), min(d, 2)].append((a, b, h, d, rank))
    assert sum(map(len, result.values())) == 729
    return dict(result)


def component_row(rule, profile, context):
    co = Counter()
    for table in D.enabled_tables(profile, context):
        for key, n in table.items():
            j = rule.index_or_none(key, context)
            if j is not None:
                co[j] += n
    return dict(co), profile[1]


def beta_key(models):
    e8 = models['E8']
    fork = next(i for i in e8 if len(e8[i]) == 3)
    ends = {}
    for neighbor in e8[fork]:
        previous, node, length = fork, neighbor, 1
        while len(e8[node]) == 2:
            nxt = next(i for i in e8[node] if i != previous)
            previous, node, length = node, nxt, length+1
        ends[length] = node
    assert set(ends) == {1, 2, 4}
    return P.encoded(P.G.code(e8, (ends[2], ends[4])))


def forest_row(rule, profiles, forest, state=None):
    context = D.exact_context(forest)
    co, baseline = Counter({Z: F(-1)}), F()
    for name in forest:
        part, base = component_row(rule, profiles[name], context)
        co.update(part); baseline += base
    rhs = -baseline
    if forest == SPECIAL['H4']:
        co[T] -= 1
    elif forest == SPECIAL['D5']:
        assert state in ('active', 'inactive')
        if state == 'active':
            rhs += 1
        else:
            graph = dict(P.G.models())['E8']
            key = beta_key({'E8': graph})
            column = rule.index_or_none(key, context)
            assert column is not None and co[column] == 2
            co[column] -= 2
    elif forest in (SPECIAL['A4'], SPECIAL['A5'], SPECIAL['B4'], SPECIAL['B5']):
        co[T] += 1
    return {j: v for j, v in co.items() if v}, rhs


def check_existing_rows(rule, profiles):
    metadata = C.H.D.forest_metadata()
    checked = 0
    for i, (literal, rhs, kind) in enumerate(C.H.rows()):
        if kind != 'forest':
            continue
        if i in metadata:
            forest = tuple(sorted(metadata[i]['types']))
            state = metadata[i].get('state')
        else:
            name = {2805: 'A5', 2806: 'B4', 2807: 'B5'}[i]
            forest, state = SPECIAL[name], None
        if i == 2805:
            literal = literal | {T: F(1)}
        row, expected_rhs = forest_row(rule, profiles, forest, state)
        assert row == literal and expected_rhs == rhs, (i, forest)
        checked += 1
    assert checked == 510
    return checked


def context_rows(rule, profiles, context):
    return {name: component_row(rule, profile, context) for name, profile in profiles.items()}


def combine(rows, counts):
    co, base = Counter(), F()
    for name, n in counts:
        row, b = rows[name]
        co.update({j: n*v for j, v in row.items()}); base += n*b
    return co, base


def extended_rows(base_rows, rule, profiles, grouped, auxiliary):
    yield from base_rows
    for context, tuples in sorted(grouped.items()):
        rows = context_rows(rule, profiles, context)
        maximum_rank = max(r[-1] for r in tuples)
        items = [name for name in profiles if name not in HUBS]
        for rank in range(2, maximum_rank+1):
            for name in items:
                size = profiles[name][0]
                remaining = rank-size
                if remaining < 0 or remaining == 1:
                    continue
                row, baseline = rows[name]
                co = Counter(row)
                co[auxiliary[context, rank]] -= 1
                if remaining:
                    co[auxiliary[context, remaining]] += 1
                yield dict(co), -baseline, 'bellman'
        for a, b, h, d, rank in tuples:
            if (a, b, h, d) in ((1, 0, 0, 0), (0, 0, 0, 0)):
                target_types = {'H4', 'E8'} if a == 1 else {'D5', 'E8'}
                fixed, baseline = combine(rows, zip(HUBS, (a, b, h, d)))
                for name in items:
                    remaining = rank-profiles[name][0]
                    if name in target_types or remaining < 0 or remaining == 1:
                        continue
                    co = Counter(fixed); co.update(rows[name][0]); co[Z] -= 1
                    if remaining:
                        co[auxiliary[context, remaining]] += 1
                    yield dict(co), -baseline-rows[name][1], 'off_target_budget'
            else:
                co, baseline = combine(rows, zip(HUBS, (a, b, h, d)))
                co[Z] -= 1
                if rank:
                    co[auxiliary[context, rank]] += 1
                yield dict(co), -baseline, 'ordinary_budget'
    # Pure target-type forests left after off-target exclusion.
    for forest in (tuple(sorted(('A1',)+('H4',)*5)),
                   tuple(sorted(('A1',)+('H4',)*3+('E8',)))):
        row, rhs = forest_row(rule, profiles, forest)
        yield row, rhs, 'pure_unchanged_budget'
    # Explicit exceptional constraints, also present in the checked old support.
    for name, forest in SPECIAL.items():
        for state in (('active', 'inactive') if name == 'D5' else (None,)):
            row, rhs = forest_row(rule, profiles, forest, state)
            yield row, rhs, 'special_budget'


def exact_forest_maximum(rule, models, profiles, grouped):
    best, attainer = F(-1), None
    def offer(value, forest, kind):
        nonlocal best, attainer
        if value > best:
            best, attainer = value, {'types': list(sorted(forest)), 'budget_kind': kind}
    for context, tuples in sorted(grouped.items()):
        values = D.component_values(rule, profiles, context)
        maximum_rank = max(t[-1] for t in tuples)
        states = {0: (F(), ())}
        for rank in range(2, maximum_rank+1):
            choices = []
            for name, profile in profiles.items():
                if name in HUBS or rank-profile[0] not in states:
                    continue
                value, forest = states[rank-profile[0]]
                choices.append((value+values[name], tuple(sorted(forest+(name,)))))
            assert choices
            states[rank] = max(choices)
        for a, b, h, d, rank in tuples:
            fixed = sum((n*values[name] for name, n in zip(HUBS, (a, b, h, d))), F())
            hubs = ('A1',)*a+('I2(5)',)*b+('H3',)*h+('A4',)*d
            if (a, b, h, d) in ((1, 0, 0, 0), (0, 0, 0, 0)):
                target_types = {'H4', 'E8'} if a else {'D5', 'E8'}
                for name, profile in profiles.items():
                    remaining = rank-profile[0]
                    if name in HUBS or name in target_types or remaining not in states:
                        continue
                    value, forest = states[remaining]
                    offer(fixed+values[name]+value, hubs+forest+(name,), 'off_target')
            else:
                value, forest = states[rank]
                offer(fixed+value, hubs+forest, 'ordinary')
    extra = [tuple(sorted(('A1',)+('H4',)*5)), tuple(sorted(('A1',)+('H4',)*3+('E8',)))]
    extra += list(SPECIAL.values())
    for forest in extra:
        context = D.exact_context(forest)
        values = D.component_values(rule, profiles, context)
        value = sum((values[name] for name in forest), F())
        if forest == SPECIAL['H4']:
            value -= rule.transfer
        elif forest == SPECIAL['D5']:
            value -= min(F(1), 2*rule.bonus(beta_key(models), context))
        elif forest in (SPECIAL['A4'], SPECIAL['A5'], SPECIAL['B4'], SPECIAL['B5']):
            value += rule.transfer
        offer(value, forest, 'explicit')
    return best, attainer


def rule_from_data(data):
    rule = D.Rule()
    rule.values = {j: F(v) for j, v in data['positive_bonus_values']}
    rule.transfer = F(data['transfer_value'])
    rule.saved_upper = F(data['exact_primal_upper'])
    return rule


def verify(data):
    rule = rule_from_data(data)
    assert all(v > 0 for v in rule.values.values()) and 0 <= rule.transfer <= 1
    base, coverage = C.build_rows()
    assert coverage == data['native_coverage']
    for row, rhs, kind in base:
        if kind != 'face':
            continue
        assert sum((a*rule.values.get(j, F()) for j, a in row.items()), F()) <= rhs
    models, profiles = D.profiles(rule)
    assert check_existing_rows(rule, profiles) == data['existing_forest_rows_reconstructed']
    maximum, attainer = exact_forest_maximum(rule, models, profiles, groups())
    assert maximum == F(data['exact_primal_upper']) and attainer == data['attaining_forest']
    return maximum


def worker():
    started, owner = time.monotonic(), int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert os.getppid() == owner
    def owner_watch():
        while True:
            if os.getppid() != owner:
                os.killpg(os.getpgrp(), signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch, daemon=True).start()
    base_rows, native = C.build_rows()
    rule = D.Rule(); models, profiles = D.profiles(rule)
    checked = check_existing_rows(rule, profiles)
    grouped = groups(); auxiliary = {}
    next_column = T+1
    for context, tuples in sorted(grouped.items()):
        for rank in range(2, max(t[-1] for t in tuples)+1):
            auxiliary[context, rank] = next_column; next_column += 1
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    rr, cc, vv, rhs_values, kinds = [], [], [], [], Counter()
    for i, (row, rhs, kind) in enumerate(extended_rows(base_rows, rule, profiles, grouped, auxiliary)):
        kinds[kind] += 1; rhs_values.append(float(rhs))
        for j, a in row.items():
            if a:
                rr.append(i); cc.append(j); vv.append(float(a))
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(rhs_values), next_column)).tocsr()
    objective = np.zeros(next_column); objective[Z] = 1
    build_seconds = time.monotonic()-started
    solve_cap = min(45., max(1., 43.-build_seconds))
    with warnings.catch_warnings(record=True) as caught:
        result = linprog(objective, A_ub=matrix, b_ub=np.asarray(rhs_values), bounds=(0, None),
                         method='highs', options={'time_limit': solve_cap, 'threads': 1, 'disp': False})
    assert result.success, f'Extended numerical proposal failed: {result.status}: {result.message}'
    grid, values = 10**10, {}
    for j in range(Z):
        if result.x[j] > 0:
            raw = F.from_float(float(result.x[j]))*grid
            values[j] = F(-(-raw.numerator//raw.denominator), grid)
    transfer = min(F(1), max(F(), F(float(result.x[T])).limit_denominator(10**8)))
    scale = F(1)
    for row, rhs, kind in base_rows:
        if kind == 'face' and rhs < 0:
            payment = -sum((a*values.get(j, F()) for j, a in row.items()), F())
            assert payment > 0
            scale = max(scale, -rhs/payment)
    values = {j: v*scale for j, v in values.items()}
    rule.values, rule.transfer = values, transfer
    maximum, attainer = exact_forest_maximum(rule, models, profiles, grouped)
    for row, rhs, kind in base_rows:
        if kind == 'face':
            assert sum((a*values.get(j, F()) for j, a in row.items()), F()) <= rhs
    data = {
        'status': 'PASS_EXACT_FULL_FACE_AND_FOREST_EPIGRAPH_VECTOR', 'native_coverage': native,
        'finite_component_types': len(models), 'exact_hub_contexts': sum(map(len, grouped.values())),
        'capped_context_groups': len(grouped), 'auxiliary_rank_states': len(auxiliary),
        'extended_rows': len(rhs_values), 'row_kinds': dict(kinds), 'nonzeros': len(vv),
        'extended_variable_columns': next_column, 'existing_forest_rows_reconstructed': checked,
        'matrix_build_seconds': round(build_seconds, 3), 'solver_time_limit_seconds': solve_cap,
        'numerical_objective': float(result.fun), 'exact_primal_upper': str(maximum),
        'transfer_value': str(transfer), 'primal_at_most_40': maximum <= 40,
        'primal_strictly_below_40': maximum < 40, 'attaining_forest': attainer,
        'positive_bonus_values': [[j, str(v)] for j, v in sorted(values.items()) if v],
        'retained_numeric_dual': [[i, float(y).hex()] for i, y in enumerate(result.ineqlin.marginals) if y],
        'numerical_dual_is_exact_proof': False, 'primal_bonus_scaling': str(scale),
        'solver_warnings': [str(w.message)[:180] for w in caught[:2]],
        'all_old_native_faces_exactly_verified': True, 'all_finite_forests_exactly_verified': True,
        'universal_coverage_composition_reviewed': False, 'new_dimension_bound': None,
        'scope': 'Exact full weight-model vector. Numerical auxiliary values are discarded; '
                 'minimal Bellman values and every forest-budget envelope are recomputed rationally. '
                 'Final dimension composition and independent verification remain separate.',
    }
    data['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    assert data['worker_peak_rss_bytes'] < C.H.LIMIT
    path = Path(os.environ['COXETER_FULL_FACE_TEMP'])
    size = C.bounded_write(path, data)
    with path.open('rb') as stream:
        digest = file_digest(stream, 'sha256').hexdigest()
    print(json.dumps({'status': data['status'], 'certificate_bytes': size,
                      'certificate_sha256': digest, 'exact_primal_upper': str(maximum)}), flush=True)


def guarded(lock_fd, temporary):
    before = C.H.S.G.resources(); C.H.S.process_tree_rss()
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[key] = '1'
    env['COXETER_SUPPORT_SCOUT_OWNER'] = str(os.getpid())
    env['COXETER_FULL_FACE_TEMP'] = str(temporary)
    child = None; output = bytearray(); started = time.monotonic(); owner = os.getppid(); peak = 0
    try:
        child = subprocess.Popen([sys.executable, '-B', str(Path(__file__).resolve()), '--worker'],
                                 cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 start_new_session=True, pass_fds=(lock_fd,))
        with selectors.DefaultSelector() as poll:
            poll.register(child.stdout, selectors.EVENT_READ)
            while poll.get_map() or child.poll() is None:
                for key, _ in poll.select(timeout=.1):
                    block = os.read(key.fd, 65536)
                    if block:
                        output.extend(block); assert len(output) < OUTPUT_LIMIT
                    else:
                        poll.unregister(key.fileobj)
                peak = max(peak, C.H.S.process_tree_rss())
                assert peak < C.H.LIMIT and time.monotonic()-started < 59
                assert os.getppid() == owner
                C.H.S.G.resources(start=False)
        assert child.wait() == 0, bytes(output[-4096:]).decode(errors='replace')
        summary = json.loads(output)
        assert temporary.stat().st_size == summary['certificate_bytes'] < ARTIFACT_LIMIT
        with temporary.open('rb') as stream:
            assert file_digest(stream, 'sha256').hexdigest() == summary['certificate_sha256']
        with temporary.open() as stream:
            data = json.load(stream)
        assert data['exact_primal_upper'] == summary['exact_primal_upper']
        data['guard'] = {'elapsed_seconds': round(time.monotonic()-started, 3),
                         'process_tree_peak_rss_bytes': peak, 'rss_limit_bytes': C.H.LIMIT,
                         'wall_limit_seconds': 59, 'output_bytes': len(output),
                         'initial_resources': before, 'thread_limit': 1,
                         'cleanup': 'signal-safe final group kill after reaping leader'}
        return data
    finally:
        if child is not None:
            SAFE.stop_group(child, grace_seconds=.3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.worker:
        worker(); return
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
        signal.signal(sig, C.H.S.G.interrupt)
    signal.alarm(59)
    target = ROOT/f'results/{TAG}.json'
    temporary, output_tmp = target.with_suffix('.candidate.tmp'), target.with_suffix('.json.tmp')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.verify:
            C.H.S.G.resources()
            assert target.stat().st_size < ARTIFACT_LIMIT
            with target.open() as stream:
                data = json.load(stream)
            for path, digest in data['sources'].items():
                with (ROOT/path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
            verify(data); C.H.S.G.resources(start=False)
        else:
            assert not target.exists() and not temporary.exists() and not output_tmp.exists()
            try:
                data = guarded(lock.fileno(), temporary)
                data['sources'] = {}
                for path in (f'scripts/{TAG}.py', f'scripts/{D.TAG}.py', f'scripts/{C.TAG}.py',
                             P.BASE, P.QUEUE, C.H.S.M.MODEL, C.CAPACITY,
                             'results/compact21_full_forest_bellman_formulation_20261018.md',
                             'results/review_compact21_full_forest_bellman_formulation_resource_20261018.md',
                             'scripts/fresh_sqrt5_dimension21_a4_source_orbits_20261016.py'):
                    with (ROOT/path).open('rb') as stream:
                        data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
                C.bounded_write(output_tmp, data); os.replace(output_tmp, target)
            finally:
                temporary.unlink(missing_ok=True); output_tmp.unlink(missing_ok=True)
    signal.alarm(0)
    print(json.dumps({k: data[k] for k in ('status', 'numerical_objective', 'exact_primal_upper',
                                          'transfer_value', 'primal_at_most_40', 'new_dimension_bound')}))


if __name__ == '__main__':
    main()
