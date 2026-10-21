#!/usr/bin/env python3
"""One reduced-column IPM min-z solve of the complete forest epigraph."""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
from hashlib import file_digest
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
import fresh_compact21_full_forest_epigraph_20261018 as E

C, D, P, ROOT, Z, T = E.C, E.D, E.P, E.ROOT, E.Z, E.T
TAG = 'fresh_compact21_epigraph_ipm_minimum_20261018'
OUTPUT_LIMIT = 512*1024


def event(kind, **fields):
    print(json.dumps({'event': kind, **fields}, separators=(',', ':')), flush=True)


def worker():
    started, owner = time.monotonic(), int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert os.getppid() == owner
    def owner_watch():
        while True:
            if os.getppid() != owner:
                os.killpg(os.getpgrp(), signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch, daemon=True).start()
    base, native = C.build_rows()
    rule = D.Rule(); models, profiles = D.profiles(rule)
    checked = E.check_existing_rows(rule, profiles)
    grouped, auxiliary, next_column = E.groups(), {}, T+1
    for context, tuples in sorted(grouped.items()):
        for rank in range(2, max(t[-1] for t in tuples)+1):
            auxiliary[context, rank] = next_column; next_column += 1
    live_bonus = set()
    for row, _, kind in base:
        if kind == 'face':
            assert all(j < Z and a <= 0 for j, a in row.items())
            live_bonus.update(j for j, a in row.items() if a < 0)
    columns = sorted(live_bonus | {Z, T} | set(auxiliary.values()))
    index = {j: i for i, j in enumerate(columns)}
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    rr, cc, vv, rhs_values, kinds = [], [], [], [], Counter()
    nonface_bonus_checks = 0
    for i, (row, rhs, kind) in enumerate(E.extended_rows(base, rule, profiles, grouped, auxiliary)):
        kinds[kind] += 1; rhs_values.append(float(rhs))
        for j, a in row.items():
            if kind != 'face' and j < Z:
                assert a >= 0, (i, kind, j, a)
                nonface_bonus_checks += 1
            if not a or (j < Z and j not in live_bonus):
                continue
            assert j in index
            rr.append(i); cc.append(index[j]); vv.append(float(a))
    assert len(rhs_values) == 79274
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(rhs_values), len(columns))).tocsr()
    bounds = [(0., None)]*len(columns)
    objective = np.zeros(len(columns)); objective[index[Z]] = 1.
    build_seconds = time.monotonic()-started
    solve_cap = min(45., max(1., 43.-build_seconds))
    event('matrix_ready', build_seconds=round(build_seconds, 3), rows=len(rhs_values),
          original_global_columns=next_column, kept_columns=len(columns),
          kept_bonus_columns=len(live_bonus), eliminated_bonus_columns=Z-len(live_bonus),
          nonzeros=len(vv), solver='highs-ipm', crossover='off',
          objective='minimize_z', solver_time_limit_seconds=solve_cap,
          nonface_bonus_coefficients_checked=nonface_bonus_checks)
    with warnings.catch_warnings(record=True) as caught:
        result = linprog(objective, A_ub=matrix, b_ub=np.asarray(rhs_values),
                         bounds=bounds, method='highs-ipm',
                         options={'time_limit': solve_cap, 'threads': 1, 'disp': False,
                                  'presolve': True, 'run_crossover': 'off'})
    event('solver_exit', success=bool(result.success), status=int(result.status),
          message=str(result.message)[:300], iterations=getattr(result, 'nit', None),
          crossover_iterations=getattr(result, 'crossover_nit', None),
          primal_available=result.x is not None)
    assert result.x is not None, f'IPM minimum proposal has no primal: {result.status}: {result.message}'
    grid, values = 10**10, {}
    for j in columns:
        if j < Z and result.x[index[j]] > 0:
            raw = F.from_float(float(result.x[index[j]]))*grid
            values[j] = F(-(-raw.numerator//raw.denominator), grid)
    transfer = min(F(1), max(F(), F(float(result.x[index[T]])).limit_denominator(10**8)))
    scale = F(1)
    for row, rhs, kind in base:
        if kind == 'face' and rhs < 0:
            payment = -sum((a*values.get(j, F()) for j, a in row.items()), F())
            assert payment > 0
            scale = max(scale, -rhs/payment)
    values = {j: v*scale for j, v in values.items()}
    rule.values, rule.transfer = values, transfer
    maximum, attainer = E.exact_forest_maximum(rule, models, profiles, grouped)
    for row, rhs, kind in base:
        if kind == 'face':
            assert sum((a*values.get(j, F()) for j, a in row.items()), F()) <= rhs
    data = {
        'status': 'PASS_EXACT_COMPLETE_EPIGRAPH_IPM_MINIMUM_VECTOR', 'native_coverage': native,
        'finite_component_types': len(models), 'exact_hub_contexts': sum(map(len, grouped.values())),
        'capped_context_groups': len(grouped), 'auxiliary_rank_states': len(auxiliary),
        'extended_rows': len(rhs_values), 'row_kinds': dict(kinds), 'nonzeros': len(vv),
        'original_global_columns': next_column, 'compressed_columns_global_indices': columns,
        'kept_bonus_columns': len(live_bonus), 'eliminated_bonus_columns': Z-len(live_bonus),
        'all_nonface_bonus_coefficients_verified_nonnegative': True,
        'existing_forest_rows_reconstructed': checked, 'matrix_build_seconds': round(build_seconds, 3),
        'solver_time_limit_seconds': solve_cap, 'solver': 'highs-ipm', 'crossover': 'off',
        'numerical_upper': float(result.x[index[Z]]), 'numerical_objective': float(result.fun),
        'numerical_solver_success': bool(result.success), 'numerical_solver_status': int(result.status),
        'exact_primal_upper': str(maximum), 'transfer_value': str(transfer),
        'primal_at_most_40': maximum <= 40, 'primal_strictly_below_40': maximum < 40,
        'attaining_forest': attainer, 'primal_bonus_scaling': str(scale),
        'positive_bonus_values': [[j, str(v)] for j, v in sorted(values.items()) if v],
        'retained_numeric_bonus_primal': [[j, float(result.x[index[j]]).hex()]
                                         for j in columns if j < Z and result.x[index[j]] > 0],
        'retained_numeric_transfer': float(result.x[index[T]]).hex(),
        'retained_numeric_dual': [[i, float(y).hex()] for i, y in enumerate(result.ineqlin.marginals) if y] if result.ineqlin.marginals is not None else [],
        'numerical_dual_is_exact_proof': False,
        'solver_warnings': [str(w.message)[:200] for w in caught[:2]],
        'all_old_native_faces_exactly_verified': True, 'all_finite_forests_exactly_verified': True,
        'universal_coverage_composition_reviewed': False, 'new_dimension_bound': None,
        'scope': 'Complete exact weight-model primal after safe zeroing of face-dead bonuses. '
                 'Auxiliary solver values are discarded and every Bellman budget is recomputed '
                 'rationally. Final universal-coverage review remains separate.',
    }
    data['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    assert data['worker_peak_rss_bytes'] < C.H.LIMIT
    path = Path(os.environ['COXETER_FULL_FACE_TEMP'])
    size = C.bounded_write(path, data)
    with path.open('rb') as stream:
        digest = file_digest(stream, 'sha256').hexdigest()
    event('certificate_ready', status=data['status'], certificate_bytes=size,
          certificate_sha256=digest, exact_primal_upper=str(maximum))


def parsed_events(output):
    result = []
    for line in output.decode(errors='replace').splitlines():
        try:
            item = json.loads(line)
            if isinstance(item, dict) and 'event' in item:
                result.append(item)
        except json.JSONDecodeError:
            pass
    return result


def guarded(lock_fd, temporary):
    before = C.H.S.G.resources(); C.H.S.process_tree_rss()
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[key] = '1'
    env['COXETER_SUPPORT_SCOUT_OWNER'] = str(os.getpid())
    env['COXETER_FULL_FACE_TEMP'] = str(temporary)
    child = None; output = bytearray(); started = time.monotonic(); owner = os.getppid(); peak = 0
    relayed = 0
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
                events = parsed_events(output)
                for item in events[relayed:]:
                    if item['event'] != 'certificate_ready':
                        print(json.dumps(item), flush=True)
                relayed = len(events)
                peak = max(peak, C.H.S.process_tree_rss())
                assert peak < C.H.LIMIT and time.monotonic()-started < 59
                assert os.getppid() == owner
                C.H.S.G.resources(start=False)
        assert child.wait() == 0, bytes(output[-4096:]).decode(errors='replace')
        events = parsed_events(output)
        summary = next(item for item in events if item['event'] == 'certificate_ready')
        assert temporary.stat().st_size == summary['certificate_bytes'] < E.ARTIFACT_LIMIT
        with temporary.open('rb') as stream:
            assert file_digest(stream, 'sha256').hexdigest() == summary['certificate_sha256']
        with temporary.open() as stream:
            data = json.load(stream)
        assert data['exact_primal_upper'] == summary['exact_primal_upper']
        data['guard'] = {'elapsed_seconds': round(time.monotonic()-started, 3),
                         'process_tree_peak_rss_bytes': peak, 'rss_limit_bytes': C.H.LIMIT,
                         'wall_limit_seconds': 59, 'output_bytes': len(output),
                         'initial_resources': before, 'thread_limit': 1, 'progress_events': events}
        return data
    except BaseException as exc:
        failure = {'status': 'FAILED_BOUNDED_EPIGRAPH_IPM_MINIMUM_ATTEMPT',
                   'error': str(exc)[-3500:], 'elapsed_seconds': round(time.monotonic()-started, 3),
                   'process_tree_peak_rss_bytes': peak, 'initial_resources': before,
                   'progress_events': parsed_events(output), 'new_dimension_bound': None}
        path = ROOT/f'results/{TAG}_failure.json'
        if not path.exists():
            C.bounded_write(path, failure)
        raise
    finally:
        if child is not None:
            E.SAFE.stop_group(child, grace_seconds=.3)


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
            assert target.stat().st_size < E.ARTIFACT_LIMIT
            with target.open() as stream:
                data = json.load(stream)
            for path, digest in data['sources'].items():
                with (ROOT/path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
            E.verify(data); C.H.S.G.resources(start=False)
        else:
            assert not target.exists() and not temporary.exists() and not output_tmp.exists()
            try:
                data = guarded(lock.fileno(), temporary)
                data['sources'] = {}
                for path in (f'scripts/{TAG}.py', f'scripts/{E.TAG}.py', f'scripts/{D.TAG}.py',
                             f'scripts/{C.TAG}.py', P.BASE, P.QUEUE, C.H.S.M.MODEL, C.CAPACITY,
                             'results/compact21_full_forest_bellman_formulation_20261018.md',
                             'results/review_compact21_full_forest_bellman_formulation_resource_20261018.md'):
                    with (ROOT/path).open('rb') as stream:
                        data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
                C.bounded_write(output_tmp, data); os.replace(output_tmp, target)
            finally:
                temporary.unlink(missing_ok=True); output_tmp.unlink(missing_ok=True)
    signal.alarm(0)
    print(json.dumps({k: data[k] for k in ('status', 'exact_primal_upper', 'transfer_value',
                                          'primal_at_most_40', 'new_dimension_bound')}))


if __name__ == '__main__':
    main()
