#!/usr/bin/env python3
"""One capped full native-face solve, streamed without the old producer caches.

All 102,917 native spectator cases are retained before exact row deduplication.
Only the existing forest support is imposed; complete forest separation follows.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import fcntl
from hashlib import file_digest
from itertools import combinations, product
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
import fresh_compact21_h4_successor_face_separation_20261018 as P

H, ROOT = P.H, P.ROOT
TAG = 'fresh_compact21_h4_full_face_successor_20261018'
CAPACITY = 'results/fresh_compact21_a5_h4_source_capacity_20261018_freeze.json'
OUTPUT_LIMIT, ARTIFACT_LIMIT = 512*1024, 8*1024**2
Z, T = H.Z, H.T


def build_rows():
    rule = P.Rule()
    unique, cases, face_count = set(), 0, 0
    for _, _, row in P.native_faces():
        face_count += 1
        edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
        corners = [P.corner(row, edges, a, b) for a, b in product((0, 2), (1, 3))]
        baseline = sum((c[0] for c in corners), F())
        for outside in P.spectators(23-row['nodes']):
            cases += 1
            co = Counter()
            for _, local, key, denominator, _ in corners:
                if key is None:
                    continue
                counts = tuple(a+b for a, b in zip(local, outside))
                factor = F(1)
                if denominator is not None:
                    k, multiplier, offset = denominator
                    divisor = multiplier*(counts[k]+offset)
                    assert divisor > 0
                    factor /= divisor
                co[rule.column(key, counts)] -= factor
            unique.add((tuple(sorted(co.items())), baseline-1))
    assert face_count == 1979 and cases == 102917
    assert len(unique) == 51443
    saved = [(dict(co), rhs, 'face') for co, rhs in sorted(unique)]
    del unique, rule
    tails = []
    for i, (row, rhs, kind) in enumerate(H.rows()):
        if kind == 'face':
            continue
        if i == 2805:
            assert kind == 'forest' and row[T] == 2
            row = row | {T: F(1)}
        tails.append((row, rhs, kind))
    assert len(tails) == 511 and tails[-1] == ({T: F(1)}, F(1), 'cap')
    saved.extend(tails)
    assert len(saved) == 51954
    ledger = {'native_faces': face_count, 'native_spectator_cases': cases,
              'exact_distinct_face_rows': 51443, 'existing_forest_rows': 510,
              'transfer_cap_rows': 1, 'A5_source_capacity': 1,
              'deduplication': 'Only exactly equal rational coefficient rows and right-hand sides.'}
    return saved, ledger


def verify(saved, data):
    values = {j: F(v) for j, v in data['positive_bonus_values']}
    assert len(values) == len(data['positive_bonus_values'])
    assert all(0 <= j < Z and v > 0 for j, v in values.items())
    values[Z], values[T] = F(data['exact_primal_upper']), F(data['transfer_value'])
    assert values[Z] >= 0 and 0 <= values[T] <= 1
    weights = {i: F(y) for i, y in data['exact_dual_rows']}
    assert len(weights) == len(data['exact_dual_rows']) and all(y < 0 for y in weights.values())
    load, bound = Counter(), F()
    minimum = None
    for i, (row, rhs, _) in enumerate(saved):
        slack = rhs-sum((a*values.get(j, F()) for j, a in row.items()), F())
        assert slack >= 0, (i, slack)
        minimum = slack if minimum is None else min(minimum, slack)
        y = weights.get(i, F())
        if y:
            bound += y*rhs
            for j, a in row.items():
                load[j] += y*a
    assert all(0 <= i < len(saved) for i in weights)
    assert all(v <= F(j == Z) for j, v in load.items())
    assert bound == F(data['exact_dual_lower']) <= values[Z]
    return minimum


def bounded_write(path, data):
    encoder, size = json.JSONEncoder(separators=(',', ':')), 0
    with path.open('x') as stream:
        for chunk in encoder.iterencode(data):
            size += len(chunk.encode())
            assert size < ARTIFACT_LIMIT
            stream.write(chunk)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
    return size+1


def worker():
    started = time.monotonic()
    owner = int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert os.getppid() == owner
    def owner_watch():
        while True:
            if os.getppid() != owner:
                os.killpg(os.getpgrp(), signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch, daemon=True).start()
    saved, ledger = build_rows()
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    columns = sorted({j for row, _, _ in saved for j in row})
    index = {j: i for i, j in enumerate(columns)}
    rr, cc, vv = [], [], []
    for i, (row, _, _) in enumerate(saved):
        for j, a in row.items():
            rr.append(i); cc.append(index[j]); vv.append(float(a))
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(saved), len(columns))).tocsr()
    objective = np.zeros(len(columns)); objective[index[Z]] = 1
    build_seconds = time.monotonic()-started
    solve_cap = min(45., max(1., 46.-build_seconds))
    with warnings.catch_warnings(record=True) as caught:
        result = linprog(objective, A_ub=matrix, b_ub=np.asarray([float(b) for _, b, _ in saved]),
                         bounds=(0, None), method='highs',
                         options={'time_limit': solve_cap, 'threads': 1, 'disp': False})
    assert result.success, f'Full-face numerical proposal failed: {result.status}: {result.message}'
    grid, values = 10**10, {}
    for j in columns:
        if j < Z and result.x[index[j]] > 0:
            raw = F.from_float(float(result.x[index[j]]))*grid
            values[j] = F(-(-raw.numerator//raw.denominator), grid)
    transfer = min(F(1), max(F(), F(float(result.x[index[T]])).limit_denominator(10**8)))
    scale = F(1)
    for row, rhs, kind in saved:
        if kind == 'face' and rhs < 0:
            payment = -sum((a*values.get(j, F()) for j, a in row.items()), F())
            assert payment > 0
            scale = max(scale, -rhs/payment)
    values = {j: v*scale for j, v in values.items()}
    values[T] = transfer
    upper = max([F()]+[sum((a*values.get(j, F()) for j, a in row.items() if j != Z), F())-rhs
                         for row, rhs, kind in saved if kind == 'forest'])
    dual, lower, clipped = H.exact_dual(saved, result.ineqlin.marginals)
    data = {
        'status': 'PASS_EXACT_FULL_NATIVE_FACES_WITH_FOREST_SUPPORT', 'coverage': ledger,
        'rows': len(saved), 'occurring_columns': len(columns), 'nonzeros': len(vv),
        'numerical_objective': float(result.fun), 'solver_time_limit_seconds': solve_cap,
        'matrix_build_seconds': round(build_seconds, 3),
        'exact_primal_upper': str(upper), 'exact_dual_lower': str(lower),
        'transfer_value': str(transfer), 'primal_at_most_40': upper <= 40,
        'dual_above_40': lower > 40, 'primal_bonus_scaling': str(scale),
        'positive_bonus_values': [[j, str(v)] for j, v in sorted(values.items()) if j < Z and v],
        'exact_dual_rows': [[i, str(y)] for i, y in sorted(dual.items())],
        'retained_numeric_dual': [[i, float(y).hex()] for i, y in enumerate(result.ineqlin.marginals) if y],
        'dual_face_rows_clipped': clipped, 'solver_warnings': [str(w.message)[:180] for w in caught[:2]],
        'all_old_native_faces_exactly_verified': True, 'all_finite_forests_verified': False,
        'new_dimension_bound': None,
        'scope': 'Complete old native-face model plus existing forest support and strengthened H4 transfer. '
                 'Full finite-forest separation is still required before any dimension conclusion.',
    }
    data['minimum_exact_row_slack'] = str(verify(saved, data))
    data['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    assert data['worker_peak_rss_bytes'] < H.LIMIT
    target = Path(os.environ['COXETER_FULL_FACE_TEMP'])
    size = bounded_write(target, data)
    with target.open('rb') as stream:
        digest = file_digest(stream, 'sha256').hexdigest()
    print(json.dumps({'status': data['status'], 'certificate_bytes': size,
                      'certificate_sha256': digest, 'exact_primal_upper': str(upper),
                      'exact_dual_lower': str(lower)}), flush=True)


def guarded(lock_fd, temporary):
    before = H.S.G.resources()
    H.S.process_tree_rss()
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
                peak = max(peak, H.S.process_tree_rss())
                assert peak < H.LIMIT and time.monotonic()-started < 59
                assert os.getppid() == owner
                H.S.G.resources(start=False)
        assert child.wait() == 0, bytes(output[-4096:]).decode(errors='replace')
        summary = json.loads(output)
        assert temporary.stat().st_size == summary['certificate_bytes'] < ARTIFACT_LIMIT
        with temporary.open('rb') as stream:
            assert file_digest(stream, 'sha256').hexdigest() == summary['certificate_sha256']
        with temporary.open() as stream:
            data = json.load(stream)
        assert data['exact_primal_upper'] == summary['exact_primal_upper']
        data['guard'] = {'elapsed_seconds': round(time.monotonic()-started, 3),
                         'process_tree_peak_rss_bytes': peak, 'rss_limit_bytes': H.LIMIT,
                         'wall_limit_seconds': 59, 'output_bytes': len(output),
                         'artifact_limit_bytes': ARTIFACT_LIMIT, 'initial_resources': before,
                         'thread_limit': 1}
        return data
    finally:
        if child is not None:
            try: os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError: pass
            try: child.wait(timeout=.3)
            except subprocess.TimeoutExpired:
                try: os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError: pass
                child.wait()
            if child.stdout is not None:
                child.stdout.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.worker:
        worker(); return
    signal.signal(signal.SIGTERM, H.S.G.interrupt)
    signal.signal(signal.SIGINT, H.S.G.interrupt)
    signal.signal(signal.SIGALRM, H.S.G.interrupt)
    signal.alarm(59)
    target = ROOT/f'results/{TAG}.json'
    temporary = target.with_suffix('.candidate.tmp')
    output_tmp = target.with_suffix('.json.tmp')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.verify:
            H.S.G.resources()
            assert target.stat().st_size < ARTIFACT_LIMIT
            with target.open() as stream:
                data = json.load(stream)
            for path, digest in data['sources'].items():
                with (ROOT/path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
            saved, ledger = build_rows()
            assert ledger == data['coverage']
            assert str(verify(saved, data)) == data['minimum_exact_row_slack']
            H.S.G.resources(start=False)
        else:
            assert not target.exists() and not temporary.exists() and not output_tmp.exists()
            frozen = json.loads((ROOT/CAPACITY).read_text())
            for path, digest in frozen['files'].items():
                with (ROOT/path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
            try:
                data = guarded(lock.fileno(), temporary)
                data['sources'] = {}
                for path in (f'scripts/{TAG}.py', f'scripts/{P.TAG}.py', f'scripts/{H.TAG}.py',
                             P.BASE, P.QUEUE, H.S.M.MODEL, CAPACITY,
                             'scripts/fresh_sqrt5_dimension21_a4_source_orbits_20261016.py'):
                    with (ROOT/path).open('rb') as stream:
                        data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
                bounded_write(output_tmp, data)
                os.replace(output_tmp, target)
            finally:
                temporary.unlink(missing_ok=True); output_tmp.unlink(missing_ok=True)
    signal.alarm(0)
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024) < H.LIMIT
    print(json.dumps({k: data[k] for k in ('status', 'exact_primal_upper', 'exact_dual_lower',
                                          'transfer_value', 'primal_at_most_40', 'dual_above_40',
                                          'new_dimension_bound')}))


if __name__ == '__main__':
    main()
