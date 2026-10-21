#!/usr/bin/env python3
"""One guarded small support scout for the proved H4 angle transfer.

Both a rational primal and a rational weak-duality certificate are checked.
Numerical dual multipliers are retained for recovery without another solve.
This does not perform full-model primal separation.
"""
import argparse
from collections import defaultdict
import fcntl
from fractions import Fraction as F
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
import fresh_compact21_priority_support_exact_dual_20261018 as D

S, ROOT = D.S, D.ROOT
TAG = 'fresh_compact21_h4_funded_successor_scout_20261018'
PROBE = 'results/fresh_compact21_h4_transition_budget_probe_20261018.json'
GEOMETRY = 'results/fresh_compact21_h4_triangle_funded_transfer_20261018_freeze.json'
Z, T = 93045, 93046
LIMIT, OUTPUT_LIMIT = 1024**3, 512*1024
SOURCE_CAPACITIES = {'A4': 1, 'A5': 2, 'B4': 1, 'B5': 1}


def rows():
    probe = json.loads((ROOT/PROBE).read_text())
    metadata = D.forest_metadata()
    matched = set()
    for i, (row, rhs, kind) in enumerate(S.rows()):
        if i in metadata:
            types = sorted(metadata[i]['types'])
            if types == ['A1', 'E8', 'E8', 'H4']:
                row = row | {T: F(-1)}
                matched.add('H4')
            elif types == ['A1', 'A4', 'E8', 'E8']:
                row = row | {T: F(1)}
                matched.add('A4')
            assert types not in (['A5', 'E8', 'E8'], ['A1', 'B4', 'E8', 'E8'],
                                 ['B5', 'E8', 'E8'])
        yield row, rhs, kind
    assert matched == {'H4', 'A4'}
    for name in ('A5', 'B4', 'B5'):
        record = probe['forest_rows'][name]
        row = {j: F(a) for j, a in record['coefficients']}
        row[Z], row[T] = F(-1), F(SOURCE_CAPACITIES[name])
        yield row, -F(record['baseline']), 'forest'
    yield {T: F(1)}, F(1), 'cap'


def verify_primal(data):
    values = {j: F(v) for j, v in data['positive_bonus_values']}
    assert len(values) == len(data['positive_bonus_values'])
    assert all(0 <= j < Z and value > 0 for j, value in values.items())
    values[Z], values[T] = F(data['exact_primal_upper']), F(data['transfer_value'])
    assert values[Z] >= 0 and 0 <= values[T] <= 1
    count = 0
    for row, rhs, _ in rows():
        assert sum((a*values.get(j, F()) for j, a in row.items()), F()) <= rhs
        count += 1
    assert count == 2809


def verify_dual(data):
    weights = {i: F(y) for i, y in data['exact_dual_rows']}
    assert len(weights) == len(data['exact_dual_rows']) and all(y < 0 for y in weights.values())
    load, value = defaultdict(F), F()
    count = 0
    for i, (row, rhs, _) in enumerate(rows()):
        y = weights.get(i, F())
        for j, a in row.items():
            load[j] += y*a
        value += y*rhs
        count += 1
    assert count == 2809 and all(0 <= i < count for i in weights)
    assert all(a <= F(j == Z) for j, a in load.items())
    assert value == F(data['exact_dual_lower'])
    assert value <= F(data['exact_primal_upper'])


def exact_dual(saved, marginals):
    grid = 10**12
    face, forest = {}, {}
    for i, raw in enumerate(marginals):
        if saved[i][2] == 'cap':
            continue
        scaled = max(F(), -F.from_float(float(raw))) * grid
        numerator = scaled.numerator // scaled.denominator
        if saved[i][2] == 'forest' and numerator*scaled.denominator < scaled.numerator:
            numerator += 1
        if numerator:
            (face if saved[i][2] == 'face' else forest)[i] = F(numerator, grid)
    capacity, demand, incident = defaultdict(F), defaultdict(F), defaultdict(list)
    for i, mass in forest.items():
        assert saved[i][0][Z] == -1
        for j, a in saved[i][0].items():
            if j not in (Z, T):
                assert a >= 0
                capacity[j] += mass*a
    for i, mass in face.items():
        for j, a in saved[i][0].items():
            assert j not in (Z, T) and a <= 0
            if a:
                demand[j] -= mass*a
                incident[j].append(i)
    clipped = set()
    for j in sorted(demand):
        if demand[j] <= capacity[j]:
            continue
        for i in sorted(incident[j], key=lambda i: (face[i], i)):
            excess = demand[j]-capacity[j]
            if excess <= 0:
                break
            take = min(face[i], excess/-saved[i][0][j])
            face[i] -= take
            if take:
                clipped.add(i)
            for k, a in saved[i][0].items():
                demand[k] += take*a
        assert demand[j] <= capacity[j]
    normalizer = max(F(1), sum(forest.values(), F()))
    weights = {i: -mass/normalizer for i, mass in face.items() if mass}
    weights.update({i: -mass/normalizer for i, mass in forest.items()})
    transfer_load = sum((y*saved[i][0].get(T, F()) for i, y in weights.items()), F())
    if transfer_load > 0:
        assert saved[-1] == ({T: F(1)}, F(1), 'cap')
        weights[len(saved)-1] = -transfer_load
    bound = sum((y*saved[i][1] for i, y in weights.items()), F())
    return weights, bound, len(clipped)


def worker():
    owner = int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert os.getppid() == owner
    def owner_watch():
        while True:
            if os.getppid() != owner:
                os.killpg(os.getpgrp(), signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch, daemon=True).start()
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix
    saved = list(rows())
    columns = sorted({j for row, _, _ in saved for j in row})
    index = {j: i for i, j in enumerate(columns)}
    assert len(saved) == 2809 and Z in index and T in index
    rr, cc, vv = [], [], []
    for i, (row, _, _) in enumerate(saved):
        for j, a in row.items():
            rr.append(i); cc.append(index[j]); vv.append(float(a))
    matrix = coo_matrix((vv, (rr, cc)), shape=(len(saved), len(columns))).tocsr()
    objective = np.zeros(len(columns)); objective[index[Z]] = 1
    with warnings.catch_warnings(record=True) as caught:
        result = linprog(objective, A_ub=matrix, b_ub=np.asarray([float(b) for _, b, _ in saved]),
                         bounds=(0, None), method='highs',
                         options={'time_limit': 45., 'threads': 1, 'disp': False})
    assert result.success, result.message
    grid = 10**10
    values = {}
    for j in columns:
        if j >= Z or result.x[index[j]] <= 0:
            continue
        raw = F.from_float(float(result.x[index[j]]))*grid
        numerator = -(-raw.numerator//raw.denominator)
        values[j] = F(numerator, grid)
    transfer = min(F(1), max(F(), F(float(result.x[index[T]])).limit_denominator(10**8)))
    values[T] = transfer
    scale = F(1)
    for row, rhs, kind in saved:
        if kind == 'face' and rhs < 0:
            payment = -sum((a*values.get(j, F()) for j, a in row.items()), F())
            assert payment > 0
            scale = max(scale, -rhs/payment)
    values = {j: a*scale for j, a in values.items() if j != T} | {T: transfer}
    upper = max([F()]+[sum((a*values.get(j, F()) for j, a in row.items() if j != Z), F())-rhs
                         for row, rhs, kind in saved if kind == 'forest'])
    weights, lower, clipped = exact_dual(saved, result.ineqlin.marginals)
    data = {
        'status': 'PASS_EXACT_H4_FUNDED_SUCCESSOR_SUPPORT_SCOUT',
        'support_rows': len(saved), 'support_columns': len(columns), 'support_nonzeros': len(vv),
        'solver_numeric_objective': float(result.fun), 'exact_primal_upper': str(upper),
        'exact_dual_lower': str(lower), 'transfer_value': str(transfer),
        'primal_at_most_40': upper <= 40, 'dual_above_40': lower > 40,
        'positive_bonus_values': [[j, str(a)] for j, a in sorted(values.items()) if j < Z and a],
        'exact_dual_rows': [[i, str(y)] for i, y in sorted(weights.items())],
        'retained_numeric_dual': [[i, float(y).hex()] for i, y in enumerate(result.ineqlin.marginals) if y],
        'dual_face_rows_clipped': clipped, 'primal_bonus_scaling': str(scale),
        'solver_warnings': [str(w.message)[:180] for w in caught[:2]],
        'full_model_primal_verified': False, 'new_dimension_bound': None,
        'scope': 'Only the saved D5 successor support plus the proved H4 transfer and three source '
                 'forests. Full-model primal separation is outstanding; no dimension bound follows.',
    }
    verify_primal(data); verify_dual(data)
    data['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    assert data['worker_peak_rss_bytes'] < LIMIT
    print(json.dumps(data, separators=(',', ':')), flush=True)


def guarded_worker(lock_fd):
    before = S.G.resources()
    S.process_tree_rss()
    with Path(S.G.__file__).open('rb') as stream:
        guard_hash = file_digest(stream, 'sha256').hexdigest()
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[key] = '1'
    env['COXETER_SUPPORT_SCOUT_OWNER'] = str(os.getpid())
    child = None; output = bytearray(); start = time.monotonic(); owner = os.getppid(); peak = 0
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
                peak = max(peak, S.process_tree_rss())
                assert peak < LIMIT and time.monotonic()-start < 59
                assert os.getppid() == owner
                S.G.resources(start=False)
        assert child.wait() == 0, bytes(output[-4096:]).decode(errors='replace')
        data = json.loads(output)
        verify_primal(data); verify_dual(data)
        data['guard'] = {'elapsed_seconds': round(time.monotonic()-start, 3),
                         'process_tree_peak_rss_bytes': peak, 'rss_limit_bytes': LIMIT,
                         'wall_limit_seconds': 59, 'output_bytes': len(output),
                         'thread_limit': 1, 'initial_resources': before}
        data['tooling_provenance'] = {'guard_module': 'scripts/run_lean_batches_20261018.py',
                                      'guard_sha256_at_start': guard_hash}
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
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(); return
    signal.signal(signal.SIGTERM, S.G.interrupt)
    signal.signal(signal.SIGINT, S.G.interrupt)
    target = ROOT/f'results/{TAG}.json'
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.verify:
            data = json.loads(target.read_text())
            verify_primal(data); verify_dual(data)
            for path, digest in data['sources'].items():
                with (ROOT/path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
        else:
            assert not target.exists()
            frozen = json.loads((ROOT/GEOMETRY).read_text())
            for path, digest in frozen['files'].items():
                with (ROOT/path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
            data = guarded_worker(lock.fileno())
            data['sources'] = {}
            for path in (f'scripts/{TAG}.py', GEOMETRY, PROBE, S.M.MODEL, S.SEED,
                         f'scripts/{S.TAG}.py', f'scripts/{D.TAG}.py'):
                with (ROOT/path).open('rb') as stream:
                    data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
            encoded = json.dumps(data, separators=(',', ':'))+'\n'
            assert len(encoded.encode()) < OUTPUT_LIMIT
            temporary = target.with_suffix('.json.tmp')
            try:
                with temporary.open('x') as stream:
                    stream.write(encoded); stream.flush(); os.fsync(stream.fileno())
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
    print(json.dumps({'status': data['status'], 'numeric_objective': data['solver_numeric_objective'],
                      'exact_primal_upper': data['exact_primal_upper'],
                      'exact_dual_lower': data['exact_dual_lower'], 'transfer_value': data['transfer_value'],
                      'primal_at_most_40': data['primal_at_most_40'], 'dual_above_40': data['dual_above_40'],
                      'new_dimension_bound': None}))


if __name__ == '__main__':
    main()
