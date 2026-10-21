#!/usr/bin/env python3
"""Recover a rational complete-model dual by acyclic flow and face clipping.

No optimizer or numerical package is imported. All global columns, including
eliminated bonus columns and Bellman potentials, are checked by weak duality.
"""
import argparse
from collections import Counter, defaultdict
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
import fresh_compact21_full_forest_epigraph_20261018 as E

C, D, P, ROOT, Z, T = E.C, E.D, E.P, E.ROOT, E.Z, E.T
TAG = 'fresh_compact21_complete_model_dual_repair_20261018'
SOURCE = 'results/fresh_compact21_epigraph_ipm_minimum_20261018.json'
OUTPUT_LIMIT, ARTIFACT_LIMIT, GRID = 512*1024, 8*1024**2, 10**10


def model():
    base, coverage = C.build_rows()
    rule = D.Rule(); models, profiles = D.profiles(rule)
    assert E.check_existing_rows(rule, profiles) == 510
    grouped, auxiliary, next_column = E.groups(), {}, T+1
    for context, tuples in sorted(grouped.items()):
        for rank in range(2, max(t[-1] for t in tuples)+1):
            auxiliary[context, rank] = next_column; next_column += 1
    rows = list(E.extended_rows(base, rule, profiles, grouped, auxiliary))
    assert len(rows) == 79274 and next_column == 94706
    return rows, coverage, models, profiles, grouped, auxiliary, next_column


def metadata(base_length, profiles, grouped):
    """Annotate final budget rows; no mathematical row is supplied by this helper."""
    result = {}
    old_metadata = C.H.D.forest_metadata()
    position = 51443
    for old, (_, _, kind) in enumerate(C.H.rows()):
        if kind == 'face':
            continue
        if kind == 'forest':
            if old in old_metadata:
                result[position] = dict(old_metadata[old])
            else:
                name = {2805: 'A5', 2806: 'B4', 2807: 'B5'}[old]
                result[position] = {'types': list(E.SPECIAL[name])}
        position += 1
    assert position == base_length
    items = [name for name in profiles if name not in E.HUBS]
    for context, tuples in sorted(grouped.items()):
        for rank in range(2, max(t[-1] for t in tuples)+1):
            for name in items:
                remaining = rank-profiles[name][0]
                if remaining >= 0 and remaining != 1:
                    position += 1
        for a, b, h, d, rank in tuples:
            if (a, b, h, d) in ((1, 0, 0, 0), (0, 0, 0, 0)):
                target = {'H4', 'E8'} if a else {'D5', 'E8'}
                for name in items:
                    remaining = rank-profiles[name][0]
                    if name in target or remaining < 0 or remaining == 1:
                        continue
                    result[position] = {'hub_tuple': [a, b, h, d],
                                        'residual_rank': rank, 'reserved_component': name}
                    position += 1
            else:
                result[position] = {'hub_tuple': [a, b, h, d], 'residual_rank': rank}
                position += 1
    for forest in (tuple(sorted(('A1',)+('H4',)*5)),
                   tuple(sorted(('A1',)+('H4',)*3+('E8',)))):
        result[position] = {'types': list(forest)}; position += 1
    for name, forest in E.SPECIAL.items():
        for state in (('active', 'inactive') if name == 'D5' else (None,)):
            result[position] = {'types': list(forest), 'state': state, 'special_name': name}
            position += 1
    assert position == 79274
    return result


def budget_witness(info, rule, profiles, cache):
    if 'types' in info:
        return info['types']
    a, b, h, d = context = tuple(info['hub_tuple'])
    if context not in cache:
        values = D.component_values(rule, profiles, context)
        states = {0: (F(), ())}
        for rank in range(2, 22):
            choices = [(states[rank-profile[0]][0]+values[name],
                        tuple(sorted(states[rank-profile[0]][1]+(name,))))
                       for name, profile in profiles.items()
                       if name not in E.HUBS and rank-profile[0] in states]
            assert choices
            states[rank] = max(choices)
        cache[context] = states
    reserved = info.get('reserved_component')
    rank = info['residual_rank']-(profiles[reserved][0] if reserved else 0)
    tail = cache[context][rank][1]
    hubs = ('A1',)*a+('I2(5)',)*b+('H3',)*h+('A4',)*d
    forest = tuple(sorted(hubs+tail+((reserved,) if reserved else ())))
    assert sum(profiles[name][0] for name in forest) == 21
    return list(forest)


def verify(data, rows, columns):
    multipliers = {i: F(y) for i, y in data['exact_dual_rows']}
    assert len(multipliers) == len(data['exact_dual_rows'])
    assert all(0 <= i < len(rows) and y < 0 for i, y in multipliers.items())
    loads, objective = defaultdict(F), F()
    for i, y in multipliers.items():
        row, rhs, _ = rows[i]
        objective += y*rhs
        for j, a in row.items():
            assert 0 <= j < columns
            loads[j] += y*a
    assert all(loads[j] <= F(j == Z) for j in range(columns))
    assert objective == F(data['exact_lower_bound'])
    assert objective <= F(data['source_exact_primal_upper'])
    return objective


def recover():
    started = time.monotonic()
    with (ROOT/SOURCE).open() as stream:
        source = json.load(stream)
    rows, coverage, models, profiles, grouped, auxiliary, columns = model()
    masses, rounded_faces, positive_noise = {}, {}, []
    for i, encoded in source['retained_numeric_dual']:
        raw = F.from_float(float.fromhex(encoded))
        assert raw <= F(1, 10**7)
        if raw > 0:
            positive_noise.append(raw)
        scaled = max(F(), -raw)*GRID
        numerator = scaled.numerator//scaled.denominator
        if rows[i][2] not in ('face', 'bellman') and numerator*scaled.denominator < scaled.numerator:
            numerator += 1
        if numerator:
            masses[i] = F(numerator, GRID)
            if rows[i][2] == 'face':
                rounded_faces[i] = masses[i]
    capacity = defaultdict(F)
    incoming = defaultdict(list)
    for i, mass in masses.items():
        row, _, kind = rows[i]
        if kind == 'face':
            continue
        for j, a in row.items():
            if j < Z:
                assert a >= 0
            capacity[j] += mass*a
            if j > T and a < 0:
                assert kind == 'bellman' and a == -1
                incoming[j].append(i)
    # The t cap may have numerical multiplier zero and hence be absent above.
    cap_row = next(i for i, (row, rhs, kind) in enumerate(rows)
                   if kind == 'cap' and row == {T: F(1)} and rhs == 1)
    flow_changes = Counter()
    for (context, rank), j in sorted(auxiliary.items(), key=lambda item: (-item[0][1], item[0][0])):
        if capacity[j] >= 0:
            continue
        for i in sorted(incoming[j], key=lambda i: (masses[i], i)):
            if capacity[j] >= 0:
                break
            take = min(masses[i], -capacity[j])
            masses[i] -= take; flow_changes[i] += take
            for k, a in rows[i][0].items():
                capacity[k] -= take*a
        assert capacity[j] >= 0
    assert all(capacity[j] >= 0 for j in auxiliary.values())
    transfer_cap_added = max(F(), -capacity[T])
    if transfer_cap_added:
        masses[cap_row] = masses.get(cap_row, F())+transfer_cap_added
        capacity[T] += transfer_cap_added
    demand, incident = defaultdict(F), defaultdict(list)
    for i, mass in rounded_faces.items():
        for j, a in rows[i][0].items():
            assert j < Z and a <= 0
            if a:
                demand[j] -= mass*a; incident[j].append(i)
    violations = sorted((j for j in demand if demand[j] > capacity[j]),
                        key=lambda j: (-(demand[j]-capacity[j]), j))
    for j in violations:
        for i in sorted(incident[j], key=lambda i: (masses[i], i)):
            excess = demand[j]-capacity[j]
            if excess <= 0:
                break
            take = min(masses[i], excess/-rows[i][0][j])
            masses[i] -= take
            for k, a in rows[i][0].items():
                demand[k] += take*a
        assert demand[j] <= capacity[j]
    assert all(demand[j] <= capacity[j] for j in demand)
    normalizer = max(F(1), -capacity[Z])
    weights = {i: -mass/normalizer for i, mass in masses.items() if mass}
    bound = sum((y*rows[i][1] for i, y in weights.items()), F())
    info = metadata(51954, profiles, grouped)
    rule = E.rule_from_data(source); cache = {}
    largest = []
    budgets = sorted(((i, -y) for i, y in weights.items() if i in info),
                     key=lambda item: (-item[1], item[0]))
    for i, mass in budgets[:32]:
        largest.append({'global_row': i, 'positive_mass': str(mass),
                        'kind': rows[i][2], **info[i],
                        'representative_forest': budget_witness(info[i], rule, profiles, cache)})
    data = {'status': 'PASS_EXACT_COMPLETE_MODEL_DUAL', 'exact_lower_bound': str(bound),
            'source_exact_primal_upper': source['exact_primal_upper'],
            'source_numerical_upper': source['numerical_upper'],
            'proves_model_infeasible_at_40': bound > 40, 'native_coverage': coverage,
            'global_rows': len(rows), 'global_columns': columns,
            'rounding_grid': GRID, 'normalizer': str(normalizer),
            'positive_numeric_multipliers_clipped': len(positive_noise),
            'maximum_positive_numeric_multiplier': str(max(positive_noise, default=F())),
            'repaired_bellman_rows': len(flow_changes),
            'bellman_mass_removed': str(sum(flow_changes.values(), F())),
            'transfer_cap_mass_added': str(transfer_cap_added),
            'initial_overloaded_bonus_columns': len(violations),
            'clipped_face_rows': sum(masses[i] != v for i, v in rounded_faces.items()),
            'exact_dual_rows': [[i, str(y)] for i, y in sorted(weights.items())],
            'largest_budget_masses': largest,
            'budget_witness_scope': 'Representative exact-primal maximizer of that final budget; '
                                    'not an asserted decomposition of its dual flow into forests.',
            'new_dimension_bound': None, 'numeric_optimizer_run': False,
            'scope': 'Exact obstruction to this complete conservative payment model only; '
                     'no nonexistence statement about polytopes or other incidence rules.',
            'algorithm': 'Round face/Bellman masses down, budget masses up; clip deficient '
                         'auxiliary capacities in descending rank; repair t by its cap; '
                         'monotonically clip face demands; normalize z mass and verify every column.'}
    verify(data, rows, columns)
    data['worker_elapsed_seconds'] = round(time.monotonic()-started, 3)
    data['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    assert data['worker_peak_rss_bytes'] < C.H.LIMIT
    return data


def worker():
    owner = int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert os.getppid() == owner
    def owner_watch():
        while True:
            if os.getppid() != owner:
                os.killpg(os.getpgrp(), signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch, daemon=True).start()
    data = recover(); path = Path(os.environ['COXETER_DUAL_TEMP'])
    size = C.bounded_write(path, data)
    with path.open('rb') as stream:
        digest = file_digest(stream, 'sha256').hexdigest()
    print(json.dumps({'certificate_bytes': size, 'certificate_sha256': digest,
                      'exact_lower_bound': data['exact_lower_bound']}), flush=True)


def guarded(lock_fd, temporary):
    before = C.H.S.G.resources(); C.H.S.process_tree_rss()
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        env[key] = '1'
    env['COXETER_SUPPORT_SCOUT_OWNER'] = str(os.getpid())
    env['COXETER_DUAL_TEMP'] = str(temporary)
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
        assert data['exact_lower_bound'] == summary['exact_lower_bound']
        data['guard'] = {'elapsed_seconds': round(time.monotonic()-started, 3),
                         'process_tree_peak_rss_bytes': peak, 'rss_limit_bytes': C.H.LIMIT,
                         'wall_limit_seconds': 59, 'output_bytes': len(output),
                         'initial_resources': before, 'thread_limit': 1}
        return data
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
            assert target.stat().st_size < ARTIFACT_LIMIT
            with target.open() as stream:
                data = json.load(stream)
            for path, digest in data['sources'].items():
                with (ROOT/path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
            rows, _, _, _, _, _, columns = model()
            verify(data, rows, columns); C.H.S.G.resources(start=False)
        else:
            assert not target.exists() and not temporary.exists() and not output_tmp.exists()
            try:
                data = guarded(lock.fileno(), temporary)
                data['sources'] = {}
                for path in (f'scripts/{TAG}.py', SOURCE, f'scripts/{E.TAG}.py',
                             f'scripts/{D.TAG}.py', f'scripts/{C.TAG}.py',
                             P.BASE, P.QUEUE, C.H.S.M.MODEL, C.CAPACITY,
                             'results/fresh_compact21_epigraph_ipm_minimum_20261018_freeze.json',
                             'results/compact21_full_forest_bellman_formulation_20261018.md'):
                    with (ROOT/path).open('rb') as stream:
                        data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
                C.bounded_write(output_tmp, data); os.replace(output_tmp, target)
            finally:
                temporary.unlink(missing_ok=True); output_tmp.unlink(missing_ok=True)
    signal.alarm(0)
    print(json.dumps({k: data[k] for k in ('status', 'exact_lower_bound',
                                          'proves_model_infeasible_at_40', 'new_dimension_bound')}))


if __name__ == '__main__':
    main()
