#!/usr/bin/env python3
"""Bounded source-level diagnosis of the exact dual's largest face consumers."""
from collections import Counter
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
from itertools import combinations, product
import json
import os
from pathlib import Path
import resource
import signal
import sys
import threading
import time
import fresh_compact21_complete_model_dual_repair_20261018 as R

C, D, P, ROOT = R.C, R.D, R.P, R.ROOT
TAG = 'fresh_compact21_dual_face_consumers_20261018'
DUAL = f'results/{R.TAG}.json'
PRIMAL = R.SOURCE
LIMIT = 1024**3


def run():
    with (ROOT/DUAL).open() as stream:
        dual = json.load(stream)
    with (ROOT/PRIMAL).open() as stream:
        primal = json.load(stream)
    rule = R.E.rule_from_data(primal)
    rows, _ = C.build_rows()
    weights = sorted(((i, -F(y)) for i, y in dual['exact_dual_rows'] if i < 51443),
                     key=lambda item: (-item[1], item[0]))[:32]
    targets = {(tuple(sorted(rows[i][0].items())), rows[i][1]):
               {'global_face_row': i, 'dual_mass': str(mass),
                'row': [[j, str(a)] for j, a in sorted(rows[i][0].items())],
                'rhs': str(rows[i][1]), 'occurrences': 0, 'witness': None}
               for i, mass in weights}
    assert len(targets) == 32
    del rows
    cases = faces = 0
    for fi, (family, index, row) in enumerate(P.native_faces()):
        faces += 1
        edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
        corners = [P.corner(row, edges, a, b) for a, b in product((0, 2), (1, 3))]
        baseline = sum((c[0] for c in corners), F())
        for outside in P.spectators(23-row['nodes']):
            cases += 1
            co, details = Counter(), []
            for marked, (base, local, key, denominator, families) in zip(
                    product((0, 2), (1, 3)), corners):
                counts = tuple(a+b for a, b in zip(local, outside))
                factor, column = F(1), None
                if key is not None:
                    if denominator is not None:
                        k, multiplier, offset = denominator
                        divisor = multiplier*(counts[k]+offset)
                        assert divisor > 0
                        factor /= divisor
                    column = rule.column(key, counts)
                    co[column] -= factor
                details.append({'marked': list(marked), 'native_forest': sorted(families),
                                'baseline': str(base), 'counts': list(counts),
                                'bonus_key': key, 'column': column,
                                'factor': str(factor),
                                'bonus_value': str(rule.values.get(column, F()))})
            identity = tuple(sorted(co.items())), baseline-1
            if identity not in targets:
                continue
            target = targets[identity]; target['occurrences'] += 1
            if target['witness'] is None:
                target['witness'] = {'native_face': fi, 'family': family, 'family_index': index,
                                     'row': row, 'outside_A1_I25_H3_A4': list(outside),
                                     'untracked_outside_rank': 23-row['nodes']-
                                     sum(n*k for n, k in zip(outside, (1, 2, 3, 4))),
                                     'corners': details}
    assert faces == 1979 and cases == 102917
    assert all(t['witness'] is not None for t in targets.values())
    return {'status': 'PASS_EXACT_DUAL_FACE_CONSUMER_DIAGNOSIS',
            'top_face_rows': sorted(targets.values(), key=lambda t: (-F(t['dual_mass']), t['global_face_row'])),
            'native_faces_checked': faces, 'spectator_cases_checked': cases,
            'geometry_exclusion_proved': False, 'new_dimension_bound': None,
            'scope': 'Each retained literal face row has a checked native source witness. '
                     'This is a diagnosis of consumers, not a proof that any case is realizable.'}


def main():
    started = time.monotonic()
    stop, failure = threading.Event(), []
    owner, pid = os.getppid(), os.getpid()
    scale = 1 if sys.platform == 'darwin' else 1024
    def monitor():
        next_resources = 0
        while not stop.wait(.1):
            try:
                assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*scale < LIMIT
                assert os.getppid() == owner
                if time.monotonic() >= next_resources:
                    C.H.S.G.resources(start=False); next_resources = time.monotonic()+1
            except BaseException as exc:
                failure.append(str(exc)); os.kill(pid, signal.SIGTERM); return
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
        signal.signal(sig, C.H.S.G.interrupt)
    signal.alarm(59)
    target = ROOT/f'results/{TAG}.json'; temporary = target.with_suffix('.json.tmp')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = C.H.S.G.resources()
        assert not target.exists() and not temporary.exists()
        watcher = threading.Thread(target=monitor, daemon=True); watcher.start()
        try:
            data = run(); assert not failure
            data['sources'] = {}
            for name in (f'scripts/{TAG}.py', DUAL, PRIMAL, f'scripts/{R.TAG}.py',
                         f'scripts/{C.TAG}.py', P.BASE, P.QUEUE, C.H.S.M.MODEL):
                with (ROOT/name).open('rb') as stream:
                    data['sources'][name] = file_digest(stream, 'sha256').hexdigest()
            data['resources'] = {'elapsed_seconds': round(time.monotonic()-started, 3),
                                 'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*scale,
                                 'initial': before, 'wall_limit_seconds': 59,
                                 'live_rss_limit_bytes': LIMIT, 'numeric_optimizer_run': False}
            assert data['resources']['peak_rss_bytes'] < LIMIT
            text = json.dumps(data, separators=(',', ':'))+'\n'
            assert len(text.encode()) < 1024**2
            with temporary.open('x') as stream:
                stream.write(text); stream.flush(); os.fsync(stream.fileno())
            os.replace(temporary, target)
        finally:
            stop.set(); watcher.join(timeout=2); temporary.unlink(missing_ok=True)
    signal.alarm(0)
    print(json.dumps({'status': data['status'], 'retained_rows': len(data['top_face_rows']),
                      'elapsed_seconds': data['resources']['elapsed_seconds'],
                      'peak_rss_bytes': data['resources']['peak_rss_bytes']}))


if __name__ == '__main__':
    main()
