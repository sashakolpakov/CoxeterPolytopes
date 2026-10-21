#!/usr/bin/env python3
"""Bounded exact native-face separation; one row and one spectator case at a time."""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as F
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
import fresh_compact21_h4_funded_successor_scout_20261018 as H
import fresh_sqrt5_forest_sensitive_weights_20261015_review as G

ROOT = H.ROOT
TAG = 'fresh_compact21_h4_successor_face_separation_20261018'
BASE = 'results/fresh_sqrt5_parameter_weight_stress_20261015.json'
QUEUE = 'results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json'
CANDIDATE = f'results/{H.TAG}.json'
OUTPUT_LIMIT = 256*1024


def objects(path, key):
    active, index, item, size = False, -1, [], 0
    with (ROOT/path).open() as stream:
        for line in stream:
            if not active:
                if line.startswith('  "'+key+'": ['):
                    active = True
                continue
            if line.startswith('  ],'):
                return
            if line.startswith('    {'):
                index += 1
                item, size = [], 0
            item.append(line); size += len(line)
            assert size < 2*1024**2
            if line.startswith('    }'):
                yield index, json.loads(''.join(item).rstrip().rstrip(','))
                item, size = [], 0
    assert active


def native_faces():
    for i, rec in objects(BASE, 'faces'):
        if rec['row']['nodes']-2 <= 21:
            yield 'prior_catalogue', i, rec['row']
    for _, rec in objects(QUEUE, 'records'):
        if rec['result']['singular'] and rec['row']['nodes']-2 <= 21:
            yield 'complete_finite_queue', rec['index'], rec['row']


def encoded(obj):
    return json.dumps(obj, separators=(',', ':'))


def context_number(a, b, h, d):
    assert 0 <= a <= 21 and 0 <= b <= 2 and 0 <= h <= 1 and 0 <= d <= 2
    return ((a*3+b)*2+h)*3+d


class Rule:
    def __init__(self):
        data = json.loads((ROOT/CANDIDATE).read_text())
        self.values = {j: F(v) for j, v in data['positive_bonus_values']}
        self.codes, self.index = {}, {}
        for j, rec in H.S.M.objects('variable_order'):
            key = rec['orbit']
            oid = self.codes.setdefault(key, len(self.codes))
            ctx = context_number(rec['A1_count'], rec['I25_context'], rec['H3_context'], rec['A4_context'])
            packed = oid*396+ctx
            assert packed not in self.index
            self.index[packed] = j
        assert len(self.index) == 93045

    def column(self, key, counts):
        assert key in self.codes, ('Missing orbit', key)
        a, b, h, d = counts
        packed = self.codes[key]*396+context_number(a, min(b, 2), min(h, 1), min(d, 2))
        assert packed in self.index, ('Missing context', key, counts)
        return self.index[packed]


def corner(row, edges, a, b):
    parts = G.pieces(G.graph(row['white']+[a, b], edges))
    left = next(g for g in parts if a in g)
    right = next(g for g in parts if b in g)
    families = [G.family(g) for g in parts]
    counts = tuple(families.count(k) for k in ('A1', 'I2(5)', 'H3', 'A4'))
    if b in left:
        baseline = G.base(left, a, b)
        return baseline, counts, (encoded(G.code(left, (a, b))) if baseline == 0 else None), None, families
    baseline, key, denominator = F(), None, None
    lf, rf = G.family(left), G.family(right)
    if lf == 'A1' and rf != 'A1':
        key, denominator = 'A1x'+rf+':'+encoded(G.code(right, (b,))), (0, 1, 0)
    elif rf == 'A1' and lf != 'A1':
        key, denominator = 'A1x'+lf+':'+encoded(G.code(left, (a,))), (0, 1, 0)
    elif lf == 'I2(5)' and rf not in ('A1', 'I2(5)'):
        key, denominator = 'I25x'+rf+':'+encoded(G.code(right, (b,))), (1, 2, 0)
    elif rf == 'I2(5)' and lf not in ('A1', 'I2(5)'):
        key, denominator = 'I25x'+lf+':'+encoded(G.code(left, (a,))), (1, 2, 0)
    elif lf == 'H3' and rf not in ('A1', 'I2(5)', 'H3'):
        key, denominator = 'H3x'+rf+':'+encoded(G.code(right, (b,))), (2, 3, 0)
    elif rf == 'H3' and lf not in ('A1', 'I2(5)', 'H3'):
        key, denominator = 'H3x'+lf+':'+encoded(G.code(left, (a,))), (2, 3, 0)
    elif lf == rf == 'A4':
        codes = sorted([encoded(G.code(left, (a,))), encoded(G.code(right, (b,)))])
        key, denominator = 'A4pair:'+encoded(codes), (3, 4, -1)
    elif lf == 'A4' and rf not in ('A1', 'I2(5)', 'H3', 'A4'):
        codes = [encoded(G.code(left, (a,))), encoded(G.code(right, (b,)))]
        key, denominator = 'A4sx'+rf+':'+encoded(codes), (3, 4, 0)
    elif rf == 'A4' and lf not in ('A1', 'I2(5)', 'H3', 'A4'):
        codes = [encoded(G.code(right, (b,))), encoded(G.code(left, (a,)))]
        key, denominator = 'A4sx'+lf+':'+encoded(codes), (3, 4, 0)
    return baseline, counts, key, denominator, families


def spectators(rank):
    for a in range(rank+1):
        for b in range((rank-a)//2+1):
            for h in range((rank-a-2*b)//3+1):
                for d in range((rank-a-2*b-3*h)//4+1):
                    if rank-a-2*b-3*h-4*d != 1:
                        yield a, b, h, d


def batch(start, count):
    rule = Rule()
    cases = violations = paid_cases = faces = 0
    minimum = None
    cuts = {}
    exhausted = True
    iterator = native_faces()
    try:
        for fi, (family, index, row) in enumerate(iterator):
            if fi < start:
                continue
            if fi >= start+count:
                exhausted = False
                break
            faces += 1
            edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
            corners = [corner(row, edges, a, b) for a, b in product((0, 2), (1, 3))]
            baseline = sum((r[0] for r in corners), F())
            for outside in spectators(23-row['nodes']):
                cases += 1
                if baseline >= 1:
                    paid_cases += 1
                    continue
                co = Counter()
                for base, local, key, denominator, _ in corners:
                    counts = tuple(x+y for x, y in zip(local, outside))
                    if key is not None:
                        factor = F(1)
                        if denominator is not None:
                            k, multiplier, offset = denominator
                            divisor = multiplier*(counts[k]+offset)
                            assert divisor > 0
                            factor /= divisor
                        co[rule.column(key, counts)] += factor
                payment = baseline+sum((v*rule.values.get(j, F()) for j, v in co.items()), F())
                minimum = payment if minimum is None else min(minimum, payment)
                if payment >= 1:
                    continue
                violations += 1
                key = tuple(sorted(co.items())), baseline-1
                if key in cuts:
                    cuts[key]['occurrences_in_batch'] += 1
                    continue
                record = {'native_face': fi, 'family': family, 'family_index': index,
                          'external_counts_A1_I25_H3_A4': list(outside),
                          'corner_forests': [sorted(c[4]) for c in corners],
                          'payment': str(payment), 'deficit': str(1-payment),
                          'coefficients': [[j, str(-a)] for j, a in sorted(co.items())],
                          'rhs': str(baseline-1), 'occurrences_in_batch': 1}
                if len(cuts) < 64:
                    cuts[key] = record
                else:
                    weakest = min(cuts, key=lambda k: (F(cuts[k]['deficit']), cuts[k]['native_face']))
                    if F(record['deficit']) > F(cuts[weakest]['deficit']):
                        del cuts[weakest]
                        cuts[key] = record
    finally:
        iterator.close()
    assert faces > 0 and faces <= count
    result = {'status': 'PASS_EXACT_NATIVE_FACE_SEPARATION_BATCH',
              'candidate': CANDIDATE, 'start_face': start, 'requested_faces': count,
              'processed_faces': faces, 'end_face_exclusive': start+faces,
              'reached_native_inventory_end': exhausted, 'spectator_cases': cases,
              'baseline_paid_cases': paid_cases, 'violating_cases': violations,
              'minimum_unpaid_case_payment': None if minimum is None else str(minimum),
              'retained_distinct_cut_count': len(cuts),
              'retained_cuts': sorted(cuts.values(), key=lambda r: (-F(r['deficit']), r['native_face'])),
              'full_model_verified': False, 'new_dimension_bound': None,
              'scope': 'All old native quadrilateral inequalities in this exact range only. '
                       'The reviewed D5/H4 geometric rules preserve their sufficiency. '
                       'Remaining faces and finite forests are separate obligations.'}
    return result


def worker(start, count):
    owner = int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert owner == os.getppid()
    def owner_watch():
        while True:
            if os.getppid() != owner:
                os.killpg(os.getpgrp(), signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch, daemon=True).start()
    data = batch(start, count)
    data['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    assert data['worker_peak_rss_bytes'] < H.LIMIT
    print(json.dumps(data, separators=(',', ':')), flush=True)


def guarded(start_face, count, lock_fd):
    before = H.S.G.resources()
    H.S.process_tree_rss()
    env = os.environ.copy(); env['COXETER_SUPPORT_SCOUT_OWNER'] = str(os.getpid())
    child = None; output = bytearray(); started = time.monotonic(); owner = os.getppid(); peak = 0
    try:
        child = subprocess.Popen([sys.executable, '-B', str(Path(__file__).resolve()), '--worker',
                                  '--start', str(start_face), '--count', str(count)], cwd=ROOT, env=env,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
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
        data = json.loads(output)
        data['guard'] = {'elapsed_seconds': round(time.monotonic()-started, 3),
                         'process_tree_peak_rss_bytes': peak, 'rss_limit_bytes': H.LIMIT,
                         'wall_limit_seconds': 59, 'output_bytes': len(output),
                         'initial_resources': before, 'numeric_optimizer_run': False}
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
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--count', type=int, default=100)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    assert args.start >= 0 and 1 <= args.count <= 100
    if args.worker:
        worker(args.start, args.count); return
    signal.signal(signal.SIGTERM, H.S.G.interrupt)
    signal.signal(signal.SIGINT, H.S.G.interrupt)
    target = ROOT/f'results/{TAG}_{args.start:04d}_{args.start+args.count:04d}.json'
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        data = guarded(args.start, args.count, lock.fileno())
        data['sources'] = {}
        for path in (f'scripts/{TAG}.py', CANDIDATE, BASE, QUEUE, H.S.M.MODEL,
                     'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
                     'scripts/fresh_sqrt5_dimension21_a4_source_orbits_20261016.py'):
            with (ROOT/path).open('rb') as stream:
                data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
        if args.verify:
            prior = json.loads(target.read_text())
            for record in (prior, data):
                record.pop('guard', None); record.pop('worker_peak_rss_bytes', None)
            assert data == prior
        else:
            assert not target.exists()
            encoded_data = json.dumps(data, separators=(',', ':'))+'\n'
            assert len(encoded_data.encode()) < OUTPUT_LIMIT
            temporary = target.with_suffix('.json.tmp')
            try:
                with temporary.open('x') as stream:
                    stream.write(encoded_data); stream.flush(); os.fsync(stream.fileno())
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
    print(json.dumps({k: data[k] for k in ('status', 'start_face', 'processed_faces', 'spectator_cases',
                                          'violating_cases', 'retained_distinct_cut_count',
                                          'minimum_unpaid_case_payment', 'new_dimension_bound')}))


if __name__ == '__main__':
    main()
