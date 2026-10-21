#!/usr/bin/env python3
"""Exact finite-forest separation of the proved shifted rank-22 seed.

No optimizer: only finite component profiles, top-three rank DP and affine
line intersections for the single proved transfer parameter.
"""
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
import fresh_compact21_h4_successor_face_separation_20261018 as P

ROOT = P.ROOT
TAG = 'fresh_compact21_shifted_seed_forest_dp_20261018'
SEED = 'results/fresh_sqrt5_dimension22_h3_cross_20261016_cap2.json'
LIMIT, OUTPUT_LIMIT = 1024**3, 1024**2


def scalars():
    active = False
    with (ROOT/SEED).open() as stream:
        for line in stream:
            if not active:
                if line.startswith('    "exact_primal": ['):
                    active = True
                continue
            if line.startswith('    ]'):
                return
            yield F(json.loads(line.strip().rstrip(',')))
    assert active


class Seed:
    def __init__(self):
        vector = list(scalars())
        assert len(vector) == 30170 and vector[-1] == F(21977, 500)
        self.values, self.keys = {}, set()
        count = 0
        for j, r in P.objects(SEED, 'variable_order'):
            count += 1
            self.keys.add(r['orbit'])
            if vector[j]:
                self.values[r['orbit'], r['A1_count'], r['I25_context'], r['H3_context']] = vector[j]
        assert count == 30169
        assert all(v >= 0 and (v*10000).denominator == 1 for v in vector[:-1])

    def bonus(self, key, a, b, h):
        # One orthogonal A1 is added when applying the proved rank-22 rule.
        return self.values.get((key, a+1, min(b, 2), min(h, 1)), F())


def component_profiles(seed):
    models = {name: g for name, g in P.G.models() if len(g) <= 21}
    assert len(models) == 68
    profiles = {}
    for name, graph in models.items():
        baseline = sum((P.G.base(graph, a, b) for a, b in combinations(graph, 2)), F())
        connected, a1, i25, h3 = Counter(), Counter(), Counter(), Counter()
        for a, b in combinations(graph, 2):
            key = P.encoded(P.G.code(graph, (a, b)))
            if key in seed.keys:
                connected[key] += 1
        family = P.G.family(graph)
        for root in graph:
            rooted = P.encoded(P.G.code(graph, (root,)))
            for prefix, excluded, target in (
                    ('A1x', ('A1',), a1),
                    ('I25x', ('A1', 'I2(5)'), i25),
                    ('H3x', ('A1', 'I2(5)', 'H3'), h3)):
                key = prefix+family+':'+rooted
                if name not in excluded and key in seed.keys:
                    target[key] += 1
        profiles[name] = len(graph), baseline, connected, a1, i25, h3
    return models, profiles


def component_values(seed, profiles, a, b, h):
    result = {}
    for name, (_, baseline, connected, a1, i25, h3) in profiles.items():
        value = baseline+sum((n*seed.bonus(k, a, b, h) for k, n in connected.items()), F())
        for enabled, table in ((a, a1), (b, i25), (h, h3)):
            if enabled:
                value += sum((n*seed.bonus(k, a, b, h) for k, n in table.items()), F())
        result[name] = value
    return result


def top_three_dp(profiles, totals, rank):
    dp = {0: [(F(), ())]}
    for r in range(1, rank+1):
        candidates = {}
        for name, (size, *_) in profiles.items():
            if name in ('A1', 'I2(5)', 'H3') or size > r:
                continue
            for old, forest in dp.get(r-size, ()):
                ff = tuple(sorted(forest+(name,)))
                candidates[ff] = old+totals[name]
        dp[r] = sorted(((v, f) for f, v in candidates.items()),
                       key=lambda pair: (-pair[0], pair[1]))[:3]
    return dp


def run():
    seed = Seed()
    models, profiles = component_profiles(seed)
    h4 = tuple(sorted(('A1', 'H4', 'E8', 'E8')))
    d5 = tuple(sorted(('D5', 'E8', 'E8')))
    target_forests = {
        'H4': h4, 'D5': d5,
        'A4': tuple(sorted(('A1', 'A4', 'E8', 'E8'))),
        'A5': tuple(sorted(('A5', 'E8', 'E8'))),
        'B4': tuple(sorted(('A1', 'B4', 'E8', 'E8'))),
        'B5': tuple(sorted(('B5', 'E8', 'E8'))),
    }
    profiles_by_context, contexts, retained = {}, 0, []
    best_unchanged, unchanged_attainer = F(-1), None
    for a in range(22):
        for bcap, hcap in product(range(3), range(2)):
            totals = component_values(seed, profiles, a, bcap, hcap)
            profiles_by_context[a, bcap, hcap] = totals
            dp = top_three_dp(profiles, totals, 21-a)
            for b in range((21-a)//2+1):
                for h in range((21-a-2*b)//3+1):
                    if (min(b, 2), min(h, 1)) != (bcap, hcap):
                        continue
                    rank = 21-a-2*b-3*h
                    if not dp[rank]:
                        continue
                    contexts += 1
                    fixed = a*totals['A1']+b*totals['I2(5)']+h*totals['H3']
                    for value, forest in dp[rank]:
                        full = tuple(sorted(forest+('A1',)*a+('I2(5)',)*b+('H3',)*h))
                        old_total = value+fixed
                        retained.append((old_total, full, (a, b, h)))
                        if full not in (h4, d5) and old_total > best_unchanged:
                            best_unchanged, unchanged_attainer = old_total, full
    # There are only two decreased physical forests. The three best distinct
    # old forests at each rank/context therefore retain a best unchanged one.
    old_targets = {}
    for name, forest in target_forests.items():
        a, b, h = (forest.count(k) for k in ('A1', 'I2(5)', 'H3'))
        totals = profiles_by_context[a, min(b, 2), min(h, 1)]
        old_targets[name] = sum((totals[k] for k in forest), F())
    e8 = models['E8']
    fork = next(i for i in e8 if len(e8[i]) == 3)
    ends = {}
    for neighbor in e8[fork]:
        previous, current, length = fork, neighbor, 1
        while len(e8[current]) == 2:
            nxt = next(i for i in e8[current] if i != previous)
            previous, current, length = current, nxt, length+1
        ends[length] = current
    assert set(ends) == {1, 2, 4}
    beta_key = P.encoded(P.G.code(e8, (ends[2], ends[4])))
    beta = seed.bonus(beta_key, 0, 0, 0)
    d5_saving = min(F(1), 2*beta)
    lines = [
        ('unchanged', best_unchanged, F(), unchanged_attainer),
        ('H4', old_targets['H4'], F(-1), h4),
        ('D5', old_targets['D5']-d5_saving, F(), d5),
    ]
    for name, capacity in (('A4', 1), ('A5', 2), ('B4', 1), ('B5', 1)):
        lines.append((name, old_targets[name], F(capacity), target_forests[name]))
    candidates = {F(), F(1), F(1, 4)}
    for (_, a, m, _), (_, b, n, _) in combinations(lines, 2):
        if m != n:
            t = (b-a)/(m-n)
            if 0 <= t <= 1:
                candidates.add(t)
    best, transfer = min((max(a+m*t for _, a, m, _ in lines), t) for t in candidates)
    fixed_quarter = max(a+m*F(1, 4) for _, a, m, _ in lines)
    largest = sorted(retained, key=lambda r: (-r[0], r[1]))[:12]
    return {
        'status': 'PASS_EXACT_SHIFTED_N22_SEED_FULL_FOREST_DP',
        'seed': SEED, 'A1_index_shift': 1,
        'face_feasibility_scope': 'Proved lift by adding one whole A1 spectator; '
                                  'normalized A1 cross payments only increase. '
                                  'Inherited completed rank-22 coverage is assumed.',
        'component_types': len(models), 'exact_hub_contexts_checked': contexts,
        'rank_dp_retained_per_state': 3, 'decreased_physical_forests': 2,
        'best_transfer': str(transfer), 'best_full_forest_maximum': str(best),
        'best_maximum_decimal': float(best), 'at_most_40': best <= 40,
        'full_forest_maximum_at_transfer_quarter': str(fixed_quarter),
        'D5_extreme_E8_bonus': str(beta), 'D5_uniform_saving': str(d5_saving),
        'old_exceptional_forest_totals': {k: str(v) for k, v in old_targets.items()},
        'envelope_lines': [{'name': name, 'intercept': str(a), 'slope': str(m),
                            'forest': list(forest)} for name, a, m, forest in lines],
        'largest_old_retained_forests': [{'old_total': str(v), 'forest': list(f),
                                         'context_A1_I25_H3': list(c)} for v, f, c in largest],
        'candidate_transfer_values_checked': [str(t) for t in sorted(candidates)],
        'numeric_optimizer_run': False, 'new_dimension_bound': None,
    }


def worker():
    owner = int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert os.getppid() == owner
    def owner_watch():
        while True:
            if os.getppid() != owner:
                os.killpg(os.getpgrp(), signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=owner_watch, daemon=True).start()
    data = run()
    data['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    assert data['worker_peak_rss_bytes'] < LIMIT
    print(json.dumps(data, separators=(',', ':')), flush=True)


def guarded(lock_fd):
    before = P.H.S.G.resources()
    P.H.S.process_tree_rss()
    env = os.environ.copy(); env['COXETER_SUPPORT_SCOUT_OWNER'] = str(os.getpid())
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
                peak = max(peak, P.H.S.process_tree_rss())
                assert peak < LIMIT and time.monotonic()-started < 59
                assert os.getppid() == owner
                P.H.S.G.resources(start=False)
        assert child.wait() == 0, bytes(output[-4096:]).decode(errors='replace')
        data = json.loads(output)
        data['guard'] = {'elapsed_seconds': round(time.monotonic()-started, 3),
                         'process_tree_peak_rss_bytes': peak, 'rss_limit_bytes': LIMIT,
                         'wall_limit_seconds': 59, 'output_bytes': len(output),
                         'initial_resources': before}
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
    signal.signal(signal.SIGTERM, P.H.S.G.interrupt)
    signal.signal(signal.SIGINT, P.H.S.G.interrupt)
    target = ROOT/f'results/{TAG}.json'
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        data = guarded(lock.fileno())
        data['sources'] = {}
        for path in (f'scripts/{TAG}.py', SEED, f'scripts/{P.TAG}.py',
                     'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
                     'scripts/fresh_sqrt5_dimension22_h3_cross_20261016_review.py',
                     'results/fresh_compact21_h4_triangle_funded_transfer_20261018_freeze.json',
                     'results/fresh_compact21_extreme_e8_quad_consumer_completeness_20261018.json'):
            with (ROOT/path).open('rb') as stream:
                data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
        if args.verify:
            prior = json.loads(target.read_text())
            for record in (prior, data):
                record.pop('guard', None); record.pop('worker_peak_rss_bytes', None)
            assert data == prior
        else:
            assert not target.exists()
            temporary = target.with_suffix('.json.tmp')
            try:
                with temporary.open('x') as stream:
                    json.dump(data, stream, indent=2); stream.write('\n')
                    stream.flush(); os.fsync(stream.fileno())
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
    print(json.dumps({k: data[k] for k in ('status', 'exact_hub_contexts_checked',
                                          'best_transfer', 'best_full_forest_maximum',
                                          'at_most_40', 'new_dimension_bound')}))


if __name__ == '__main__':
    main()
