#!/usr/bin/env python3
"""Complete four-hub finite-forest separation of the exact full-face primal."""
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
import fresh_compact21_h4_full_face_successor_20261018 as C
import fresh_compact21_shifted_seed_forest_dp_20261018 as K

P, ROOT, Z, T = C.P, C.ROOT, C.Z, C.T
TAG = 'fresh_compact21_full_candidate_forest_separation_20261018'
CANDIDATE = f'results/{C.TAG}.json'
OUTPUT_LIMIT = 1024**2
HUBS = ('A1', 'I2(5)', 'H3', 'A4')


class Rule(P.Rule):
    def __init__(self):
        super().__init__()
        with (ROOT/CANDIDATE).open() as stream:
            data = json.load(stream)
        self.values = {j: F(v) for j, v in data['positive_bonus_values']}
        self.transfer = F(data['transfer_value'])
        self.saved_upper = F(data['exact_primal_upper'])

    def index_or_none(self, key, context):
        if key not in self.codes:
            return None
        a, b, h, d = context
        packed = self.codes[key]*396+P.context_number(a, min(b, 2), min(h, 1), min(d, 2))
        return self.index.get(packed)

    def bonus(self, key, context):
        return self.values.get(self.index_or_none(key, context), F())


def profiles(rule):
    models = {name: g for name, g in P.G.models() if len(g) <= 21}
    assert len(models) == 68
    source = models['A4']
    source_codes = sorted({P.encoded(P.G.code(source, (r,))) for r in source})
    assert len(source_codes) == 2
    result = {}
    for name, graph in models.items():
        baseline = sum((P.G.base(graph, u, v) for u, v in combinations(graph, 2)), F())
        connected, a1, i25, h3, a4, pair = (Counter() for _ in range(6))
        for u, v in combinations(graph, 2):
            key = P.encoded(P.G.code(graph, (u, v)))
            if key in rule.codes:
                connected[key] += 1
        family = P.G.family(graph)
        for root in graph:
            rooted = P.encoded(P.G.code(graph, (root,)))
            for prefix, excluded, target in (
                    ('A1x', ('A1',), a1),
                    ('I25x', ('A1', 'I2(5)'), i25),
                    ('H3x', ('A1', 'I2(5)', 'H3'), h3)):
                key = prefix+family+':'+rooted
                if name not in excluded and key in rule.codes:
                    target[key] += 1
            if name not in HUBS:
                for sc in source_codes:
                    key = 'A4sx'+family+':'+P.encoded([sc, rooted])
                    if key in rule.codes:
                        # Each source-root orbit has two roots; normalization 4d.
                        a4[key] += F(1, 2)
        if name == 'A4':
            for u, v in product(graph, repeat=2):
                codes = sorted([P.encoded(P.G.code(graph, (u,))),
                                P.encoded(P.G.code(graph, (v,)))])
                key = 'A4pair:'+P.encoded(codes)
                if key in rule.codes:
                    pair[key] += F(1, 8)
        result[name] = len(graph), baseline, connected, a1, i25, h3, a4, pair
    return models, result


def enabled_tables(profile, context):
    _, _, connected, a1, i25, h3, a4, pair = profile
    a, b, h, d = context
    yield connected
    for enabled, table in ((a, a1), (b, i25), (h, h3), (d, a4), (d >= 2, pair)):
        if enabled:
            yield table


def component_values(rule, data, context):
    return {name: profile[1]+sum((n*rule.bonus(key, context)
                                  for table in enabled_tables(profile, context)
                                  for key, n in table.items()), F())
            for name, profile in data.items()}


def exact_context(forest):
    return tuple(forest.count(name) for name in HUBS)


def run():
    rule = Rule()
    models, data = profiles(rule)
    ordinary = {name: row for name, row in data.items() if name != 'A4'}
    targets = {
        'H4': tuple(sorted(('A1', 'H4', 'E8', 'E8'))),
        'D5': tuple(sorted(('D5', 'E8', 'E8'))),
        'A4': tuple(sorted(('A1', 'A4', 'E8', 'E8'))),
        'A5': tuple(sorted(('A5', 'E8', 'E8'))),
        'B4': tuple(sorted(('A1', 'B4', 'E8', 'E8'))),
        'B5': tuple(sorted(('B5', 'E8', 'E8'))),
    }
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
    beta_key = P.encoded(P.G.code(e8, (ends[2], ends[4])))
    beta = rule.bonus(beta_key, (0, 0, 0, 0))
    d5_saving = min(F(1), 2*beta)
    evaluated, component_cache = {}, {}
    contexts = cap_contexts = 0

    def adjusted(old, forest):
        if forest == targets['H4']:
            return old-rule.transfer
        if forest == targets['D5']:
            return old-d5_saving
        if forest in (targets['A4'], targets['A5'], targets['B4'], targets['B5']):
            return old+rule.transfer
        return old

    for a in range(22):
        for bcap, hcap, dcap in product(range(3), range(2), range(3)):
            remaining_max = 21-a-2*bcap-3*hcap-4*dcap
            if remaining_max < 0:
                continue
            context = a, bcap, hcap, dcap
            totals = component_values(rule, data, context)
            component_cache[context] = totals
            cap_contexts += 1
            dp = K.top_three_dp(ordinary, totals, remaining_max)
            for b in range((21-a)//2+1):
                for h in range((21-a-2*b)//3+1):
                    for d in range((21-a-2*b-3*h)//4+1):
                        if (min(b, 2), min(h, 1), min(d, 2)) != (bcap, hcap, dcap):
                            continue
                        rank = 21-a-2*b-3*h-4*d
                        if not dp[rank]:
                            continue
                        contexts += 1
                        fixed = sum((count*totals[name] for count, name in zip((a, b, h, d), HUBS)), F())
                        for value, forest in dp[rank]:
                            full = tuple(sorted(forest+('A1',)*a+('I2(5)',)*b+('H3',)*h+('A4',)*d))
                            evaluated[full] = adjusted(value+fixed, full)
    # Positively adjusted source forests may not be among the old top three.
    # Evaluate all six exceptions explicitly. Only two physical forests decrease.
    for forest in targets.values():
        a, b, h, d = exact_context(forest)
        totals = component_cache[a, min(b, 2), min(h, 1), min(d, 2)]
        evaluated[forest] = adjusted(sum((totals[name] for name in forest), F()), forest)
    maximum = max(evaluated.values())
    largest = sorted(evaluated, key=lambda f: (-evaluated[f], f))

    def literal_row(forest):
        context = exact_context(forest)
        co, baseline = Counter({Z: F(-1)}), F()
        for name in forest:
            baseline += data[name][1]
            for table in enabled_tables(data[name], context):
                for key, n in table.items():
                    j = rule.index_or_none(key, context)
                    if j is not None:
                        co[j] += n
        rhs, state = -baseline, None
        if forest == targets['H4']:
            co[T] -= 1
        elif forest == targets['D5']:
            if 2*beta <= 1:
                column = rule.index_or_none(beta_key, context)
                assert column is not None and co[column] == 2
                co[column] -= 2
                state = 'inactive'
            else:
                rhs += 1
                state = 'active'
        elif forest in (targets['A4'], targets['A5'], targets['B4'], targets['B5']):
            co[T] += 1
        row = {j: n for j, n in co.items() if n}
        value = sum((n*rule.values.get(j, F()) for j, n in row.items() if j not in (Z, T)), F())
        value += row.get(T, F())*rule.transfer-rhs
        assert value == evaluated[forest]
        return {'types': list(forest), 'context_A1_I25_H3_A4': list(context),
                'vertex_budget': str(value), 'coefficients': [[j, str(n)] for j, n in sorted(row.items())],
                'rhs': str(rhs), 'D5_state': state}

    violated = [f for f in largest if evaluated[f] > 40]
    retained = [literal_row(f) for f in violated[:64]]
    attainer = literal_row(largest[0])
    return {
        'status': 'PASS_EXACT_FULL_FINITE_FOREST_SEPARATION',
        'candidate': CANDIDATE, 'finite_component_types': len(models),
        'capped_contexts_checked': cap_contexts, 'exact_hub_contexts_checked': contexts,
        'rank_dp_retained_distinct_forests': 3, 'decreased_physical_forest_types': 2,
        'exact_full_forest_maximum': str(maximum), 'maximum_decimal': float(maximum),
        'at_most_40': maximum <= 40, 'strictly_below_40': maximum < 40,
        'saved_support_upper': str(rule.saved_upper), 'transfer_value': str(rule.transfer),
        'D5_extreme_bonus': str(beta), 'D5_uniform_saving': str(d5_saving),
        'attaining_forest': attainer, 'violating_retained_maximizers': retained,
        'retained_violation_count': len(retained),
        'distinct_top_three_forests_violating_40': len(violated),
        'all_connected_and_four_hub_cross_coefficients_included': True,
        'A4_pair_component_coefficient': 'one eighth of all sixteen ordered source-root pairs',
        'full_native_face_vector_already_verified': True, 'numeric_optimizer_run': False,
        'new_dimension_bound': None,
        'scope': 'Exact full finite-forest model separation. A dimension conclusion additionally '
                 'requires final independent review and universal arithmetic/face-coverage composition.',
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
    assert data['worker_peak_rss_bytes'] < C.H.LIMIT
    print(json.dumps(data, separators=(',', ':')), flush=True)


def guarded(lock_fd):
    before = C.H.S.G.resources(); C.H.S.process_tree_rss()
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
                peak = max(peak, C.H.S.process_tree_rss())
                assert peak < C.H.LIMIT and time.monotonic()-started < 59
                assert os.getppid() == owner
                C.H.S.G.resources(start=False)
        assert child.wait() == 0, bytes(output[-4096:]).decode(errors='replace')
        data = json.loads(output)
        data['guard'] = {'elapsed_seconds': round(time.monotonic()-started, 3),
                         'process_tree_peak_rss_bytes': peak, 'rss_limit_bytes': C.H.LIMIT,
                         'wall_limit_seconds': 59, 'output_bytes': len(output), 'initial_resources': before}
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
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
        signal.signal(sig, C.H.S.G.interrupt)
    signal.alarm(59)
    target = ROOT/f'results/{TAG}.json'
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        data = guarded(lock.fileno())
        data['sources'] = {}
        for path in (f'scripts/{TAG}.py', CANDIDATE, f'scripts/{K.TAG}.py',
                     f'scripts/{P.TAG}.py', C.H.S.M.MODEL,
                     'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
                     'scripts/fresh_sqrt5_dimension21_a4_source_orbits_20261016.py', C.CAPACITY):
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
    signal.alarm(0)
    print(json.dumps({k: data[k] for k in ('status', 'exact_hub_contexts_checked',
                                          'exact_full_forest_maximum', 'at_most_40',
                                          'retained_violation_count', 'new_dimension_bound')}))


if __name__ == '__main__':
    main()
