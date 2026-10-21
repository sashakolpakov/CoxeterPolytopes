#!/usr/bin/env python3
"""Compare streamed row reconstruction with frozen literal producer witnesses."""
import argparse
from collections import Counter, defaultdict
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
from itertools import combinations, product
import json
import os
import resource
import signal
import sys
import fresh_compact21_h4_successor_face_separation_20261018 as P

TAG = 'check_compact21_streamed_face_rows_20261018'


def run(start, count):
    expected = defaultdict(list)
    for ti, term in P.H.S.M.objects('terms', 4):
        if term['metadata']['kind'] != 'face':
            continue
        for w in term['metadata']['witnesses']:
            if start <= w['face'] < start+count:
                expected[w['face']].append((ti, w, term['coefficients'], term['rhs']))
    assert expected
    rule = P.Rule()
    checked, terms, faces = 0, set(), set()
    native = P.native_faces()
    try:
        for fi, (_, _, row) in enumerate(native):
            if fi >= start+count:
                break
            if fi not in expected:
                continue
            edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
            corners = [P.corner(row, edges, a, b) for a, b in product((0, 2), (1, 3))]
            baseline = sum((r[0] for r in corners), F())
            for ti, w, literal, rhs in expected[fi]:
                outside = tuple(w[k] for k in ('external_A1', 'external_I25', 'external_H3', 'external_A4'))
                co = Counter()
                for _, local, key, denominator, _ in corners:
                    counts = tuple(a+b for a, b in zip(local, outside))
                    if key is not None:
                        factor = F(1)
                        if denominator is not None:
                            k, multiplier, offset = denominator
                            factor /= multiplier*(counts[k]+offset)
                        co[rule.column(key, counts)] -= factor
                assert dict(co) == {j: F(v) for j, v in literal}, (ti, w)
                assert baseline-1 == F(rhs), (ti, w)
                checked += 1; terms.add(ti); faces.add(fi)
    finally:
        native.close()
    assert checked == sum(map(len, expected.values()))
    return {'status': 'PASS_STREAMED_ROWS_EQUAL_FROZEN_LITERAL_WITNESSES',
            'start_face': start, 'count_faces': count, 'literal_witness_cases_checked': checked,
            'distinct_literal_terms_checked': len(terms), 'distinct_native_faces_checked': len(faces),
            'new_dimension_bound': None}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--count', type=int, default=100)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    assert args.start >= 0 and 1 <= args.count <= 100
    signal.signal(signal.SIGALRM, P.H.S.G.interrupt)
    signal.alarm(59)
    with (P.ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        P.H.S.G.resources()
        data = run(args.start, args.count)
        data['sources'] = {}
        for path in (f'scripts/{TAG}.py', f'scripts/{P.TAG}.py', P.H.S.M.MODEL, P.BASE, P.QUEUE):
            with (P.ROOT/path).open('rb') as stream:
                data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
        target = P.ROOT/f'results/{TAG}_{args.start:04d}_{args.start+args.count:04d}.json'
        if args.verify:
            assert json.loads(target.read_text()) == data
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
        P.H.S.G.resources(start=False)
    signal.alarm(0)
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024) < 1024**3
    print(json.dumps({k: v for k, v in data.items() if k != 'sources'}))
