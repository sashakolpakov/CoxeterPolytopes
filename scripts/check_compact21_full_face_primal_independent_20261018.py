#!/usr/bin/env python3
"""Independent exact evaluation of every native rank21 old face inequality.

Does not import the full-face producer, its Rule/corner helpers or its solver.
Only the frozen finite-diagram graph primitives and streaming input reader
are shared. The coefficient evaluator and spectator loops are reconstructed.
"""
from collections import Counter
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
from itertools import combinations
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_sqrt5_forest_sensitive_weights_20261015_review as G
import fresh_compact21_extreme_e8_priority_saved_dual_20261018 as M
from run_lean_batches_signal_safe_20261018 import resources

ROOT = Path(__file__).resolve().parents[1]
TAG = 'check_compact21_full_face_primal_independent_20261018'
CANDIDATE = 'results/fresh_compact21_h4_full_face_successor_20261018.json'
BASE = 'results/fresh_sqrt5_parameter_weight_stress_20261015.json'
QUEUE = 'results/fresh_sqrt5_nonsingular_rank29_fast_20261015.json'
LIMIT = 1024**3


def memory_gate():
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert peak*(1 if sys.platform == 'darwin' else 1024) < LIMIT


def objects(path, key):
    active, lines, size = False, None, 0
    with (ROOT/path).open() as stream:
        for line in stream:
            if not active:
                active = line.startswith('  "'+key+'": [')
                continue
            if line.startswith('  ],'):
                return
            if line.startswith('    {'):
                assert lines is None
                lines, size = [], 0
            assert lines is not None
            lines.append(line); size += len(line)
            assert size < 2*1024**2
            if line.startswith('    }'):
                yield json.loads(''.join(lines).rstrip().rstrip(','))
                lines = None
    assert active and lines is None


def encode(value):
    return json.dumps(value, separators=(',', ':'))


def angle(row, edges, marks):
    blocks = G.pieces(G.graph(row['white']+list(marks), edges))
    left = next(block for block in blocks if marks[0] in block)
    right = next(block for block in blocks if marks[1] in block)
    names = [G.family(block) for block in blocks]
    local = tuple(names.count(k) for k in ('A1', 'I2(5)', 'H3', 'A4'))
    if left is right:
        base = F(G.base(left, *marks))
        return base, local, encode(G.code(left, marks)) if base == 0 else None, None
    lf, rf = G.family(left), G.family(right)
    for hub, prefix, slot, multiplier in (
            ('A1', 'A1x', 0, 1), ('I2(5)', 'I25x', 1, 2), ('H3', 'H3x', 2, 3)):
        prior = ('A1', 'I2(5)', 'H3')[:slot+1]
        if lf == hub and rf not in prior:
            return F(), local, prefix+rf+':'+encode(G.code(right, (marks[1],))), (slot, multiplier, 0)
        if rf == hub and lf not in prior:
            return F(), local, prefix+lf+':'+encode(G.code(left, (marks[0],))), (slot, multiplier, 0)
    if lf == rf == 'A4':
        codes = sorted((encode(G.code(left, (marks[0],))), encode(G.code(right, (marks[1],)))))
        return F(), local, 'A4pair:'+encode(codes), (3, 4, -1)
    if lf == 'A4' or rf == 'A4':
        source, sr, target, tr = (left, marks[0], right, marks[1]) if lf == 'A4' else (right, marks[1], left, marks[0])
        family = G.family(target)
        assert family not in ('A1', 'I2(5)', 'H3', 'A4')
        key = 'A4sx'+family+':'+encode([encode(G.code(source, (sr,))), encode(G.code(target, (tr,)))])
        return F(), local, key, (3, 4, 0)
    return F(), local, None, None


def main():
    with (ROOT/CANDIDATE).open() as stream:
        saved = json.load(stream)
    for path, want in saved['sources'].items():
        with (ROOT/path).open('rb') as stream:
            assert file_digest(stream, 'sha256').hexdigest() == want
    indexed = {j: F(v) for j, v in saved['positive_bonus_values']}
    assert len(indexed) == len(saved['positive_bonus_values'])
    assert all(0 <= j < 93045 and v > 0 for j, v in indexed.items())
    values, records, used = {}, 0, set()
    for j, rec in M.objects('variable_order'):
        records += 1
        if j in indexed:
            key = (rec['orbit'], rec['A1_count'], rec['I25_context'], rec['H3_context'], rec['A4_context'])
            assert key not in values
            values[key] = indexed[j]; used.add(j)
        if records % 4096 == 0:
            memory_gate()
    assert records == 93045 and used == set(indexed)
    def faces():
        for rec in objects(BASE, 'faces'):
            if rec['row']['nodes'] <= 23:
                yield rec['row']
        for rec in objects(QUEUE, 'records'):
            if rec['result']['singular'] and rec['row']['nodes'] <= 23:
                yield rec['row']
    cases, count, minimum, witness = 0, 0, None, None
    for face, row in enumerate(faces()):
        count += 1
        edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
        angles = [angle(row, edges, marks) for marks in ((0, 1), (0, 3), (2, 1), (2, 3))]
        baseline = sum((a[0] for a in angles), F())
        rank = 23-row['nodes']
        for a in range(rank+1):
            for b in range((rank-a)//2+1):
                for h in range((rank-a-2*b)//3+1):
                    for d in range((rank-a-2*b-3*h)//4+1):
                        if rank-a-2*b-3*h-4*d == 1:
                            continue
                        outside = (a, b, h, d)
                        payment = baseline
                        for _, local, orbit, denominator in angles:
                            if orbit is None:
                                continue
                            aa, bb, hh, dd = (x+y for x, y in zip(local, outside))
                            v = values.get((orbit, aa, min(bb, 2), min(hh, 1), min(dd, 2)), F())
                            if denominator is not None:
                                slot, multiplier, offset = denominator
                                divisor = multiplier*((aa, bb, hh, dd)[slot]+offset)
                                assert divisor > 0
                                v /= divisor
                            payment += v
                        assert payment >= 1, (face, outside, str(payment))
                        cases += 1
                        if minimum is None or payment < minimum:
                            minimum, witness = payment, (face, outside)
        if face % 32 == 0:
            memory_gate(); resources(start=False)
    assert count == 1979 and cases == 102917
    memory_gate()
    print('PASS_INDEPENDENT_FULL_NATIVE_FACE_PRIMAL', count, cases, records,
          len(values), 'minimum_payment='+str(minimum), 'witness='+str(witness), flush=True)


if __name__ == '__main__':
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        signal.alarm(59); resources(); main(); resources(start=False)
