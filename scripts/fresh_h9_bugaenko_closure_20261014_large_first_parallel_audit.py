#!/usr/bin/env python3
"""Replay a saved batch with separate finite and modular arithmetic.

The DFS, explicit canonical-state permutation, and Lanner tables are reused.
Workers only replay saved inputs; this does not enlarge the search scope.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_closure_20261014_large_first_parallel_audit'
INPUT = 'results/fresh_h9_bugaenko_closure_20261014_large_first_0_505.json'
INC = 'results/fresh_h9_bugaenko_surgery_20261014_incidence.json'


def replay(task):
    import fresh_h9_bugaenko_closure_20261014 as C
    import fresh_h9_bugaenko_closure_20261014_large_first as P
    import fresh_h9_bugaenko_closure_20261014_large_first_audit as A
    row, expected, limit = task
    C.finite = P.finite = A.finite
    C.residue = A.determinant
    actual = json.loads(json.dumps(P.search(row, limit)))
    assert actual == expected, expected['case']
    return dict(case=actual['case'], verdict=actual['verdict'],
                evaluated_states=actual['nodes'], visits=actual['visits'])


def run(workers):
    saved = json.loads((ROOT / INPUT).read_text())
    for p, h in saved['sources'].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    incidence = json.loads((ROOT / INC).read_text())
    tasks = [(incidence['cases'][r['case']], r, saved['node_limit'])
             for r in saved['cases']]
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for row in pool.map(replay, tasks, chunksize=1):
            rows.append(row)
            if len(rows) % 50 == 0:
                print('replayed', len(rows), flush=True)
    files = sorted(set(saved['sources']) | {INPUT, INC,
        str(Path(__file__).resolve().relative_to(ROOT)),
        'scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py',
        'scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py'})
    return dict(status='PASS_SEPARATE_ARITHMETIC_SAVED_BATCH_REPLAY',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        input=INPUT, cases_replayed=len(rows),
        exact_exclusions=sum(r['verdict'] == 'EXACT_COMPLETE_OBSTRUCTION' for r in rows),
        limited_cases=[r['case'] for r in rows if r['verdict'] == 'NODE_LIMIT_UNRESOLVED'],
        records=rows,
        scope='Separate finite Coxeter classifier and SymPy finite-field determinants; '
              'the DFS, canonical-state construction, incidence catalogue and Lanner tables '
              'are reused and source-bound. Every saved record is reproduced, including '
              'the three node-limited records. Limited records are not exclusions.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    assert 1 <= args.workers <= 4
    print('PID', os.getpid(), flush=True)
    result = run(args.workers)
    target = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['cases_replayed'], result['exact_exclusions'], flush=True)
