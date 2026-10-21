#!/usr/bin/env python3
"""Separate arithmetic replay of newly needed nonstructural exclusions."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_search_audit'
INC = 'results/fresh_h9_bugaenko_surgery_20261014_incidence.json'
INPUTS = ['results/fresh_h9_bugaenko_remaining_20261014_' + x + '.json'
          for x in ['650_1150', '1150_1950', 'selected',
                    'propagated_1506_1507', 'propagated_2421_2422']]
GATE = 'results/fresh_h9_bugaenko_remaining_20261014_block_gate.json'
OLD = 'results/fresh_h9_bugaenko_closure_20261014_checkpoint.json'


def replay(task):
    import fresh_h9_bugaenko_closure_20261014 as C
    import fresh_h9_bugaenko_closure_20261014_large_first as P
    import fresh_h9_bugaenko_remaining_20261014_propagated as Q
    import fresh_h9_bugaenko_closure_20261014_large_first_audit as A
    row, expected, limit, propagated = task
    C.finite = P.finite = Q.finite = A.finite
    C.residue = A.determinant
    actual = json.loads(json.dumps((Q if propagated else P).search(row, limit)))
    assert actual == expected, expected['case']
    return dict(case=actual['case'], propagated=propagated,
                nodes=actual['nodes'], visits=actual['visits'], verdict=actual['verdict'])


def run():
    files = INPUTS + [GATE, OLD, INC]
    data = {p: json.loads((ROOT / p).read_text()) for p in files}
    for d in data.values():
        for p, h in d['sources'].items():
            assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    gate = {r['case'] for r in data[GATE]['applications']}
    old = set(map(int, data[OLD]['exact_exclusion_sources']))
    todo = {}
    for p in INPUTS:
        d = data[p]
        for r in d['cases']:
            if r['case'] in gate or r['case'] in old or r['verdict'] != 'EXACT_COMPLETE_OBSTRUCTION':
                continue
            todo[r['case']] = (data[INC]['cases'][r['case']], r, d['node_limit'], 'propagated' in p)
    assert len(todo) == 476
    assert gate | old | todo.keys() == set(range(2953))
    records = []
    with ProcessPoolExecutor(max_workers=4) as pool:
        for row in pool.map(replay, [todo[i] for i in sorted(todo)], chunksize=1):
            records.append(row)
            if len(records) % 50 == 0:
                print('replayed', len(records), flush=True)
    files += [str(Path(__file__).resolve().relative_to(ROOT)),
        'scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py',
        'scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
    return dict(status='PASS_SEPARATE_ARITHMETIC_REPLAY_476_NEEDED_EXCLUSIONS',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        records=records, newly_checked=476, inherited_nonstructural_exclusions=len(old - gate),
        structural_exclusions=len(gate), catalogue_coverage=2953,
        scope='Separate finite-tree classification and SymPy finite-field determinants, '
              'with the corresponding saved DFS/canonicalizer and relation propagation '
              'reused. Only the 476 newly necessary nonstructural records are rerun; '
              '112 prior nonstructural exclusions retain their frozen separate replay. '
              'All saved proof records and counters must match exactly.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    a = ap.parse_args()
    print('PID', os.getpid(), flush=True)
    data = run()
    target = ROOT / f'results/{TAG}.json'
    if a.verify:
        assert data == json.loads(target.read_text())
    else:
        target.write_text(json.dumps(data, indent=2) + '\n')
    print(data['status'], flush=True)
