#!/usr/bin/env python3
"""Bounded continuation on previously untouched surgery types."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import json
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014'
INC = 'results/fresh_h9_bugaenko_surgery_20261014_incidence.json'


def worker(task):
    from fresh_h9_bugaenko_closure_20261014_large_first import search
    row, bound = task
    return search(row, bound)


def run(ids, bound, workers, output):
    inc = json.loads((ROOT / INC).read_text())
    for p, h in inc['sources'].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h
    files = [INC, str(Path(__file__).resolve().relative_to(ROOT)),
        'scripts/fresh_h9_bugaenko_closure_20261014_large_first.py',
        'scripts/fresh_h9_bugaenko_closure_20261014.py',
        'scripts/fresh_h9_highbase_construction_20261014_cell.py',
        'scripts/independent_h9_14_record30_a1_labels.py']
    manifest = {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files}
    rows = []
    start = time.monotonic()
    partial = output.with_name(output.stem + '_partial.json')
    def snapshot(complete):
        return dict(status='BOUNDED_LARGEST_FIRST_CONTINUATION' if complete else 'RUNNING',
            sources=manifest, selected_case_ids=ids, node_limit=bound,
            all_selected_cases_processed=complete,
            verdict_counts=dict(Counter(r['verdict'] for r in rows)),
            cases=sorted(rows, key=lambda r: r['case']),
            scope='The frozen largest-first necessary-label search is reused with its exact '
                  'finite, negative-index and fixed-minor gates. Completed exclusions cover '
                  'all finite labels and arbitrary dotted distances for the selected compact '
                  'incidence. Node-limited cases remain unresolved. A full-category survivor '
                  'requires actual Lorentzian metric solving.')
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(worker, (inc['cases'][i], bound)): i for i in ids}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            if row['verdict'] == 'FULL_CATEGORY_METRIC_REQUIRED':
                print('METRIC CANDIDATE', row['case'], json.dumps(row), flush=True)
            if len(rows) % 25 == 0 or len(rows) == len(ids):
                partial.write_text(json.dumps(snapshot(False), indent=2) + '\n')
                print(len(rows), '/', len(ids), dict(Counter(r['verdict'] for r in rows)),
                      'seconds', round(time.monotonic() - start, 2), flush=True)
    data = snapshot(True)
    output.write_text(json.dumps(data, indent=2) + '\n')
    partial.unlink(missing_ok=True)
    return data


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--start', type=int, default=650)
    ap.add_argument('--end', type=int, default=1150)
    ap.add_argument('--nodes', type=int, default=500)
    ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args()
    assert 0 <= a.start < a.end <= 2953 and a.nodes > 0 and 1 <= a.workers <= 4
    print('PID', os.getpid(), flush=True)
    target = ROOT / f'results/{TAG}_{a.start}_{a.end}.json'
    result = run(list(range(a.start, a.end)), a.nodes, a.workers, target)
    print(result['status'], result['verdict_counts'], flush=True)
