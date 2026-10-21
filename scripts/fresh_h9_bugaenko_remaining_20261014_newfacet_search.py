#!/usr/bin/env python3
"""Prioritize a small changed-incidence bank and fit every full survivor."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_newfacet_search'
INPUT = 'results/fresh_h9_bugaenko_remaining_20261014_newfacet.json'


def score(row):
    c = Counter(m.bit_count() for m in row['canonical_minimal_nonfaces'])
    return (-c[2], c[4] + c[5], len(row['canonical_minimal_nonfaces']), len(row['canonical_vertices']))


def worker(row):
    from fresh_h9_bugaenko_closure_20261014_large_first import search
    return search(row, 500)


def run():
    data = json.loads((ROOT / INPUT).read_text())
    for p, h in data['sources'].items():
        assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    eligible = [data['cases'][i] for i in data['genuinely_new_unexcluded_types']]
    best = sorted(eligible, key=lambda r: (score(r), r['case']))[:24]
    chosen = {r['case'] for r in best}
    other = sorted(eligible, key=lambda r: (score(r)[1:], score(r)[0], r['case']))
    for row in other:
        if row['case'] not in chosen:
            best.append(row)
            chosen.add(row['case'])
            if len(best) == 32:
                break
    records = []
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(worker, row): row for row in best}
        for future in as_completed(futures):
            row = future.result()
            records.append(row)
            print(row['case'], row['verdict'], row['nodes'], flush=True)
    metrics = []
    for row in records:
        if row['verdict'] != 'FULL_CATEGORY_METRIC_REQUIRED':
            continue
        from fresh_h9_bugaenko_surgery_20261014_metric_scout import fit
        source = data['cases'][row['case']]
        vertices = [tuple(i for i in range(14) if v >> i & 1) for v in source['canonical_vertices']]
        result = fit(row['label_vector'], vertices, row['dotted_pairs'], 2000)
        metrics.append(dict(case=row['case'], **result))
        print('METRIC', row['case'], result['residual'], flush=True)
    files = [INPUT, str(Path(__file__).resolve().relative_to(ROOT)),
        'scripts/fresh_h9_bugaenko_closure_20261014_large_first.py',
        'scripts/fresh_h9_bugaenko_closure_20261014.py',
        'scripts/fresh_h9_highbase_construction_20261014_cell.py',
        'scripts/independent_h9_14_record30_a1_labels.py',
        'scripts/fresh_h9_bugaenko_surgery_20261014_metric_scout.py']
    return dict(status='BOUNDED_CHANGED_INCIDENCE_METRIC_FRONTIER',
        sources={p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
        selected_case_ids=[r['case'] for r in best], node_limit=500,
        verdict_counts=dict(Counter(r['verdict'] for r in records)),
        cases=sorted(records, key=lambda r: r['case']), metric_scouts=metrics,
        scope='32 explicitly selected changed-incidence types, favouring more dotted '
              'parameters and fewer large Lanner constraints. All are absent from the '
              'closed old catalogue and survive the two-block incidence gate. Completed '
              'label searches are exact exclusions; limits are unresolved. Every full '
              'necessary-label survivor receives a numerical Schur fit with three starts, '
              'dotted distances<=30 and2000 evaluations; numerical failure is not a theorem.')


if __name__ == '__main__':
    print('PID', os.getpid(), flush=True)
    result = run()
    (ROOT / f'results/{TAG}.json').write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], result['verdict_counts'], flush=True)
