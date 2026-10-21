#!/usr/bin/env python3
"""Search only the still untouched types outside the proved block gate."""
from hashlib import sha256
import json
import os
from pathlib import Path
import fresh_h9_bugaenko_remaining_20261014 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_bugaenko_remaining_20261014_selected'

if __name__ == '__main__':
    files = ['results/fresh_h9_bugaenko_remaining_20261014_block_gate.json',
             'results/fresh_h9_bugaenko_remaining_20261014_650_1150.json',
             'results/fresh_h9_bugaenko_remaining_20261014_1150_1950.json',
             'results/fresh_h9_bugaenko_closure_20261014_checkpoint.json']
    data = [json.loads((ROOT / p).read_text()) for p in files]
    for d in data:
        for p, h in d['sources'].items():
            assert sha256((ROOT / p).read_bytes()).hexdigest() == h, p
    assert data[0]['status'] == 'PASS_EXACT_TWO_LANNER_BLOCK_OBSTRUCTION'
    gate = {r['case'] for r in data[0]['applications']}
    ids = [i for i in range(1950, 2953) if i not in gate]
    assert len(ids) == 86
    print('PID', os.getpid(), 'selected', len(ids), flush=True)
    output = ROOT / f'results/{TAG}.json'
    result = P.run(ids, 500, 2, output)
    files.append(str(Path(__file__).resolve().relative_to(ROOT)))
    result['sources'].update({p: sha256((ROOT / p).read_bytes()).hexdigest() for p in files})
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(result['verdict_counts'], flush=True)
