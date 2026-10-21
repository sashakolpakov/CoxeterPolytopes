#!/usr/bin/env python3
"""Thin, source-bound wrapper resuming T31's saved cursor with155 proved predicates."""
import argparse
import json
import fresh_h9_tnext_20261015_search as engine

engine.TAG = 'fresh_h9_tnext_20261015_resume'
engine.GATE = 'results/fresh_h9_tnext_20261015_case31_gates.json'
PRIOR = 'results/fresh_h9_tnext_20261015_search_case_31_batch_001.json'

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--batch', type=int, default=2)
    p.add_argument('--budget', type=int, default=100000)
    p.add_argument('--resume', default=PRIOR)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    if a.verify:
        path = engine.ROOT / f'results/{engine.TAG}_case_31_batch_{a.batch:03d}.json'
        saved = json.loads(path.read_text())
        a.budget = saved['settings']['node_budget']
        a.resume = saved['settings']['resume_input']
    result = engine.run(31, a.batch, a.budget, a.resume, a.verify)
    if a.verify:
        assert result == saved
    print(result['status'], result['new_nodes'], len(result['pending_frames']))
