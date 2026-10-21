#!/usr/bin/env python3
"""Resume the saved T10 cursor only after binding its independent metric review."""
import argparse
import json
import fresh_h9_tnext_20261015_search as engine

engine.TAG='fresh_h9_tremaining_20261015_case10_metric_search'
engine.GATE='results/fresh_h9_tremaining_20261015_case10_metric_gates.json'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--batch',type=int,default=2)
    p.add_argument('--budget',type=int,default=100000)
    p.add_argument('--resume',default='results/fresh_h9_tremaining_20261015_search_case_10_batch_001.json')
    p.add_argument('--verify',action='store_true');a=p.parse_args()
    if a.verify:
        path=engine.ROOT/f'results/{engine.TAG}_case_10_batch_{a.batch:03d}.json'
        saved=json.loads(path.read_text());a.budget=saved['settings']['node_budget'];a.resume=saved['settings']['resume_input']
    d=engine.run(10,a.batch,a.budget,a.resume,a.verify)
    if a.verify:assert d==saved
    print(d['status'],d['new_nodes'],len(d['pending_frames']))
