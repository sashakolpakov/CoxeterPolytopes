#!/usr/bin/env python3
"""Source-bound wrapper for resumable searches on newly bound remaining T cases."""
import argparse
import json
import fresh_h9_tnext_20261015_search as engine

engine.TAG = 'fresh_h9_tremaining_20261015_search'
engine.GATE = 'results/fresh_h9_tremaining_20261015_gates_case_{case:02d}.json'

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',type=int,required=True)
    p.add_argument('--batch',type=int,required=True);p.add_argument('--budget',type=int,default=100000)
    p.add_argument('--resume');p.add_argument('--verify',action='store_true');a=p.parse_args()
    if a.verify:
        path=engine.ROOT/f'results/{engine.TAG}_case_{a.case:02d}_batch_{a.batch:03d}.json'
        saved=json.loads(path.read_text());a.budget=saved['settings']['node_budget'];a.resume=saved['settings']['resume_input']
    d=engine.run(a.case,a.batch,a.budget,a.resume,a.verify)
    if a.verify:
        assert d==saved
    print(d['status'],d['new_nodes'],len(d['pending_frames']))
