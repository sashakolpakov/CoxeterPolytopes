#!/usr/bin/env python3
"""Fresh temporary-database replay of the bounded rank9 pilot."""
import argparse,json,tempfile
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_triangle_extension_frontier_20261016 as P
ROOT=P.ROOT;TAG='fresh_sqrt5_triangle_extension_frontier_20261016_replay'
def run():
    saved=ROOT/'results/fresh_sqrt5_triangle_extension_frontier_20261016_rank9.json'
    with tempfile.TemporaryDirectory(prefix='sqrt5_triangle_rank9_replay_')as d:
        out=Path(d)/'replay.json';P.run(Path(d)/'replay.sqlite',out,9,100000)
        result=json.loads(out.read_text());assert result==json.loads(saved.read_text())
    src=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_triangle_extension_frontier_20261016.py','results/fresh_sqrt5_triangle_extension_frontier_20261016_rank9.json']
    return dict(status='PASS_FRESH_BOUNDED_RANK9_DATABASE_REPLAY',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in src},connected_states=result['connected_states'],counts=result['counts'],scope='Same producer, fresh temporary database and exact full-export equality. This is a reproducibility replay, not an independent mathematical review or full rank22 coverage.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
