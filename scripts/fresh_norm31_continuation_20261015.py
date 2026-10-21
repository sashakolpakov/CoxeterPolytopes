#!/usr/bin/env python3
"""Complete height100 continuation, preserving the frozen norm31 pilot."""
import argparse,json
from pathlib import Path
from hashlib import sha256
import fresh_frontier_construction_20261015 as A
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_norm31_continuation_20261015'
PILOT='results/fresh_frontier_construction_20261015.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old=json.loads((ROOT/PILOT).read_text())
    for f,h in old['sources'].items():assert sha(f)==h,f
    d=A.run(100)
    assert d['levels'][:len(old['levels'])]==old['levels']
    assert d['roots'][:len(old['roots'])]==old['roots']
    d['sources'].update({f'scripts/{TAG}.py':sha(f'scripts/{TAG}.py'),PILOT:sha(PILOT)})
    d['status']='PASS_COMPLETE_DISCRIMINANT31_HEIGHT100_PREFIX'
    d['pilot_levels_reproduced']=len(old['levels'])
    return d
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS height100 saved replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['counts'])
