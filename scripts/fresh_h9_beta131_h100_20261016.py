#!/usr/bin/env python3
"""Complete height100 continuation, preserving the frozen norm131 pilot."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta131_pilot_20261016 as P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_beta131_h100_20261016'
OLD='results/fresh_h9_beta131_pilot_20261016.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old=json.loads((ROOT/OLD).read_text())
    for p,h in old['sources'].items():assert sha(p)==h,p
    d=P.run(F(100));n=len(old['levels'])
    assert d['levels'][:n]==old['levels'] and d['roots'][:len(old['roots'])]==old['roots']
    d['sources'].update({p:sha(p)for p in[f'scripts/{TAG}.py',OLD]})
    d['preserved_prefix']=dict(height=20,levels=n,walls=len(old['roots']))
    return d
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved norm131 height100 replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
