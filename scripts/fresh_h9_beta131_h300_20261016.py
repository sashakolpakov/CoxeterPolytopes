#!/usr/bin/env python3
"""Complete height-300 continuation of the frozen norm-131 prefix."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta131_pilot_20261016 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_beta131_h300_20261016'
OLD='results/fresh_h9_beta131_h100_20261016.json'

def sha(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()

def run():
    old=json.loads((ROOT/OLD).read_text())
    for path,digest in old['sources'].items():
        assert sha(path)==digest,path
    data=P.run(F(300))
    count=len(old['levels'])
    assert data['levels'][:count]==old['levels']
    assert data['roots'][:len(old['roots'])]==old['roots']
    data['sources'].update({path:sha(path) for path in [f'scripts/{TAG}.py',OLD]})
    data['preserved_prefix']=dict(height=100,levels=count,walls=len(old['roots']))
    return data

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(path.read_text())
        print('PASS saved norm-131 height-300 replay')
    else:
        path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
