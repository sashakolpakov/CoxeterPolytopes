#!/usr/bin/env python3
"""Fresh higher-height continuation, preserving the frozen height100 prefix."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta11_vinberg_20261014 as V
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_beta11_vinberg_h300_20261014'
BASE='results/fresh_h9_beta11_vinberg_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    old=json.loads((ROOT/BASE).read_text())
    for p,h in old['sources'].items():assert sha(p)==h,p
    result=V.run(F(300))
    count=old['counts']['positive_levels']
    assert result['levels'][:count]==old['levels']
    assert result['roots'][:old['counts']['total_walls']]==old['roots']
    result['sources'].update({BASE:sha(BASE),f'scripts/{TAG}.py':sha(f'scripts/{TAG}.py')})
    result['previous_complete_height']=100
    result['previous_levels_preserved']=count
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
