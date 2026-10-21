#!/usr/bin/env python3
"""Fresh complete height100 wrapper, preserving frozen height20 artifacts."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta55_pilot_20261014 as P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_beta55_vinberg_h100_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 previous=json.loads((ROOT/'results/fresh_h9_beta55_pilot_20261014.json').read_text())
 for f,h in previous['sources'].items():assert sha(f)==h,f
 d=P.run(F(100));n=len(previous['levels'])
 assert d['levels'][:n]==previous['levels']
 assert d['roots'][:len(previous['roots'])]==previous['roots']
 d['sources'].update({p:sha(p)for p in [f'scripts/{TAG}.py','results/fresh_h9_beta55_pilot_20261014.json']})
 d['previous_complete_height']='20';d['previous_levels_preserved']=n
 return d
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text());print('PASS saved beta55 height100 replay')
 else:out.write_text(json.dumps(d,indent=2)+'\n');print('PASS complete beta55 height100',d['counts'])
