#!/usr/bin/env python3
"""Append the complete last-internal4/5 family to frozen star0 gates."""
import argparse,json
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_lastblock_20261014_gates'
OLD='results/fresh_h9_star_zero_20261014_opposite_gates.json';PROOF='results/fresh_h9_star_lastblock_20261014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p not in seen:
   seen.add(p)
   for q,h in d.get('sources',{}).items():
    assert sha(q)==h,q
    if q.endswith('.json'):read(q)
  return d
 old=read(OLD);new=read(PROOF);assert new['status']=='PASS_COMPLETE_LAST_BLOCK_4_5_OWN_PROFILE_EXCLUSION'
 fs=old['families'][:]
 for m in new['last_internal_values']:
  fs.append(dict(id=f'last_block_{m}_arbitrary_own',origin='complete_last_block_profile_gate',fixed=new['fixed_first_three_own']+[[12,13,m]],bounds=new['internal_lower_bounds']))
 paths=[OLD,PROOF,f'scripts/{TAG}.py']
 return dict(status='PASS_STAR_ZERO_WITH_LAST_BLOCK_PROFILE_FAMILIES',case=0,sources={p:sha(p)for p in paths},families=fs,scope='Adds two complete last-internal4/5 families with arbitrary own finite labels. First three blocks retain all-own3. No whole-incidence claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
