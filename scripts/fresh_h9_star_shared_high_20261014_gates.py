#!/usr/bin/env python3
"""Bind the two literal shared-own-root orientations for continuation."""
import argparse,json
from pathlib import Path
from hashlib import sha256
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_shared_high_20261014_gates'
OLD='results/fresh_h9_star_lastblock_20261014_gates.json';PROOF='results/fresh_h9_star_shared_high_20261014.json'
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
 old=read(OLD);new=read(PROOF);assert new['status']=='PASS_SHARED_SINGLE_OWN_ROOT_ALL_HIGH_LABELS_EXCLUDED'
 fs=old['families'][:]
 for swap in(False,True):
  perm=new['literal_twin_swap']if swap else list(range(18))
  def mapped(rows):return sorted(sorted((perm[a],perm[b]))+[m]for a,b,m in rows)
  fs.append(dict(id='shared_single_high'+('_swapped'if swap else''),origin='continuous_shared_root_high_gate',fixed=mapped(new['fixed_own_contacts']),bounds=mapped(new['internal_lower_bounds'])))
 paths=[OLD,PROOF,f'scripts/{TAG}.py']
 return dict(status='PASS_STAR_ZERO_WITH_SHARED_HIGH_FAMILIES',case=0,sources={p:sha(p)for p in paths},families=fs,scope='Adds two literal same-root orientations, with the whole finite internal tail>=7. No whole-incidence claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
