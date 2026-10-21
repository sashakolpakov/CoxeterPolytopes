#!/usr/bin/env python3
"""Bind the exact all-own-three family to the unchanged prior149 conditions."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_zero_20261014_family_gates'
OLD='results/fresh_h9_star_zero_20261014_gates.json';PROOF='results/fresh_h9_star_zero_20261014_own_three.json'
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
 old,proof=read(OLD),read(PROOF)
 assert old['status']=='PASS_LITERAL_STAR_ZERO_CONDITIONS_BOUND'
 assert proof['status']=='PASS_ALL_OWN_THREE_ARBITRARY_BLOCK_LABELS_EXCLUDED' and not proof['unresolved']
 assert proof['residual_cases']==21952 and len(proof['arbitrary_high_block_reduction'])==4
 f=dict(id='all_own_three_arbitrary_internal',origin='exact_rank_family',fixed=proof['all_own_contacts'],bounds=[])
 assert len(f['fixed'])==16 and all(m==3 for a,b,m in f['fixed'])
 paths=[OLD,PROOF,f'scripts/{TAG}.py']
 return dict(status='PASS_STAR_ZERO_WITH_ALL_OWN_THREE_FAMILY',case=0,
  sources={p:sha(p)for p in paths},families=old['families']+[f],
  scope='Inherits149 conditional exclusions and adds the all-eight-leaves own3/3 theorem for arbitrary finite internal labels. No other label restriction imposed.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
