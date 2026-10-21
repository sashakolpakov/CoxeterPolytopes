#!/usr/bin/env python3
"""Bind both opposite-single internal regimes and their literal twin images."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_zero_20261014_opposite_gates'
OLD='results/fresh_h9_star_zero_20261014_family_gates.json'
PROOFS=['results/fresh_h9_star_zero_20261014_'+s+'.json'for s in('opposite_single','opposite_three')]
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
 old=read(OLD);large,three=[read(p)for p in PROOFS]
 assert large['status']=='PASS_ALL_LABEL_OPPOSITE_SINGLE_FAMILY_EXCLUDED'
 assert three['status']=='PASS_OPPOSITE_SINGLE_INTERNAL_THREE_NEGATIVE_INDEX'
 families=old['families'][:]
 for name,proof,fixed,bounds in [('internal_ge4',large,large['fixed_other_own_contacts'],large['nonright_own_contacts']+large['internal_label_lower_bounds']),('internal3',three,three['fixed'],three['own_bounds'])]:
  for swap in(False,True):
   perm=proof['literal_twin_swap']if swap else list(range(18))
   def mapped(rows):return sorted(sorted((perm[a],perm[b]))+[m]for a,b,m in rows)
   families.append(dict(id='opposite_single_'+name+('_swapped'if swap else''),origin='reviewed_opposite_single',fixed=mapped(fixed),bounds=mapped(bounds)))
 paths=[OLD]+PROOFS+[f'scripts/{TAG}.py']
 return dict(status='PASS_STAR_ZERO_WITH_OPPOSITE_SINGLE_FAMILIES',case=0,sources={p:sha(p)for p in paths},families=families,
  scope='The prior150 conditions plus both literal orientations of two separately certified opposite-single regimes. No unrestricted whole-incidence claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
