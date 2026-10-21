#!/usr/bin/env python3
"""Require three independent reviews before adding local-core predicates."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_core_transfers_20261014_gates'
OLD='results/fresh_h9_star_mixed_singles_20261014_gates.json'
NEW='results/fresh_h9_star_core_transfers_20261014.json'
REVIEW='results/fresh_h9_star_core_transfers_20261014_review.json'
DEPENDENCIES=[('results/fresh_h9_star_mixed_singles_20261014.json','results/fresh_h9_star_mixed_singles_20261014_review.json'),
 ('results/fresh_h9_star_opposite_core_20261014.json','results/fresh_h9_star_opposite_core_20261014_review.json')]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 seen={}
 def read(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for q,h in d.get('sources',{}).items():
   assert sha(q)==h,q
   if q.endswith('.json'):read(q)
  return d
 old,new,review=[read(p)for p in (OLD,NEW,REVIEW)]
 assert old['status']=='PASS_SHARED_SINGLE_CORE_AND_REVERSAL_COMPOSITION'
 assert new['status']=='PASS_TWELVE_LITERAL_CORE_TRANSFERS'
 assert review['status'].startswith('PASS_')and review['sources'][NEW]==sha(NEW)
 for primary,independent in DEPENDENCIES:
  p,r=read(primary),read(independent)
  assert p['status'].startswith('PASS_')and r['status'].startswith('PASS_')
  assert r['sources'][primary]==sha(primary)and new['sources'][primary]==sha(primary)
 assert old['case']==new['case']==0
 def key(f):return tuple(tuple(sorted(tuple(row)for row in f[k]))for k in ('fixed','bounds'))
 families=old['families'][:];seen_keys={key(f)for f in families};added=[];duplicate=[]
 assert len(families)==209 and len(new['families'])==48
 for f in new['families']:
  if key(f)in seen_keys:duplicate.append(f['id']);continue
  seen_keys.add(key(f));families.append(f);added.append(f['id'])
 assert len(families)==253 and len(added)==44 and len(duplicate)==4
 paths=[OLD,NEW,REVIEW,*[x for pair in DEPENDENCIES for x in pair],f'scripts/{TAG}.py']
 return dict(status='PASS_REVIEWED_LOCAL_CORE_GATE_COMPOSITION',case=0,sources={p:sha(p)for p in paths},
  families=families,inherited_predicates=209,new_predicates=44,duplicate_images=duplicate,added_ids=added,
  scope='253 reviewed predicates for star0, including48 literal local-core transfers (four already present). Both arithmetic theorem reviews and the independent transfer review are mandatory dependencies. No whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
