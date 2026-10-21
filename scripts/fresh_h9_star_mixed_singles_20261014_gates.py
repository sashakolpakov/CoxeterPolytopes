#!/usr/bin/env python3
"""Bind the independently reviewed shared-single theorem and literal reversal."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_mixed_singles_20261014_gates'
OLD='results/fresh_h9_star_single_right_20261014_finite_core_gates.json'
NEW='results/fresh_h9_star_mixed_singles_20261014.json'
REVIEW='results/fresh_h9_star_mixed_singles_20261014_review.json'
SYM='results/fresh_h9_star_single_right_20261014_symmetry.json'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
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
 old,new,review,sym,geometry=[read(p)for p in (OLD,NEW,REVIEW,SYM,GEOM)]
 assert old['status']=='PASS_CORE_ONLY_AND_REVERSAL_GATE_COMPOSITION'
 assert new['status']=='PASS_CORE_SHARED_SINGLE_ALL_INTERNAL_LABELS'
 assert review['status'].startswith('PASS_')and review['sources'][NEW]==sha(NEW)
 assert old['case']==new['case']==sym['case']==0
 p=sym['permutation'];assert sorted(p)==list(range(18))and all(p[p[i]]==i for i in range(18))
 g=geometry['representatives'][0]['geometry']
 vertices={tuple(sorted(v['incident']))for v in g['vertices']}
 missing={tuple(sorted(m))for m in g['minimal_nonfaces']}
 assert {tuple(sorted(p[i]for i in v))for v in vertices}==vertices
 assert {tuple(sorted(p[i]for i in m))for m in missing}==missing
 def key(f):return tuple(tuple(sorted(tuple(row)for row in f[k]))for k in ('fixed','bounds'))
 families=old['families'][:];seen_keys={key(f)for f in families};transfers=[]
 assert len(families)==205 and len(new['families'])==2
 for f in new['families']:
  for reverse in (False,True):
   image={k:sorted([sorted((p[a],p[b]))+[m]for a,b,m in f[k]])if reverse else f[k]for k in ('fixed','bounds')}
   assert key(image)not in seen_keys;seen_keys.add(key(image))
   ident=('reversal_'if reverse else'')+f['id']
   families.append(dict(id=ident,origin='reviewed_core_shared_single_all_internal_labels',**image))
   transfers.append(dict(source=f['id'],reverse=reverse,target=ident,**image))
 assert len(families)==209
 paths=[OLD,NEW,REVIEW,SYM,GEOM,f'scripts/{TAG}.py']
 return dict(status='PASS_SHARED_SINGLE_CORE_AND_REVERSAL_COMPOSITION',case=0,
  sources={p:sha(p)for p in paths},families=families,transfers=transfers,
  inherited_predicates=205,new_predicates=4,
  scope='Two reviewed shared-single families and their exact literal reversal images, added to205 previous predicates. Every new predicate retains its stated core contacts and internal lower bounds; no whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
