#!/usr/bin/env python3
"""Compose reviewed finite-core exclusions with the literal star reversal."""
import argparse,json
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_single_right_20261014_finite_core_gates'
OLD='results/fresh_h9_star_single_right_20261014_symmetry.json'
NEW='results/fresh_h9_star_single_right_20261014_finite_core.json'
REVIEW='results/fresh_h9_star_single_right_20261014_finite_core_review.json'
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
 assert old['status']=='PASS_LITERAL_STAR_REVERSAL_GATE_TRANSFER'
 assert new['status']=='PASS_CORE_ONLY_FINITE_INTERNAL_LABEL_FAMILY'
 assert review['status'].startswith('PASS_')and review['sources'][NEW]==sha(NEW)
 assert old['case']==new['case']==0
 p=old['permutation'];assert sorted(p)==list(range(18))and all(p[p[i]]==i for i in range(18))
 def key(f):return tuple(tuple(sorted(tuple(row)for row in f[k]))for k in ('fixed','bounds'))
 families=old['families'][:];seen_keys={key(f)for f in families};transfers=[]
 assert len(families)==173 and len(new['families'])==16
 for f in new['families']:
  for reverse in (False,True):
   image={k:sorted([sorted((p[a],p[b]))+[m]for a,b,m in f[k]])if reverse else f[k]for k in ('fixed','bounds')}
   assert key(image)not in seen_keys;seen_keys.add(key(image))
   ident=('reversal_'if reverse else'')+f['id']
   families.append(dict(id=ident,origin='reviewed_core_only_finite_block_rank',**image))
   transfers.append(dict(source=f['id'],reverse=reverse,target=ident,**image))
 assert len(families)==205
 paths=[OLD,NEW,REVIEW,f'scripts/{TAG}.py']
 return dict(status='PASS_CORE_ONLY_AND_REVERSAL_GATE_COMPOSITION',case=0,
  sources={p:sha(p)for p in paths},families=families,transfers=transfers,
  inherited_predicates=173,new_predicates=32,
  scope='The16 reviewed finite-core label families and their exact literal reversal images, added to173 previous predicates. No additional arithmetic or whole-incidence classification.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['families']))
