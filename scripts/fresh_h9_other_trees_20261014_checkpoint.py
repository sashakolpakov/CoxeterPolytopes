#!/usr/bin/env python3
"""Quick provenance/result consistency check after four full independent DFS runs."""
import argparse,json
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_other_trees_20261014_checkpoint'
TRANSFERS='results/fresh_h9_other_trees_20261014_transfers.json'
REVIEW='results/fresh_h9_other_trees_20261014_transfer_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 seen=set();refs={}
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p in seen:return d
  seen.add(p)
  for q,h in d.get('sources',{}).items():
   assert sha(q)==h,q
   assert q not in refs or refs[q]==h
   refs[q]=h
   if q.endswith('.json'):read(q)
  return d
 transfers,review=read(TRANSFERS),read(REVIEW)
 assert 'PASS' in review['status']
 rows=[];paths=[TRANSFERS,REVIEW,f'scripts/{TAG}.py']
 for case,expected in[(17,48022),(22,47372),(25,53274),(26,62497)]:
  p=f'results/fresh_h9_other_trees_20261014_search_case_{case:02d}.json'
  q=f'results/fresh_h9_other_trees_20261014_independent_dfs_case_{case:02d}.json'
  primary,independent=read(p),read(q);paths +=[p,q]
  assert primary['status']=='EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
  assert independent['status']=='PASS_INDEPENDENT_COMPLETE_OTHER_TREE_DFS_EXCLUSION'
  assert primary['solution']is None and primary['unexplored_frontier']is None
  assert independent['survivor']is None and independent['unexplored_frontier']is None
  assert independent['cutoff']is None
  assert primary['nodes']==independent['nodes']==expected<primary['settings']['node_limit']
  assert primary['prunes']['spherical_vertex']==independent['prunes']['nonpositive_vertex']
  assert primary['prunes']['orthogonal_negative_blocks']==independent['prunes']['two_negative_blocks']
  for i,family in enumerate(primary['proved_family_conditions']):
   assert primary['prunes'].get('certified_'+family['id']+'_partial',0)==independent['prunes'].get(f'theorem_{i}',0)
  rows.append(dict(case=case,nodes=expected,vertices=independent['vertex_count'],
   missing_triples=independent['triangle_count'],bridge_clauses=independent['bridge_clauses'],
   primary=p,independent=q,complete=True,finite_masks=independent['finite_masks_checked']))
 assert sum(r['nodes']for r in rows)==211165
 return dict(status='PASS_FOUR_OTHER_TREE_EXCLUSIONS_AND_INDEPENDENT_TRAVERSALS',
  sources={p:sha(p)for p in paths},checked_json_manifests=len(seen),checked_unique_sources=len(refs),
  cases=rows,total_nodes=211165,
  scope='Four exact compact-Coxeter incidence exclusions. Every primary traversal is independently repeated with a separate exact forest-Schur engine and no cutoff. This fast checker verifies provenance/count consistency only; it does not rerun the full DFS. Other catalogue cases, including separate case5 work, are outside this certificate.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
 d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['total_nodes'],d['checked_json_manifests'],d['checked_unique_sources'])
