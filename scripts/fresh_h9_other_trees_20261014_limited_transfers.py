#!/usr/bin/env python3
"""Literal local-premise transfers to the previously limited case5.

No incidence isomorphism is assumed. Only maps for which every proof
support is an actual face and every required missing triple is literal
are retained. Algebra is inherited from source-bound independent proofs.
"""
import argparse
import json
from hashlib import sha256
from itertools import combinations, permutations
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_other_trees_20261014_limited_transfers'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
HIGH='results/fresh_h9_t8_sparse_search_20261014_two_high_blocks.json'
HIGH_REVIEW='results/fresh_h9_t8_sparse_search_20261014_two_high_blocks_review.json'
ALL='results/fresh_h9_t8_sparse_search_20261014_all_blocks.json'
ALL_REVIEW='results/fresh_h9_t8_all_blocks_20261014_review.json'
SPARSE='results/fresh_h9_t8_sparse_frontier_20261014.json'
SPARSE_REVIEW='results/fresh_h9_t8_sparse_frontier_20261014_review.json'
DOMAINS='results/fresh_h9_t8_asymmetric_20261014.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def run():
 seen=set()
 def read(p):
  d=json.loads((ROOT/p).read_text())
  if p in seen:return d
  seen.add(p)
  for q,h in d.get('sources',{}).items():
   assert sha(q)==h,q
   if q.endswith('.json'):read(q)
  return d
 geom,high,allblocks,sparse,domains=[read(p)for p in[GEOM,HIGH,ALL,SPARSE,DOMAINS]]
 for p in[HIGH_REVIEW,ALL_REVIEW,SPARSE_REVIEW]:assert 'PASS' in read(p)['status']
 high_supports=[('positive_support',high['positive_support'])]
 high_supports += [('meeting_pair',r['pair'])for r in high['all45_core_meeting_pairs']]
 for k in['cross_right_guards','anchor_high_right_budgets','third_other_blocks_right_budgets','third_anchor_upper_stars']:
  high_supports += [(k,r['support'])for r in high[k]]
 sparse_supports=[('positive_support',domains['positive_support'])]
 sparse_supports += [('meeting_pair',list(p))for p in combinations(sparse['core_roots'],2)]
 for k in['cross_right_guards','anchor_budgets','third_budgets']:
  sparse_supports += [(k,r['support'])for r in domains[k]]
 for pattern in domains['own_support_patterns']:
  sparse_supports.append(('anchor_pair_guard',pattern['anchor_pair_guard'].get('support',[11,17])))
 sparse_supports += [('third_upper_guard',r['support'])for r in sparse['replacement_third_contact_guards']]
 cases=[]
 for cid in[5]:
  g=geom['representatives'][cid]['geometry'];vertices=[set(v['incident'])for v in g['vertices']]
  missing={frozenset(m)for m in g['minimal_nonfaces']}
  blocks=[b for b in g['blocks']if len(b)==2];assert len(blocks)==4
  outside=sorted(set(range(18))-set(sum(blocks,[])))
  def actual(s):return next((i for i,v in enumerate(vertices)if set(s)<=v),None)
  def transfer(kind,mp,supports,triples,proof):
   checks=[]
   for key,s in supports:
    target=[mp[x]for x in s];vi=actual(target)
    if vi is None:return None
    checks.append(dict(kind=key,source_support=s,target_support=target,actual_vertex=vi))
   target_triples=[sorted(mp[x]for x in t)for t in triples]
   if any(frozenset(t)not in missing for t in target_triples):return None
   def translate(rows):return sorted([sorted((mp[a],mp[b]))+[m]for a,b,m in rows])
   return dict(kind=kind,mapping=sorted([[a,b]for a,b in mp.items()]),
    fixed_labels=translate(proof['fixed_labels']),finite_label_lower_bounds=translate(proof['finite_label_lower_bounds']),
    actual_support_checks=checks,required_missing_triples=target_triples)
  rows=[]
  for hi in combinations(range(4),2):
   low=[i for i in range(4)if i not in hi]
   for own in low:
    order=[*hi,next(i for i in low if i!=own),own]
    for bridge in outside:
     if not all(frozenset(blocks[i]+[bridge])in missing for i in hi):continue
     for anchor in outside:
      if anchor==bridge or frozenset(blocks[own]+[anchor])not in missing:continue
      mp=dict(zip(high['positive_support'],sum([blocks[i]for i in order],[])))
      mp.update({14:bridge,11:anchor})
      row=transfer('two_high_negative_plane',mp,high_supports,high['missing_triples'],high)
      if row:rows.append(row)
  for order in permutations(range(4)):
   a,b,c,d=[blocks[i]for i in order]
   for bridge in outside:
    if not all(frozenset(block+[bridge])in missing for block in[a,b]):continue
    anchors=[x for x in outside if x!=bridge and frozenset(d+[x])in missing]
    for x,y in combinations(anchors,2):
     for reverse in[False,True]:
      mp=dict(zip(domains['positive_support'],(a[::-1]if reverse else a)+b+c+d))
      mp.update({11:x,14:bridge,17:y})
      row=transfer('sparse_all_blocks',mp,sparse_supports,[[11,12,13],[12,13,17]],allblocks)
      if row:rows.append(row)
  families=[];keys=set()
  for j,row in enumerate(rows):
   key=(tuple(map(tuple,row['fixed_labels'])),tuple(map(tuple,row['finite_label_lower_bounds'])))
   if key in keys:continue
   keys.add(key);families.append(dict(witness_map=j,kind=row['kind'],fixed=row['fixed_labels'],bounds=row['finite_label_lower_bounds']))
  saved_path=f'results/fresh_h9_tree_cube_wedges_20261014_case_{cid:02d}.json';saved=read(saved_path)
  lab={(a,b):m for a,b,m in saved['solution']['pair_labels']} if saved['solution'] else {}
  hits=[j for j,f in enumerate(families)if lab and all(lab[a,b]==m for a,b,m in f['fixed'])and all(lab[a,b]>=m for a,b,m in f['bounds'])]
  assert saved['status']=='NODE_LIMIT_UNRESOLVED'
  cases.append(dict(case=cid,vertices=len(vertices),blocks=blocks,maps=rows,families=families,
    saved_survivor_source=saved_path,saved_survivor_excluded_by=hits,
    map_counts={kind:sum(r['kind']==kind for r in rows)for kind in['two_high_negative_plane','sparse_all_blocks']}))
 paths=[GEOM,HIGH,HIGH_REVIEW,ALL,ALL_REVIEW,SPARSE,SPARSE_REVIEW,DOMAINS,f'scripts/{TAG}.py']
 paths +=[r['saved_survivor_source']for r in cases]
 return dict(status='PASS_LITERAL_REVIEWED_THEOREM_TRANSFERS_TO_LIMITED_CASE5',
  sources={p:sha(p)for p in paths},cases=cases,
  scope='Sufficient literal local-premise maps, not an isomorphism assertion or exhaustive map catalogue. Every selected proof face and missing triple is checked anew. Analytic arithmetic is inherited from separate frozen reviews. The old case5 search was node-limited; no old survivor exclusion is asserted. The incidence remains open until a separate continuation.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 d=run();out=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],[(r['case'],r['map_counts'],len(r['families']))for r in d['cases']])
