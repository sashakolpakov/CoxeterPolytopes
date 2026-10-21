#!/usr/bin/env python3
"""Bounded numerical finite-graph screen of exact-corner endpoint proposals."""
import argparse,gzip,json
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
import fresh_h9_joint_corner_20261015_rescore as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_joint_corner_20261015_graph_screen'
BASE='results/fresh_h9_joint_pair_20261015_two_extension.json';SCOUT='results/fresh_h9_joint_corner_20261015_scout.json'
def run():
 K.mp.mp.dps=100;d=json.loads((ROOT/BASE).read_text());meta=json.loads((ROOT/SCOUT).read_text());bp=meta['candidate_bank']['path'];blob=(ROOT/bp).read_bytes();assert sha256(blob).hexdigest()==meta['candidate_bank']['sha256'];bank=json.loads(gzip.decompress(blob))
 scale=np.array([np.sqrt(2/float(K.scalar(x)))for x in d['norms']]);G=np.array([[float(K.scalar(x))for x in row]for row in d['gram']])*scale[:,None]*scale[None,:]
 old=[set(v['incident'])for v in d['finite_vertices']];oldedges=[e['vertices']for e in d['finite_edges']];joint=next(i for i,x in enumerate(old)if{30,31}<=x)
 def graph(supports):
  parent=list(range(len(supports)))
  def find(i):
   while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
   return i
  def union(i,j):parent[find(i)]=find(j)
  for i,j in oldedges:union(i,j)
  for j in range(len(old),len(supports)):
   for i in range(j):
    if len(supports[i]&supports[j])==8:union(i,j)
  return [find(i)for i in range(len(supports))]
 base=graph(old);rows=[];tested=0
 for k,c in enumerate(bank):
  new=[];v=np.array(c['normalized_wall_pairings'])
  for S in combinations(c['finite_neighbors'],8):
   tested+=1;M=np.zeros((9,9));M[:8,:8]=G[np.ix_(S,S)];M[:8,8]=M[8,:8]=v[list(S)];M[8,8]=2
   if np.linalg.eigvalsh(M)[0]>1e-7:new.append(set(S)|{32})
  assert new
  final=graph(old+new);merged=[i for i in range(len(old))if final[i]==final[joint]and base[i]!=base[joint]]
  rows.append(dict(bank_index=k,new_vertices=[sorted(x)for x in new],merged_old_components=len(set(base))-len(set(final)),
     joint_corner_newly_connected_old_vertices=merged,new_edges=sum(len(s&t)==8 for j,s in enumerate(new)for t in old+new[:j]),escape_pairing=c['escape_pairing']))
 ranked=sorted(rows,key=lambda r:(-len(r['joint_corner_newly_connected_old_vertices']),-r['merged_old_components'],-len(r['new_vertices']),-r['new_edges']))
 paths=[BASE,SCOUT,bp,f'scripts/{TAG}.py','scripts/fresh_h9_joint_corner_20261015_rescore.py']
 return dict(status='BOUNDED_NUMERICAL_FINITE_EDGE_COMPONENT_SCREEN',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
    bank_count=len(bank),candidate_supports_tested=tested,old_component_sizes=sorted(Counter(base).values(),reverse=True),joint_corner_index=joint,
    merge_histogram={str(k):v for k,v in Counter(r['merged_old_components']for r in rows).items()},
    new_vertex_histogram={str(k):v for k,v in Counter(len(r['new_vertices'])for r in rows).items()},best_candidates=ranked[:24],
    scope='Numericalpositive Gram tests with margin1e-7, on bounded savedendpointbank. Counts may include duplicategeometric proposals. Only selected bridge is separately exactified; this screen excludes no arbitrarycompletion and proves no newwall.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['merge_histogram'])
