#!/usr/bin/env python3
"""Exhaustive small overlapping-Lanner relations used by the new frontier."""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h9_highbase_construction_20261014_cell import codes,finite
from fresh_h9_newfacet_frontier_20261014_patterns import bits,maximal_faces,local_search,inventory,INC,FRONT
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_newfacet_frontier_20261014_relations'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def complete_six_relation():
 n=6;missing=(61,62);pairs=list(combinations(range(n),2));pi={p:i for i,p in enumerate(pairs)}
 vv=[(f.bit_count(),tuple(pi[p]for p in combinations(bits(f),2)))for f in maximal_faces(n,missing)]
 tt=[(tuple(pi[p]for p in combinations(bits(m),2)),codes(m.bit_count()))for m in missing]
 rows=[];nodes=0
 def rec(left,labels,known):
  nonlocal nodes
  nodes+=1
  if any(not finite(k,tuple(labels[p]for p in pp))for k,pp in vv):return
  if left:
   pp,pool=tt[left-1]
   for code in pool:
    if any(p in known and labels[p]!=m for p,m in zip(pp,code)):continue
    child=labels[:]
    for p,m in zip(pp,code):child[p]=m
    rec(left-1,child,known|set(pp))
  else:
   unknown=set(range(len(pairs)))-known
   if not unknown:rows.append(labels[:]);return
   p=min(unknown)
   for m in range(2,8):
    child=labels[:];child[p]=m;rec(0,child,known|{p})
 rec(2,[2]*15,set());assert len(rows)==len({tuple(c)for c in rows})==96
 return dict(order=6,minimal_nonfaces=missing,maximal_faces=maximal_faces(n,missing),pair_order=pairs,nodes=nodes,complete_labels=sorted(rows))
def run():
 data=json.loads((ROOT/INC).read_text());front=json.loads((ROOT/FRONT).read_text())
 for d in(data,front):
  for p,h in d['sources'].items():assert sha(p)==h
 bank=inventory(data,front['changed_order']['open_case_ids']);key=(7,(31,99));r=local_search(*key,100000)
 assert r['verdict']=='EXACT_LOCAL_OBSTRUCTION'and r['nodes']==16237
 applications=[]
 for i,perm in sorted(bank[key].items()):
  row=data['cases'][i];union=sum(1<<j for j in perm)
  assert {sum(1<<perm[j]for j in bits(m))for m in key[1]}=={m for m in row['canonical_minimal_nonfaces']if m&union==m}
  for face in r['maximal_faces']:
   lift=sum(1<<perm[j]for j in bits(face));assert any(v&lift==lift for v in row['canonical_vertices'])
  applications.append([i,perm])
 primary='results/fresh_h9_newfacet_frontier_20261014_patterns.json';p=json.loads((ROOT/primary).read_text())
 for f,h in p['sources'].items():assert sha(f)==h
 combined=set(i for i,k in p['excluded_cases'])|{i for i,_ in applications};remaining=sorted(set(front['changed_order']['open_case_ids'])-combined)
 paths=[INC,FRONT,primary,f'scripts/{TAG}.py','scripts/fresh_h9_newfacet_frontier_20261014_patterns.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='PASS_COMPLETE_OVERLAP_RELATION_AND_SEVEN_NODE_GATE',sources={f:sha(f)for f in paths},six_node_relation=complete_six_relation(),seven_node_obstruction=dict(order=7,minimal_nonfaces=key[1],**r,applications=applications),combined_excluded_cases=sorted(combined),remaining_case_ids=remaining,scope='Two order5 Lanner blocks sharing4 nodes admit exactly96 finite-category completions on their6-node union under the actual-face constraints. This relation is necessary, not a metric realization. An order5/order4 pair sharing2 nodes, with no other induced missing face, admits no finite-category completion; 16237 states exhaust it. Both conclusions cover all integer finite labels via category7, with no dotted parameters in these induced configurations.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],len(out['combined_excluded_cases']),len(out['remaining_case_ids']))
