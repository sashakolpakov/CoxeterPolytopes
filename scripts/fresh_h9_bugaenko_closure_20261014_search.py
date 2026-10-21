#!/usr/bin/env python3
"""Early necessary spherical/signature/rank gates during Lanner assignment."""
import argparse,json,os,time
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import pynauty
from fresh_h9_bugaenko_closure_20261014 import Gates,ROOT,OLD,PAIRS,PI,sha
from fresh_h9_highbase_construction_20261014_cell import finite,codes
TAG='fresh_h9_bugaenko_closure_20261014_search'
def search(row,limit):
 g=Gates(row);constraints=[(tuple(PI[p]for p in combinations(a,2)),codes(len(a)))for a in g.M if len(a)>=3]
 affected=[tuple(i for i,v in enumerate(g.vp)if p in v)for p in range(91)]
 nodes=visits=0;stats=Counter();proofs=[];candidate=None;dead=set()
 twins=defaultdict(set)
 for i in range(14):twins[tuple(j for j,m in enumerate(g.M)if i in m)].add(i)
 twin_blocks=list(twins.values());adj={i:[]for i in range(105)}
 for k,(i,j)in enumerate(PAIRS,14):adj[i].append(k);adj[j].append(k);adj[k]=[i,j]
 def canonical(left,labels,known):
  values=[-1 if p in g.dotted else(labels[p]if p in known else 0)for p in range(91)];ec=defaultdict(set)
  for p,c in enumerate(values,14):ec[c].add(p)
  graph=pynauty.Graph(105,adjacency_dict=adj,vertex_coloring=twin_blocks+[ec[c]for c in sorted(ec)])
  order=pynauty.canon_label(graph)[:14];assert set(order)==set(range(14))
  pos=0
  for b in twin_blocks:assert set(order[pos:pos+len(b)])==b;pos+=len(b)
  # The key is the explicitly relabelled state, not only a graph digest.
  return left,tuple(values[PI[tuple(sorted((order[i],order[j])))]]for i,j in PAIRS)
 def force_right(labels,known):
  out=set(known)
  for p in range(91):
   if p in known or p in g.dotted:continue
   labels[p]=3
   if any(not finite(9,tuple(labels[k]for k in g.vp[i]))for i in affected[p]):out.add(p)
   labels[p]=2
  stats['forced_right_occurrences']+=len(out)-len(known)
  return out
 def walk(left,labels,known):
  nonlocal nodes,visits,candidate
  visits+=1;key=canonical(left,labels,known)
  if key in dead:stats['cache_hits']+=1;return False
  nodes+=1
  if nodes>limit:return None
  bad=g.bad_vertex(labels)
  if bad is not None:stats['NONSPHERICAL_VERTEX']+=1;dead.add(key);return False
  known=force_right(labels,known)
  pr=g.proof(labels,known)
  if pr:
   stats[pr['kind']]+=1
   if len(proofs)<100:proofs.append(dict(labels=labels[:],known=sorted(known),proof=pr))
   dead.add(key);return False
  pending=[i for i in range(len(constraints))if left>>i&1]
  if pending:
   # Shared constraints are checked against forced right entries as well.
   options={}
   for i in pending:
    pp,pool=constraints[i];keep=[code for code in pool if all(p not in known or code[j]==labels[p]for j,p in enumerate(pp))]
    if not keep:stats['NO_LANNER_EXTENSION']+=1;dead.add(key);return False
    options[i]=keep
   k=min(pending,key=lambda i:(len(options[i]),-sum(p in known for p in constraints[i][0]),-len(constraints[i][0]),i))
   pp=constraints[k][0]
   for code in options[k]:
    child=labels[:]
    for p,m in zip(pp,code):child[p]=m
    r=walk(left&~(1<<k),child,known|set(pp))
    if r is not False:return r
  else:
   unknown=set(range(91))-known-g.dotted
   if not unknown:candidate=labels[:];return True
   # Complete the nearest fixed minor rather than deferring rank checks.
   targets=[(sum(p not in known for p in pp),pp)for a,pp in g.cliques if all(labels[p]<=6 for p in pp if p in known)]
   target=min(targets,key=lambda t:t[0])[1]if targets else []
   p=min(unknown,key=lambda p:(p not in target,-len(affected[p]),p))
   for m in range(2,8):
    child=labels[:];child[p]=m;r=walk(left,child,known|{p})
    if r is not False:return r
  dead.add(key);return False
 result=walk((1<<len(constraints))-1,[2]*91,set())
 out=dict(case=row['case'],verdict='FULL_CATEGORY_METRIC_REQUIRED'if result is True else('NODE_LIMIT_UNRESOLVED'if result is None else'EXACT_COMPLETE_OBSTRUCTION'),nodes=nodes,visits=visits,twin_blocks=[sorted(b)for b in twin_blocks],counters=dict(stats),sample_proofs=proofs,distinct_modular_minor_codes=len(g.det_cache))
 if candidate is not None:out.update(label_vector=candidate,dotted_pairs=[PAIRS[p]for p in sorted(g.dotted)])
 return out
def run(start,end,limit):
 src=f'results/{OLD}_incidence.json';inc=json.loads((ROOT/src).read_text());end=min(end,len(inc['cases']));rows=[];t=time.time()
 for row in inc['cases'][start:end]:
  r=search(row,limit);rows.append(r);print(row['case'],r['verdict'],r['nodes'],'seconds',round(time.time()-t,2),flush=True)
 files=[src,f'scripts/{TAG}.py','scripts/fresh_h9_bugaenko_closure_20261014.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='EXACT_EARLY_GATE_SEARCH_CHECKPOINT',sources={p:sha(p)for p in files},range=[start,end],node_limit=limit,verdict_counts=dict(Counter(r['verdict']for r in rows)),cases=rows,scope='Necessary finite-category Lanner assignments and all actual vertex diagrams; unknown meeting pairs relaxed to2 only for pruning. Label3 failure forces a pair right by cosine-matrix monotonicity. Orthogonal disjoint minimal nonfaces force two negative directions even before their internal labels are assigned. Fixed11minors avoid all unknown/dotted/>=7 entries. Canonical caching uses explicit partial-state relabelings within twin blocks, each preserving every missing face. Limits count evaluated uncached states. Untruncated exhaustive results exclude incidence; limits remain unresolved and surviving full categories require actual metrics.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,default=0);ap.add_argument('--end',type=int,default=25);ap.add_argument('--nodes',type=int,default=1000);ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('PID',os.getpid(),flush=True);d=json.loads(json.dumps(run(a.start,a.end,a.nodes)));p=ROOT/f'results/{TAG}_{a.start}_{d["range"][1]}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['verdict_counts'])
