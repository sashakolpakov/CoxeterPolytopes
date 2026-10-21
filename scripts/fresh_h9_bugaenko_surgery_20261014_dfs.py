#!/usr/bin/env python3
"""Depth-first necessary-label search, with explicit unresolved node limits."""
import argparse,json,time,os
from collections import Counter
from itertools import combinations
from hashlib import sha256
from pathlib import Path
from fresh_h9_highbase_construction_20261014_cell import codes,finite

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_surgery_20261014';SOURCE=f'results/{TAG}_incidence.json'
def tup(m):return tuple(i for i in range(m.bit_length())if m>>i&1)
def search(row,limit):
 pairs=list(combinations(range(14),2));pi={p:i for i,p in enumerate(pairs)};M=list(map(tup,row['canonical_minimal_nonfaces']));V=list(map(tup,row['canonical_vertices']));supports=[a for a in M if len(a)>=3]
 con=[(tuple(pi[p]for p in combinations(a,2)),codes(len(a)))for a in supports]
 vpairs=[tuple(pi[p]for p in combinations(v,2))for v in V]
 affected=[[j for j,pv in enumerate(vpairs)if len(set(a)&set(pv))]for a,_ in con]
 labels=[2]*91;assigned=0;visits=0;terminal=None;dead=set();checks=0
 def walk(left,bits):
  nonlocal visits,terminal,checks
  visits+=1
  if visits>limit:return None
  if not left:terminal=labels.copy();return True
  key=(left,tuple((i,labels[i])for i in range(91)if bits>>i&1))
  if key in dead:return False
  candidates=[k for k in range(len(con))if left>>k&1]
  k=max(candidates,key=lambda k:(sum(bits>>p&1 for p in con[k][0]),-sum(not(bits>>p&1)for p in con[k][0]),len(supports[k]),-k))
  pp,options=con[k];new=[p for p in pp if not(bits>>p&1)];nbits=bits|sum(1<<p for p in new);rest=left&~(1<<k)
  known=[(i,p)for i,p in enumerate(pp)if bits>>p&1]
  for code in options:
   if any(code[i]!=labels[p]for i,p in known):continue
   for p,m in zip(pp,code):labels[p]=m
   valid=True
   for vi in affected[k]:
    checks+=1
    if not finite(9,tuple(labels[p]for p in vpairs[vi])):valid=False;break
   if valid:
    # Every pending Lanner table must retain at least one shared-pair match.
    for j in candidates:
     if j==k:continue
     jj,oo=con[j];known_j=[(i,p)for i,p in enumerate(jj)if nbits>>p&1]
     if known_j and not any(all(c[i]==labels[p]for i,p in known_j)for c in oo):valid=False;break
   if valid:
    result=walk(rest,nbits)
    if result is not False:
     for p in new:labels[p]=2
     return result
   for p in new:labels[p]=2
  dead.add(key);return False
 result=walk((1<<len(con))-1,0)
 out=dict(case=row['case'],nodes=visits,vertex_checks=checks,node_limit=limit,verdict='NECESSARY_LABELING_FOUND'if result is True else('EXACT_DFS_LABEL_OBSTRUCTION'if result is False else'NODE_LIMIT_UNRESOLVED'))
 if terminal is not None:
  assert all(finite(9,tuple(terminal[p]for p in a))for a in vpairs)
  assert all(tuple(terminal[p]for p in a)in oo for a,oo in con)
  covered={p for a,_ in con for p in a};missingpairs={pi[a]for a in M if len(a)==2}
  out.update(label_vector=terminal,pair_order=[list(p)for p in pairs],unconstrained_meeting_pairs=[list(pairs[p])for p in range(91)if p not in covered|missingpairs],dotted_pairs=[list(a)for a in M if len(a)==2])
 return out

def run(start,end,limit):
 src=json.loads((ROOT/SOURCE).read_text());end=min(end,len(src['cases']));rows=[];counts=Counter();t=time.time();print('PID',os.getpid(),'cases',start,end,flush=True)
 for r in src['cases'][start:end]:
  s=search(r,limit);rows.append(s);counts[s['verdict']]+=1
  if len(rows)%10==0 or end-start<=25:print('done',len(rows),'case',r['case'],s['verdict'],'nodes',s['nodes'],'seconds',round(time.time()-t,2),flush=True)
 paths=[SOURCE,f'scripts/{TAG}_dfs.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='EXACT_DEPTH_FIRST_NECESSARY_LABEL_CHECKPOINT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},range=[start,end],total_types=len(src['cases']),node_limit=limit,verdict_counts=dict(counts),cases=rows,scope='Necessary exact Lanner and all-vertex finite categories. Unassigned meeting pairs are relaxed to2, dotted distances unassigned. A found labeling is a metric proposal only; node limits are unresolved.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--start',type=int,default=0);p.add_argument('--end',type=int,default=25);p.add_argument('--nodes',type=int,default=50000);p.add_argument('--verify',action='store_true');a=p.parse_args();out=run(a.start,a.end,a.nodes);path=ROOT/f'results/{TAG}_dfs_{a.start}_{out["range"][1]}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['verdict_counts'],flush=True)
