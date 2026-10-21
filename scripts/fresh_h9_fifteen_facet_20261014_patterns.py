#!/usr/bin/env python3
"""New bounded local-pattern trials, allowing fifteen-facet missing triangles."""
import argparse,json,os,time
from collections import Counter,defaultdict
from itertools import combinations
from pathlib import Path
from fresh_h9_fifteen_facet_20261014_incidence import ROOT,TAG,sha,load,bits
from fresh_h9_bugaenko_surgery_20261014_incidence import canonical
from fresh_h9_newfacet_frontier_20261014_patterns import local_search
INC=f'results/{TAG}_incidence.json'
def run(trials,limit):
 data=load(INC);todo_cases=[r for r in data['cases']if r['new_obstruction']is None]
 bank=defaultdict(dict);cache={};start=time.monotonic()
 for row in todo_cases:
  mm=row['canonical_minimal_nonfaces']
  unions={a|b for a,b in combinations(mm,2)if min(a.bit_count(),b.bit_count())>=3 and(a|b).bit_count()<=9}
  for u in unions:
   inside=[m for m in mm if m&u==m]
   if any(m.bit_count()==2 for m in inside):continue
   vv=bits(u);red=tuple(sorted(sum(1<<k for k,j in enumerate(vv)if m>>j&1)for m in inside));key=(len(vv),red)
   if key not in cache:cache[key]=canonical(red,len(vv))
   cc,order=cache[key];bank[len(vv),cc].setdefault(row['case'],[vv[j]for j in order])
 print('inventory',len(bank),'cases',len(todo_cases),'seconds',round(time.monotonic()-start,1),flush=True)
 open_ids={r['case']for r in todo_cases};todo=set(bank);records=[];closed={}
 for k in range(trials):
  useful=[p for p in todo if set(bank[p])&open_ids]
  if not useful:break
  key=min(useful,key=lambda p:(-len(set(bank[p])&open_ids),p[0],-len(p[1]),p[1]));todo.remove(key)
  n,mm=key;out=local_search(n,mm,limit)
  record=dict(pattern=len(records),order=n,minimal_nonfaces=mm,applications=sorted(bank[key].items()),**out)
  if out['verdict']=='EXACT_LOCAL_OBSTRUCTION':
   for case,perm in record['applications']:closed.setdefault(case,len(records));open_ids.discard(case)
  records.append(record)
  print('trial',k,'n',n,'missing',[m.bit_count()for m in mm],out['verdict'],out['nodes'],'open',len(open_ids),flush=True)
 paths=[INC,f'scripts/{TAG}_patterns.py',f'scripts/{TAG}_incidence.py','scripts/fresh_h9_newfacet_frontier_20261014_patterns.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='BOUNDED_FIFTEEN_FACET_LOCAL_PATTERN_TRIALS',sources={p:sha(p)for p in paths},trial_limit=trials,node_limit=limit,
  starting_frontier=len(todo_cases),pattern_inventory=len(bank),verdict_counts=dict(Counter(r['verdict']for r in records)),patterns=records,
  excluded_cases=sorted(closed.items()),remaining_cases=sorted(open_ids),scope='All recorded exact local failures are exhaustive necessary finite/Lanner constraints; local label survivors and node limits make no exclusion. Missing triangles are allowed. Finite categories2..6 plus7 represent every finite integer label. Only the selected bounded pattern trials are claimed; remaining nerves still need global labels and metric checks.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=30);ap.add_argument('--nodes',type=int,default=10000);ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('PID',os.getpid(),flush=True);out=json.loads(json.dumps(run(a.trials,a.nodes)));p=ROOT/f'results/{TAG}_patterns.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
 print(out['status'],len(out['excluded_cases']),len(out['remaining_cases']),flush=True)
