#!/usr/bin/env python3
"""Bounded diverse-parent full-label search beyond the structural gates."""
import argparse,json,os,time
from collections import defaultdict,Counter
from fresh_h9_fifteen_facet_20261014_incidence import ROOT,TAG,SOURCE,load,sha,expand_case
from fresh_h9_fifteen_facet_20261014_engine import search
INC=f'results/{TAG}_incidence.json';PAT=f'results/{TAG}_patterns.json';INITIAL=f'results/{TAG}_initial_searches.json'
def run(cases,limit,verify=False):
 inc=load(INC);source=load(SOURCE);p=load(PAT);first=json.loads((ROOT/INITIAL).read_text());old={r['case']for r in first}
 assert len(first)==50 and all(r['verdict']=='EXACT_COMPLETE_OBSTRUCTION'for r in first)
 # The initial batch deliberately tested the50 simplest necessary nerves.
 # These literal saved outcomes are independently rerun by the separate checker.
 groups=defaultdict(list)
 for i in p['remaining_cases']:
  if i in old:continue
  row=inc['cases'][i];groups[row['parent']].append(row)
 def priority(row):
  mm=row['canonical_minimal_nonfaces']
  return(-sum(m.bit_count()==2 for m in mm),sum(m.bit_count()>=4 for m in mm),row['case'])
 for rows in groups.values():rows.sort(key=priority)
 todo=[]
 while groups and len(todo)<cases:
  for parent in sorted(list(groups),key=lambda k:priority(groups[k][0])):
   todo.append(groups[parent].pop(0))
   if not groups[parent]:del groups[parent]
   if len(todo)>=cases:break
 records=[];start=time.monotonic()
 for row in todo:
  r=search(expand_case(row,source),limit);records.append(r)
  print('case',row['case'],'parent',row['parent'],r['verdict'],r['nodes'],'seconds',round(time.monotonic()-start,1),flush=True)
  if not verify:(ROOT/f'results/{TAG}_search_partial.json').write_text(json.dumps(dict(status='RUNNING',records=records),separators=(',',':'))+'\n')
  if r['verdict']=='FULL_CATEGORY_METRIC_REQUIRED':break
 paths=[INC,PAT,SOURCE,INITIAL,f'scripts/{TAG}_search.py',f'scripts/{TAG}_engine.py',f'scripts/{TAG}_incidence.py',
  'scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/fresh_h9_bugaenko_closure_20261014.py','scripts/fresh_h9_newfacet_frontier_20261014_relations.py','results/fresh_h9_newfacet_frontier_20261014_relations.json']
 combined=first+records
 return dict(status='BOUNDED_FIFTEEN_FACET_COMPLETE_LABEL_SEARCHES',sources={f:sha(f)for f in paths},
  initial_node_limit=3000,new_node_limit=limit,new_case_limit=cases,initial_case_count=50,starting_frontier=len(p['remaining_cases']),
  selected_new_case_ids=[r['case']for r in todo],cases=combined,verdict_counts=dict(Counter(r['verdict']for r in combined)),
  states=sum(r['nodes']for r in combined),scope='Initial50 simplest nerves followed by a bounded round-robin of distinct parent nerves, preferring more dotted pairs. Every EXCLUDED result exhausts all necessary finite categories; NODE_LIMIT is unresolved. A full category survivor still requires an actual Lorentzian metric, and category>=7 must then be resolved to an integer label. No unsearched case is excluded.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--cases',type=int,default=100);ap.add_argument('--nodes',type=int,default=2000);ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('PID',os.getpid(),flush=True);out=json.loads(json.dumps(run(a.cases,a.nodes,a.verify)));p=ROOT/f'results/{TAG}_search.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n');(ROOT/f'results/{TAG}_search_partial.json').unlink(missing_ok=True)
 print(out['status'],out['verdict_counts'],out['states'],flush=True)
