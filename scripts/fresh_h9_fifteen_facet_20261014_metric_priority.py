#!/usr/bin/env python3
"""Small final search of untested nerves with more metric freedom.

Select one representative per minimal-nonface count signature, preferring
more dotted pairs and missing triangles. Stop at a full category survivor.
"""
import argparse,json,time,os
from collections import Counter
from fresh_h9_fifteen_facet_20261014_incidence import ROOT,TAG,SOURCE,load,sha,expand_case
from fresh_h9_fifteen_facet_20261014_engine import search
INC=f'results/{TAG}_incidence.json';PAT=f'results/{TAG}_patterns.json';FIRST=f'results/{TAG}_search.json'
def run(count,limit):
 inc=load(INC);pat=load(PAT);first=load(FIRST);source=load(SOURCE);done={r['case']for r in first['cases']}
 def signature(i):
  mm=inc['cases'][i]['canonical_minimal_nonfaces']
  return tuple(sum(m.bit_count()==k for m in mm)for k in range(2,6))
 def priority(i):
  d,t,f,p=signature(i)
  return(-d,-t,p,f,i)
 pool=sorted(set(pat['remaining_cases'])-done,key=priority);seen=set();chosen=[]
 for i in pool:
  s=signature(i)
  if s in seen:continue
  seen.add(s);chosen.append(i)
  if len(chosen)==count:break
 records=[];start=time.monotonic()
 for i in chosen:
  r=search(expand_case(inc['cases'][i],source),limit);records.append(r)
  print(i,signature(i),r['verdict'],r['nodes'],'seconds',round(time.monotonic()-start,1),flush=True)
  if r['verdict']=='FULL_CATEGORY_METRIC_REQUIRED':break
 paths=[INC,PAT,FIRST,SOURCE,f'scripts/{TAG}_metric_priority.py',f'scripts/{TAG}_engine.py',f'scripts/{TAG}_incidence.py']
 return dict(status='BOUNDED_METRIC_PRIORITY_FIFTEEN_FACET_SEARCH',sources={p:sha(p)for p in paths},
  initial_case_count=0,initial_node_limit=limit,new_node_limit=limit,new_case_limit=count,starting_frontier=len(pool),
  chosen_cases=chosen,cases=records,verdict_counts=dict(Counter(r['verdict']for r in records)),states=sum(r['nodes']for r in records),
  scope='One untested representative per missing-face count signature, prioritizing more dotted distances and missing triangles, is tested with a finite node limit. A survivor is only a necessary full-label category assignment and would require an actual metric. This bounded selection does not cover the remaining fifteen-facet family.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--cases',type=int,default=24);ap.add_argument('--nodes',type=int,default=5000);ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('PID',os.getpid(),flush=True)
 out=json.loads(json.dumps(run(a.cases,a.nodes)));p=ROOT/f'results/{TAG}_metric_priority.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['verdict_counts'],out['states'],flush=True)
