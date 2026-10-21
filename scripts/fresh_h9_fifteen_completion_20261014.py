#!/usr/bin/env python3
"""Resumable necessary-label completion of the1,334 frozen open nerves.

Each immutable case record binds this driver and the unchanged exact engine.
At most two cases are dispatched concurrently. The first category survivor
stops dispatch; cases that were already running finish and are archived too.
"""
import argparse,json,os,time
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from collections import Counter
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_fifteen_completion_20261014'
OLD='fresh_h9_fifteen_facet_20261014'
LEDGER=f'results/{OLD}_checkpoint.json';INC=f'results/{OLD}_incidence.json'
SOURCE='results/fresh_h9_gale_frontier_20261014_geometry.json'
ENGINE=f'scripts/{OLD}_engine.py';FOLDER=ROOT/f'results/{TAG}_cases'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def manifest():
 paths=[LEDGER,INC,SOURCE,ENGINE,f'scripts/{TAG}.py',f'scripts/{OLD}_incidence.py',
  'scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/fresh_h9_bugaenko_closure_20261014.py',
  'scripts/fresh_h9_newfacet_frontier_20261014_relations.py','results/fresh_h9_newfacet_frontier_20261014_relations.json']
 return {p:sha(p)for p in paths}
def initialize():
 global INCIDENCE,GEOMETRY
 INCIDENCE=load(INC);GEOMETRY=load(SOURCE)
def worker(i,limit):
 from fresh_h9_fifteen_facet_20261014_incidence import expand_case
 from fresh_h9_fifteen_facet_20261014_engine import search
 row=expand_case(INCIDENCE['cases'][i],GEOMETRY)
 return json.loads(json.dumps(search(row,limit)))
def atomic(path,d):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(d,separators=(',',':'))+'\n');tmp.replace(path)
def run(limit,jobs,max_cases):
 ledger=load(LEDGER);inc=load(INC);ids=ledger['open_case_ids'];assert len(ids)==len(set(ids))==1334
 source=manifest();FOLDER.mkdir(exist_ok=True);records={};refs={}
 for p in sorted(FOLDER.glob('case_*.json')):
  d=load(str(p.relative_to(ROOT)));assert d['sources']==source and d['node_limit']==limit
  i=d['result']['case'];assert i in ids and i not in records
  records[i]=d['result'];refs[i]=str(p.relative_to(ROOT))
 def priority(i):
  mm=inc['cases'][i]['canonical_minimal_nonfaces']
  return(-sum(m.bit_count()==2 for m in mm),-sum(m.bit_count()==3 for m in mm),sum(m.bit_count()==5 for m in mm),i)
 todo=sorted(set(ids)-records.keys(),key=priority)
 if max_cases is not None:todo=todo[:max_cases]
 survivor=next((i for i,r in records.items()if r['verdict']=='FULL_CATEGORY_METRIC_REQUIRED'),None)
 summary_path=ROOT/f'results/{TAG}.json';start=time.monotonic();new_count=0
 def save(status):
  remaining=sorted(set(ids)-records.keys());limited=sorted(i for i,r in records.items()if r['verdict']=='NODE_LIMIT_UNRESOLVED')
  data=dict(status=status,sources=source,node_limit_per_case=limit,workers=jobs,starting_open_cases=1334,
   completed_cases=len(records),new_cases_this_invocation=new_count,verdict_counts=dict(Counter(r['verdict']for r in records.values())),
   total_states=sum(r['nodes']for r in records.values()),first_survivor=survivor,remaining_unsearched=remaining,node_limited_cases=limited,
   records=[dict(case=i,verdict=records[i]['verdict'],nodes=records[i]['nodes'],path=refs[i],sha256=sha(refs[i]))for i in sorted(records)],
   scope='Only exact-complete searches exclude cases. Each case has the stated state limit; node-limited and undispatched cases remain open. Category>=7 means every finite integer label>=7, never a literal m7 restriction. No fourteen-facet special restrictions are used. A category survivor requires an actual Lorentzian metric; no metric is asserted by this search.')
  atomic(summary_path,data);return data
 if survivor is not None:return save('STOPPED_AT_FULL_CATEGORY_SURVIVOR')
 iterator=iter(todo);pending={}
 print('PID',os.getpid(),'resumed',len(records),'queued',len(todo),'limit',limit,'workers',jobs,flush=True)
 save('RUNNING_RESUMABLE_COMPLETE_LABEL_PASS')
 with ProcessPoolExecutor(max_workers=jobs,initializer=initialize)as pool:
  def dispatch():
   try:i=next(iterator)
   except StopIteration:return False
   pending[pool.submit(worker,i,limit)]=i;return True
  for k in range(jobs):dispatch()
  while pending:
   done,_=wait(pending,return_when=FIRST_COMPLETED)
   for future in done:
    i=pending.pop(future);r=future.result();assert r['case']==i
    file=FOLDER/f'case_{i:05d}.json';atomic(file,dict(status='EXACT_CASE_RESULT',sources=source,node_limit=limit,result=r))
    records[i]=r;refs[i]=str(file.relative_to(ROOT));new_count+=1
    print('case',i,r['verdict'],'states',r['nodes'],'done',len(records),'seconds',round(time.monotonic()-start,1),flush=True)
    if r['verdict']=='FULL_CATEGORY_METRIC_REQUIRED'and survivor is None:survivor=i
   save('RUNNING_RESUMABLE_COMPLETE_LABEL_PASS'if survivor is None else'STOPPED_DISPATCH_AT_FULL_CATEGORY_SURVIVOR')
   if survivor is None:
    for k in range(jobs-len(pending)):
     if not dispatch():break
 status='STOPPED_AT_FULL_CATEGORY_SURVIVOR'if survivor is not None else('FINISHED_ALL_1334_BOUNDED_CASES'if len(records)==1334 else'CHECKPOINT_WITH_UNSEARCHED_CASES')
 return save(status)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--nodes',type=int,default=20000);ap.add_argument('--jobs',type=int,default=2);ap.add_argument('--cases',type=int);a=ap.parse_args()
 assert a.nodes>0 and 1<=a.jobs<=2
 out=run(a.nodes,a.jobs,a.cases);print(out['status'],out['verdict_counts'],out['total_states'],flush=True)
