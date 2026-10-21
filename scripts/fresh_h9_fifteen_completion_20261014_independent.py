#!/usr/bin/env python3
"""Streaming separate-arithmetic replay of immutable completion records.

The primary traversal is reused. A separate finite-diagram classifier and
SymPy GF241 determinants are substituted, and every output is compared.
"""
import argparse,json,os,time
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_fifteen_completion_20261014'
PRIMARY=f'results/{TAG}.json';FOLDER=ROOT/f'results/{TAG}_audit_cases'
INC='results/fresh_h9_fifteen_facet_20261014_incidence.json'
SRC='results/fresh_h9_gale_frontier_20261014_geometry.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def checked(p):
 d=read(p)
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def atomic(path,d):
 tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(d,separators=(',',':'))+'\n');tmp.replace(path)
def initialize():
 global ENGINE,INCIDENCE,GEOMETRY,RECONSTRUCT,MOVE
 import fresh_h9_fifteen_facet_20261014_engine as E
 from fresh_h9_bugaenko_closure_20261014_large_first_audit import finite,determinant
 from fresh_h9_fifteen_facet_20261014_incidence_independent import reconstruct,move
 E.finite=finite;E.residue=determinant;ENGINE=E;RECONSTRUCT=reconstruct;MOVE=move
 INCIDENCE=checked(INC);GEOMETRY=checked(SRC)
def replay(path,expected_hash):
 assert sha(path)==expected_hash;d=read(path);r=d['result'];i=r['case'];old=INCIDENCE['cases'][i]
 mm,vv=RECONSTRUCT(GEOMETRY['cases'][old['parent']],old['cut_face']);order=old['canonical_permutation'];inv=[order.index(j)for j in range(15)]
 assert sorted(MOVE(m,inv)for m in mm)==old['canonical_minimal_nonfaces']
 row=dict(old,canonical_vertices=sorted(MOVE(v,inv)for v in vv))
 answer=json.loads(json.dumps(ENGINE.search(row,d['node_limit'])));assert answer==r,(i,answer,r)
 return dict(case=i,status='PASS_SEPARATE_ARITHMETIC_LITERAL_CASE',input=path,input_sha256=expected_hash,
  verdict=r['verdict'],states=r['nodes'])
def run(jobs,follow):
 primary=checked(PRIMARY)
 paths=[INC,SRC,f'scripts/{TAG}_independent.py',f'scripts/{TAG}.py',
  'scripts/fresh_h9_fifteen_facet_20261014_engine.py','scripts/fresh_h9_fifteen_facet_20261014_incidence_independent.py',
  'scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py','scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
 sources={p:sha(p)for p in paths};source_manifest=primary['sources'];FOLDER.mkdir(exist_ok=True);records={}
 for path in sorted(FOLDER.glob('case_*.json')):
  d=checked(str(path.relative_to(ROOT)));assert d['sources']==sources
  r=d['result'];assert sha(r['input'])==r['input_sha256'];records[r['case']]=r
 pending={};start=time.monotonic();terminal=False
 def save(final=False):
  d=read(PRIMARY);available={r['case']for r in d['records']};assert set(records)<=available
  refs={r['case']:r for r in d['records']}
  for i,r in records.items():assert refs[i]['sha256']==r['input_sha256']
  manifest=dict(sources)
  if final:manifest[PRIMARY]=sha(PRIMARY)
  out=dict(status='PASS_SEPARATE_ARITHMETIC_ALL_SAVED_COMPLETION_CASES'if final and set(records)==available else'RUNNING_STREAMING_SEPARATE_ARITHMETIC',
   sources=manifest,primary_status=d['status'],cases_checked=len(records),primary_cases=len(available),states=sum(r['states']for r in records.values()),
   complete_exclusions_checked=sum(r['verdict']=='EXACT_COMPLETE_OBSTRUCTION'for r in records.values()),
   remaining_saved_cases=sorted(available-set(records)),records=[records[i]for i in sorted(records)],
   scope='Every listed case is freshly replayed with a separately implemented finite-tree classifier and SymPy GF241 determinants, comparing the entire result. Literal maximal faces are reconstructed separately. DFS traversal, propagation and classical Lanner tables are reused. Node-limited cases, if any, remain unresolved despite a matching replay.')
  atomic(ROOT/f'results/{TAG}_independent.json',out);return out
 print('PID',os.getpid(),'resumed',len(records),'workers',jobs,flush=True)
 with ProcessPoolExecutor(max_workers=jobs,initializer=initialize)as pool:
  while True:
   d=read(PRIMARY);assert d['sources']==source_manifest
   terminal=not d['status'].startswith(('RUNNING','STOPPED_DISPATCH'))
   waiting={i for i,p,h in pending.values()}
   todo=[r for r in d['records']if r['case']not in records and r['case']not in waiting]
   for r in todo[:jobs-len(pending)]:
    case=read(r['path']);assert case['sources']==source_manifest
    pending[pool.submit(replay,r['path'],r['sha256'])]=(r['case'],r['path'],r['sha256'])
   if not pending:
    if terminal or not follow:break
    time.sleep(2);continue
   done,_=wait(pending,timeout=3,return_when=FIRST_COMPLETED)
   for f in done:
    i,p,h=pending.pop(f);r=f.result();assert r['case']==i
    atomic(FOLDER/f'case_{i:05d}.json',dict(sources=sources,result=r));records[i]=r
    print('audit',i,r['verdict'],'states',r['states'],'done',len(records),'seconds',round(time.monotonic()-start,1),flush=True)
   if done:save()
 return save(final=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=2);ap.add_argument('--follow',action='store_true');a=ap.parse_args();assert 1<=a.jobs<=2
 out=run(a.jobs,a.follow);print(out['status'],out['cases_checked'],out['states'],flush=True)
