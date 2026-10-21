#!/usr/bin/env python3
"""Re-arithmetize the132 complete global searches; traversal is shared."""
import argparse,json,os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_newfacet_frontier_20261014_search_independent'
INPUT='results/fresh_h9_newfacet_frontier_20261014_search.json';INC='results/fresh_h9_bugaenko_remaining_20261014_newfacet.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def worker(payload):
 row,saved,limit=payload
 import fresh_h9_bugaenko_closure_20261014 as C
 import fresh_h9_newfacet_frontier_20261014_engine as E
 from fresh_h9_bugaenko_closure_20261014_large_first_audit import finite,determinant
 C.finite=E.finite=finite;C.residue=determinant
 result=json.loads(json.dumps(E.search(row,limit)))
 assert result==saved and result['verdict']=='EXACT_COMPLETE_OBSTRUCTION'
 return dict(case=row['case'],states=result['nodes'],visits=result['visits'],six_node_relations=result['six_node_relations'])
def run():
 d=json.loads((ROOT/INPUT).read_text());inc=json.loads((ROOT/INC).read_text())
 local='results/fresh_h9_newfacet_frontier_20261014_local_independent.json';audit=json.loads((ROOT/local).read_text())
 for doc in(d,inc,audit):
  for p,h in doc['sources'].items():assert sha(p)==h,p
 assert audit['six_node_relation_rows']==96 and audit['obstruction_count']==44
 rows=[]
 with ProcessPoolExecutor(max_workers=2)as pool:
  jobs=[pool.submit(worker,(inc['cases'][r['case']],r,d['node_limit']))for r in d['cases']]
  for f in as_completed(jobs):
   row=f.result();rows.append(row)
   if len(rows)%20==0:print('replayed',len(rows),flush=True)
 assert len(rows)==132
 paths=[INPUT,INC,local,f'scripts/{TAG}.py','scripts/fresh_h9_newfacet_frontier_20261014_engine.py','scripts/fresh_h9_bugaenko_closure_20261014.py','scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py','scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
 return dict(status='PASS_SEPARATE_FINITE_ARITHMETIC_ALL132_SEARCHES',sources={p:sha(p)for p in paths},complete_cases=len(rows),total_states=sum(r['states']for r in rows),maximum_states=max(r['states']for r in rows),records=sorted(rows,key=lambda r:r['case']),scope='All132 saved complete searches are reproduced using a separate finite-tree classifier and SymPy finite-field determinant implementation. The DFS, twin canonicalization, exact relation consistency and variable ordering are shared, not independently reimplemented. The complete96-row relation and44 local obstructions have a separate-traversal replay. No limits, full necessary label survivors or metric attempts occur in these132 cases.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('PID',os.getpid(),flush=True);out=run();path=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['total_states'],flush=True)
