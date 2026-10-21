#!/usr/bin/env python3
"""Metric-directed search on the132 types surviving induced-nerve gates."""
import json,os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor,as_completed
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_newfacet_frontier_20261014_search'
INC='results/fresh_h9_bugaenko_remaining_20261014_newfacet.json';REL='results/fresh_h9_newfacet_frontier_20261014_relations.json'
def score(r):
 c=Counter(m.bit_count()for m in r['canonical_minimal_nonfaces']);return(-c[2],c[4]+c[5],len(r['canonical_minimal_nonfaces']),len(r['canonical_vertices']),r['case'])
def worker(row):
 from fresh_h9_newfacet_frontier_20261014_engine import search
 return search(row,3000)
def run():
 inc=json.loads((ROOT/INC).read_text());rel=json.loads((ROOT/REL).read_text())
 for d in(inc,rel):
  for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 cases=sorted((inc['cases'][i]for i in rel['remaining_case_ids']),key=score);records=[];metrics=[]
 with ProcessPoolExecutor(max_workers=2)as pool:
  jobs={pool.submit(worker,r):r['case']for r in cases}
  for f in as_completed(jobs):
   r=f.result();records.append(r);print(r['case'],r['verdict'],r['nodes'],flush=True)
 for r in records:
  if r['verdict']!='FULL_CATEGORY_METRIC_REQUIRED':continue
  from fresh_h9_bugaenko_surgery_20261014_metric_scout import fit
  row=inc['cases'][r['case']];vertices=[tuple(i for i in range(14)if v>>i&1)for v in row['canonical_vertices']]
  out=fit(r['label_vector'],vertices,r['dotted_pairs'],4000)
  metrics.append(dict(case=r['case'],**out));print('METRIC',r['case'],out['residual'],flush=True)
 paths=[INC,REL,f'scripts/{TAG}.py','scripts/fresh_h9_newfacet_frontier_20261014_engine.py','scripts/fresh_h9_bugaenko_closure_20261014.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py','scripts/fresh_h9_bugaenko_surgery_20261014_incidence.py','scripts/fresh_h9_bugaenko_surgery_20261014_metric_scout.py']
 return dict(status='BOUNDED_RELATION_PROPAGATED_METRIC_FRONTIER',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},node_limit=3000,case_order=[r['case']for r in cases],verdict_counts=dict(Counter(r['verdict']for r in records)),cases=sorted(records,key=lambda r:r['case']),metric_scouts=metrics,scope='Every132 type remaining after the induced-nerve gates is tested with the complete96-row overlapping-Lanner relation, shared-label propagation, all actual finite vertices, orthogonal-negative-block and fixed11-minor gates. Complete searches are exact incidence exclusions;3000-state limits remain unresolved. Every full necessary-label survivor receives a numerical Schur fit; metric failure alone is not an exclusion. All inherited field/finite classification dependencies remain explicit.')
if __name__=='__main__':
 print('PID',os.getpid(),flush=True);out=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['verdict_counts'],flush=True)
