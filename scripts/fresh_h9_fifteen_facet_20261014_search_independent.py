#!/usr/bin/env python3
"""Separate arithmetic replay of every completed fifteen-facet label search.

The search traversal and classical Lanner tables are reused. Spherical
classification and finite-field determinants have separate implementations.
The incidence is rebuilt from the old maximal faces, not producer cut code.
"""
import argparse,json,time
from hashlib import sha256
from pathlib import Path
import fresh_h9_fifteen_facet_20261014_engine as E
from fresh_h9_bugaenko_closure_20261014_large_first_audit import finite,determinant
from fresh_h9_fifteen_facet_20261014_incidence_independent import reconstruct,move
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_fifteen_facet_20261014'
SRC='results/fresh_h9_gale_frontier_20261014_geometry.json'
INC=f'results/{TAG}_incidence.json';SEARCH=f'results/{TAG}_search.json';AUDIT=f'results/{TAG}_incidence_independent.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def run():
 source=load(SRC);inc=load(INC);search=load(SEARCH);load(AUDIT);start=time.monotonic()
 E.finite=finite;E.residue=determinant;records=[]
 for j,r in enumerate(search['cases']):
  if r['verdict']!='EXACT_COMPLETE_OBSTRUCTION':continue
  old=inc['cases'][r['case']];mm,vv=reconstruct(source['cases'][old['parent']],old['cut_face'])
  order=old['canonical_permutation'];inverse=[order.index(i)for i in range(15)]
  assert sorted(move(m,inverse)for m in mm)==old['canonical_minimal_nonfaces']
  row=dict(old,canonical_vertices=sorted(move(v,inverse)for v in vv))
  limit=search['initial_node_limit']if j<search['initial_case_count']else search['new_node_limit']
  replay=json.loads(json.dumps(E.search(row,limit)));assert replay==r,(r['case'],replay,r)
  records.append(dict(case=r['case'],states=r['nodes'],verdict=replay['verdict']))
  print('separate arithmetic',len(records),r['case'],r['nodes'],'seconds',round(time.monotonic()-start,1),flush=True)
 paths=[SRC,INC,SEARCH,AUDIT,f'scripts/{TAG}_search_independent.py',f'scripts/{TAG}_engine.py',f'scripts/{TAG}_incidence_independent.py',
  'scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py','scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
 return dict(status='PASS_SEPARATE_ARITHMETIC_ALL_COMPLETED_FIFTEEN_FACET_SEARCHES',sources={p:sha(p)for p in paths},
  complete_cases_checked=len(records),states=sum(r['states']for r in records),records=records,
  scope='Every completed case, including the50 initially saved trials, is freshly replayed with a separate finite-tree classifier and SymPy GF241 determinants. Maximal faces are rebuilt by a separate set-face truncation routine. DFS propagation/canonical-state traversal and classical Lanner tables are reused and source-bound. No node-limited or unsearched case is excluded, and this is not an independently written full search.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_search_independent.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['complete_cases_checked'],d['states'],flush=True)
