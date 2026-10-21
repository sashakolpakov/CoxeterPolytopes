#!/usr/bin/env python3
"""Separate-arithmetic replay of the small final metric-priority batch."""
import argparse,json,time
from fresh_h9_fifteen_facet_20261014_search_independent import ROOT,TAG,SRC,INC,load,sha,E,finite,determinant,reconstruct,move
INPUT=f'results/{TAG}_metric_priority.json'
def run():
 d=load(INPUT);inc=load(INC);source=load(SRC);E.finite=finite;E.residue=determinant;rows=[];start=time.monotonic()
 for r in d['cases']:
  assert r['verdict']=='EXACT_COMPLETE_OBSTRUCTION'
  old=inc['cases'][r['case']];mm,vv=reconstruct(source['cases'][old['parent']],old['cut_face'])
  order=old['canonical_permutation'];inverse=[order.index(i)for i in range(15)]
  assert sorted(move(m,inverse)for m in mm)==old['canonical_minimal_nonfaces']
  row=dict(old,canonical_vertices=sorted(move(v,inverse)for v in vv))
  out=json.loads(json.dumps(E.search(row,d['new_node_limit'])));assert out==r
  rows.append(dict(case=r['case'],states=r['nodes']))
  print('metric-priority replay',len(rows),r['case'],'seconds',round(time.monotonic()-start,1),flush=True)
 paths=[INPUT,INC,SRC,f'scripts/{TAG}_metric_priority_independent.py',f'scripts/{TAG}_search_independent.py',f'scripts/{TAG}_engine.py',
  f'scripts/{TAG}_incidence_independent.py','scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py','scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
 return dict(status='PASS_SEPARATE_ARITHMETIC_METRIC_PRIORITY_BATCH',sources={p:sha(p)for p in paths},complete_cases_checked=len(rows),states=sum(r['states']for r in rows),records=rows,
  scope='Separate finite-tree classifier and SymPy finite-field determinants; shared DFS traversal and Lanner tables. Literal incidence reconstructed separately. All24 completed cases are replayed; unsearched nerves remain open.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}_metric_priority_independent.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['complete_cases_checked'],d['states'],flush=True)
