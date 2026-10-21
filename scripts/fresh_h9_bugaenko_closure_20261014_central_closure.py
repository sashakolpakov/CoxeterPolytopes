#!/usr/bin/env python3
"""Bind the complete central-carrier incidence census to exact exclusions."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_closure_20261014_central_closure'
def read(p):return json.loads((ROOT/p).read_text())
def bind(p):
 d=read(p)
 for source,h in d['sources'].items():assert sha256((ROOT/source).read_bytes()).hexdigest()==h,source
 return d
def run():
 cp='results/fresh_h9_bugaenko_closure_20261014_central.json';central=bind(cp)
 lp='results/fresh_h9_bugaenko_closure_20261014_large_first_20_21.json';last=bind(lp)
 ap='results/fresh_h9_bugaenko_closure_20261014_large_first_audit_20_21.json';audit=bind(ap)
 assert audit['complete_cases_checked']==1 and audit['records'][0]['case']==20
 providers={int(i):p for i,p in central['prior_exact_exclusions'].items()}
 for r in central['cases']:
  if r['verdict']=='EXACT_COMPLETE_OBSTRUCTION':providers[r['case']]=cp
 for r in last['cases']:
  if r['verdict']=='EXACT_COMPLETE_OBSTRUCTION':providers[r['case']]=lp
 ids=sorted({r['case']for r in central['eligible_recipes']});assert len(ids)==128 and len(central['eligible_recipes'])==208
 used={i:providers[i]for i in ids};cache={};records=[]
 for i,p in used.items():
  if p not in cache:cache[p]=bind(p)
  data=cache[p];r=next(r for r in data['cases']if r['case']==i)
  if r['verdict']=='EXACT_LANNER_VERTEX_OBSTRUCTION':assert not r['state_limit_reached']and r['stages'][-1]['passing_states']==0
  else:assert r['verdict']=='EXACT_COMPLETE_OBSTRUCTION'and r['nodes']<=data['node_limit']
  records.append(dict(case=i,exact_source=p,verdict=r['verdict']))
 assert central['rejected_large_missing_face_recipes']==15582
 sources={cp,lp,ap,*used.values(),str(Path(__file__).resolve().relative_to(ROOT))}
 return dict(status='PASS_COMPLETE_CENTRAL_SINGLETON_CARRIER_EXCLUSION',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(sources)},dimension=9,facets=14,base='Bugaenko compact H8',wedge_carrier=4,ordered_truncation_codimensions=[2,3,4,5],ordered_operation_orbits=15790,excluded_by_missing_order=15582,eligible_operation_orbits=208,distinct_eligible_types=128,exclusion_records=records,source_counts=dict(Counter(used.values())),scope='Every central-singleton wedge followed by two ordered face truncations of codimensions2..5 is excluded as a compact H9 Coxeter incidence. Includes second truncations through the first new facet and all overlaps with types first introduced by another carrier. No metric-field or distance cutoff. Completeness of the incidence and Lanner tables remains in the bound upstream certificates; this aggregation does not rerun all searches.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==read(str(p.relative_to(ROOT)))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['distinct_eligible_types'])
