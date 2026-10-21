#!/usr/bin/env python3
"""Separate arithmetic replay of completed searches; traversal is reused."""
import argparse,json
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from sympy import GF
from sympy.polys.matrices import DomainMatrix
import fresh_h9_bugaenko_closure_20261014 as C
import fresh_h9_bugaenko_closure_20261014_search as P
from fresh_h9_bugaenko_surgery_20261014_rank_independent import spherical

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_closure_20261014_search_audit'
@lru_cache(None)
def finite(n,code):assert n==9;return spherical(tuple(code))
@lru_cache(None)
def determinant(code):
 p=241;field=GF(p);r={d:next(i for i in range(p)if i*i%p==d)for d in[2,3,5]}
 inv2=pow(2,-1,p);vals={2:0,3:inv2,4:r[2]*inv2%p,5:(1+r[5])*pow(4,-1,p)%p,6:r[3]*inv2%p}
 A=[[field(int(i==j))for j in range(11)]for i in range(11)];k=0
 for i in range(11):
  for j in range(i+1,11):A[i][j]=A[j][i]=field(-vals[code[k]]);k+=1
 return int(DomainMatrix(A,(11,11),field).det())%p
def run(path):
 d=json.loads((ROOT/path).read_text())
 for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 ip='results/fresh_h9_bugaenko_surgery_20261014_incidence.json';inc=json.loads((ROOT/ip).read_text())
 C.finite=finite;P.finite=finite;C.residue=determinant
 rows=[]
 for r in d['cases']:
  if r['verdict']!='EXACT_COMPLETE_OBSTRUCTION':continue
  q=P.search(inc['cases'][r['case']],d['node_limit']);assert json.loads(json.dumps(q))==r
  rows.append(dict(case=r['case'],evaluated_states=r['nodes'],visits=r['visits'],verdict=q['verdict']))
 files=[path,ip,str(Path(__file__).resolve().relative_to(ROOT)),
  'scripts/fresh_h9_bugaenko_closure_20261014_search.py','scripts/fresh_h9_bugaenko_closure_20261014.py',
  'scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
 return dict(status='PASS_SEPARATE_ARITHMETIC_REPLAY_OF_COMPLETE_CASES',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},input=path,range=d['range'],complete_cases_checked=len(rows),records=rows,finite_classifier_cache=list(finite.cache_info()),modular_determinant_cache=list(determinant.cache_info()),scope='The finite Coxeter diagram classifier is separately implemented and the modular determinant uses SymPy finite-field arithmetic. The DFS traversal, canonical-state construction, incidence catalog and Lanner tables are reused and source-bound. Limited cases are deliberately not counted as exclusions; this is not a second independently written search.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--input',default='results/fresh_h9_bugaenko_closure_20261014_search_0_25.json');ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run(a.input);lo,hi=d['range'];p=ROOT/f'results/{TAG}_{lo}_{hi}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['complete_cases_checked'])
