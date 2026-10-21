#!/usr/bin/env python3
"""Reduce proposal-graph escape clauses, then classify every core pair exactly."""
import argparse,json
from itertools import combinations
import numpy as np
import z3
from scripts.fresh_h8_cap_completion_20261013 import load,ROOT,read,gate,real

def run(source,output,exact=True):
 data=read(source);G,I,center,nodes,polys,cn,qn,hn,oldgraph=load();N=len(nodes);active=data['active_nodes'];aset=set(active);allowed=[float((lo+hi)/8)for _,_,lo,hi in polys];gn=np.array([[real(x)for x in r]for r in G]);p=cn@gn@cn.T;norms=np.array([real(n['q'])for n in nodes]);rat=p*p/(norms[:,None]*norms[None,:]);solver=z3.Solver();xs=[z3.Bool(f'x{i}')for i in range(N)]
 if data['polar_required']:solver.add(xs[348])
 for i in set(range(N))-aset:solver.add(z3.Not(xs[i]))
 for i,j in combinations(active,2):
  if i<493 and j<493:possible=f'{i},{j}'in oldgraph['edges']
  else:possible=p[i,j]<=1e-6 and(rat[i,j]>1+1e-7 or min(abs(rat[i,j]-v)for v in allowed)<1e-6)
  if f'{i},{j}'in data['exact_pair_checks']:possible=data['exact_pair_checks'][f'{i},{j}']['compatible']
  if not possible:solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
 cs=[z3.Bool(f'e{i}')for i in range(len(data['escape_clauses']))]
 for c,row in zip(cs,data['escape_clauses']):solver.add(z3.Implies(c,z3.Or(*(xs[i]for i in row['complete_dictionary_cutters']))))
 assert solver.check(*cs)==z3.unsat
 core=sorted(int(str(c)[1:])for c in solver.unsat_core())
 for i in core.copy():
  trial=[j for j in core if j!=i]
  if solver.check(*(cs[j]for j in trial))==z3.unsat:core=trial
 union=set().union(*(set(data['escape_clauses'][i]['complete_dictionary_cutters'])for i in core))
 if data['polar_required']:union.add(348)
 union=sorted(union);out={'source':source,'status':'PROPOSAL_GRAPH_CORE','clause_indices':core,'clauses':[data['escape_clauses'][i]for i in core],'core_nodes':union,'polar_required':data['polar_required'],'scope':'Initial reduction uses the numerical proposal graph; only a completed exact angle classification plus exact Boolean proof can certify a finite-bank obstruction.'};(ROOT/output).write_text(json.dumps(out,indent=2)+'\n');print({'core':core,'clause_sizes':[len(r['complete_dictionary_cutters'])for r in out['clauses']],'union_nodes':len(union)},flush=True)
 if not exact:return
 checks=[];counts={};s=z3.Solver();x={i:z3.Bool(f'u{i}')for i in union}
 if data['polar_required']:s.add(x[348])
 for row in out['clauses']:s.add(z3.Or(*(x[i]for i in row['complete_dictionary_cutters'])))
 for i,j in combinations(union,2):
  row=gate(nodes[i],nodes[j],polys);checks.append({'i':i,'j':j,**row});counts[row['kind']]=counts.get(row['kind'],0)+1
  # Include parallel pairs: this is the finite-volume graph, not just compact.
  if not row['compatible']and row['kind']!='parallel':s.add(z3.Or(z3.Not(x[i]),z3.Not(x[j])))
 status=s.check();out['exact_pair_checks']=checks;out['exact_pair_kind_counts']=counts;out['finite_volume_core_status']=str(status);out['complete_safe_labels']=[r[0]for r in polys]
 if status==z3.sat:out['surviving_exact_core_indices']=[i for i in union if z3.is_true(s.model().eval(x[i],model_completion=True))];out['status']='NUMERIC_GRAPH_OMITTED_A_FINITE_VOLUME_CORE_SOLUTION'
 elif status==z3.unsat:out['status']='EXACT_FINITE_VOLUME_CORE_UNSAT'
 else:out['status']='EXACT_CORE_SOLVER_UNDECIDED'
 for row in checks:
  if row['kind']=='parallel':s.add(z3.Or(z3.Not(x[row['i']]),z3.Not(x[row['j']])))
 compact=s.check();out['compact_core_status']=str(compact)
 if compact==z3.sat:out['surviving_compact_core_indices']=[i for i in union if z3.is_true(s.model().eval(x[i],model_completion=True))]
 elif compact==z3.unsat and status==z3.sat:out['status']='EXACT_COMPACT_CORE_UNSAT_FINITE_VOLUME_CORE_SAT'
 (ROOT/output).write_text(json.dumps(out,indent=2)+'\n');print({'exact_status':out['status'],'pair_counts':counts,'survivor':out.get('surviving_exact_core_indices')},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',default='data/enumerations/fresh_h8_cap_completion_20261013.json');p.add_argument('--output',default='data/enumerations/fresh_h8_cap_completion_20261013_core.json');p.add_argument('--numeric-only',action='store_true');a=p.parse_args();run(a.source,a.output,not a.numeric_only)
