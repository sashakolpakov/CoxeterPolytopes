#!/usr/bin/env python3
"""Reduce the numerical proposal graph's escape clauses for adaptive guidance."""
import json
from fractions import Fraction as F
from itertools import combinations
import numpy as np
import z3
from scripts.fresh_h8_alternate_selection_20261013 import load
from scripts.fresh_h8_alternate_caps_20261013 import *

ROOT=ROOT
OUTPUT=ROOT/'data/enumerations/fresh_h8_threefree_core_20261013.json'

def run():
 data=json.loads((ROOT/'data/enumerations/fresh_h8_threefree_selection_20261013.json').read_text());new=json.loads((ROOT/'data/enumerations/fresh_h8_threefree_caps_20261013.json').read_text());oldgraph=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_selection_20261013.json').read_text());angledata=json.loads((ROOT/'results/fresh_h8_alternate_core_angles_independent_20261013.json').read_text());base,pol,raw,G,I,center,nodes=load()
 for j,c in enumerate(new['caps']):nodes.append({'source':'threefree','index':j,'d':parse(c['discriminant']),'a':tuple(map(parse,c['coefficients_rational'])),'b':tuple(map(parse,c['coefficients_radical'])),'pa':tuple(map(parse,c['pairings_rational'])),'pb':tuple(map(parse,c['pairings_radical'])),'q':scale(O,2)})
 active=data['active_nodes'];N=len(nodes);cn=np.array([[real(x)+real(r['d'])**.5*real(y)for x,y in zip(r['a'],r['b'])]for r in nodes]);gn=np.array([[real(x)for x in r]for r in G]);qs=np.array([real(r['q'])for r in nodes]);pn=cn@gn@cn.T;ratio=pn**2/(qs[:,None]*qs[None,:]);allowed=[float((F(r['largest_root_interval'][0])+F(r['largest_root_interval'][1]))/8)for r in angledata['real_cyclotomic_polynomials']]
 solver=z3.Solver();xs=[z3.Bool(f'w{i}')for i in range(N)];solver.add(xs[348]);conflicts=[]
 for i,j in combinations(active,2):
  if i<493 and j<493:possible=f'{i},{j}'in oldgraph['edges']
  else:possible=pn[i,j]<=1e-6 and (ratio[i,j]>1+1e-7 or min(abs(ratio[i,j]-v)for v in allowed)<1e-6)
  if f'{i},{j}'in data['exact_pair_checks']:possible=data['exact_pair_checks'][f'{i},{j}']['compatible']
  if not possible:conflicts.append([i,j]);solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
 cs=[z3.Bool(f'c{i}')for i in range(len(data['escape_clauses']))]
 for c,row in zip(cs,data['escape_clauses']):solver.add(z3.Implies(c,z3.Or(*(xs[i]for i in row['complete_cutters_among_active_nodes']))))
 assert solver.check(*cs)==z3.unsat
 core=sorted(int(str(c)[1:])for c in solver.unsat_core())
 for j in list(core):
  trial=[i for i in core if i!=j]
  if solver.check(*(cs[i]for i in trial))==z3.unsat:core=trial
 union=sorted(set().union(*(set(data['escape_clauses'][i]['complete_cutters_among_active_nodes'])for i in core)))
 out={'status':'GUIDANCE_CORE_OF_NUMERICAL_PROPOSAL_GRAPH','clause_indices':core,'clauses':[data['escape_clauses'][i]for i in core],'core_nodes':union,'numeric_proposal_conflicts_within_core':[p for p in conflicts if p[0]in union and p[1]in union],'scope':'Adaptive search guidance only: omitted proposal-graph edges have not been exactly classified.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print({'core':core,'clause_sizes':[len(r['complete_cutters_among_active_nodes'])for r in out['clauses']],'union':len(union)},flush=True)
if __name__=='__main__':run()
