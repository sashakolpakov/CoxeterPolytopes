#!/usr/bin/env python3
"""Counterexample-guided subset search with exact positive-pair conflicts.

The SAT constraints are necessary and deliberately overinclusive: only
positive outward pairings are initially forbidden. Thus an UNSAT result
is stronger than exclusion for any chosen finite-angle alphabet.
"""
import json
from itertools import combinations
import z3
from scripts.fresh_h8_alternate_selection_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import esign
OUTPUT=ROOT/'data/enumerations/fresh_h8_alternate_sat_20261013.json'

def run():
 base,pol,raw,G,I,center,nodes=load();N=len(nodes);solver=z3.Solver();xs=[z3.Bool(f'w{i}')for i in range(N)];conflicts={};escapes=[];clauses=set();checks=[]
 def add_escape(row,source):
  y=tuple(map(parse,row['y']));assert sign(pair(y,I,y))>0 and sign(pair(y,I,center))<0 and all(sign(x)>0 for x in y)
  cutters=[i for i,n in enumerate(nodes)if esign((dot(n['a'],y),dot(n['b'],y)),n['d'])<0]
  if tuple(cutters)in clauses:return
  clauses.add(tuple(cutters));solver.add(z3.Or(*(xs[i]for i in cutters)));escapes.append({'source':source,'certificate':row,'cutters':cutters});print({'escape':len(escapes),'cutters':len(cutters)},flush=True)
 add_escape(base['volume_gate'],'base');add_escape(pol['volume_gate'],'base-plus-polar')
 for j,r in enumerate(raw['rounds']):
  if r['volume_gate']['status']=='EXACT_STRICT_SPACELIKE_ESCAPE':add_escape(r['volume_gate'],f'coupled-round-{j}')
 selection=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_selection_20261013.json').read_text())
 for j,c in enumerate(selection['cases']):add_escape(c['volume_gate'],f'greedy-{j}')
 def save(status,chosen=None):
  out={'schema':'alternate-H8-lazy-positive-conflict-SAT-v1','status':status,'node_count':N,'source_dictionary':'data/enumerations/fresh_h8_alternate_selection_20261013.json','positive_pair_conflicts':[{'indices':list(k),'pairing':vec(v)}for k,v in conflicts.items()],'escape_clauses':escapes,'checks':checks,'last_selected_indices':chosen,'scope':'Necessary positive-pair incompatibilities plus exact spacelike escape covering clauses. UNSAT would exclude all Coxeter subsets of this finite dictionary, independent of angle alphabet; satisfiable or undecided is not a polytope certificate.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 for iteration in range(200):
  status=solver.check()
  if status==z3.unsat:save('UNSAT');print('UNSAT',flush=True);return
  assert status==z3.sat;model=solver.model();chosen=[i for i,x in enumerate(xs)if z3.is_true(model.eval(x,model_completion=True))]
  added=0
  for i,j in combinations(chosen,2):
   p=pairing(nodes[i],nodes[j])
   if bsign(p,nodes[i]['d'],nodes[j]['d'])>0:
    assert(i,j)not in conflicts;conflicts[i,j]=p;solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])));added+=1
  if added:print({'iteration':iteration,'selected':len(chosen),'new_positive_conflicts':added},flush=True);continue
  old=[nodes[i]['a']for i in chosen if nodes[i]['source']!='coupled'];caps=[raw['caps'][nodes[i]['index']]for i in chosen if nodes[i]['source']=='coupled'];gate=escape(I,old,caps,center,tries=500);checks.append({'selected_indices':chosen,'volume_gate':gate});print({'iteration':iteration,'facets':10+len(chosen),'gate':gate['status']},flush=True)
  if gate['status']!='EXACT_STRICT_SPACELIKE_ESCAPE':save('UNRESOLVED_NONPOSITIVE_SUBSET',chosen);return
  add_escape(gate,f'sat-{len(checks)-1}');save('SEARCHING',chosen)
  if len(checks)>=20:save('BOUNDED_SEARCH_INCOMPLETE',chosen);return
 save('BOUNDED_SEARCH_INCOMPLETE',chosen)
if __name__=='__main__':run()
