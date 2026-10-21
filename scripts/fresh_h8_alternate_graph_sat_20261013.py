#!/usr/bin/env python3
"""Exact accepted-edge graph SAT search; absent edges are a search restriction."""
import json
from itertools import combinations
import z3
from scripts.fresh_h8_alternate_selection_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import esign
OUTPUT=ROOT/'data/enumerations/fresh_h8_alternate_graph_sat_20261013.json'

def run():
 base,pol,raw,G,I,center,nodes=load();selection=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_selection_20261013.json').read_text());previous=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_sat_20261013.json').read_text());N=len(nodes);neighbors=[set()for _ in nodes]
 for key in selection['edges']:
  i,j=map(int,key.split(','));neighbors[i].add(j);neighbors[j].add(i)
 solver=z3.Solver();xs=[z3.Bool(f'w{i}')for i in range(N)]
 for i,j in combinations(range(N),2):
  if j not in neighbors[i]:solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
 rows=list(previous['escape_clauses']);checks=[]
 for r in rows:solver.add(z3.Or(*(xs[i]for i in r['cutters'])))
 def save(status):
  out={'schema':'alternate-H8-accepted-edge-SAT-v1','status':status,'graph_source':'data/enumerations/fresh_h8_alternate_selection_20261013.json','inherited_escapes':'data/enumerations/fresh_h8_alternate_sat_20261013.json','new_escape_clauses':rows[len(previous['escape_clauses']):],'checks':checks,'scope':'Search restricted to the exact accepted-edge graph. An UNSAT result applies to that graph; omitted angles and numeric proposal-filter omissions are not ruled out without a separate graph-completeness certificate.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 for iteration in range(20):
  status=solver.check()
  if status==z3.unsat:save('UNSAT_FOR_ACCEPTED_EDGE_GRAPH');print('UNSAT_FOR_ACCEPTED_EDGE_GRAPH',flush=True);return
  assert status==z3.sat;model=solver.model();chosen=[i for i,x in enumerate(xs)if z3.is_true(model.eval(x,model_completion=True))];pool=set(range(N))-set(chosen)
  for i in chosen:pool&=neighbors[i]
  while pool:
   i=max(pool,key=lambda j:len(neighbors[j]&pool));chosen.append(i);pool&=neighbors[i]
  chosen.sort();old=[nodes[i]['a']for i in chosen if nodes[i]['source']!='coupled'];caps=[raw['caps'][nodes[i]['index']]for i in chosen if nodes[i]['source']=='coupled'];gate=escape(I,old,caps,center,tries=1000);checks.append({'selected_indices':chosen,'volume_gate':gate});print({'case':iteration,'facets':10+len(chosen),'gate':gate['status']},flush=True)
  if gate['status']!='EXACT_STRICT_SPACELIKE_ESCAPE':save('UNRESOLVED_COMPATIBLE_SUBSET');return
  y=tuple(map(parse,gate['y']));cutters=[i for i,n in enumerate(nodes)if esign((dot(n['a'],y),dot(n['b'],y)),n['d'])<0];rows.append({'source':f'graph-sat-{iteration}','certificate':gate,'cutters':cutters});solver.add(z3.Or(*(xs[i]for i in cutters)));save('SEARCHING')
 save('BOUNDED_SEARCH_INCOMPLETE')
if __name__=='__main__':run()
