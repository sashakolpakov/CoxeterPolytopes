#!/usr/bin/env python3
"""Coherent replacements forced by an uncovered gap, then multi-gap completion."""
import json
from itertools import combinations
import numpy as np
import z3
from scripts.fresh_h9_joint_completion_20261013 import setup,ROOT,read
from scripts.fresh_h8_cap_completion_20261013 import node,gate,vertex_escapes
from scripts.fresh_h8_alternate_caps_20261013 import parse,real,dot
from scripts.fresh_h8_release_h9_free_caps_20261013 import esign
OUTPUT=ROOT/'data/enumerations/fresh_h9_joint_completion_20261013_replacements.json'

def run():
 G,I,center,old,polys,cn,qn,hn,original=setup();data=read('data/enumerations/fresh_h9_joint_completion_20261013.json');mirrors=read('data/enumerations/fresh_h9_joint_completion_20261013_mirrors.json');targets=read(data['target_source']);ys=[tuple(map(parse,r['certificate']['y']))for r in targets['targets']];new=[node(c,G)for c in data['caps']];nodes=old+new;offset=len(old);cn=np.array([[real(a)+np.sqrt(real(n['d']))*real(b)for a,b in zip(n['a'],n['b'])]for n in nodes]);scenarios=[];cases=[];crosschecks={};pairchecks={}
 def save():
  out={'schema':'joint-gap-forced-coherent-replacements-v1','source':'data/enumerations/fresh_h9_joint_completion_20261013.json','target_source':data['target_source'],'old_dictionary_size':offset,'original_configuration':original,'candidate_representatives':mirrors['new_representative_indices'],'scenarios':scenarios,'exact_cross_checks':crosschecks,'exact_new_pair_checks':pairchecks,'cases':cases,'scope':'Each scenario deliberately replaces every original cap incompatible with one exact cutter of gap8. Surviving joint caps must preserve the full remaining coherent configuration. Exact angle gates and strict escapes support only the archived configurations; unfinished solver or numeric vertex outcomes are not no-go claims.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 def compatible(i,j):
  key=f'{min(i,j)},{max(i,j)}'
  if key not in pairchecks:pairchecks[key]=gate(nodes[i],nodes[j],polys)
  return pairchecks[key]['compatible']
 choices=sorted((r for r in read('data/enumerations/fresh_h9_joint_completion_20261013_gap8.json')['checks']if len(r['dropped'])<=6),key=lambda r:len(r['dropped']))
 for choice in choices:
  cutter=choice['node'];kept=choice['kept'];fixed=sorted(kept+[cutter]);eligible=[]
  for j in mirrors['new_representative_indices']:
   if not set(kept)<=set(data['caps'][j]['exact_retained_nodes']):continue
   row=gate(new[j],old[cutter],polys);crosschecks[f'{j},{cutter}']=row
   if row['compatible']:eligible.append(j)
  already=[t for t,y in enumerate(ys)if esign((dot(old[cutter]['a'],y),dot(old[cutter]['b'],y)),old[cutter]['d'])<0];clauses=[{'target':t,'cutters':[j for j in eligible if t in data['caps'][j]['cut_target_indices']]}for t in range(len(ys))if t not in already];scenario={'cutter':cutter,'dropped_original_caps':choice['dropped'],'fixed_cap_nodes':fixed,'eligible_new_cap_indices':eligible,'targets_cut_by_cutter':already,'necessary_cover_clauses':clauses};scenarios.append(scenario)
  empty=[r['target']for r in clauses if not r['cutters']]
  if empty:scenario['status']='EXACT_UNCOVERED_TARGET';scenario['uncut_target_indices']=empty;save();print({'cutter':cutter,'dropped':len(choice['dropped']),'eligible':len(eligible),'uncut_targets':empty},flush=True);continue
  solver=z3.Solver();xs={j:z3.Bool(f's{cutter}_{j}')for j in eligible};adj={j:set()for j in eligible}
  for i,j in combinations(eligible,2):
   if compatible(offset+i,offset+j):adj[i].add(j);adj[j].add(i)
   else:solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
  for row in clauses:solver.add(z3.Or(*(xs[j]for j in row['cutters'])))
  scenario['status']='READY';save();print({'cutter':cutter,'dropped':len(choice['dropped']),'eligible':len(eligible),'clause_sizes':[len(r['cutters'])for r in clauses]},flush=True)
  for step in range(4):
   status=solver.check()
   if status==z3.unsat:scenario['status']='EXACT_COMPATIBILITY_COVER_UNSAT';break
   assert status==z3.sat;model=solver.model();selected=[j for j in eligible if z3.is_true(model.eval(xs[j],model_completion=True))];pool=set(eligible)-set(selected)
   for j in selected:pool&=adj[j]
   while pool:
    j=max(pool,key=lambda j:(len(data['caps'][j]['cut_target_indices']),len(adj[j]&pool)));selected.append(j);pool&=adj[j]
   chosen=sorted(fixed+[offset+j for j in selected]);summary,rows=vertex_escapes(chosen,I,center,nodes,cn,qn,hn,limit=8);case={'scenario_cutter':cutter,'step':step,'selected_node_indices':chosen,'new_cap_indices':selected,'facets':10+len(chosen),'vertex_gate':summary,'escapes':rows};cases.append(case);print({'cutter':cutter,'step':step,'facets':10+len(chosen),**summary},flush=True)
   if not rows:scenario['status']='UNRESOLVED_VOLUME_CANDIDATE';save();return
   for r in rows:
    cuts=[j for j in eligible if offset+j in r['complete_dictionary_cutters']];solver.add(z3.Or(*(xs[j]for j in cuts)))
   scenario['status']='BOUNDED_EXACT_NEGATIVE_CASES';save()
  save()
 print({'scenarios':len(scenarios),'volume_cases':len(cases)},flush=True)
if __name__=='__main__':run()
