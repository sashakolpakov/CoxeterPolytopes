#!/usr/bin/env python3
"""Compatible subset search with exact pair gates and witness-driven SAT.

A numerical angle graph proposes cliques; every tested subset is replayed
exactly before its volume gate. An UNSAT result for the proposal graph is
not an exhaustive angle/dictionary theorem.
"""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
import numpy as np
import z3
from scripts.fresh_h8_alternate_selection_20261013 import load,pairing,bsign,bmul
from scripts.fresh_h8_alternate_caps_20261013 import *
from scripts.fresh_h8_alternate_coupled_20261013 import escape
from scripts.fresh_h8_release_h9_free_caps_20261013 import esign
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import finite_label
SOURCE=ROOT/'data/enumerations/fresh_h8_threefree_caps_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_threefree_selection_20261013.json'

def run(max_cases=30,resume=False,inherited=None,check_only=False):
 base,pol,raw,G,I,center,nodes=load();oldgraph=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_selection_20261013.json').read_text());new=json.loads(SOURCE.read_text());core=json.loads((ROOT/'results/fresh_h8_alternate_core_cover_independent_20261013.json').read_text());angledata=json.loads((ROOT/'results/fresh_h8_alternate_core_angles_independent_20261013.json').read_text());polys=[(r['label'],tuple(map(F,r['coefficients_descending'])),F(r['largest_root_interval'][0]),F(r['largest_root_interval'][1]))for r in angledata['real_cyclotomic_polynomials']]
 for j,c in enumerate(new['caps']):nodes.append({'source':'threefree','index':j,'d':parse(c['discriminant']),'a':tuple(map(parse,c['coefficients_rational'])),'b':tuple(map(parse,c['coefficients_radical'])),'pa':tuple(map(parse,c['pairings_rational'])),'pb':tuple(map(parse,c['pairings_radical'])),'q':scale(O,2)})
 N=len(nodes);polar=348;old_neighbors=set()
 for key in oldgraph['edges']:
  i,j=map(int,key.split(','))
  if i==polar:old_neighbors.add(j)
  if j==polar:old_neighbors.add(i)
 active=sorted(old_neighbors|{polar}|set(range(493,N)));nn=[set()for _ in nodes];cn=np.array([[real(x)+real(r['d'])**.5*real(y)for x,y in zip(r['a'],r['b'])]for r in nodes]);gn=np.array([[real(x)for x in r]for r in G]);qs=np.array([real(r['q'])for r in nodes]);pn=cn@gn@cn.T;ratio=pn**2/(qs[:,None]*qs[None,:]);allowed=[float((lo+hi)/8)for _,_,lo,hi in polys]
 for i,j in combinations(active,2):
  if i<493 and j<493:possible=f'{i},{j}'in oldgraph['edges']
  else:possible=pn[i,j]<=1e-6 and (ratio[i,j]>1+1e-7 or min(abs(ratio[i,j]-v)for v in allowed)<1e-6)
  if possible:nn[i].add(j);nn[j].add(i)
 solver=z3.Solver();xs=[z3.Bool(f'w{i}')for i in range(N)];solver.add(xs[polar])
 for i in set(range(N))-set(active):solver.add(z3.Not(xs[i]))
 for i,j in combinations(active,2):
  if j not in nn[i]:solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
 clauses=[];pair_checks={};cases=[];rng=np.random.default_rng(20261013)
 def add_escape(row,source):
  y=tuple(map(parse,row['y']));assert sign(pair(y,I,y))>0 and sign(pair(y,I,center))<0
  cutters=[i for i in active if esign((dot(nodes[i]['a'],y),dot(nodes[i]['b'],y)),nodes[i]['d'])<0]
  solver.add(z3.Or(*(xs[i]for i in cutters)));clauses.append({'source':source,'certificate':row,'complete_cutters_among_active_nodes':cutters});return len(cutters)
 for j,row in enumerate(core['five_complete_escape_clauses']):add_escape({'status':'EXACT_STRICT_SPACELIKE_ESCAPE','y':row['escape_y'],'norm':row['escape_norm'],'center_pairing':row['escape_center_pairing']},f'old-core-{j}')
 if inherited:
  previous=json.loads((ROOT/inherited).read_text())
  for row in previous['escape_clauses'][5:]:add_escape(row['certificate'],row['source'])
 if resume and OUTPUT.exists():
  previous=json.loads(OUTPUT.read_text());assert previous['node_count']==N and previous['active_nodes']==active
  pair_checks.update(previous['exact_pair_checks']);cases.extend(previous['cases'])
  for row in previous['escape_clauses'][5:]:add_escape(row['certificate'],row['source'])
  for key,row in pair_checks.items():
   if not row['compatible']:
    i,j=map(int,key.split(','));nn[i].discard(j);nn[j].discard(i);solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
 def pairgate(i,j):
  key=f'{i},{j}'
  if key in pair_checks:return pair_checks[key]['compatible']
  a,b=nodes[i],nodes[j];d,e=a['d'],b['d'];p=pairing(a,b);ps=bsign(p,d,e);label=None
  if ps>0:kind='positive';ok=False
  else:
   rat=tuple(mul(div(scale(O,4),mul(a['q'],b['q'])),x)for x in bmul(p,p,d,e));cmp=bsign((sub(rat[0],scale(O,4)),*rat[1:]),d,e)
   if cmp>0:kind='dotted';ok=True
   elif cmp==0:kind='parallel';ok=False
   else:label=finite_label(rat,d,e,polys);ok=label is not None;kind='finite'if ok else'noncoxeter'
  pair_checks[key]={'compatible':ok,'kind':kind,'label':label};return ok
 def save(status,last=None):
  out={'schema':'three-free-compatible-subset-search-v1','status':status,'new_cap_source':str(SOURCE.relative_to(ROOT)),'old_dictionary_source':'data/enumerations/fresh_h8_alternate_selection_20261013.json','inherited_escape_source':inherited,'node_count':N,'active_nodes':active,'forced_polar_node':polar,'exact_pair_checks':pair_checks,'escape_clauses':clauses,'cases':cases,'last_proposed_indices':last,'scope':'Numerical angle graph proposes cliques. Every archived tested subset has exact pair compatibility and strict escape verification. UNSAT of this proposal graph is not a complete angle or unrestricted dictionary obstruction.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 print({'all_nodes':N,'active_nodes':len(active),'new_caps':len(new['caps']),'initial_clause_sizes':[len(r['complete_cutters_among_active_nodes'])for r in clauses]},flush=True)
 if check_only:
  assert resume;status=solver.check();assert status in(z3.sat,z3.unsat);certificate={'status':str(status)}
  if status==z3.sat:
   model=solver.model();chosen=[i for i in active if z3.is_true(model.eval(xs[i],model_completion=True))];assert all(j in nn[i]for i,j in combinations(chosen,2));assert all(set(chosen)&set(r['complete_cutters_among_active_nodes'])for r in clauses);certificate['chosen_indices']=chosen
  previous=json.loads(OUTPUT.read_text());previous['final_proposal_clause_gate']=certificate;OUTPUT.write_text(json.dumps(previous,indent=2)+'\n');print(certificate,flush=True);return
 for iteration in range(200):
  status=solver.check()
  if status==z3.unsat:save('UNSAT_FOR_NUMERIC_PROPOSAL_GRAPH');print('UNSAT_FOR_NUMERIC_PROPOSAL_GRAPH',flush=True);return
  assert status==z3.sat;model=solver.model();chosen=[i for i in active if z3.is_true(model.eval(xs[i],model_completion=True))];pool=set(active)-set(chosen)
  for i in chosen:pool&=nn[i]
  noise=rng.random(N)
  while pool:
   i=max(pool,key=lambda j:(len(nn[j]&pool),noise[j]));chosen.append(i);pool&=nn[i]
  chosen.sort();bad=[]
  for i,j in combinations(chosen,2):
   if not pairgate(i,j):bad.append((i,j))
  if bad:
   for i,j in bad:nn[i].discard(j);nn[j].discard(i);solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
   print({'iteration':iteration,'exact_rejected_pairs':len(bad)},flush=True);save('SEARCHING',chosen);continue
  fixed=[nodes[i]['a']for i in chosen if nodes[i]['d']==Z];fields=[]
  for i in chosen:
   n=nodes[i]
   if n['d']==Z:continue
   fields.append(raw['caps'][n['index']]if n['source']=='coupled'else new['caps'][n['index']])
  gate=escape(I,fixed,fields,center,tries=1000);row={'case':len(cases),'selected_node_indices':chosen,'total_facets':10+len(chosen),'new_caps':sum(i>=493 for i in chosen),'volume_gate':gate};cases.append(row);print({k:v for k,v in row.items()if k not in('selected_node_indices','volume_gate')},gate['status'],flush=True)
  if gate['status']!='EXACT_STRICT_SPACELIKE_ESCAPE':save('UNRESOLVED_EXACT_COMPATIBLE_SUBSET',chosen);return
  nc=add_escape(gate,f'new-subset-{len(cases)-1}');print({'new_clause_cutters':nc},flush=True);save('SEARCHING',chosen)
  if len(cases)>=max_cases:save('BOUNDED_SEARCH_INCOMPLETE',chosen);return
 save('BOUNDED_SEARCH_INCOMPLETE',chosen)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cases',type=int,default=30);p.add_argument('--resume',action='store_true');p.add_argument('--source');p.add_argument('--output');p.add_argument('--inherited');p.add_argument('--check-only',action='store_true');args=p.parse_args()
 if args.source:SOURCE=ROOT/args.source
 if args.output:OUTPUT=ROOT/args.output
 run(args.cases,args.resume,args.inherited,args.check_only)
