#!/usr/bin/env python3
"""Diversified compatible completion using batches of extremal escape clauses.

HalfspaceIntersection only proposes vertices. Every negative volume certificate
is reconstructed rationally and checked exactly; no numeric closure is a proof.
"""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import numpy as np
from scipy.spatial import HalfspaceIntersection,QhullError
import z3
from scripts.fresh_h8_threefree_verify_20261013 import (
 old_load,node,gate,pair,parse,enc,read,matrix,mul,add,div,scale,sign,dot,qsign,Z,O)
from scripts.fresh_h8_alternate_caps_20261013 import real

ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_h8_threefree_combined_caps_20261013.json'
PREVIOUS='data/enumerations/fresh_h8_threefree_combined_selection_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_cap_completion_20261013.json'

def load():
 caps=read(SOURCE);base=read(caps['base_source']);G=matrix(base['seed_gram']);I=matrix(base['seed_gram_inverse']);center=tuple(map(parse,base['strict_center_y']));oldgraph,_,nodes=old_load();nodes+=list(node(c,G)for c in caps['caps']);ad=read('results/fresh_h8_alternate_core_angles_independent_20261013.json');polys=[(r['label'],tuple(map(F,r['coefficients_descending'])),F(r['largest_root_interval'][0]),F(r['largest_root_interval'][1]))for r in ad['real_cyclotomic_polynomials']]
 cn=np.array([[real(a)+np.sqrt(real(n['d']))*real(b)for a,b in zip(n['a'],n['b'])]for n in nodes]);qn=np.array([[real(x)for x in row]for row in I]);hn=np.array([real(x)for x in center]);hn/=hn.sum()
 return G,I,center,nodes,polys,cn,qn,hn,oldgraph

def certify(y0,chosen,I,center,nodes,hn):
 for eps in (1e-3,1e-5,1e-7):
  target=y0+eps*hn
  for den in (10**9,10**12,10**15):
   y=tuple(scale(O,F(round(float(x)*den),den))for x in target)
   if any(sign(t)<=0 for t in y)or sign(pair(y,I,y))<=0:continue
   if any(qsign((dot(nodes[i]['a'],y),dot(nodes[i]['b'],y)),nodes[i]['d'])<=0 for i in chosen):continue
   cp=pair(y,I,center)
   if sign(cp)>0:
    lam=div(scale(cp,-2),pair(center,I,center));y=tuple(add(x,mul(lam,c))for x,c in zip(y,center))
   elif sign(cp)==0:continue
   q=pair(y,I,y);cp=pair(y,I,center);assert sign(q)>0 and sign(cp)<0
   cuts=[i for i,n in enumerate(nodes)if qsign((dot(n['a'],y),dot(n['b'],y)),n['d'])<0]
   assert not set(cuts)&set(chosen)
   return {'status':'EXACT_STRICT_SPACELIKE_ESCAPE','y':[enc(x)for x in y],'norm':enc(q),'center_pairing':enc(cp),'complete_dictionary_cutters':cuts}
 return None

def vertex_escapes(chosen,I,center,nodes,cn,qn,hn,limit=24):
 A=np.vstack((np.eye(10),cn[chosen]));half=np.c_[A[:,9,None]-A[:,:9],-A[:,9]];half/=np.linalg.norm(half[:,:9],axis=1)[:,None]
 try:vertices=HalfspaceIntersection(half,hn[:9],qhull_options='Qx').intersections
 except QhullError as e:return {'status':'NUMERICAL_VERTEX_ENUMERATION_FAILED','detail':str(e).splitlines()[0]},[]
 ys=np.c_[vertices,1-vertices.sum(axis=1)];qs=np.einsum('bi,ij,bj->b',ys,qn,ys);space=ys[qs>1e-8];norms=qs[qs>1e-8];cut_counts=(cn@space.T< -1e-8).sum(axis=0)if len(space)else np.array([]);order=sorted(range(len(space)),key=lambda i:(cut_counts[i],-norms[i]));rows=[];seen=set();attempted=0
 for j in order:
  # Use complete cutter sets to distinguish open ends; do not retain many
  # barely different rounded witnesses that impose the same condition.
  nk=tuple(np.flatnonzero(cn@space[j]< -1e-8))
  if nk in seen:continue
  attempted+=1;row=certify(space[j],chosen,I,center,nodes,hn)
  if not row:continue
  key=tuple(row['complete_dictionary_cutters']);seen.add(nk)
  if any(set(old['complete_dictionary_cutters'])<=set(key)for old in rows):continue
  rows=[old for old in rows if not set(key)<set(old['complete_dictionary_cutters'])];rows.append(row)
  if len(rows)>=limit or not key:break
 summary={'status':'EXACT_ESCAPES_FOUND'if rows else'UNRESOLVED_NO_CERTIFIED_ESCAPE','numeric_vertices':len(ys),'numeric_spacelike_vertices':len(space),'numeric_maximum_norm':float(qs.max()),'exact_escape_count':len(rows),'rationalization_attempts':attempted}
 return summary,rows

def run(rounds=12,optional_polar=False,saved_only=False,extra_clauses=None):
 G,I,center,nodes,polys,cn,qn,hn,oldgraph=load();previous=read(PREVIOUS);N=len(nodes);allowed=[float((lo+hi)/8)for _,_,lo,hi in polys];gn=np.array([[real(x)for x in r]for r in G]);p=cn@gn@cn.T;norms=np.array([real(n['q'])for n in nodes]);rat=p*p/(norms[:,None]*norms[None,:]);active=list(range(N))if optional_polar else previous['active_nodes'];neighbors=[set()for _ in nodes];pair_checks={};cases=[];clauses=[];solver=z3.Solver();xs=[z3.Bool(f'c{i}')for i in range(N)];rng=np.random.default_rng(20261013)
 for i,j in combinations(active,2):
  if i<493 and j<493:possible=f'{i},{j}'in oldgraph['edges']
  else:possible=p[i,j]<=1e-6 and(rat[i,j]>1+1e-7 or min(abs(rat[i,j]-v)for v in allowed)<1e-6)
  if possible:neighbors[i].add(j);neighbors[j].add(i)
  else:solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
 if not optional_polar:solver.add(xs[348])
 for i in set(range(N))-set(active):solver.add(z3.Not(xs[i]))
 def add_clause(row,source):
  cuts=row['complete_dictionary_cutters'];solver.add(z3.Or(*(xs[i]for i in cuts)));clauses.append({'source':source,**row})
 for j,row in enumerate(previous['escape_clauses']):
  r=row['certificate'];y=tuple(map(parse,r['y']));assert sign(pair(y,I,y))>0 and sign(pair(y,I,center))<0;cuts=[i for i,n in enumerate(nodes)if qsign((dot(n['a'],y),dot(n['b'],y)),n['d'])<0];add_clause({'status':r['status'],'y':r['y'],'norm':r['norm'],'center_pairing':r['center_pairing'],'complete_dictionary_cutters':cuts},f'inherited-{j}')
 if extra_clauses:
  for j,row in enumerate(read(extra_clauses)['escape_clauses'][len(clauses):]):add_clause({k:v for k,v in row.items()if k!='source'},f'extra-{j}')
 def exact_pairs(chosen):
  bad=[]
  for i,j in combinations(sorted(chosen),2):
   key=f'{i},{j}'
   if key not in pair_checks:pair_checks[key]=gate(nodes[i],nodes[j],polys)
   if not pair_checks[key]['compatible']:
    bad.append([i,j]);neighbors[i].discard(j);neighbors[j].discard(i);solver.add(z3.Or(z3.Not(xs[i]),z3.Not(xs[j])))
  return bad
 def save(status):
  out={'schema':'H9-cap-completion-multi-escape-v1','status':status,'new_cap_source':SOURCE,'inherited_selection_source':PREVIOUS,'node_count':N,'active_nodes':active,'polar_required':not optional_polar,'cases':cases,'exact_pair_checks':pair_checks,'escape_clauses':clauses,'scope':'Numerical graph and numeric vertex enumeration propose candidates only. Every tested configuration has exact pair checks; every archived escape and full-dictionary cutter list is exact. No absence-of-escape or incomplete-graph nonexistence claim.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 saved=previous['final_proposal_clause_gate']['chosen_indices'];bad=exact_pairs(saved)
 if bad:cases.append({'source':'saved-model','selected_node_indices':saved,'rejected_pairs':bad});print('Saved model has rejected pairs',bad,flush=True)
 else:
  summary,rows=vertex_escapes(saved,I,center,nodes,cn,qn,hn);cases.append({'source':'saved-model','selected_node_indices':saved,'facets':10+len(saved),'vertex_gate':summary,'escape_clause_indices':list(range(len(clauses),len(clauses)+len(rows)))})
  for row in rows:add_clause(row,'saved-model')
  print('Saved model',summary,flush=True)
 save('SAVED_MODEL_CHECKED')
 if saved_only:return
 for iteration in range(rounds):
  if solver.check()==z3.unsat:save('UNSAT_FOR_NUMERICAL_PROPOSAL_GRAPH');print('UNSAT_FOR_NUMERICAL_PROPOSAL_GRAPH',flush=True);return
  # Seek redundant coverage of the hardest open ends, then select a maximal
  # clique favoring coverage rarity instead of only graph degree.
  order=sorted(range(len(clauses)),key=lambda j:len(set(clauses[j]['complete_dictionary_cutters'])&set(active)));soft=[]
  for j in order[:12]:
   ids=[i for i in clauses[j]['complete_dictionary_cutters']if i in active]
   if len(ids)>1:soft.append(z3.PbGe([(xs[i],1)for i in ids],2))
  solver.push()
  for c in soft:
   solver.push();solver.add(c)
   if solver.check()!=z3.sat:solver.pop()
  assert solver.check()==z3.sat;model=solver.model();chosen=[i for i in active if z3.is_true(model.eval(xs[i],model_completion=True))]
  # Restore all trial coverage preferences; learned escape clauses stay hard.
  while solver.num_scopes():solver.pop()
  pool=set(active)-set(chosen)
  for i in chosen:pool&=neighbors[i]
  weights=np.zeros(N)
  for row in clauses:
   ids=[i for i in row['complete_dictionary_cutters']if i in active]
   if ids:weights[ids]+=1/len(ids)
  noise=rng.random(N)
  while pool:
   i=max(pool,key=lambda j:(weights[j],len(neighbors[j]&pool),noise[j]));chosen.append(i);pool&=neighbors[i]
  chosen.sort();bad=exact_pairs(chosen)
  if bad:print({'iteration':iteration,'rejected_pairs':bad},flush=True);save('SEARCHING');continue
  summary,rows=vertex_escapes(chosen,I,center,nodes,cn,qn,hn);cases.append({'source':f'batch-{iteration}','selected_node_indices':chosen,'facets':10+len(chosen),'vertex_gate':summary,'escape_clause_indices':list(range(len(clauses),len(clauses)+len(rows)))})
  for row in rows:add_clause(row,f'batch-{iteration}')
  print({'iteration':iteration,'facets':10+len(chosen),**summary,'new_cutter_counts':[len(row['complete_dictionary_cutters'])for row in rows]},flush=True)
  if not rows:save('UNRESOLVED_COMPATIBLE_CONFIGURATION');return
  save('SEARCHING')
 save('BOUNDED_MULTI_ESCAPE_SEARCH_INCOMPLETE')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--rounds',type=int,default=12);p.add_argument('--optional-polar',action='store_true');p.add_argument('--saved-only',action='store_true');p.add_argument('--output');p.add_argument('--extra-clauses');a=p.parse_args()
 if a.output:OUTPUT=ROOT/a.output
 run(a.rounds,a.optional_polar,a.saved_only,a.extra_clauses)
