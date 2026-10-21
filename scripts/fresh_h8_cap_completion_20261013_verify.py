#!/usr/bin/env python3
"""Producer-free exact replay of the finite dictionary's compact obstruction.

No NumPy, LP, floating point, Qhull or SAT solver is used. Exact field and
cyclotomic arithmetic and a Boolean proof-tree replay come from prior audits.
"""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from scripts.fresh_h8_threefree_verify_20261013 import (
 old_load,node,gate,read,pair,parse,enc,matrix,sign,dot,qsign,add,sub,mul,Z)
from scripts.fresh_h8_alternate_core_cover_independent_20261013 import dpll
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'data/enumerations/fresh_h8_cap_completion_20261013_verify.json'
def run(verify=False):
 source='data/enumerations/fresh_h8_threefree_combined_caps_20261013.json';caps=read(source);base=read(caps['base_source']);G=matrix(base['seed_gram']);I=matrix(base['seed_gram_inverse']);center=tuple(map(parse,base['strict_center_y']));_,_,nodes=old_load();nodes+=list(node(c,G)for c in caps['caps']);assert len(nodes)==1035
 assert sign(pair(center,I,center))<0 and all(sign(x)>0 for x in center)
 assert all(qsign((dot(n['a'],center),dot(n['b'],center)),n['d'])>0 for n in nodes)
 for n in nodes:
  assert sign(n['q'])>0
  square=(add(dot(n['pa'],n['a']),mul(n['d'],dot(n['pb'],n['b']))),add(dot(n['pa'],n['b']),dot(n['pb'],n['a'])))
  assert qsign((sub(square[0],n['q']),square[1]),n['d'])==0
  assert all(qsign((a,b),n['d'])<=0 for a,b in zip(n['pa'],n['pb']))
 ad=read('results/fresh_h8_alternate_core_angles_independent_20261013.json');polys=[(r['label'],tuple(map(F,r['coefficients_descending'])),F(r['largest_root_interval'][0]),F(r['largest_root_interval'][1]))for r in ad['real_cyclotomic_polynomials']]
 core=read('data/enumerations/fresh_h8_cap_completion_20261013_optional_core.json');raw=read(core['source']);assert not core['polar_required'];assert core['clauses']==[raw['escape_clauses'][i]for i in core['clause_indices']]
 def check_escape(row,chosen=None):
  y=tuple(map(parse,row['y']));q=pair(y,I,y);cp=pair(y,I,center);assert all(sign(x)>0 for x in y)and sign(q)>0 and sign(cp)<0;assert q==parse(row['norm'])and cp==parse(row['center_pairing'])
  if chosen is not None:assert all(qsign((dot(nodes[i]['a'],y),dot(nodes[i]['b'],y)),nodes[i]['d'])>0 for i in chosen)
  cuts=[i for i,n in enumerate(nodes)if qsign((dot(n['a'],y),dot(n['b'],y)),n['d'])<0];assert cuts==row['complete_dictionary_cutters'];return cuts
 clauses=[check_escape(row)for row in core['clauses']];union=sorted(set().union(*map(set,clauses)));assert union==core['core_nodes'];local={n:i for i,n in enumerate(union)};adj=[0]*len(union);counts={};assert {(r['i'],r['j'])for r in core['exact_pair_checks']}==set(combinations(union,2))
 for row in core['exact_pair_checks']:
  i,j=row['i'],row['j'];got=gate(nodes[i],nodes[j],polys);assert got=={k:row[k]for k in('compatible','kind','label')};counts[got['kind']]=counts.get(got['kind'],0)+1
  if got['compatible']:
   a,b=local[i],local[j];adj[a]|=1<<b;adj[b]|=1<<a
 assert counts==core['exact_pair_kind_counts'];masks=[sum(1<<local[i]for i in row)for row in clauses];proof=dpll(masks,adj,False)
 assert core['compact_core_status']=='unsat'
 # Including parallel pairs leaves an exact core solution, so the theorem
 # must not be upgraded from compact to finite-volume completion.
 fv=core['surviving_exact_core_indices'];assert all(set(fv)&set(c)for c in clauses)
 for i,j in combinations(fv,2):
  r=gate(nodes[i],nodes[j],polys);assert r['compatible']or r['kind']=='parallel'
 checked_cases=0;new_escape_count=0
 search=read('data/enumerations/fresh_h8_cap_completion_20261013.json')
 for case in search['cases']:
  chosen=case['selected_node_indices']
  for i,j in combinations(chosen,2):assert gate(nodes[i],nodes[j],polys)['compatible']
  for ci in case['escape_clause_indices']:check_escape(search['escape_clauses'][ci],chosen);new_escape_count+=1
  checked_cases+=1
 polar=read('data/enumerations/fresh_h8_cap_completion_20261013_core.json');assert len(polar['clauses'])==1;cuts=check_escape(polar['clauses'][0]);assert len(cuts)==13
 for i in cuts:
  r=gate(nodes[348],nodes[i],polys);assert not r['compatible']and r['kind']!='parallel'
 para=read('data/enumerations/fresh_h8_cap_completion_20261013_parabolic_candidate.json');chosen=para['selected_node_indices'];assert chosen==fv
 for i,j in combinations(chosen,2):
  r=gate(nodes[i],nodes[j],polys);assert r['compatible']or r['kind']=='parallel'
 for r in para['escapes']:check_escape(r,chosen)
 result={'status':'EXACT_FINITE_DICTIONARY_COMPACT_OBSTRUCTION','dictionary_nodes':len(nodes),'core_clause_indices':core['clause_indices'],'core_nodes':union,'core_cutter_counts':list(map(len,clauses)),'complete_pair_counts':counts,'common_center_y':base['strict_center_y'],'common_center_norm':enc(pair(center,I,center)),'boolean_compact_unsat_proof':proof,'finite_volume_core_solution':fv,'selected_configurations_checked':checked_cases,'new_exact_escapes_checked':new_escape_count,'parallel_candidate_exact_escapes_checked':len(para['escapes']),'polar_facet_bottleneck_cutters':cuts,'scope':'No subset of this explicit1035-entry dictionary completes the ten fixed seed walls to a compact Coxeter polytope; polar optional. Parallel edges leave a core solution, so this is not a finite-volume dictionary theorem. All negative volume claims use strict exact future spacelike witnesses.'}
 if verify:assert read(str(OUTPUT.relative_to(ROOT)))==result
 else:OUTPUT.write_text(json.dumps(result,indent=2)+'\n')
 print({k:v for k,v in result.items()if k not in('core_nodes','common_center_y','boolean_compact_unsat_proof')},flush=True);print({'boolean_proof_states':proof['distinct_states']},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
