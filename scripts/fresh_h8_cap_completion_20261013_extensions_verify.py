#!/usr/bin/env python3
"""Exact independent replay of the49 targeted caps and all coherent extensions."""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from scripts.fresh_h8_threefree_verify_20261013 import (
 run as replay_caps,old_load,node,gate,read,pair,parse,matrix,sign,dot,qsign)
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'data/enumerations/fresh_h8_cap_completion_20261013_extensions_verify.json'
def run(verify=False):
 source='data/enumerations/fresh_h8_cap_completion_20261013_bottleneck_caps.json';replay_caps(source,None,'data/enumerations/fresh_h8_cap_completion_20261013_bottleneck_caps_verify.json',True)
 oldcaps=read('data/enumerations/fresh_h8_threefree_combined_caps_20261013.json');base=read(oldcaps['base_source']);G=matrix(base['seed_gram']);I=matrix(base['seed_gram_inverse']);center=tuple(map(parse,base['strict_center_y']));_,_,nodes=old_load();nodes+=list(node(c,G)for c in oldcaps['caps']);assert len(nodes)==1035;nodes+=list(node(c,G)for c in read(source)['caps']);assert len(nodes)==1084
 ad=read('results/fresh_h8_alternate_core_angles_independent_20261013.json');polys=[(r['label'],tuple(map(F,r['coefficients_descending'])),F(r['largest_root_interval'][0]),F(r['largest_root_interval'][1]))for r in ad['real_cyclotomic_polynomials']];data=read('data/enumerations/fresh_h8_cap_completion_20261013_extensions.json');selectedbase=data['base_selected_nodes'];cache={}
 def compatible(i,j):
  key=tuple(sorted((i,j)))
  if key not in cache:cache[key]=gate(nodes[key[0]],nodes[key[1]],polys)
  return cache[key]['compatible']
 eligible=[i for i in range(49)if all(compatible(j,1035+i)for j in selectedbase)];assert eligible==data['eligible_new_cap_indices']
 cliques=[]
 for mask in range(1<<len(eligible)):
  ids=[eligible[j]for j in range(len(eligible))if(mask>>j)&1]
  if not all(compatible(1035+i,1035+j)for i,j in combinations(ids,2)):continue
  if any(all(compatible(1035+i,1035+j)for j in ids)for i in eligible if i not in ids):continue
  cliques.append(ids)
 assert sorted(cliques)==sorted(data['maximal_extension_cliques'])
 for row in data['exact_pair_checks']:
  i,j=row['i'],row['j'];compatible(i,j);assert cache[tuple(sorted((i,j)))]=={k:row[k]for k in('compatible','kind','label')}
 # Full1035 cutter completeness at the target proves that all49 new caps,
 # which strictly cut it, are outside the previously certified dictionary.
 target=read('data/enumerations/fresh_h8_cap_completion_20261013_core.json')['clauses'][0];ty=tuple(map(parse,target['y']));oldcuts=[i for i,n in enumerate(nodes[:1035])if qsign((dot(n['a'],ty),dot(n['b'],ty)),n['d'])<0];assert oldcuts==target['complete_dictionary_cutters']
 assert all(qsign((dot(n['a'],ty),dot(n['b'],ty)),n['d'])<0 for n in nodes[1035:])
 # None of the13 old cutters is polar compatible; every new cap is.
 assert all(not compatible(348,i)for i in oldcuts)
 assert all(compatible(348,i)for i in range(1035,1084))
 escapes=0
 for case in data['cases']:
  chosen=case['selected_node_indices'];assert chosen==sorted(selectedbase+[1035+i for i in case['added_cap_indices']]);assert case['added_cap_indices']in cliques
  assert all(compatible(i,j)for i,j in combinations(chosen,2))
  for row in case['escapes']:
   y=tuple(map(parse,row['y']));assert all(sign(x)>0 for x in y);q=pair(y,I,y);cp=pair(y,I,center);assert sign(q)>0 and sign(cp)<0;assert q==parse(row['norm'])and cp==parse(row['center_pairing'])
   assert all(qsign((dot(nodes[i]['a'],y),dot(nodes[i]['b'],y)),nodes[i]['d'])>0 for i in chosen)
   cuts=[i for i,n in enumerate(nodes)if qsign((dot(n['a'],y),dot(n['b'],y)),n['d'])<0];assert cuts==row['complete_dictionary_cutters'];escapes+=1
 assert len(data['cases'])==len(cliques)
 out={'status':'EXACT_TARGETED_CAP_AND_ALL_MAXIMAL_EXTENSION_REPLAY_PASS','new_caps':49,'all_new_caps_outside_old1035':True,'base_facet_count':10+len(selectedbase),'eligible_new_cap_indices':eligible,'maximal_extension_cliques':cliques,'exact_pair_checks':len(cache),'exact_escaped_extensions':len(cliques),'exact_escapes_and_complete1084_cutter_lists':escapes,'scope':'All maximal pairwise compact-compatible extensions from these49 caps that preserve every facet of the fixed30-facet configuration have exact strict escapes. This does not exclude dropping/replacing base facets or unrestricted new caps.'}
 if verify:assert read(str(OUTPUT.relative_to(ROOT)))==out
 else:OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
