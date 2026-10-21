#!/usr/bin/env python3
"""Test all maximal compatible extensions of one coherent30-facet configuration."""
import json
from itertools import combinations
import numpy as np
from scripts.fresh_h8_cap_completion_20261013 import (
 load,read,node,gate,vertex_escapes,ROOT,real)
OUTPUT=ROOT/'data/enumerations/fresh_h8_cap_completion_20261013_extensions.json'
def run():
 G,I,center,nodes,polys,cn,qn,hn,_=load();source='data/enumerations/fresh_h8_cap_completion_20261013_bottleneck_caps.json';new=read(source);oldN=len(nodes);nodes+=list(node(c,G)for c in new['caps']);cn=np.array([[real(a)+np.sqrt(real(n['d']))*real(b)for a,b in zip(n['a'],n['b'])]for n in nodes]);prior=read('data/enumerations/fresh_h8_cap_completion_20261013_target_compatibility.json');base=prior['base_selected_nodes'];eligible=[];pairchecks=[]
 for j in range(len(new['caps'])):
  ok=True
  for i in base:
   row=gate(nodes[i],nodes[oldN+j],polys);pairchecks.append({'i':i,'j':oldN+j,**row});ok&=row['compatible']
  if ok:eligible.append(oldN+j)
 adjacency={i:set()for i in eligible}
 for i,j in combinations(eligible,2):
  row=gate(nodes[i],nodes[j],polys);pairchecks.append({'i':i,'j':j,**row})
  if row['compatible']:adjacency[i].add(j);adjacency[j].add(i)
 cliques=[]
 def bk(r,p,x):
  if not p and not x:cliques.append(sorted(r));return
  for v in sorted(p.copy()):bk(r|{v},p&adjacency[v],x&adjacency[v]);p.remove(v);x.add(v)
 bk(set(),set(eligible),set());cases=[]
 print({'eligible_new_caps':[i-oldN for i in eligible],'maximal_clique_sizes':[len(c)for c in cliques]},flush=True)
 for k,clique in enumerate(cliques):
  chosen=sorted(base+clique);summary,escapes=vertex_escapes(chosen,I,center,nodes,cn,qn,hn,limit=12);cases.append({'selected_node_indices':chosen,'added_cap_indices':[i-oldN for i in clique],'facets':10+len(chosen),'vertex_gate':summary,'escapes':escapes});out={'new_cap_source':source,'base_selected_nodes':base,'old_dictionary_nodes':oldN,'eligible_new_cap_indices':[i-oldN for i in eligible],'exact_pair_checks':pairchecks,'maximal_extension_cliques':[[i-oldN for i in c]for c in cliques],'cases':cases,'scope':'All maximal compact-compatible subsets of the9 exact caps retaining the full30-facet configuration. Numerical vertices only propose exact escapes; absence of a witness is unresolved.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print({'case':k,'added':[i-oldN for i in clique],**summary},flush=True)
if __name__=='__main__':run()
