#!/usr/bin/env python3
"""Exact selected traditional-angle cliques; rational strict escape gates."""
from itertools import combinations
import json
import numpy as np
from scripts.fresh_h8_alternate_caps_20261013 import *
from scripts.fresh_h8_release_h9_free_caps_20261013 import bsign,bmul
from scripts.fresh_h8_alternate_coupled_20261013 import escape
BASE=ROOT/'data/enumerations/fresh_h8_alternate_3_6_20261013.json'
POLAR=ROOT/'data/enumerations/fresh_h8_alternate_3_6_polar_20261013.json'
CAPS=ROOT/'data/enumerations/fresh_h8_alternate_polar_coupled_20261013.json'
OUTPUT=ROOT/'data/enumerations/fresh_h8_alternate_selection_20261013.json'

def load():
 base=json.loads(BASE.read_text());pol=json.loads(POLAR.read_text());raw=json.loads(CAPS.read_text());G=matrix(base['seed_gram']);I=matrix(base['seed_gram_inverse']);center=tuple(map(parse,base['strict_center_y']));nodes=[]
 for j,c in enumerate(base['caps']):nodes.append({'source':'base','index':j,'d':Z,'a':tuple(map(parse,c['seed_coefficients'])),'b':(Z,)*10,'pa':tuple(map(parse,c['seed_pairings'])),'pb':(Z,)*10,'q':parse(c['norm'])})
 for j,c in enumerate(pol['polar_caps']):
  a=tuple(map(parse,c['seed_coefficients']));nodes.append({'source':'polar','index':j,'d':Z,'a':a,'b':(Z,)*10,'pa':mv(G,a),'pb':(Z,)*10,'q':parse(c['norm'])})
 for j,c in enumerate(raw['caps']):nodes.append({'source':'coupled','index':j,'d':parse(c['discriminant']),'a':tuple(map(parse,c['coefficients_rational'])),'b':tuple(map(parse,c['coefficients_radical'])),'pa':tuple(map(parse,c['pairings_rational'])),'pb':tuple(map(parse,c['pairings_radical'])),'q':scale(O,2)})
 return base,pol,raw,G,I,center,nodes

def pairing(a,b):return (dot(a['pa'],b['a']),dot(a['pb'],b['a']),dot(a['pa'],b['b']),dot(a['pb'],b['b']))
def subtract(p,c):return (sub(p[0],c),*p[1:])
def test(a,b):
 p=pairing(a,b);d,e=a['d'],b['d']
 if bsign(p,d,e)>0:return False,'positive'
 sq=bmul(p,p,d,e);ratio=tuple(mul(div(scale(O,4),mul(a['q'],b['q'])),x)for x in sq)
 if bsign(subtract(ratio,scale(O,4)),d,e)>0:return True,'dotted'
 for m,c in ((2,Z),(3,O),(4,scale(O,2)),(5,add(O,P)),(6,scale(O,3)),(10,add(scale(O,2),P))):
  if bsign(subtract(ratio,c),d,e)==0:return True,str(m)
 return False,'outside_selected_alphabet_or_parallel'

def run():
 base,pol,raw,G,I,center,nodes=load();N=len(nodes);cn=np.array([[real(x)+real(r['d'])**.5*real(y)for x,y in zip(r['a'],r['b'])]for r in nodes]);gn=np.array([[real(x)for x in r]for r in G]);Q=np.array([real(r['q'])for r in nodes]);pn=cn@gn@cn.T;ratio=pn**2/(Q[:,None]*Q[None,:]);allowed=[0,.25,.5,(1+real(P))/4,.75,(2+real(P))/4];neighbors=[set()for _ in nodes];edge_types={};exact_pairs=0
 for i,j in combinations(range(N),2):
  if pn[i,j]>1e-7:continue
  if ratio[i,j]<=1+1e-7 and min(abs(ratio[i,j]-v)for v in allowed)>1e-7:continue
  ok,why=test(nodes[i],nodes[j]);exact_pairs+=1
  if ok:neighbors[i].add(j);neighbors[j].add(i);edge_types[f'{i},{j}']=why
 print({'nodes':N,'selected_alphabet_edges':len(edge_types),'exact_checked_pairs':exact_pairs},flush=True)
 cases=[];seen=set();rng=np.random.default_rng(20261013);polarid=len(base['caps'])
 for trial in range(20):
  forced=[polarid]
  if trial>0:
   newer=[i for i in neighbors[polarid]if nodes[i]['source']=='coupled']
   if newer:forced.append(newer[(trial-1)%len(newer)])
  chosen=list(forced);pool=set(range(N))-set(chosen)
  for i in chosen:pool&=neighbors[i]
  noise=rng.random(N)
  while pool:
   if trial%3==0:i=max(pool,key=lambda j:(len(neighbors[j]&pool),noise[j]))
   elif trial%3==1:i=max(pool,key=lambda j:(nodes[j]['source']=='coupled',len(neighbors[j]&pool),noise[j]))
   else:i=max(pool,key=lambda j:noise[j])
   chosen.append(i);pool&=neighbors[i]
  chosen=sorted(chosen);key=tuple(chosen)
  if key in seen:continue
  seen.add(key);old=[nodes[i]['a']for i in chosen if nodes[i]['source']!='coupled'];caps=[raw['caps'][nodes[i]['index']]for i in chosen if nodes[i]['source']=='coupled'];gate=escape(I,old,caps,center,tries=500)
  row={'index':len(cases),'selected_node_indices':chosen,'total_facet_inequalities':10+len(chosen),'coupled_caps':len(caps),'volume_gate':gate};cases.append(row)
  out={'schema':'alternate-H8-compatible-subsets-v1','base_source':str(BASE.relative_to(ROOT)),'polar_source':str(POLAR.relative_to(ROOT)),'coupled_source':str(CAPS.relative_to(ROOT)),'nodes':[{'source':r['source'],'index':r['index']}for r in nodes],'selected_finite_labels':[2,3,4,5,6,10],'edges':edge_types,'cases':cases,'scope':'Exact compatible pairs in a restricted finite-angle search graph, with numeric proposal prefilter. Absent edges do not constitute an exhaustive incompatibility proof; subset escapes certify only the archived subsets.'};OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print({'case':row['index'],'facets':row['total_facet_inequalities'],'new_caps':row['coupled_caps'],'gate':gate['status']},flush=True)
if __name__=='__main__':run()
