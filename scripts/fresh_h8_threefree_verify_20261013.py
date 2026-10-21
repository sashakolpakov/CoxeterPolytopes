#!/usr/bin/env python3
"""Exact replay of new cap equations and escaped compatible configurations.

Uses the previously independent field arithmetic, not the cap or SAT producers.
Floating point, LP and Z3 are absent. An unresolved proposal graph is not audited
as a complete dictionary graph.
"""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from scripts.fresh_h8_alternate_core_angles_independent_20261013 import (
 load as old_load,qa,qs,qm,qsign,bsign,bmul,pairing,finite_label,
 Z,O,P,parse,enc,add,sub,neg,scale,mul,div,sign,total,dot,mv,matrix)

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'data/enumerations/fresh_h8_threefree_verify_20261013.json'
def pair(a,G,b):return dot(a,mv(G,b))
def read(p):return json.loads((ROOT/p).read_text())
def node(c,G):
 a=tuple(map(parse,c['coefficients_rational']));b=tuple(map(parse,c['coefficients_radical']));d=parse(c['discriminant'])
 return {'a':a,'b':b,'pa':mv(G,a),'pb':mv(G,b),'d':d,'q':scale(O,2)}
def sq(a,b,d):return qm((a,b),(a,b),d)
def gate(a,b,polys):
 d,e=a['d'],b['d'];p=pairing(a,b)
 if bsign(p,d,e)>0:return {'compatible':False,'kind':'positive','label':None}
 ratio=tuple(mul(div(scale(O,4),mul(a['q'],b['q'])),x)for x in bmul(p,p,d,e));cmp=bsign((sub(ratio[0],scale(O,4)),*ratio[1:]),d,e)
 if cmp>0:return {'compatible':True,'kind':'dotted','label':None}
 if cmp==0:return {'compatible':False,'kind':'parallel','label':None}
 label=finite_label(ratio,d,e,polys)
 return {'compatible':label is not None,'kind':'finite'if label is not None else'noncoxeter','label':label}
def run(source,selection=None,output=None,verify=False):
 data=read(source);base=read(data['base_source']);G=matrix(base['seed_gram']);I=matrix(base['seed_gram_inverse']);center=tuple(map(parse,base['strict_center_y']));_,_,old=old_load();assert len(old)==493
 assert all(dot(G[i],tuple(row[j]for row in I))==(O if i==j else Z)for i in range(10)for j in range(10));assert sign(pair(center,I,center))<0 and all(sign(x)>0 for x in center)
 ad=read('results/fresh_h8_alternate_core_angles_independent_20261013.json');polys=[(r['label'],tuple(map(F,r['coefficients_descending'])),F(r['largest_root_interval'][0]),F(r['largest_root_interval'][1]))for r in ad['real_cyclotomic_polynomials']]
 proof=read(data['witness_source']);ys=[tuple(map(parse,r['escape_y']))for r in proof['five_complete_escape_clauses']]if'five_complete_escape_clauses'in proof else[tuple(map(parse,r['certificate']['y']))for r in proof['clauses']]
 nodes=old.copy();checks=[];old_dotted=[]
 for n in old:
  assert qsign((dot(n['a'],center),dot(n['b'],center)),n['d'])>0
  old_dotted.append(tuple(j for j in range(10)if qsign(qs(tuple(mul(div(scale(O,2),n['q']),x)for x in sq(n['pa'][j],n['pb'][j],n['d'])),(scale(O,4),Z)),n['d'])>0))
 for ci,c in enumerate(data['caps']):
  n=node(c,G);d=n['d'];assert sign(d)>0;assert tuple(map(parse,c['pairings_rational']))==n['pa'] and tuple(map(parse,c['pairings_radical']))==n['pb']
  norm=(add(dot(n['pa'],n['a']),mul(d,dot(n['pb'],n['b']))),add(dot(n['pa'],n['b']),dot(n['pb'],n['a'])));assert qsign(qs(norm,(scale(O,2),Z)),d)==0
  w0=tuple(map(parse,c['line_constant']));w1=tuple(map(parse,c['line_direction']));A=pair(w1,I,w1);B=scale(pair(w0,I,w1),2);C=sub(pair(w0,I,w0),scale(O,2));assert [A,B,C]==list(map(parse,c['norm_polynomial_descending']));assert sub(mul(B,B),scale(mul(A,C),4))==d
  x=(parse(c['parameter']['rational']),parse(c['parameter']['radical']));assert qsign(qa(qa(tuple(mul(A,t)for t in qm(x,x,d)),tuple(mul(B,t)for t in x)),(C,Z)),d)==0
  for j in range(10):assert qsign(qa((n['pa'][j],n['pb'][j]),qa((w0[j],Z),tuple(mul(w1[j],t)for t in x))),d)==0
  dotted=[]
  for j,(a,b)in enumerate(zip(n['pa'],n['pb'])):
   assert qsign((a,b),d)<=0;value=sq(a,b,d)
   if qsign(qs(value,(scale(O,4),Z)),d)>0:dotted.append(j)
   else:assert finite_label((value[0],value[1],Z,Z),d,Z,polys)is not None
  assert dotted==sorted(c['free_seed_indices']) and len(dotted)==3
  assert c['free_original_walls']==[base['retained_original_walls'][j]for j in c['free_seed_indices']]
  assert qsign((dot(n['a'],center),dot(n['b'],center)),d)>0
  for ai,mag in zip(c['anchor_node_indices'],map(parse,c['anchor_pairing_magnitudes'])):
   assert old[ai]['d']==Z and bsign((add(pairing(n,old[ai])[0],mag),*pairing(n,old[ai])[1:]),d,Z)==0
   assert gate(n,old[ai],polys)['compatible']
  cuts=[]
  for j,y in enumerate(ys):
   val=dot(n['a'],y),dot(n['b'],y);saved=c['core_constraint_values'][j];assert val==(parse(saved['rational']),parse(saved['radical']))
   if qsign(val,d)<0:cuts.append(j)
  assert cuts==c['cut_core_witness_indices'] and cuts
  for oi,m in enumerate(old):
   if old_dotted[oi]!=tuple(dotted):continue
   factor=div(scale(O,2),m['q']);different=False
   for j in range(10):
    p=sq(n['pa'][j],n['pb'][j],d);q=sq(m['pa'][j],m['pb'][j],m['d'])
    if bsign((sub(p[0],mul(factor,q[0])),p[1],neg(mul(factor,q[1])),Z),d,m['d'])!=0:different=True;break
   assert different,(ci,oi)
  nodes.append(n);checks.append({'cap':ci,'three_dotted_seed_indices':dotted,'cut_witnesses':cuts})
 paircount=casecount=clausecount=0
 if selection:
  sel=read(selection);assert sel['node_count']==len(nodes);assert sel['new_cap_source']==source
  for key,expected in sel['exact_pair_checks'].items():
   i,j=map(int,key.split(','));assert gate(nodes[i],nodes[j],polys)==expected;paircount+=1
  for case in sel['cases']:
   chosen=case['selected_node_indices'];assert 348 in chosen
   for i,j in combinations(chosen,2):assert sel['exact_pair_checks'][f'{i},{j}']['compatible']
   r=case['volume_gate'];y=tuple(map(parse,r['y']));assert r['status']=='EXACT_STRICT_SPACELIKE_ESCAPE';assert all(sign(t)>0 for t in y);q=pair(y,I,y);cp=pair(y,I,center);assert sign(q)>0 and sign(cp)<0;assert q==parse(r['norm']) and cp==parse(r['center_pairing'])
   assert all(qsign((dot(nodes[i]['a'],y),dot(nodes[i]['b'],y)),nodes[i]['d'])>0 for i in chosen);casecount+=1
  for row in sel['escape_clauses']:
   r=row['certificate'];y=tuple(map(parse,r['y']));assert sign(pair(y,I,y))>0 and sign(pair(y,I,center))<0 and all(sign(t)>0 for t in y)
   cuts=[i for i in sel['active_nodes']if qsign((dot(nodes[i]['a'],y),dot(nodes[i]['b'],y)),nodes[i]['d'])<0];assert cuts==row['complete_cutters_among_active_nodes'];clausecount+=1
 result={'status':'EXACT_CAP_AND_ESCAPED_SUBSET_REPLAY_PASS','cap_source':source,'selection_source':selection,'cap_count':len(checks),'exact_pair_classifications':paircount,'escaped_compatible_subsets':casecount,'complete_active_node_escape_clauses':clausecount,'cap_checks':checks,'scope':'Exact roots, seed/anchor angles, strict common center, distinction from old493, selected pair classifications, and strict future spacelike escapes. This does not certify proposal-graph completeness or a global no-go.'}
 path=ROOT/output if output else OUTPUT
 if verify:assert read(str(path.relative_to(ROOT)))==result
 else:path.write_text(json.dumps(result,indent=2)+'\n')
 print({k:v for k,v in result.items()if k!='cap_checks'},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',default='data/enumerations/fresh_h8_threefree_caps_20261013.json');p.add_argument('--selection',default='data/enumerations/fresh_h8_threefree_selection_20261013.json');p.add_argument('--output');p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.source,a.selection or None,a.output,a.verify)
