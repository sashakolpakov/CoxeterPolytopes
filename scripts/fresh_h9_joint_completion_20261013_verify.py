#!/usr/bin/env python3
"""Producer-free Fraction replay of joint anchors, cut masks and replacements."""
import argparse,json
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
from scripts.fresh_h8_threefree_verify_20261013 import (
 old_load,node,gate,read,pair,parse,enc,matrix,sign,dot,qsign,qa,qs,qm,
 Z,O,P,add,sub,neg,mul,scale,mv,div,bsign)
ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'data/enumerations/fresh_h9_joint_completion_20261013_verify.json'
def inertia(matrix):
 a=[list(row)for row in matrix];counts=[0,0,0]
 while a:
  pivot=next((i for i in range(len(a))if a[i][i]!=Z),None)
  if pivot is not None:
   order=[pivot]+[i for i in range(len(a))if i!=pivot]
   a=[[a[i][j]for j in order]for i in order];p=a[0][0]
   counts[0 if sign(p)>0 else 1]+=1
   a=[[sub(a[i][j],div(mul(a[i][0],a[0][j]),p))for j in range(1,len(a))]for i in range(1,len(a))]
   continue
  edge=next(((i,j)for i in range(len(a))for j in range(i+1,len(a))if a[i][j]!=Z),None)
  if edge is None:counts[2]+=len(a);break
  i,j=edge;order=[i,j]+[k for k in range(len(a))if k not in edge]
  a=[[a[i][j]for j in order]for i in order];p=a[0][1]
  counts[0]+=1;counts[1]+=1
  a=[[sub(a[i][j],div(add(mul(a[i][0],a[1][j]),mul(a[i][1],a[0][j])),p))for j in range(2,len(a))]for i in range(2,len(a))]
 return counts
def mirror_census(old,new):
 """Independent exact normalized-pairing-square comparison.

 All ten seed pairings are nonpositive. Their normalized squares therefore
 determine the oriented normal, since the seed normals are independent.
 The key only discards pairs whose equality is already exactly impossible.
 """
 groups=defaultdict(list);aliases=[];novel=[];tests=0
 def pattern(n):
  f=div(scale(O,2),n['q'])
  values=[tuple(mul(f,t)for t in qm((a,b),(a,b),n['d']))for a,b in zip(n['pa'],n['pb'])]
  dotted=tuple(i for i,s in enumerate(values)if qsign(qs(s,(scale(O,4),Z)),n['d'])>0)
  if len(dotted)!=3:return None,values
  fixed=[]
  for i,s in enumerate(values):
   if i in dotted:continue
   v=next((j for j,t in enumerate((Z,O,add(O,P)))if qsign(qs(s,(t,Z)),n['d'])==0),None)
   if v is None:return None,values
   fixed.append(v)
  return(dotted,tuple(fixed)),values
 for i,n in enumerate(old):
  k,values=pattern(n)
  if k is not None:groups[k].append((n['d'],values,{'dictionary':'old','index':i}))
 for i,n in enumerate(new):
  k,values=pattern(n);assert k is not None;alias=None
  for d,other,origin in groups[k]:
   tests+=1
   if all(bsign((sub(a[0],b[0]),a[1],neg(b[1]),Z),n['d'],d)==0 for a,b in zip(values,other)):
    alias=origin;break
  if alias is None:
   novel.append(i);groups[k].append((n['d'],values,{'dictionary':'joint','index':i}))
  else:aliases.append({'candidate':i,**alias})
 archive=read('data/enumerations/fresh_h9_joint_completion_20261013_mirrors.json')
 assert novel==archive['new_representative_indices']and aliases==archive['aliases']
 assert tests==archive['exact_same_pattern_pair_checks']and len(novel)==archive['distinct_new_mirrors']
 return novel,{'distinct_new_mirrors':len(novel),'alias_entries':len(aliases),'exact_same_pattern_pair_checks':tests}
def run(verify=False):
 data=read('data/enumerations/fresh_h9_joint_completion_20261013.json');base=read(data['base_source']);G=matrix(base['seed_gram']);I=matrix(base['seed_gram_inverse']);center=tuple(map(parse,base['strict_center_y']));_,_,old=old_load()
 for source in('data/enumerations/fresh_h8_threefree_combined_caps_20261013.json','data/enumerations/fresh_h8_cap_completion_20261013_bottleneck_caps.json'):old+=list(node(c,G)for c in read(source)['caps'])
 assert len(old)==1084 and sign(pair(center,I,center))<0 and all(sign(x)>0 for x in center)
 assert inertia(G)==[9,1,0]
 assert all(sum_k== (O if i==j else Z)for i,row in enumerate(G)for j,sum_k in enumerate(mv(I,row)))
 for n in old:
  assert sign(n['q'])>0
  norm=(add(dot(n['pa'],n['a']),mul(n['d'],dot(n['pb'],n['b']))),add(dot(n['pa'],n['b']),dot(n['pb'],n['a'])))
  assert qsign(qs(norm,(n['q'],Z)),n['d'])==0
  assert all(qsign((a,b),n['d'])<=0 for a,b in zip(n['pa'],n['pb']))
  assert qsign((dot(n['a'],center),dot(n['b'],center)),n['d'])>0
 ad=read('results/fresh_h8_alternate_core_angles_independent_20261013.json');polys=[(r['label'],tuple(map(F,r['coefficients_descending'])),F(r['largest_root_interval'][0]),F(r['largest_root_interval'][1]))for r in ad['real_cyclotomic_polynomials']];target=read(data['target_source']);chosen=target['selected_node_indices'];assert chosen==data['selected_configuration'];ys=[]
 for row in target['targets']:
  c=row['certificate'];y=tuple(map(parse,c['y']));q=pair(y,I,y);cp=pair(y,I,center);assert all(sign(x)>0 for x in y)and sign(q)>0 and sign(cp)<0;assert q==parse(c['norm'])and cp==parse(c['center_pairing']);assert all(qsign((dot(old[i]['a'],y),dot(old[i]['b'],y)),old[i]['d'])>0 for i in chosen);cuts=[i for i,n in enumerate(old)if qsign((dot(n['a'],y),dot(n['b'],y)),n['d'])<0];assert cuts==c['complete_dictionary_cutters'];ys.append(y)
 assert all(gate(old[i],old[j],polys)['compatible']for h,i in enumerate(chosen)for j in chosen[h+1:])
 new=[];paircount=0;cutcount=0
 for ci,c in enumerate(data['caps']):
  n=node(c,G);d=n['d'];assert sign(d)>0 and n['pa']==tuple(map(parse,c['pairings_rational']))and n['pb']==tuple(map(parse,c['pairings_radical']));w0=tuple(map(parse,c['line_constant']));w1=tuple(map(parse,c['line_direction']));A=pair(w1,I,w1);B=scale(pair(w0,I,w1),2);C=sub(pair(w0,I,w0),scale(O,2));assert [A,B,C]==list(map(parse,c['norm_polynomial_descending']));assert sub(mul(B,B),scale(mul(A,C),4))==d;x=parse(c['parameter']['rational']),parse(c['parameter']['radical']);assert qsign(qa(qa(tuple(mul(A,t)for t in qm(x,x,d)),tuple(mul(B,t)for t in x)),(C,Z)),d)==0
  for j in range(10):
   assert qsign(qa((n['pa'][j],n['pb'][j]),qa((w0[j],Z),tuple(mul(w1[j],t)for t in x))),d)==0;assert qsign((n['pa'][j],n['pb'][j]),d)<=0;sq=qm((n['pa'][j],n['pb'][j]),(n['pa'][j],n['pb'][j]),d)
   if j in c['free_seed_indices']:assert qsign(qs(sq,(scale(O,4),Z)),d)>0
   else:assert any(qsign(qs(sq,(t,Z)),d)==0 for t in(Z,O,add(O,P)))
  norm=(add(dot(n['pa'],n['a']),mul(d,dot(n['pb'],n['b']))),add(dot(n['pa'],n['b']),dot(n['pb'],n['a'])));assert qsign(qs(norm,(scale(O,2),Z)),d)==0;assert qsign((dot(n['a'],center),dot(n['b'],center)),d)>0
  for i,M in zip(c['anchor_node_indices'],map(parse,c['anchor_pairing_magnitudes'])):
   assert i in chosen and old[i]['d']==Z;assert dot(old[i]['a'],w0)==M and dot(old[i]['a'],w1)==Z
  kept=[];dropped=[]
  assert [r['old_node']for r in c['exact_pair_checks']]==chosen
  for row in c['exact_pair_checks']:
   i=row['old_node'];got=gate(n,old[i],polys);assert got=={k:row[k]for k in('compatible','kind','label')};(kept if got['compatible']else dropped).append(i);paircount+=1
  assert kept==c['exact_retained_nodes']and dropped==c['exact_dropped_nodes']and len(kept)>=data['minimum_retained_caps'];assert set(c['anchor_node_indices'])<=set(kept);cuts=[i for i,y in enumerate(ys)if qsign((dot(n['a'],y),dot(n['b'],y)),d)<0];assert cuts==c['cut_target_indices']and len(cuts)>=2;cutcount+=len(ys);new.append(n)
 replacements=read('data/enumerations/fresh_h9_joint_completion_20261013_replacements.json');reps,mirror_result=mirror_census(old,new);scenario_results=[]
 for scenario in replacements['scenarios']:
  cutter=scenario['cutter'];retained=[];dropped=[]
  for j in chosen:
   (retained if gate(old[cutter],old[j],polys)['compatible']else dropped).append(j)
  assert dropped==scenario['dropped_original_caps'];assert sorted(retained+[cutter])==scenario['fixed_cap_nodes'];eligible=[]
  for j in reps:
   if not set(retained)<=set(data['caps'][j]['exact_retained_nodes']):continue
   got=gate(new[j],old[cutter],polys);assert got==replacements['exact_cross_checks'][f'{j},{cutter}']
   if got['compatible']:eligible.append(j)
  assert eligible==scenario['eligible_new_cap_indices'];already=[j for j,y in enumerate(ys)if qsign((dot(old[cutter]['a'],y),dot(old[cutter]['b'],y)),old[cutter]['d'])<0];assert already==scenario['targets_cut_by_cutter'];clauses=[{'target':t,'cutters':[j for j in eligible if t in data['caps'][j]['cut_target_indices']]}for t in range(len(ys))if t not in already];assert clauses==scenario['necessary_cover_clauses'];missing=[r['target']for r in clauses if not r['cutters']];assert missing==scenario['uncut_target_indices']and missing;assert scenario['status']=='EXACT_UNCOVERED_TARGET';scenario_results.append({'cutter':cutter,'dropped':dropped,'eligible_count':len(eligible),'exact_uncovered_targets':missing})
 common=[i for i in range(len(ys))if all(i not in c['cut_target_indices']for c in data['caps'])];assert common==[8]
 anchor_aq=pair(old[1039]['a'],G,old[1039]['a']);anchor_bq=pair(old[1039]['b'],G,old[1039]['b']);anchor_ab=pair(old[1039]['a'],G,old[1039]['b'])
 assert anchor_aq==add(scale(O,F(7,2)),scale(P,4))and anchor_bq==add(scale(O,F(-9,4)),scale(P,F(-29,8)))and anchor_ab==Z
 out={'status':'EXACT_JOINT_ROOT_TARGET_AND_REPLACEMENT_REPLAY_PASS','seed_inertia':[9,1,0],'old_norm_and_orientation_checks':len(old),'starting_configuration_pair_checks':len(chosen)*(len(chosen)-1)//2,'accepted_cap_entries':len(new),'exact_target_witnesses':len(ys),'old_dictionary_target_sign_checks':len(old)*len(ys),'exact_retention_pair_checks':paircount,'exact_new_cap_target_sign_checks':cutcount,'joint_cut_count_range':[min(len(c['cut_target_indices'])for c in data['caps']),max(len(c['cut_target_indices'])for c in data['caps'])],'mirror_identification':mirror_result,'full_new_bank_surviving_target_indices':common,'quartic_anchor_decomposition':{'node':1039,'rational_part_norm':enc(anchor_aq),'radical_part_norm':enc(anchor_bq),'cross_pairing':enc(anchor_ab)},'coherent_replacement_checks':scenario_results,'scope':'Exact accepted roots and multiple cuts, all retained/dropped wall relations, distinct-mirror comparisons, and complete eligible-cap coverage checks for six specified coherent replacement schemes. This is not an exhaustive cap-family or global H9 obstruction.'}
 if verify:assert read(str(OUTPUT.relative_to(ROOT)))==out
 else:OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 print(out,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
