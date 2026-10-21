#!/usr/bin/env python3
"""Core-only opposite-single high-label obstruction; inherit polynomial bank."""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
from sympy import QQ,sqrt,symbols,expand
from fresh_h9_star_zero_20261014_own_three import K,z,o,h,co,sign,enc
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_opposite_core_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
OLD='results/fresh_h9_star_opposite_high_20261014.json'
REVIEW='results/fresh_h9_star_opposite_high_20261014_review.json'
S=[3,4,6,7,9,10,12,13];CORE=sorted(S+[11,14,17]);R=K.from_sympy(sqrt(2))
THETA=K.from_sympy(sqrt(2)+sqrt(5))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def decode(xs):
 v=z
 for x in xs:v=v*THETA+K.convert(QQ(x))
 return v
def run():
 seen={}
 def read(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for q,hh in d.get('sources',{}).items():
   assert sha(q)==hh,q
   if q.endswith('.json'):read(q)
  return d
 geo,old,review=[read(p)for p in (GEOM,OLD,REVIEW)]
 assert old['status'].startswith('PASS_')and review['status'].startswith('PASS_')
 assert review['sources'][OLD]==sha(OLD)
 g=geo['representatives'][0]['geometry'];V=[sorted(v['incident'])for v in g['vertices']];VS=list(map(set,V))
 missing={frozenset(m)for m in g['minimal_nonfaces']}
 assert all(any(set(p)<=v for v in VS)for p in combinations(CORE,2))
 assert any(set(S)<=v for v in VS)
 allfaces={tuple(sorted(set(v)&set(CORE)))for v in V}
 maximal=sorted([f for f in allfaces if not any(set(f)<set(g)for g in allfaces)],key=lambda f:(len(f),f))
 assert len(maximal)==9
 def face(nodes):
  nodes=sorted(nodes);vi=next(i for i,v in enumerate(VS)if set(nodes)<=v)
  return dict(support=nodes,actual_vertex=vi)
 guards=[]
 def guard(path,labels,w,kind,pair,test):
  assert len(set(path))==len(path)and set(path)<=set(CORE)
  n=len(path);B=[[o if i==j else z for j in range(n)]for i in range(n)]
  for i,m in enumerate(labels):B[i][i+1]=B[i+1][i]=-co[m]
  assert all(sign(t)>0 for t in w)
  q=sum((w[i]*B[i][j]*w[j]for i in range(n)for j in range(n)),z)
  assert sign(q)<=0
  if q==z:assert all(sum((B[i][j]*w[j]for j in range(n)),z)==z for i in range(n))
  guards.append(dict(kind=kind,ordered_support=path,labels=labels,
    test_pair=sorted(pair),test_label=test,positive_witness=[enc(x)for x in w],
    quadratic_value=enc(q),gram=[[enc(x)for x in row]for row in B],**face(path)))
 for i,j in combinations(range(4),2):
  for a,b in product(S[2*i:2*i+2],S[2*j:2*j+2]):
   aa=next(x for x in S[2*i:2*i+2]if x!=a);bb=next(x for x in S[2*j:2*j+2]if x!=b)
   guard([aa,a,b,bb],[4,3,4],[o,R,R,o],'cross_block_right',[a,b],3)
 guard([3,14,11,13],[3,5,3],[h,o,o,h],'third_contact_upper',[11,14],5)
 guard([3,14,17,12],[3,5,3],[h,o,o,h],'third_contact_upper',[14,17],5)
 guard([11,14,17],[4,4],[o,R,o],'two_label_four_contacts_impossible',[11,14],4)
 assert len(guards)==27
 mutual=dict(ordered_support=[13,11,17],test_label=6,known_label=3,
   edge_cosine_squares=['1/4','3/4'],determinant='0',**face([13,11,17]))
 assert QQ(1)-QQ(1,4)-QQ(3,4)==0
 budgets=[dict(leaf=leaf,removed_own_root=drop,retained_own_root=keep,
   **face((set(S)-{drop})|{leaf}))for leaf,drop,keep in [(11,12,13),(17,13,12),(14,3,4)]]
 assert all(len(row['support'])==8 for row in budgets)
 assert QQ(1,4)+2*QQ(1,4)*2>1 and QQ(1,4)+QQ(1,2)*2>1
 assert QQ(1,4)+4*QQ(1,4)>1
 high_third=face([12,13,14]);assert QQ(1)-QQ(3,4)-QQ(1,4)==0
 assert all(frozenset([12,13,a])in missing for a in (11,17))
 # A local face-list automorphism suffices for these core-only arguments.
 local=list(range(18));local[11],local[17]=17,11;local[12],local[13]=13,12
 assert {tuple(sorted(local[i]for i in f))for f in maximal}==set(maximal)
 a,b=symbols('a b',nonnegative=True)
 assert expand((2+b*b)*(1+a*a)-(sqrt(2)+a*b)**2-(sqrt(2)*a-b)**2)==0
 assert (1/(2*(1-sqrt(3)/2))-(2+sqrt(3))).simplify()==0
 plane=dict(ambient_vectors=['r14','r11+r17'],first_own_cost='X=1/[2(1-c_A)]>=2+sqrt3',
   last_sum_own_cost='Y=1/[2(1-c_D)]>2+sqrt3',ambient_sum_norm='2+2<r11,r17><=2',
   residual_norms=['-(K+H)','-(L+F)'],parameters='K=X-1>2,L=Y-(r11+r17)^2>1,F>=0,H>=0',
   mutual_absolute_bound='sqrt2+sqrt(F*H)',
   determinant_lower_identity='(2+H)*(1+F)-(sqrt2+sqrt(F*H))^2=(sqrt(2F)-sqrt(H))^2',
   conclusion='Strictly positive determinant and negative diagonals give a negative plane. Strictness follows from K>2 and L>1 even if the displayed square is zero.')
 base={(3,14):3,(4,14):3,(11,12):2,(11,13):3,(12,17):3,(13,17):2,(12,13):4}
 cache={}
 def component_failure(nodes,lab):
  nodes=sorted(nodes);key=(len(nodes),tuple(lab.get((a,b),2)for a,b in combinations(nodes,2)))
  if key not in cache:
   B=[[o if a==b else-co[lab.get(tuple(sorted((a,b))),2)]for b in nodes]for a in nodes]
   piv=[];answer=None
   for k in range(len(nodes)):
    q=B[k][k];piv.append(enc(q))
    if sign(q)<=0:answer=dict(prefix_order=k+1,pivots=piv);break
    for i in range(k+1,len(nodes)):
     for j in range(k+1,len(nodes)):B[i][j]-=B[i][k]*B[k][j]/q
   cache[key]=answer
  ans=cache[key]
  if ans is None:return None
  return dict(support=nodes[:ans['prefix_order']],component_order=nodes,component_pair_labels=list(key[1]),**ans)
 def failure(lab):
  for f in maximal:
   left=set(f)
   while left:
    start=min(left);left.remove(start);component={start};queue=[start]
    for a in queue:
     nxt={b for b in left if lab.get(tuple(sorted((a,b))),2)>2}
     left-=nxt;component|=nxt;queue+=sorted(nxt)
    if len(component)<3:continue
    ans=component_failure(component,lab)
    if ans is not None:return dict(actual_vertex=next(i for i,v in enumerate(VS)if set(f)<=v),containing_core_face=list(f),**ans)
  return None
 expected={tuple(block)+(ia,ib,ic,mu,ma,mb)for block in product((4,5),repeat=3)
   for ia,ib,ic in product([None]+list(range(6)),[None]+list(range(6)),[None,2,3,4,5])
   for mu,ma,mb in product((2,3,4,5),(2,3),(2,3,4))}
 keys=set();new_fail=[];counts=Counter()
 for idx,row in enumerate(old['cases']):
  block=tuple(row['block_labels']);ia,ib,ic=row['foreign_contacts'];mu=row['mutual_label'];ma,mb=row['third_labels']
  key=block+(ia,ib,ic,mu,ma,mb);assert key not in keys;keys.add(key);counts[row['status']]+=1
  if row['status']=='STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE':
   coeff=list(map(decode,row['bernstein_6_7']));assert all(sign(x)>=0 for x in coeff)and any(sign(x)>0 for x in coeff);continue
  if row['status']=='STRICT_BERNSTEIN_9_10_TO_ONE':
   assert all(sign(decode(x))>0 for x in row['bernstein_9_10']);continue
  assert row['status']=='NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR'
  lab=base|{tuple(S[2*i:2*i+2]):block[i]for i in range(3)}|{(11,17):mu,(11,14):ma,(14,17):mb}
  for leaf,extra in zip((11,17,14),(ia,ib,ic)):
   if extra is not None:lab[tuple(sorted((leaf,S[extra])))]=3
  witness=failure(lab);assert witness is not None,(idx,key)
  assert set(witness['support'])<=set(CORE)
  mapped=sorted(local[i]for i in witness['support']);assert any(set(mapped)<=v for v in VS)
  new_fail.append(dict(old_case_index=idx,block_labels=list(block),foreign_contacts=[ia,ib,ic],mutual_label=mu,
    third_labels=[ma,mb],relaxed_core_labels=sorted([list(p)+[m]for p,m in lab.items()]),witness=witness))
 assert keys==expected and len(keys)==47040 and len(new_fail)==3612
 assert counts=={'STRICT_INTERIOR_BERNSTEIN_6_7_TO_ONE':43412,'NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR':3612,'STRICT_BERNSTEIN_9_10_TO_ONE':16}
 swap=list(range(18));swap[12],swap[13]=13,12
 assert {tuple(sorted(swap[i]for i in v))for v in V}==set(map(tuple,V))
 assert {frozenset(swap[i]for i in m)for m in missing}==missing
 fixed=sorted([list(p)+[m]for p,m in base.items()if p!=(12,13)])
 bounds=[S[2*i:2*i+2]+[4 if i<3 else 7]for i in range(4)]
 families=[]
 for reverse in (False,True):
  f=sorted([sorted((swap[a],swap[b]))+[m]for a,b,m in fixed])if reverse else fixed
  families.append(dict(id='core_opposite_single_tail'+('_swapped'if reverse else''),origin='core_only_opposite_single_all_internal_labels',fixed=f,bounds=bounds))
 paths=[GEOM,OLD,REVIEW,'scripts/fresh_h9_star_zero_20261014_own_three.py',f'scripts/{TAG}.py']
 selected={(len(x['witness']['component_order']),tuple(x['witness']['component_pair_labels']))for x in new_fail}
 return dict(status='PASS_CORE_OPPOSITE_SINGLE_ALL_INTERNAL_LABELS',case=0,sources={p:sha(p)for p in paths},rank_core=CORE,
  core_positive_quadratic_guards=guards,mutual_upper_six=mutual,projection_budgets=budgets,high_last_third_right=high_third,
  complete_maximal_core_faces=[list(f)for f in maximal],local_core_swap=local,
  local_swap_scope='Preserves maximal induced core faces only; no full-incidence automorphism claimed. It exchanges the two third-contact labels, preserving the opposite-single own pattern. The4/4 possibility is excluded by an actual affine path; otherwise one third label is<=3.',
  high_first_negative_plane=plane,inherited_cases=len(keys),inherited_status_counts=dict(counts),
  new_core_fallback_count=len(new_fail),new_core_fallbacks=new_fail,ldl_matrix_types_tried=len(cache),selected_full_matrix_types=len(selected),
  literal_last_root_swap=swap,families=families,
  high_middle_reduction='Any foreign nonright contact to an internal>=6 block costs at least1 plus retained own1/4. Thus all core contacts to a high middle block vanish, and replacing its internal label by4 leaves the residual unchanged.',
  inherited_tail='The frozen independently reviewed polynomial bank covers every finite last m>=7: cos(pi/m)>6/7;16 cases use the independently proved sharper cos(pi/m)>9/10. No finite-label cutoff.',
  verification_scope='Fresh core domain/negative-plane proof and3612 exact core fallback LDL witnesses; inherited separately reviewed47040-case determinant/Bernstein bank. Local core-face symmetry extends asymmetric third-contact domains.',
  scope='Star0 first3internal labels>=4,last>=7;14own3/3;11and17 opposite single ownpi/3 contacts. Every noncore own contact and all remaining finite contacts/dotted distances free. Both last-root orientations. Not a whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_core_fallback_count'],d['ldl_matrix_types_tried'],d['selected_full_matrix_types'])
