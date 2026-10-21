#!/usr/bin/env python3
"""Core-only shared-single tail obstruction, with arbitrary internal labels.

The frozen affine determinant bank is inherited. Fresh exact arithmetic
replaces all448 fallback witnesses and verifies the core-only domain guards.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
from sympy import QQ,sqrt,symbols,expand
from fresh_h9_star_zero_20261014_own_three import K,z,o,h,co,sign,enc

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_mixed_singles_20261014'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
OLD='results/fresh_h9_star_shared_high_20261014.json'
REVIEW='results/fresh_h9_star_shared_high_20261014_review.json'
S=[3,4,6,7,9,10,12,13];CORE=sorted(S+[11,14,17])
THETA=K.from_sympy(sqrt(2)+sqrt(5));R=K.from_sympy(sqrt(2))
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
 assert old['status']=='PASS_SHARED_SINGLE_OWN_ROOT_ALL_HIGH_LABELS_EXCLUDED'
 assert review['status'].startswith('PASS_')and review['sources'][OLD]==sha(OLD)
 g=geo['representatives'][0]['geometry'];V=[sorted(v['incident'])for v in g['vertices']];VS=list(map(set,V))
 missing={frozenset(m)for m in g['minimal_nonfaces']}
 assert all(any(set(p)<=v for v in VS)for p in combinations(CORE,2))
 assert any(set(S)<=v for v in VS)
 corefaces={tuple(sorted(set(v)&set(CORE)))for v in V}
 maximal=sorted([f for f in corefaces if not any(set(f)<set(g)for g in corefaces)],key=lambda f:(len(f),f))
 assert len(maximal)==9
 def face(nodes):
  nodes=sorted(nodes);vi=next(i for i,v in enumerate(VS)if set(nodes)<=v)
  return dict(support=nodes,actual_vertex=vi)
 # Fresh positive null kernels require only rank-core contacts.
 guards=[]
 def guard(path,labels,kernel,kind,pair,test,triangle=False):
  assert len(set(path))==len(path)and set(path)<=set(CORE)
  n=len(path);B=[[o if i==j else z for j in range(n)]for i in range(n)]
  for i,m in enumerate(labels):B[i][i+1]=B[i+1][i]=-co[m]
  if triangle:B[0][-1]=B[-1][0]=-h
  assert all(sign(x)>0 for x in kernel)
  assert all(sum((B[i][j]*kernel[j]for j in range(n)),z)==z for i in range(n))
  guards.append(dict(kind=kind,ordered_support=path,labels=labels,complete_triangle=triangle,
    test_pair=sorted(pair),test_label=test,positive_kernel=[enc(x)for x in kernel],
    gram=[[enc(x)for x in row]for row in B],**face(path)))
 for i,j in combinations(range(4),2):
  for a,b in product(S[2*i:2*i+2],S[2*j:2*j+2]):
   aa=next(x for x in S[2*i:2*i+2]if x!=a);bb=next(x for x in S[2*j:2*j+2]if x!=b)
   guard([aa,a,b,bb],[4,3,4],[o,R,R,o],'cross_block_right',[a,b],3)
 guard([11,13,17],[3,3],[o,o,o],'mutual_anchor_right',[11,17],3,True)
 kernel=[2*o,4*o,3*R,2*R,R]
 guard([3,14,11,13,17],[3,4,3,3],kernel,'third_contact_upper',[11,14],4)
 guard([3,14,17,13,11],[3,4,3,3],kernel,'third_contact_upper',[14,17],4)
 assert len(guards)==27
 budgets=[dict(leaf=leaf,removed_own_root=drop,retained_own_root=keep,
   **face((set(S)-{drop})|{leaf}))for leaf,drop,keep in [(11,12,13),(17,12,13),(14,3,4)]]
 assert all(len(row['support'])==8 for row in budgets)
 assert QQ(1,4)+2*QQ(1,4)*2>1 and QQ(1,4)+QQ(1,2)*2>1
 # The same retained-own cost rules out any contact to an internal>=6 block.
 assert QQ(1,4)+4*QQ(1,4)>1
 high_third=face([12,13,14]);assert o-K.convert(QQ(3,4))-h*h==z
 assert all(frozenset([12,13,a])in missing for a in (11,17))
 # Exact all-parameter negative-plane identity for first internal>=6.
 a,b=symbols('a b',nonnegative=True)
 lhs=(1+b*b)*(2+a*a)-(1+a*b)**2
 assert expand(lhs-(1+(a-b)**2+b*b))==0
 plane=dict(first_own_projection_cost='X=D_A*(1+c_A)/2>=2',
   shared_anchor_sum_projection_cost='D_last>4',
   residual_norms=['-(K+H)','-(L+F)'],parameters='K=X-1>=1,L=D_last-2>2,F>=0,H>=0',
   mutual_absolute_bound='1+sqrt(F*H)',
   determinant_lower_identity='(1+H)*(2+F)-(1+sqrt(F*H))^2=1+(sqrt(F)-sqrt(H))^2+H',
   determinant_strict_lower_bound='1',
   ambient_vectors=['r14','r11+r17'],
   justification='Both vectors are right to the other own block. Remaining projections lie in the positive middle-block space; Cauchy bounds their product. The mutual ambient product has magnitude<=1 because each third contact is right orpi/3. Negative diagonals and positive determinant give two negative directions, impossible in a Lorentz space.')
 # Build NEW fallback matrices with every noncore contact completely absent.
 base={(3,14):3,(4,14):3,(11,12):2,(11,13):3,(12,17):2,(13,17):3,(12,13):4}
 cache={};arithmetic_checks=0
 def component_failure(nodes,lab):
  nonlocal arithmetic_checks
  nodes=sorted(nodes);key=(len(nodes),tuple(lab.get((a,b),2)for a,b in combinations(nodes,2)))
  if key not in cache:
   B=[[o if a==b else-co[lab.get(tuple(sorted((a,b))),2)]for b in nodes]for a in nodes]
   piv=[];answer=None
   for k in range(len(nodes)):
    q=B[k][k];piv.append(enc(q))
    if sign(q)<=0:answer=dict(prefix_order=k+1,pivots=piv);break
    for i in range(k+1,len(nodes)):
     for j in range(k+1,len(nodes)):B[i][j]-=B[i][k]*B[k][j]/q
   cache[key]=answer;arithmetic_checks+=1
  answer=cache[key]
  if answer is None:return None
  return dict(support=nodes[:answer['prefix_order']],component_order=nodes,
    component_pair_labels=list(key[1]),**answer)
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
    if ans is not None:
     vi=next(i for i,v in enumerate(VS)if set(f)<=v)
     return dict(actual_vertex=vi,containing_core_face=list(f),**ans)
  return None
 expected={tuple(block)+(ia,ib,ic,ma,mb)for block in product((4,5),repeat=3)
   for ia,ib,ic in product([None]+list(range(6)),[None]+list(range(6)),[None,2,3,4,5])
   for ma,mb in product((2,3),repeat=2)}
 keys=set();new_fail=[];positive=0
 for idx,row in enumerate(old['cases']):
  block=tuple(row['block_labels']);ia,ib,ic=row['foreign_contacts'];ma,mb=row['third_labels']
  key=block+(ia,ib,ic,ma,mb);assert key not in keys;keys.add(key)
  slope=decode(row['slope']);bound=decode(row['boundary_D4']);intercept=decode(row['intercept'])
  assert sign(slope)>0 and bound==intercept+4*slope
  if row['status']=='STRICT_POSITIVE_FOR_ALL_D_GREATER_THAN_FOUR':
   assert sign(bound)>0;positive+=1;continue
  assert row['status']=='NONPOSITIVE_ACTUAL_VERTEX_AT_RELAXED_LAST_INTERNAL_FOUR'
  lab=base|{tuple(S[2*i:2*i+2]):block[i]for i in range(3)}|{(11,14):ma,(14,17):mb}
  for leaf,extra in zip((11,17,14),(ia,ib,ic)):
   if extra is not None:lab[tuple(sorted((leaf,S[extra])))]=3
  assert all(set(p)<=set(CORE)for p in lab)
  witness=failure(lab);assert witness is not None,(idx,key)
  assert set(witness['support'])<=set(CORE)
  new_fail.append(dict(old_case_index=idx,block_labels=list(block),foreign_contacts=[ia,ib,ic],
    third_labels=[ma,mb],relaxed_core_labels=sorted([list(p)+[m]for p,m in lab.items()]),witness=witness))
 assert keys==expected and len(keys)==7840 and positive==7392 and len(new_fail)==448
 swap=list(range(18));swap[12],swap[13]=13,12
 assert {tuple(sorted(swap[i]for i in v))for v in V}==set(map(tuple,V))
 assert {frozenset(swap[i]for i in m)for m in missing}==missing
 fixed=sorted([list(p)+[m]for p,m in base.items()if p!=(12,13)])
 bounds=[S[2*i:2*i+2]+[4 if i<3 else 7]for i in range(4)]
 families=[]
 for reverse in (False,True):
  f=sorted([sorted((swap[a],swap[b]))+[m]for a,b,m in fixed])if reverse else fixed
  families.append(dict(id='core_shared_single_tail'+('_swapped'if reverse else''),
    origin='core_only_shared_single_all_internal_labels',fixed=f,bounds=bounds))
 paths=[GEOM,OLD,REVIEW,'scripts/fresh_h9_star_zero_20261014_own_three.py',f'scripts/{TAG}.py']
 return dict(status='PASS_CORE_SHARED_SINGLE_ALL_INTERNAL_LABELS',case=0,
  sources={p:sha(p)for p in paths},rank_core=CORE,core_positive_kernel_guards=guards,
  projection_budgets=budgets,high_last_third_right=high_third,
  high_foreign_reduction='An internal>=6 block has inverse diagonal>=4. Any nonright foreign contact costs>=1 plus the retained own cost1/4, so all core leaves are right to a high middle block; its internal label drops out of the residual matrix. A high first block similarly forces both anchors right to it.',
  high_first_negative_plane=plane,complete_maximal_core_faces=[list(f)for f in maximal],
  inherited_affine_cases=len(keys),inherited_strict_positive_cases=positive,
  new_core_fallbacks=new_fail,new_core_fallback_count=len(new_fail),
  fresh_ldl_matrix_types=arithmetic_checks,ldl_nonpositive_support_sizes=dict(Counter(len(r['witness']['support'])for r in new_fail)),
  literal_last_root_swap=swap,families=families,
  low_first_reduction='If first internal4or5, replace each middle internal>=6 by4 after proving all core contacts to that block right. The residual equals an inherited finite-block matrix. New fallback witnesses hold at last internal4 and any replaced middle4; increasing those internal labels and adding other acute contacts preserves nonpositivity.',
  tail='Last labelm>=7 gives D=1/sin^2(pi/m)>4, with no upper cutoff or coefficient-field assumption.',
  verification_scope='Fresh exact core-only domain and448 fallback witnesses; inherited separately reviewed affine determinant bank. New all-parameter high-first negative-plane proof uses only Cauchy and the displayed exact square identity.',
  scope='Star0: first3internal labels>=4,last>=7; leaf14own3/3; leaves11and17 share one ownpi/3 contact and are right to the other last root. Every noncore own contact and all remaining finite contacts/dotted distances free. Both last-root orientations. Not a whole-incidence exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_core_fallback_count'],d['fresh_ldl_matrix_types'])
