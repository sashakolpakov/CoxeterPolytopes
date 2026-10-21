#!/usr/bin/env python3
"""Independent exact induced-subcore search; no producer canonicalization.
Uses the frozen independent Fraction field only for the458 full inertias.
"""
import argparse,json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
import fresh_sqrt5_quad_contacts_20261015_review as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_sealed_subcore_closure_20261015_review'
OPP={frozenset((0,2)),frozenset((1,3))}
BLACKS=[p for p in permutations(range(4))if{frozenset((p[0],p[2])),frozenset((p[1],p[3]))}==OPP]
assert len(BLACKS)==8

def matrix(row):
 n=row['nodes'];M=[[2]*n for _ in range(n)]
 for(i,j),m in zip(combinations(range(n),2),row['labels']):M[i][j]=M[j][i]=m
 return M

def invariant(M):
 n=len(M)
 def sig(i):return(tuple(sorted(M[i][j]for j in range(4)if j!=i)),tuple(sorted(M[i][j]for j in range(4,n)if j!=i)))
 return(n,tuple(sorted(sig(i)for i in range(4))),tuple(sorted(sig(i)for i in range(4,n))))

def isomorphism(M,T):
 n=len(M)
 if len(T)!=n:return None
 for p in BLACKS:
  if any(M[i][j]!=T[p[i]][p[j]]for i,j in combinations(range(4),2)):continue
  choices={i:[j for j in range(4,n)if all(M[i][b]==T[j][p[b]]for b in range(4))and sorted(M[i][z]for z in range(4,n)if z!=i)==sorted(T[j][z]for z in range(4,n)if z!=j)]for i in range(4,n)}
  if any(not v for v in choices.values()):continue
  order=sorted(choices,key=lambda i:(len(choices[i]),i));mp=dict(enumerate(p))
  def rec(k,used):
   if k==len(order):return [mp[i]for i in range(n)]
   i=order[k]
   for j in choices[i]:
    if j in used or any(M[i][z]!=T[j][v]for z,v in mp.items()):continue
    mp[i]=j;ans=rec(k+1,used|{j})
    if ans is not None:return ans
    del mp[i]
   return None
  ans=rec(0,set(p))
  if ans is not None:return ans
 return None

def check_map(M,T,mp):
 n=len(M);assert sorted(mp)==list(range(n))
 assert all((i<4)==(mp[i]<4)for i in range(n))
 assert{frozenset((mp[0],mp[2])),frozenset((mp[1],mp[3]))}==OPP
 assert all(M[i][j]==T[mp[i]][mp[j]]for i,j in combinations(range(n),2))

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def run():
 primary=read('results/fresh_sqrt5_sealed_subcore_closure_20261015.json');old=read('results/fresh_sqrt5_quad_contacts_20261015.json');base=read('results/fresh_sqrt5_quad_exhaustion_20261015.json')
 pool=defaultdict(list);sizes=set();targets={}
 for r in old['records']:
  if r['status']=='UNPAID_FINITE_CONTACT_DOMAIN'and r['singular']and not r['accepted']:
   row=base['labelled_union_supports'][r['union_index']];M=matrix(row);pool[invariant(M)].append((r['class_index'],M));targets[r['class_index']]=M;sizes.add(len(row['white']))
 assert len(targets)==112
 subsets=0;counts=Counter();records=[];trials=0;maps_checked=0
 for r in old['records']:
  if r['status']!='UNPAID_FINITE_CONTACT_DOMAIN'or r['singular']or not r['accepted']:continue
  source=base['labelled_union_supports'][r['union_index']]
  for ti,t in enumerate(r['trials']):
   if t['status']!='PASS_EXTRA_WHITE_NECESSARY_TESTS':continue
   trials+=1;ext=K.extension(source,t['contacts']);M=matrix(ext);found=[]
   inertia=K.inertia(K.tree_scaled(ext['nodes'],K.edges(ext)));assert inertia[1]==1
   singular=inertia[2]>0
   for size in sorted(sizes):
    for white in combinations(ext['white'],size):
     ids=list(range(4))+list(white);S=[[M[i][j]for j in ids]for i in ids];subsets+=1
     for ci,T in pool.get(invariant(S),[]):
      mp=isomorphism(S,T)
      if mp is not None:check_map(S,T,mp);found.append((ci,tuple(ids)));maps_checked+=1
   ex=next(z for z in primary['records']if z['source_class_index']==r['class_index']and z['trial_index']==ti)
   expected=[]
   for z in ex['sealed_subcores']:
    ids=z['nodes_in_extension'];S=[[M[i][j]for j in ids]for i in ids];check_map(S,targets[z['class_index']],z['to_sealed_representative']);expected.append((z['class_index'],tuple(ids)))
   assert sorted(found)==sorted(expected)and singular==ex['singular']
   assert not found or singular
   status='PAID_BY_REVIEWED_SEALED_SUBCORE'if found else'OPEN_SINGULAR'if singular else'OPEN_NONSINGULAR'
   assert status==ex['status'];counts[status]+=1
   records.append(dict(source_class_index=r['class_index'],trial_index=ti,full_inertia=inertia,status=status,contained=[dict(class_index=ci,nodes=list(ids))for ci,ids in sorted(found)]))
 assert trials==458 and subsets==15922 and dict(counts)==primary['counts']
 assert counts=={'PAID_BY_REVIEWED_SEALED_SUBCORE':86,'OPEN_SINGULAR':61,'OPEN_NONSINGULAR':311}
 paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_sealed_subcore_closure_20261015.py','results/fresh_sqrt5_sealed_subcore_closure_20261015.json','scripts/fresh_sqrt5_quad_contacts_20261015_review.py','results/fresh_sqrt5_quad_contacts_20261015.json','results/fresh_sqrt5_quad_exhaustion_20261015.json','results/fresh_sqrt5_forest_sensitive_weights_20261015.json','results/fresh_sqrt5_forest_sensitive_weights_20261015_review.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_INDEPENDENT_SEALED_SUBCORE_CONTAINMENT',sources={p:sha(p)for p in paths},extensions=trials,white_subsets=subsets,maps_checked=maps_checked,counts=dict(sorted(counts.items())),records=records,method='Independent colored complete-graph isomorphism backtracking, all8 opposite-pair-preserving black permutations and every compatible white bijection; no producer canonicalization. All458 physical full-Gram inertias independently recomputed with the frozen independent Fraction field.',scope='Finite first-extension containment and inherited arbitrary-extension payment of found sealed subcores only. Unmatched61/311 and24 dotted fronts remain open; no arbitrary-depth closure or dimension theorem.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'],out['white_subsets'])
