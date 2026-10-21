#!/usr/bin/env python3
"""Separate finite-diagram classification and modular-rank replay."""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
from sympy import GF
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_surgery_20261014'
PAIRS=list(combinations(range(14),2));INDEX={p:i for i,p in enumerate(PAIRS)}
def support(m):return tuple(i for i in range(14)if m>>i&1)
@lru_cache(None)
def spherical(code):
 A=[{}for _ in range(9)]
 for (i,j),m in zip(combinations(range(9),2),code):
  if m>2:A[i][j]=A[j][i]=m
 left=set(range(9))
 while left:
  start=min(left);cc={start};todo=[start]
  for i in todo:
   for j in A[i]:
    if j not in cc:cc.add(j);todo.append(j)
  left-=cc;n=len(cc)
  if n<=2:continue
  if sum(len(A[i])for i in cc)!=2*(n-1):return False
  degrees=Counter(len(A[i])for i in cc)
  if max(degrees)<=2:
   start=next(i for i in cc if len(A[i])==1);seq=[];prev=None
   while True:
    nxt=[j for j in A[start]if j!=prev]
    if not nxt:break
    j=nxt[0];seq.append(A[start][j]);prev,start=start,j
   seq=min(tuple(seq),tuple(seq[::-1]))
   if not(all(m==3 for m in seq)or seq==tuple([3]*(n-2)+[4])or(n==3 and seq==(3,5))or(n==4 and seq in[(3,3,5),(3,4,3)])):return False
  else:
   if degrees[3]!=1 or max(degrees)>3 or any(m!=3 for i in cc for m in A[i].values()):return False
   branch=next(i for i in cc if len(A[i])==3);arms=[]
   for j in A[branch]:
    prev,cur,length=branch,j,1
    while len(A[cur])==2:
     nxt=next(k for k in A[cur]if k!=prev);prev,cur=cur,nxt;length+=1
    arms.append(length)
   if tuple(sorted(arms))not in[(1,1,n-3),(1,2,2),(1,2,3),(1,2,4)]:return False
 return True
def run():
 paths={s:f'results/{TAG}_{s}.json'for s in['incidence','dfs_0_5','metric_scout','rank_search']};D={s:json.loads((ROOT/p).read_text())for s,p in paths.items()}
 for d in D.values():
  for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 p=241;roots={d:next(a for a in range(p)if a*a%p==d)for d in[2,3,5]};cos={2:0,3:pow(2,-1,p),4:roots[2]*pow(2,-1,p)%p,5:(1+roots[5])*pow(4,-1,p)%p,6:roots[3]*pow(2,-1,p)%p};field=GF(p)
 records=[];det_cache={};total=0
 for base in D['dfs_0_5']['cases']:
  if 'label_vector'not in base:continue
  row=D['incidence']['cases'][base['case']];V=[support(m)for m in row['canonical_vertices']];vv=[tuple(INDEX[a]for a in combinations(v,2))for v in V];l=base['label_vector'][:];options=[]
  for pair in base['unconstrained_meeting_pairs']:
   k=INDEX[tuple(pair)];aa=[]
   for m in range(2,8):
    l[k]=m
    if all(spherical(tuple(l[j]for j in a))for a in vv):aa.append(m)
   l[k]=2;assert aa and aa[0]==2
   if len(aa)>1:options.append((k,aa))
  expected={tuple(r['labels'])for r in next(r for r in D['metric_scout']['cases']if r['case']==base['case'])['fits']};actual=set();witnesses=[]
  dotted={tuple(a)for a in base['dotted_pairs']}
  for choices in product(*(aa for k,aa in options)):
   for(k,aa),m in zip(options,choices):l[k]=m
   if not all(spherical(tuple(l[j]for j in a))for a in vv):continue
   actual.add(tuple(l));total+=1
   for a in combinations(range(14),11):
    pp=list(combinations(a,2))
    if any(pair in dotted or l[INDEX[pair]]>=7 for pair in pp):continue
    code=tuple(l[INDEX[pair]]for pair in pp)
    if code not in det_cache:
     mat=[[field(1 if i==j else -cos[l[INDEX[tuple(sorted((i,j)))]]])for j in a]for i in a]
     det_cache[code]=int(DomainMatrix(mat,(11,11),field).det())%p
    residue=det_cache[code]
    if residue:
     witnesses.append(dict(labels=l[:],support=a,residue=residue));break
   else:raise AssertionError('An extension has no modular rank witness')
  assert actual==expected
  records.append(dict(case=base['case'],free_pair_options=[dict(pair=PAIRS[k],labels=aa)for k,aa in options],finite_extensions=len(actual),rank_witnesses=witnesses))
 # Replay every archived rank witness from the bounded DFS, independently
 # of its Gaussian-elimination implementation.
 search_witnesses=search_multiplicity=0
 for row in D['rank_search']['cases']:
  missing=[support(m)for m in D['incidence']['cases'][row['case']]['canonical_minimal_nonfaces']if m.bit_count()==2]
  for w in row['rank_witnesses']:
   a=w['support'];code=tuple(w['labels']);assert len(a)==11 and len(set(a))==11
   assert not any(set(b)<=set(a)for b in missing)and len(code)==55 and set(code)<=set(cos)
   if code not in det_cache:
    entries=dict(zip(combinations(range(11),2),code));mat=[[field(1 if i==j else -cos[entries[tuple(sorted((i,j)))]])for j in range(11)]for i in range(11)]
    det_cache[code]=int(DomainMatrix(mat,(11,11),field).det())%p
   assert det_cache[code]==w['residue']!=0
   search_witnesses+=1;search_multiplicity+=w['multiplicity']
 assert search_multiplicity==sum(r['rank_rejections']for r in D['rank_search']['cases'])
 # The separate search's only rank survivor has two orthogonal negative blocks.
 r=next(r for r in D['rank_search']['cases']if 'label_vector'in r);assert r['case']==318
 M=[support(m)for m in D['incidence']['cases'][318]['canonical_minimal_nonfaces']];ds={tuple(a)for a in r['dotted_pairs']};l=r['label_vector']
 a,b=next((a,b)for a,b in combinations(M,2)if not set(a)&set(b)and all(tuple(sorted((i,j)))not in ds and l[INDEX[tuple(sorted((i,j)))]]==2 for i in a for j in b))
 assert len(a)==len(b)==3
 for v in[a,b]:
  labels=[l[INDEX[pair]]for pair in combinations(v,2)];assert sorted(labels)==[2,3,7]
 files=list(paths.values())+[str(Path(__file__).resolve().relative_to(ROOT))]
 return dict(status='PASS_SEPARATE_FINITE_CATEGORY_AND_RANK_REPLAY',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},prime=p,square_roots=roots,cosine_residues=cos,finite_extensions=total,distinct_minor_codes=len(det_cache),bounded_search_witnesses_checked=search_witnesses,bounded_search_witness_multiplicity=search_multiplicity,cases=records,additional_signature_obstruction=dict(case=318,orthogonal_Lanner_blocks=[a,b],each_labels=[2,3,'every integer >=7'],negative_index_at_least=2),scope='Excludes the three saved partial Lanner labelings in every meeting-edge category completion and all dotted distances. Also excludes the one saved rank-surviving full category on case318. This does not exclude their incidence types, whose other Lanner assignments remain open. Finite classification and modular determinants are independently implemented; source incidence and Lanner-table provenance remain upstream. All2358 bounded-search minor witnesses are replayed, but the DFS coverage itself is not independently rerun here.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();o=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}_rank_independent.json'
 if a.verify:assert o==json.loads(path.read_text())
 else:path.write_text(json.dumps(o,indent=2)+'\n')
 print(o['status'],o['finite_extensions'],o['additional_signature_obstruction'])
