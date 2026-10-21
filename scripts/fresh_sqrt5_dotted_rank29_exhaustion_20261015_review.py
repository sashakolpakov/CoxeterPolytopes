#!/usr/bin/env python3
"""Independent complete 438-state dotted queue, with fresh finite-tree Schur tests.
Only independent tuple-Fraction arithmetic/parameter review is imported.
"""
import argparse,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_dotted_parameters_20261015_review as P
R=P.R;F=P.F;Z,ONE,TWO=P.Z,P.ONE,P.TWO
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_sqrt5_dotted_rank29_exhaustion_20261015';TAG=BASE+'_review'
CORNERS=((0,1),(0,3),(2,1),(2,3));ALPHABET=(10,6,5,4,3,2)
SQUARES={m:R.div(v,P.q(4))for m,v in R.SQUARES.items()}
def pair(i,j):return tuple(sorted((i,j)))
def extension(row,c):
 n=row['nodes'];g=P.edges(row);g.update({(i,n):m for i,m in enumerate(c)})
 return dict(nodes=n+1,white=list(range(4,n+1)),labels=[g[p]for p in combinations(range(n+1),2)])

@lru_cache(maxsize=200000)
def finite(n,word):
 # Every positive definite Coxeter diagram is a forest: a cycle with all
 # cosine magnitudes >=1/2 has a nonpositive all-ones quadratic witness.
 adj=[{}for _ in range(n)]
 for (i,j),m in zip(combinations(range(n),2),word):
  if m==0:return False
  if m!=2:adj[i][j]=m;adj[j][i]=m
 d=[ONE]*n;todo=[i for i in range(n)if len(adj[i])<=1];alive=set(range(n));done=0
 while todo:
  i=todo.pop()
  if i not in alive:continue
  if R.sign(d[i])<=0:return False
  alive.remove(i);done+=1
  if adj[i]:
   j,m=next(iter(adj[i].items()));d[j]=R.sub(d[j],R.div(SQUARES[m],d[i]));del adj[j][i]
   if len(adj[j])<=1:todo.append(j)
 return done==n

def domains(row):
 n=row['nodes'];g=P.edges(row);c=[2]*n;out=[];visits=0;white=list(range(4,n+1));order=list(reversed(range(n)))
 corner_ids=[sorted(white+[a,b])for a,b in CORNERS]
 def good():
  for (a,b),ids in zip(CORNERS,corner_ids):
   if g.get(pair(a,b),2)!=2 or any(g.get(pair(a,w),2)!=2 and g.get(pair(b,w),2)!=2 for w in white):return False
   if not finite(len(ids),tuple(g.get(pair(i,j),2)for i,j in combinations(ids,2))):return False
  return True
 def walk(k):
  nonlocal visits
  visits+=1
  if k==n:
   if any(m!=2 for m in c):out.append(tuple(c))
   return
  i=order[k]
  for m in ALPHABET:
   c[i]=m;g[i,n]=m
   if good():walk(k+1)
  del g[i,n];c[i]=2
 walk(0);assert len(out)==len(set(out));return out,visits

def h4_paid(row):
 g=P.edges(row);white=row['white'];out=[]
 for a,b in CORNERS:
  ids=white+[a,b];adj={i:{j:g[pair(i,j)]for j in ids if i!=j and g[pair(i,j)]!=2}for i in ids};seen={a};todo=[a]
  for i in todo:
   for j in adj[i]:
    if j not in seen:seen.add(j);todo.append(j)
  if b not in seen or len(seen)!=4:continue
  if len(adj[a])!=1 or len(adj[b])!=1 or sorted(len(adj[i])for i in seen)!=[1,1,2,2]:continue
  labels=sorted(g[pair(i,j)]for i,j in combinations(seen,2)if g[pair(i,j)]!=2)
  if labels==[3,3,5]:out.append(tuple(sorted((a,b))))
 return out

def sparse_inertia(M,white,conjugate=False):
 # Eliminate the white forest first; all its pivots must be positive in
 # either embedding. Sparse congruences avoid the producer's dense engine.
 A=[[R.conj(x)if conjugate else x for x in r]for r in M];alive=set(range(len(A)));W=set(white);positive=0
 while W:
  i=min(W,key=lambda i:sum(A[i][j]!=Z for j in W if j!=i));assert sum(A[i][j]!=Z for j in W if j!=i)<=1
  d=A[i][i];assert R.sign(d)>0;positive+=1;I=[j for j in alive if j!=i and A[i][j]!=Z]
  for k,j in enumerate(I):
   v=R.div(A[i][j],d)
   for l in I[k:]:A[j][l]=R.sub(A[j][l],R.mul(v,A[i][l]));A[l][j]=A[j][l]
  alive.remove(i);W.remove(i)
 B=[[A[i][j]for j in sorted(alive)]for i in sorted(alive)];ans=R.inertia(B);ans[0]+=positive;return ans

def check_map(row,target,mp):
 n=row['nodes'];assert target['nodes']==n and sorted(mp)==list(range(n))
 assert set(mp[:4])==set(range(4))and {mp[0],mp[2]}=={0,2}and {mp[1],mp[3]}=={1,3}
 assert set(mp[i]for i in row['white'])==set(target['white'])
 g=P.edges(row);h=P.edges(target)
 for (i,j),m in g.items():assert h[pair(mp[i],mp[j])]==m
 return n*(n-1)//2

def run():
 data=P.read('results/'+BASE+'.json');params=P.read('results/'+P.BASE+'.json');assert data['status']=='COMPLETE_DOTTED_ARITHMETIC_ATTACHMENT_QUEUE_THROUGH_DIMENSION29'
 assert data['dimension_bound']==29 and data['maximum_star_nodes']==31 and data['processed']==data['marked_states']==len(data['records'])==438
 reached=set();maps=0;comparisons=0;rootkeys=set()
 for seed in data['seeds']:
  ri,ci=seed['record'],seed['accepted'];assert(ri,ci)not in rootkeys;rootkeys.add((ri,ci));rec=params['records'][ri];accepted=rec['accepted'][ci];target=data['records'][seed['index']]
  trial=rec['trials'][accepted['trial_index']]
  assert P.dec(accepted['ratio'])==R.div(P.dec(trial['twice_gram_cycle']),P.FOUR)
  assert accepted['singular']==(trial['status']=='ACCEPT_SINGULAR')
  assert accepted['ratio']==target['ratio'];comparisons+=check_map(rec['row'],target['row'],seed['map']);maps+=1;reached.add(seed['index'])
 assert rootkeys=={(i,j)for i,r in enumerate(params['records'])for j in range(len(r['accepted']))}and len(rootkeys)==32
 counts=Counter();tc=Counter();records=[];obligations=[];dfs=0;cycle_products=0;totaldomains=0;inertias=0
 for index,rec in enumerate(data['records']):
  assert rec['index']==index and index in reached;row=rec['row'];n=row['nodes'];assert row['white']==list(range(4,n))and n<=31
  assert [p for p,m in P.edges(row).items()if m==0]==[(0,2)]
  g=P.edges(row)
  for a,b in CORNERS:
   ids=sorted(row['white']+[a,b]);assert finite(len(ids),tuple(g[pair(i,j)]for i,j in combinations(ids,2)))
   assert g[pair(a,b)]==2 and not any(g[pair(a,w)]!=2 and g[pair(b,w)]!=2 for w in row['white'])
  assert not h4_paid(row)
  ratio=P.dec(rec['ratio']);M,error=P.matrix(row,ratio);assert error is None
  physical=sparse_inertia(M,row['white']);conjugate=sparse_inertia(M,row['white'],True);inertias+=2;saved=rec['result']
  assert physical==saved['physical']and conjugate==saved['conjugate'];assert physical[1]==1 and physical[2]in(0,1)and conjugate==[n-physical[2],0,physical[2]]
  singular=bool(physical[2]);assert singular==saved['singular'];cc=P.cycles(M);assert all(ok for v,ok in cc.values());cycle_products+=len(cc)
  if singular:obligations.append(index)
  local=Counter();accept=[]
  if n==31:
   assert saved['trials']==saved['accepted']==saved['children']==[]
   status='DIMENSION29_SINGULAR_WEIGHT_OBLIGATION'if singular else'REJECT_RANK_EXCEEDS_DIMENSION29';visits=0;dom=[]
  else:
   dom,visits=domains(row);dfs+=visits;totaldomains+=len(dom);archive={tuple(t['contacts']):t for t in saved['trials']}
   assert len(archive)==len(saved['trials'])and set(dom)==set(archive),(index,len(dom),len(archive))
   for c in sorted(dom):
    ext=extension(row,c);t=archive[c]
    if h4_paid(ext):ts='PAID_H4'
    else:
     G,error=P.matrix(ext,ratio)
     if error:ts=error
     else:
      pi=sparse_inertia(G,ext['white']);co=sparse_inertia(G,ext['white'],True);inertias+=2
      assert pi==t['physical']and co==t['conjugate'],(index,c,pi,co,t)
      if pi[1]>1:ts='REJECT_TWO_NEGATIVE'
      elif co[1]:ts='REJECT_CONJUGATE'
      else:
       assert pi[1]==1 and pi[2]in(0,1)and co==[n+1-pi[2],0,pi[2]]
       if singular:assert pi[2]==1
       cy=P.cycles(G);cycle_products+=len(cy)
       assert all(ok for v,ok in cy.values());ts='ACCEPT_UNPAID';accept.append(c)
    assert ts==t['status'],(index,c,ts,t['status']);local[ts]+=1;tc[ts]+=1
   assert set(accept)==set(map(tuple,saved['accepted']))and len(accept)==len(saved['accepted'])
   status=('EXPANDED_SINGULAR'if singular else'EXPANDED_NONSINGULAR')if accept else('SEALED_SINGULAR_WEIGHT_OBLIGATION'if singular else'REJECT_SEALED_NONSINGULAR_RANK')
   children={tuple(c['contacts']):c for c in saved['children']};assert len(children)==len(saved['children'])and set(children)==set(accept)
   for c in accept:
    child=children[c];target=data['records'][child['index']];assert target['ratio']==rec['ratio'];comparisons+=check_map(extension(row,c),target['row'],child['map']);maps+=1;reached.add(child['index'])
  assert status==saved['status'];counts[status]+=1
  records.append(dict(index=index,nodes=n,singular=singular,status=status,domain_size=len(dom),independent_DFS_nodes=visits,domain_sha256=sha256(json.dumps(sorted(dom),separators=(',',':')).encode()).hexdigest(),trial_counts=dict(sorted(local.items())),accepted_count=len(accept)))
  if index%25==0:print('independent dotted queue',index,'profiles',totaldomains,'maps',maps,flush=True)
 assert reached==set(range(438))and obligations==data['weight_obligations']and len(obligations)==264
 assert dict(counts)==data['counts']and dict(tc)==data['trial_counts']and maps==454
 paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',f'results/{BASE}.md',f'scripts/{P.TAG}.py',f'results/{P.TAG}.json','results/fresh_sqrt5_quad_contacts_20261015_review.json','results/fresh_sqrt5_triangle_weights_20261015_review.json','results/arithmetic_dimension_bounds_and_effective_high_dimension_search.md']
 return dict(status='PASS_INDEPENDENT_COMPLETE_DOTTED_ATTACHMENT_QUEUE_DIMENSION29',sources=P.manifest(paths),marked_states=438,parameter_seeds=32,complete_domains=423,profiles=totaldomains,independent_DFS_nodes=dfs,exact_inertia_checks=inertias,cyclic_products=cycle_products,seed_and_child_maps=maps,label_comparisons=comparisons,weight_obligations=len(obligations),counts=dict(counts),trial_counts=dict(tc),records=records,method='No producer imports. Reversed complete contact DFS; independent exact finite-forest leaf Schur classifier; minimal radical gauge; sparse full-matrix white elimination plus remaining4x4 inertia; elementary cycles; literal seed/child full-label bijections.',scope='Complete necessary attachment queue of the32 parameter seeds for ambient dimensions at most29. Every singular state is retained as a zero-further-touch payment obligation; no payment, geometric realization, all-rank sealing, or improved dimension theorem is asserted. Short paths and entireH4 payments remain inherited dependencies.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==P.read(str(p.relative_to(ROOT)))
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'],out['trial_counts'])
