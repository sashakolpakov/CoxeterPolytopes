#!/usr/bin/env python3
"""Independent51 singular domains and extended forest-weight review.
Physical positive-definiteness uses the separate tuple/Fraction radical engine,
not the producer's finite-type predicate or field class. No producer imports.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_quad_filters_20261015_review as A
import fresh_sqrt5_forest_sensitive_weights_20261015_review as W
import fresh_sqrt5_sealed_subcore_closure_20261015_review as S
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_extended_forest_weights_20261015_review'

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def edges(r):return dict(zip(combinations(range(r['nodes']),2),r['labels']))
def extend(r,c):
 n=r['nodes'];g=edges(r);g.update({(i,n):m for i,m in enumerate(c)});return dict(nodes=n+1,white=list(range(4,n+1)),labels=[g[p]for p in combinations(range(n+1),2)])
def conj(M):return[[A.conj(x)for x in r]for r in M]
def gram(r):
 c=A.coloring(r);return None if c is None else A.matrix(r,c)

@lru_cache(maxsize=100000)
def pd(n,word):
 # Every Coxeter cycle is non-PD: its norm-two cycle submatrix has
 # all-one quadratic value <=0. Otherwise use exact matrix inertia.
 parent=list(range(n))
 def find(i):
  while i!=parent[i]:i=parent[i]
  return i
 for(i,j),m in zip(combinations(range(n),2),word):
  if m==2:continue
  a,b=find(i),find(j)
  if a==b:return False
  parent[a]=b
 M=gram(dict(nodes=n,labels=word));assert M is not None
 return A.inertia(M)==[n,0,0]
def short(r):
 g=edges(r)
 for a,b in product((0,2),(1,3)):
  H=W.graph(r['white']+[a,b],g);d=W.distances(H,a).get(b)
  if d is not None and d<=2:return True
 return False

def domains(r):
 n=r['nodes'];g=edges(r);c=[2]*n;out=[]
 I=[tuple(sorted(r['white']+[n,a,b]))for a,b in product((0,2),(1,3))]
 def valid():
  for ids in I:
   a,b=[i for i in ids if i<4]
   if g.get(tuple(sorted((a,b))),2)!=2 or any(g.get(tuple(sorted((a,w))),2)!=2 and g.get(tuple(sorted((b,w))),2)!=2 for w in r['white']+[n]):return False
   if not pd(len(ids),tuple(g.get(tuple(sorted((i,j))),2)for i,j in combinations(ids,2))):return False
  return True
 def rec(k):
  if k==n:
   if any(m!=2 for m in c):out.append(tuple(c))
   return
  i=n-1-k
  for m in(10,6,5,4,3,2):
   c[i]=m;g[i,n]=m
   if valid():rec(k+1)
  del g[i,n];c[i]=2
 assert not short(r);rec(0);return out

def positive_kernel(M):
 n=len(M);B=[r[:]for r in M];pivs=[];k=0
 for j in range(n):
  p=next((i for i in range(k,n)if B[i][j]!=A.Z),None)
  if p is None:continue
  B[k],B[p]=B[p],B[k];z=B[k][j];B[k]=[A.div(x,z)for x in B[k]]
  for i in range(n):
   if i!=k:
    z=B[i][j];B[i]=[A.sub(x,A.mul(z,y))for x,y in zip(B[i],B[k])]
  pivs.append(j);k+=1
 free=[j for j in range(n)if j not in pivs];assert len(free)==1
 v=[A.Z]*n;v[free[0]]=A.ONE
 for i,j in enumerate(pivs):v[j]=A.neg(B[i][free[0]])
 signs=[A.sign(x)for x in v];assert all(s>0 for s in signs)or all(s<0 for s in signs)
 assert all(sum((A.mul(a,b)[0]for a,b in zip(r,v)),F())==0 and sum((A.mul(a,b)[1]for a,b in zip(r,v)),F())==0 for r in M)
 return[list(map(str,x))for x in v]
def affine(r,w):
 g=edges(r);ids=w['nodes'];word=[g[tuple(sorted((i,j)))]for i,j in combinations(ids,2)];assert word==w['labels'];M=gram(dict(nodes=len(ids),labels=word));assert M is not None and A.inertia(M)==[len(ids)-1,0,1];return positive_kernel(M)
def h4paid(r):
 out=[];g=edges(r)
 for a,b in product((0,2),(1,3)):
  cs=W.pieces(W.graph(r['white']+[a,b],g));C=next(c for c in cs if a in c)
  if b in C and W.family(C)=='H4'and len(C[a])==len(C[b])==1:out.append((a,b,tuple(sorted(C))))
 return out

def audit_domain(rec):
 r=rec['row'];n=r['nodes'];M=gram(r);assert M is not None
 assert A.inertia(M)==[n-2,1,1]and A.inertia(conj(M))==[n-1,0,1]
 dom=domains(r);assert set(dom)=={tuple(t['contacts'])for t in rec['trials']}and len(dom)==rec['profiles']
 out=[];accepted=[];counts=Counter()
 for t in rec['trials']:
  ext=extend(r,t['contacts']);status=t['status'];counts[status]+=1
  if status=='BECOMES_PAID_ENTIRE_H4':
   actual=h4paid(ext);assert actual
   assert{(tuple(z['marked']),tuple(z['entire_H4']))for z in t['witness']}=={((a,b),ids)for a,b,ids in actual}
   check=dict(payments=len(actual))
  elif status=='REJECT_AFFINE':check=dict(positive_kernel=affine(ext,t['witness']))
  else:
   G=gram(ext)
   if G is None:assert status=='REJECT_FINITE_CYCLE_FIELD';check={}
   else:
    physical=A.inertia(G);other=A.inertia(conj(G));check=dict(physical_inertia=physical,conjugate_inertia=other)
    if status=='REJECT_PHYSICAL_RELATION_AND_CONJUGATE_KERNEL':assert physical[1]>1 and other[1]>0
    elif status=='REJECT_TWO_NEGATIVE':assert physical[1]>1
    elif status=='REJECT_CONJUGATE_NOT_PSD':assert other[1]>0
    else:assert status=='PASS_UNPAID_EXTRA_WHITE'and physical==[n-1,1,1]and other==[n,0,1];accepted.append(tuple(t['contacts']))
  out.append(dict(contacts=t['contacts'],status=status,**check))
 assert set(accepted)==set(map(tuple,rec['accepted']));return dict(index=rec['index'],profiles=len(dom),trials=out,accepted=[list(c)for c in accepted]),counts

def old_weights(data):
 weights={}
 for name,g in W.models():
  r=next(z for z in data['component_rows']if z['type']==name)
  for p in r['pair_weights']:
   key=W.code(g,p['marked']);value=tuple(map(F,p['weights']))
   if key in weights:assert weights[key]==value
   weights[key]=value
 return weights

def payments(r,weights):
 out=[];g=edges(r)
 for a,b in product((0,2),(1,3)):
  cs=W.pieces(W.graph(r['white']+[a,b],g));C=next(c for c in cs if a in c);flag=any(len(c)==1 for c in cs)
  key=W.code(C,(a,b))if b in C else None;vals=[F(0),F(0)]if key is None else[weights[key][int(flag)],weights[key][1]]
  out.append(dict(marked=[a,b],core_has_a1=flag,weights=list(map(str,vals))))
 return out,[sum(F(c['weights'][i])for c in out)for i in(0,1)]

def run():
 frontier=read('results/fresh_sqrt5_singular_extension_frontier_20261015.json');primary=read('results/fresh_sqrt5_extended_forest_weights_20261015.json');prior=read('results/fresh_sqrt5_sealed_subcore_closure_20261015.json');base=read('results/fresh_sqrt5_quad_exhaustion_20261015.json')['labelled_union_supports'];old=read('results/fresh_sqrt5_forest_sensitive_weights_20261015.json')
 selected=[r for r in prior['records']if r['status']=='OPEN_SINGULAR'];assert len(selected)==61 and len(frontier['records'])==51
 # Independent exact quotient check: every origin is covered, and the51
 # representatives are pairwise nonisomorphic, preserving allblack opposites.
 buckets=defaultdict(list);seen=set()
 for rec in frontier['records']:
  T=S.matrix(rec['row']);iv=S.invariant(T)
  assert all(S.isomorphism(T,U)is None for U in buckets[iv]);buckets[iv].append(T)
  for origin in rec['origins']:
   i=origin['open_index'];assert i not in seen;seen.add(i);r=selected[i];row=extend(base[r['source_union_index']],r['contacts']);M=S.matrix(row)
   assert S.isomorphism(M,T)is not None;S.check_map(M,T,origin['to_representative'])
   assert all(origin[k]==r[k]for k in('source_class_index','source_union_index','trial_index'))
 assert seen==set(range(61))
 weights=old_weights(old);records=[];counts=Counter();trialcounts=Counter();sealed=[];unpaid=[]
 for rec in frontier['records']:
  row,cnt=audit_domain(rec);trialcounts.update(cnt);corners,sums=payments(rec['row'],weights);saved=rec['payment']
  assert sums==list(map(F,saved['sums_without_with_external_a1']))and(min(sums)>=1)==saved['paid']
  for c in corners:
   z=next(z for z in saved['corners']if z['marked']==c['marked']);assert c['weights']==z['weights']and c['core_has_a1']==z['core_has_a1']
  status='OPEN_SINGULAR_WITH_CONTACTS'if rec['accepted']else'SEALED_PAID_BY_FOREST_RULE'if min(sums)>=1 else'SEALED_UNPAID_BY_CURRENT_FOREST_RULE'
  assert status==rec['status'];counts[status]+=1;row.update(status=status,old_payment=list(map(str,sums)));records.append(row)
  if not rec['accepted']:
   sealed.append(rec)
   if min(sums)<1:unpaid.append(rec)
 assert dict(counts)==frontier['counts']and dict(trialcounts)==frontier['trial_counts']
 assert len(sealed)==46 and len(unpaid)==14
 upgrades=set()
 for rec in unpaid:
  g=edges(rec['row'])
  for a,b in product((0,2),(1,3)):
   C=next(c for c in W.pieces(W.graph(rec['row']['white']+[a,b],g))if a in c)
   if b in C:upgrades.add(W.code(C,(a,b)))
 assert len(upgrades)==15 and upgrades=={W.decode(z['key'])for z in primary['upgraded_orbits']}
 for k in upgrades:weights[k]=tuple(max(v,F(1))for v in weights[k])
 face_checks=[]
 inputs=[('old',r['class_index'],base[r['union_index']])for r in old['sealed_quad_rows']]+[('new',r['index'],r['row'])for r in sealed]
 for fam,i,row in inputs:
  corners,sums=payments(row,weights);saved=next(z for z in primary['face_rows']if z['family']==fam and z['index']==i)
  assert min(sums)>=1 and sums==list(map(F,saved['sums_without_with_external_a1']))and corners==saved['corners']
  face_checks.append(dict(family=fam,index=i,sums=list(map(str,sums))))
 components=[]
 for name,g in W.models():
  totals=[sum((weights[W.code(g,p)][flag]for p in combinations(g,2)),F())for flag in(0,1)];saved=next(z for z in primary['component_rows']if z['type']==name)
  assert totals==list(map(F,saved['totals']));components.append(dict(type=name,rank=len(g),totals=list(map(str,totals))))
 dps=[]
 for context in(0,1):
  dp={(0,False):(F(0),())}
  for r in components:
   for n in range(r['rank'],30):
    for flag in(False,True):
     prev=dp.get((n-r['rank'],flag))
     if prev is None:continue
     key=(n,flag or r['type']=='A1');value=(prev[0]+F(r['totals'][context]),tuple(sorted(prev[1]+(r['type'],))))
     if key not in dp or value>dp[key]:dp[key]=value
  dps.append(dp)
 no=dps[0][29,False];yes=dps[1][29,True];assert no[0]==F(163,2)and yes[0]==84
 assert no[0]==F(primary['rank29_without_a1']['maximum'])and yes[0]==F(primary['rank29_with_a1']['maximum'])
 assert primary['no_a1_rank_maxima']==[None if(n,False)not in dps[0]else str(dps[0][n,False][0])for n in range(30)]
 assert primary['with_a1_tail_rank_maxima']==[str(max(v[0]for(k,_),v in dps[1].items()if k==n))for n in range(29)]
 paths=[f'scripts/{TAG}.py','scripts/fresh_sqrt5_singular_extension_frontier_20261015.py','results/fresh_sqrt5_singular_extension_frontier_20261015.json','scripts/fresh_sqrt5_extended_forest_weights_20261015.py','results/fresh_sqrt5_extended_forest_weights_20261015.json','scripts/fresh_sqrt5_quad_filters_20261015_review.py','scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py','scripts/fresh_sqrt5_sealed_subcore_closure_20261015_review.py','results/fresh_sqrt5_sealed_subcore_closure_20261015_review.json','results/fresh_sqrt5_forest_sensitive_weights_20261015_review.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_INDEPENDENT_51_SINGULAR_DOMAINS_AND_EXTENDED_FOREST_RULE',sources={p:sha(p)for p in paths},input_origins=61,distinct_marked_classes=51,domain_profiles=502,domain_counts=dict(sorted(counts.items())),trial_counts=dict(sorted(trialcounts.items())),singular_records=records,upgraded_orbits=15,sealed_face_cases=316,face_checks=face_checks,component_totals=components,rank29_without_A1=dict(maximum=str(no[0]),forest=list(no[1])),rank29_with_A1=dict(maximum=str(yes[0]),forest=list(yes[1])),method='Separate tuple-Fraction exact radical matrices; terminal and pruning cornerPD checked by exact inertia, not the producer finite-diagram classifier. Independent marked graph isomorphisms, center-rooted marked-tree weights, and type-ordered two-state rankDP.',scope='Complete51 next-white necessary domains conditional on the61-input frontier; all46 newly sealed stars and112old ones receive arbitrary-extension payments under the reviewed sealing premise. Five classes with7accepted profiles, nonsingular and dotted frontiers still remain. No dimension theorem.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['domain_counts'],out['rank29_without_A1'],out['rank29_with_A1'])
