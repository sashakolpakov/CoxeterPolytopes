#!/usr/bin/env python3
"""Separate literal-graph/Fraction review of enlarged three-context weights.
No project imports, no optimizer and no producer canonical-key parsing.
"""
import argparse,json
from collections import Counter,defaultdict,deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_sqrt5_enlarged_a1_weights_20261015';TAG=PRIMARY+'_review'
OPPOSITES={frozenset((0,2)),frozenset((1,3))}

def load(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def label_matrix(row):
 n=row['nodes'];assert row['white']==list(range(4,n));assert len(row['labels'])==n*(n-1)//2
 M=[[2]*n for _ in range(n)]
 for(i,j),v in zip(combinations(range(n),2),row['labels']):M[i][j]=M[j][i]=v
 return M
def check_bijection(A,B,p):
 assert len(A)==len(B)==len(p)and sorted(p)==list(range(len(p)))
 assert set(p[:4])==set(range(4))and{frozenset((p[0],p[2])),frozenset((p[1],p[3]))}==OPPOSITES
 assert all(A[i][j]==B[p[i]][p[j]]for i,j in combinations(range(len(p)),2))
def face_invariant(M):
 return len(M),tuple(sorted(tuple(sorted(M[i][j]for j in range(len(M))if i!=j))for i in range(4))),tuple(sorted(tuple(sorted(M[i][j]for j in range(len(M))if i!=j))for i in range(4,len(M))))
def face_iso(A,B):
 if face_invariant(A)!=face_invariant(B):return False
 n=len(A)
 for p in permutations(range(4)):
  if{frozenset((p[0],p[2])),frozenset((p[1],p[3]))}!=OPPOSITES:continue
  if any(A[i][j]!=B[p[i]][p[j]]for i,j in combinations(range(4),2)):continue
  mp=dict(enumerate(p));available=set(range(4,n))
  def recurse():
   if len(mp)==n:return True
   choices=[]
   for u in range(4,n):
    if u in mp:continue
    opts=[v for v in available if all(A[u][a]==B[v][b]for a,b in mp.items())and Counter(A[u][w]for w in range(4,n)if w not in mp and w!=u)==Counter(B[v][w]for w in available if w!=v)]
    if not opts:return False
    choices.append((len(opts),u,opts))
   _,u,opts=min(choices)
   for v in opts:
    mp[u]=v;available.remove(v)
    if recurse():return True
    available.add(v);del mp[u]
   return False
  if recurse():return True
 return False

def graph(vertices,edges):
 g={v:{}for v in vertices}
 for u,v,m in edges:
  if m!=2 and u in g and v in g:g[u][v]=m;g[v][u]=m
 return g
def components(g):
 unseen=set(g);answer=[]
 while unseen:
  root=min(unseen);queue=deque([root]);seen={root}
  while queue:
   for v in g[queue.popleft()]:
    if v not in seen:seen.add(v);queue.append(v)
  unseen-=seen;answer.append({v:dict(g[v])for v in seen})
 return answer
def tree_key(g,marks=()):
 assert len(components(g))==1 and sum(len(v)for v in g.values())==2*len(g)-2
 marked=set(marks)
 def rooted(v,parent):return(int(v in marked),tuple(sorted((m,rooted(w,v))for w,m in g[v].items()if w!=parent)))
 # Minimum over every root, rather than the producer's center-rooted encoding.
 return min(rooted(v,None)for v in g)
def model_graphs():
 out={}
 def path(name,labels):out[name]=graph(range(len(labels)+1),[(i,i+1,m)for i,m in enumerate(labels)])
 def arms(name,lengths):
  es=[];nextnode=1
  for length in lengths:
   previous=0
   for _ in range(length):es.append((previous,nextnode,3));previous=nextnode;nextnode+=1
  out[name]=graph(range(nextnode),es)
 for r in range(1,30):path('A'+str(r),[3]*(r-1))
 for r in range(2,30):path('B'+str(r),[3]*(r-2)+[4])
 for r in range(4,30):arms('D'+str(r),[r-3,1,1])
 for r in(6,7,8):arms('E'+str(r),[r-4,2,1])
 path('F4',[3,4,3]);path('H3',[3,5]);path('H4',[3,3,5])
 for m in(5,6,10):path('I2(%d)'%m,[m])
 assert len(out)==92 and len({tree_key(g)for g in out.values()})==92
 return out
MODELS=model_graphs();FAMILIES={tree_key(g):name for name,g in MODELS.items()}
def type_of(g):return FAMILIES[tree_key(g)]
def alias(name):return{'I2(3)':'A2','I2(4)':'B2'}.get(name,name)
def distance(g,a,b):
 queue=deque([(a,0)]);seen={a}
 while queue:
  u,d=queue.popleft()
  if u==b:return d
  for v in g[u]:
   if v not in seen:seen.add(v);queue.append((v,d+1))
 return None
def baseline(g,a,b):
 d=distance(g,a,b)
 return F(int(d is not None and d<=2)+int(type_of(g)=='H4'and(d==3 or g[a].get(b)==5)))
def corner_data(row):
 M=label_matrix(row);es=[(i,j,M[i][j])for i,j in combinations(range(len(M)),2)];out=[]
 for a,b in product((0,2),(1,3)):
  parts=components(graph(row['white']+[a,b],es));g=next(h for h in parts if a in h)
  out.append(dict(marked=[a,b],types=sorted(type_of(h)for h in parts),A1=sum(len(h)==1 for h in parts),baseline=baseline(g,a,b)if b in g else F(0),key=tree_key(g,(a,b))if b in g else None))
 return out

def literal_inputs():
 unions=load('results/fresh_sqrt5_quad_exhaustion_20261015.json')['labelled_union_supports'];old=load('results/fresh_sqrt5_forest_sensitive_weights_20261015.json')
 first=load('results/fresh_sqrt5_singular_extension_frontier_20261015.json');end=load('results/fresh_sqrt5_singular_attachment_closure_20261015.json')
 rows=[('old112',r['class_index'],unions[r['union_index']])for r in old['sealed_quad_rows']]
 rows.extend(('new46',r['index'],r['row'])for r in first['records']if not r['accepted'])
 rows.extend(('closure11',r['index'],r['row'])for r in end['records']if r['status']=='NEW_SEALED_PAID')
 rows.extend(('zero_touch5',r['index'],r['row'])for r in first['records']if r['accepted'])
 bfs=load('results/fresh_sqrt5_nonsingular_attachment_search_20261015_batch610.json')
 rows.extend(('bfs20',r['index'],r['row'])for r in bfs['records'])
 export=load('results/fresh_sqrt5_two_chain_star_closure_20261015_export.json')
 for r in export['records']:
  if r.get('necessary_weight_obligation'):rows.append(('two_chain3',r['source_index'],r['row']))
 for r in export['records']:
  if not r.get('necessary_weight_obligation'):rows.append(('two_chain_'+r['source_kind'],r['source_index'],r['row']))
 chain=load('results/fresh_sqrt5_two_chain_star_closure_20261015.json')
 rows.extend(('two_chain_zero_touch',r['index'],r['row'])for r in chain['records']if r['status']=='OPEN_WITH_CONTACTS')
 assert len(rows)==227;return rows

def forest_dp(totals):
 # Type-first unbounded coin DP, unlike the producer's rank-first recurrence.
 answers=[]
 for context in range(3):
  dp={(0,0):F(0)}
  for name in sorted(MODELS,reverse=True):
   rank=len(MODELS[name]);inc=int(name=='A1');value=totals[name][context]
   for n in range(rank,30):
    for count in range(3):
     if(n-rank,count)not in dp:continue
     key=n,min(2,count+inc);v=dp[n-rank,count]+value
     if key not in dp or v>dp[key]:dp[key]=v
  answers.append(dp)
 return answers

def run():
 P=load('results/'+PRIMARY+'.json');x=list(map(F,P['solution']['exact_primal']));zi=len(x)-1
 assert zi==282 and min(x)>=0
 keymap={};variable_index={};orbit_labels={}
 for r in P['orbits']:
  g=graph(r['nodes'],[(a,b,m)for(a,b),m in zip(combinations(r['nodes'],2),r['labels'])]);k=tree_key(g,r['marked'])
  assert k not in keymap and type_of(g)==alias(r['type'])and len(g)==r['rank'];assert baseline(g,*r['marked'])==0
  keymap[k]=r['key'];orbit_labels[r['key']]=k
 assert len(keymap)==94
 for j,r in enumerate(P['variable_order']):
  key=orbit_labels[r['orbit']],r['A1_count_capped2'];assert key not in variable_index and key[1]in range(3);variable_index[key]=j
 assert set(variable_index)=={(k,c)for k in keymap for c in range(3)}
 raw=literal_inputs();literal_by_origin={(fam,i):row for fam,i,row in raw};assert len(literal_by_origin)==227
 assert len(P['literal_maps'])==227 and len(P['faces'])==224
 face_matrices=[label_matrix(r['row'])for r in P['faces']];origin_count=Counter()
 for (fam,i,row),mp in zip(raw,P['literal_maps']):
  assert(mp['family'],mp['index'])==(fam,i);j=mp['representative'];check_bijection(label_matrix(row),face_matrices[j],mp['to_representative']);origin_count[j]+=1
  assert dict(family=fam,index=i,to_representative=mp['to_representative'])in P['faces'][j]['origins']
 assert len(origin_count)==224 and sum(origin_count.values())==227
 assert sum(len(r['origins'])for r in P['faces'])==227
 tested=0;pool=defaultdict(list)
 for M in face_matrices:
  iv=face_invariant(M)
  for N in pool[iv]:tested+=1;assert not face_iso(M,N)
  pool[iv].append(M)
 rows_unique={};face_payments=[];all_seen_orbits=set();case=0
 for fi,saved in enumerate(P['faces']):
  cs=corner_data(saved['row'])
  for c,s in zip(cs,saved['corners']):
   assert c['marked']==s['marked']and c['types']==sorted(alias(v)for v in s['core_forest'])and c['baseline']==F(s['baseline'])
   if s['orbit']is None:assert c['key']is None or c['baseline']>0
   else:assert keymap[c['key']]==s['orbit'];all_seen_orbits.add(c['key'])
  sums=[]
  for external in range(3):
   coef=Counter();base=sum(c['baseline']for c in cs)
   for c in cs:
    if c['key']in keymap and c['baseline']==0:coef[variable_index[c['key'],min(2,c['A1']+external)]]+=1
   lhs=base+sum(v*x[j]for j,v in coef.items());assert lhs>=1;sums.append(str(lhs))
   expected=dict(face=fi,external_A1_count_capped2=external,baseline=str(base),coefficients=[[j,v]for j,v in sorted(coef.items())]);assert expected==P['face_constraints'][case];case+=1
   key=tuple(sorted(coef.items())),base-1;rows_unique.setdefault(key,[]).append(dict(face=fi,origins=saved['origins'],external_A1_count_capped2=external))
  face_payments.append(dict(face=fi,payments=sums))
 assert case==672 and len(rows_unique)==385 and all_seen_orbits==set(keymap)
 rows=[];rhs=[];metadata=[]
 for(co,b),witnesses in sorted(rows_unique.items()):rows.append({j:-v for j,v in co});rhs.append(b);metadata.append(dict(kind='face',witnesses=witnesses))
 totals={};component_coefficients={};baselines={}
 for name,g in MODELS.items():
  base=F(0);coef=[Counter()for _ in range(3)]
  for a,b in combinations(g,2):
   base+=baseline(g,a,b);key=tree_key(g,(a,b))
   if key in keymap:
    assert baseline(g,a,b)==0
    for c in range(3):coef[c][variable_index[key,c]]+=1
  totals[name]=[base+sum(v*x[j]for j,v in co.items())for co in coef];component_coefficients[name]=coef;baselines[name]=base
 # Component labels in the saved table serve only as index names for LP cuts.
 table=P['exact_certificate']['component_totals'];names=[r['type']for r in table]
 assert len(names)==len(set(names))==89 and set(names)=={k for k,g in MODELS.items()if len(g)<=28}
 for r in table:assert r['rank']==len(MODELS[r['type']])and F(r['baseline'])==baselines[r['type']]and list(map(F,r['totals']))==totals[r['type']]
 for cut in P['solution']['forest_constraints']:
  forest=[names[i]for i in cut['components']];context=cut['context'];assert sum(len(MODELS[t])for t in forest)==28 and min(2,forest.count('A1'))==context
  co=Counter({zi:-1});b=F(0)
  for name in forest:co.update(component_coefficients[name][context]);b-=baselines[name]
  rows.append(dict(co));rhs.append(b);metadata.append(cut)
 slacks=[b-sum(v*x[j]for j,v in row.items())for row,b in zip(rows,rhs)];assert min(slacks)>=0
 dual=[F(0)]*len(rows)
 for record in P['solution']['nonzero_dual']:
  i=record['row'];assert dual[i]==0 and record['constraint']==metadata[i];dual[i]=F(record['multiplier']);assert dual[i]<0
 columns=[F(0)]*(zi+1)
 for y,row in zip(dual,rows):
  for j,a in row.items():columns[j]+=y*a
 assert all(v<=int(j==zi)for j,v in enumerate(columns))
 objective=sum(y*b for y,b in zip(dual,rhs));assert objective==x[-1]==F(P['solution']['optimum'])==F(153,2)
 dp=forest_dp(totals);forests=[]
 for n,expected in [(28,P['exact_certificate']['fronts']),(29,P['rank29_same_rule']['fronts'])]:
  for r in expected:
   c=r['context'];value=dp[c][n,c];assert value==F(r['maximum'])
   forest=r['maximizing_forest'];assert sum(len(MODELS[t])for t in forest)==n and min(2,forest.count('A1'))==c and sum(totals[t][c]for t in forest)==value
   if n==28:assert r['rank_maxima']==[str(dp[c][k,c])if(k,c)in dp[c]else None for k in range(29)]
   forests.append(dict(rank=n,context=c,maximum=str(value),witness=forest))
 assert max(dp[c][28,c]for c in range(3))==objective and max(dp[c][29,c]for c in range(3))==84
 assert[(name,str(baselines[name]))for name in('A29','B29','D29')]==[('A29','55'),('B29','55'),('D29','56')]
 # Directly replay the 76 baseline triangle markings. All added bonuses are nonnegative.
 triangle=load('results/fresh_sqrt5_triangle_weights_20261015.json');triangle_sums=[]
 for r in triangle['markings']:
  es=[(i,j,m)for(i,j),m in zip(combinations(range(r['order']),2),r['labels'])];v=F(0)
  for corner in r['corners']:
   a,b=corner['marked'];parts=components(graph(r['white']+[a,b],es));g=next(h for h in parts if a in h)
   v+=baseline(g,a,b)if b in g else F(0)
  assert v>=2;triangle_sums.append(str(v))
 assert len(triangle_sums)==76
 paths=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json','results/fresh_sqrt5_a1_count_weights_20261015.json','results/fresh_sqrt5_nonsingular_attachment_search_20261015_batch610.json','results/fresh_sqrt5_two_chain_star_closure_20261015_export.json','results/fresh_sqrt5_two_chain_star_closure_20261015.json','results/fresh_sqrt5_quad_exhaustion_20261015.json','results/fresh_sqrt5_forest_sensitive_weights_20261015.json','results/fresh_sqrt5_singular_extension_frontier_20261015.json','results/fresh_sqrt5_singular_attachment_closure_20261015.json','results/fresh_sqrt5_triangle_weights_20261015.json','results/fresh_sqrt5_triangle_weights_20261015_review.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in load(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_SEPARATE_LITERAL_GRAPH_RATIONAL_WEIGHT_REPLAY',sources={p:sha(p)for p in paths},literal_graphs=227,distinct_marked_graphs=224,dedup_comparisons=tested,face_cases=672,unique_face_constraints=385,orbits=94,variables=282,finite_component_types=92,primal=str(x[-1]),dual=str(objective),dual_inequalities=zi+1,primal_inequalities=len(rows),minimum_primal_slack=str(min(slacks)),face_payments=face_payments,forests=forests,triangles=76,triangle_baseline_payments=triangle_sums,
  method='No project imports or optimizer. All-root marked-tree canonicalization, arm-built finite diagrams, dynamic full-marked graph isomorphism, exact Fraction primal/dual arithmetic and reverse-type unbounded coin DP with actual capped A1 count.',
  scope='Independent implementation replay of the224-face weight certificate and both forest budgets. It inherits input sealing and attachment completeness from the bound source packages; no independent arithmetic review of those seals, no complete remaining-face coverage and no dimension theorem. Separate code is written in the same construction lane, not a separate mathematical author.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();p=ROOT/'results'/f'{TAG}.json';out=json.loads(json.dumps(run()))
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['primal'],out['forests'])
