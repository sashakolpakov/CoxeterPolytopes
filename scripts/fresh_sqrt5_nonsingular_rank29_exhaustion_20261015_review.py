#!/usr/bin/env python3
"""Independent finite rank29 attachment queue. No producer imports.
An immutable completed-batch JSON may be reviewed before the full queue closes.
"""
import argparse,json,sqlite3
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_dotted_rank29_exhaustion_20261015_review as D
P,R=D.P,D.R
Z,ONE=P.Z,P.ONE
ROOT=Path(__file__).resolve().parents[1];BASE='fresh_sqrt5_nonsingular_rank29_exhaustion_20261015';TAG=BASE+'_review'

def finite_matrix(row):
 assert 0 not in row['labels'];colors=R.coloring(row)
 return None if colors is None else R.matrix(row,colors)

def kernel(M):
 A=[r[:]for r in M];n=len(A);pivots=[];r=0
 for j in range(n):
  i=next((i for i in range(r,n)if A[i][j]!=Z),None)
  if i is None:continue
  A[i],A[r]=A[r],A[i];v=A[r][j];A[r]=[R.div(x,v)for x in A[r]]
  for i in range(n):
   if i!=r and A[i][j]!=Z:
    t=A[i][j];A[i]=[R.sub(x,R.mul(t,y))for x,y in zip(A[i],A[r])]
  pivots.append(j);r+=1
 out=[]
 for j in range(n):
  if j in pivots:continue
  x=[ONE if i==j else Z for i in range(n)]
  for k,i in enumerate(pivots):x[i]=R.neg(A[k][j])
  assert all(P.dot(z,x)==Z for z in M);out.append(x)
 return out

def graph(row,ids):
 g=P.edges(row);return {i:{j:g[D.pair(i,j)]for j in ids if j!=i and g[D.pair(i,j)]!=2}for i in ids}

def components(H):
 remaining=set(H);out=[]
 while remaining:
  seen={min(remaining)};todo=list(seen)
  for i in todo:
   for j in H[i]:
    if j not in seen:seen.add(j);todo.append(j)
  remaining-=seen;out.append({i:H[i]for i in sorted(seen)})
 return out

def typename(H):
 n=len(H);word=sorted(m for i in H for j,m in H[i].items()if i<j)
 if n==1:return 'A1'
 if n==2:return 'I2('+str(word[0])+')'
 degree={i:len(H[i])for i in H};branch=[i for i,d in degree.items()if d==3]
 if not branch:
  assert max(degree.values())<=2
  if 5 in word:assert n in(3,4);return 'H'+str(n)
  if 4 in word:
   special=next((i,j)for i in H for j,m in H[i].items()if m==4)
   return 'B'+str(n)if any(degree[i]==1 for i in special)else'F4'
  assert set(word)=={3};return 'A'+str(n)
 assert len(branch)==1 and max(degree.values())==3 and set(word)=={3}
 c=branch[0];arms=[]
 for x in H[c]:
  length=1;prev=c
  while degree[x]==2:
   nxt=next(j for j in H[x]if j!=prev);prev,x=x,nxt;length+=1
  arms.append(length)
 return ('D'if sorted(arms)[:2]==[1,1]else'E')+str(n)

def marked_key(H,a,b):
 # Full weighted rooted-tree codes, minimized over every possible root.
 def code(i,parent):return(int(i in(a,b)),tuple(sorted((m,code(j,i))for j,m in H[i].items()if j!=parent)))
 return json.dumps(min(code(i,None)for i in H),separators=(',',':'))

def corner_patterns(row):
 out=[]
 for a,b in D.CORNERS:
  cc=components(graph(row,row['white']+[a,b]));marked=next(H for H in cc if a in H);connected=b in marked
  out.append(dict(marked=[a,b],forest=[typename(H)for H in cc],a1_count=sum(len(H)==1 for H in cc),marked_component=typename(marked)if connected else None,pair_orbit=marked_key(marked,a,b)if connected else None))
 return out

def check_map(row,target,mp):
 n=row['nodes'];assert target['nodes']==n and sorted(mp)==list(range(n))and set(mp[:4])==set(range(4))
 assert {frozenset((mp[0],mp[2])),frozenset((mp[1],mp[3]))}=={frozenset((0,2)),frozenset((1,3))}
 assert set(mp[i]for i in row['white'])==set(target['white'])
 a,b=P.edges(row),P.edges(target)
 for(i,j),m in a.items():assert b[D.pair(mp[i],mp[j])]==m
 return n*(n-1)//2

def verify_database(data):
 # Only a complete immutable export binds the terminal database.
 path=ROOT/data['database'];assert P.sha(data['database'])==data['database_sha256']
 con=sqlite3.connect('file:'+str(path)+'?mode=ro',uri=True)
 assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
 assert json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])==data['state']
 assert [json.loads(x[0])for x in con.execute('SELECT payload FROM roots ORDER BY literal_index')]==data['root_maps']
 stored=list(con.execute('SELECT id,row,result FROM cores ORDER BY id'));assert len(stored)==len(data['records'])
 for (i,row,result),rec in zip(stored,data['records']):
  assert i==rec['index']and json.loads(row)==rec['row']and(json.loads(result)if result else None)==rec['result']
 con.close()

def run(input_path,check_database=False):
 data=P.read(input_path);complete=data['complete_queue'];state=data['state'];processed=state['next_core'];rows=data['records'];N=len(rows)
 assert data['dimension_bound']==29 and data['maximum_star_nodes']==31 and data['marked_states']==N
 assert complete==(processed==N)and 0<=processed<=N
 assert state['initial_literal_roots']==311 and state['initial_marked_roots']==199
 assert data['sources']==state['sources']
 for p,h in data['sources'].items():assert P.sha(p)==h,(p,h)
 if check_database:assert complete;verify_database(data)
 base=P.read('results/fresh_sqrt5_quad_exhaustion_20261015.json')['labelled_union_supports']
 origins=[r for r in P.read('results/fresh_sqrt5_sealed_subcore_closure_20261015.json')['records']if r['status']=='OPEN_NONSINGULAR'];assert len(origins)==len(data['root_maps'])==311
 reached=set();maps=0;comparisons=0
 for i,(origin,mp)in enumerate(zip(origins,data['root_maps'])):
  assert mp['literal_index']==i
  for k in('source_union_index','source_class_index','trial_index','contacts'):assert mp[k]==origin[k]
  row=D.extension(base[origin['source_union_index']],origin['contacts']);comparisons+=check_map(row,rows[mp['target']]['row'],mp['map']);maps+=1;reached.add(mp['target'])
 assert reached==set(range(199))
 counts=Counter();trials=Counter();records=[];obligations=[];nodes=0;profiles=0;inertias=0;cycle_products=0;relations=0;corners=0
 for i,rec in enumerate(rows):
  assert rec['index']==i and i in reached
  row=rec['row'];n=row['nodes'];assert n<=31 and row['white']==list(range(4,n))
  s=rec['result'];assert (s is not None)==(i<processed)
  if s is None:counts['PENDING']+=1;records.append(dict(index=i,nodes=n,status='PENDING'));continue
  g=P.edges(row);assert 0 not in row['labels']
  for a,b in D.CORNERS:
   ids=sorted(row['white']+[a,b]);assert D.finite(len(ids),tuple(g[D.pair(x,y)]for x,y in combinations(ids,2)))
   assert g[D.pair(a,b)]==2 and not any(g[D.pair(a,w)]!=2 and g[D.pair(b,w)]!=2 for w in row['white'])
  assert not D.h4_paid(row);M=finite_matrix(row);assert M is not None
  pi=D.sparse_inertia(M,row['white']);co=D.sparse_inertia(M,row['white'],True);inertias+=2
  assert pi==s['physical']and co==s['conjugate']and pi[1]==1 and pi[2]in(0,1)and co==[n-pi[2],0,pi[2]]
  singular=bool(pi[2]);assert singular==s['singular'];ker=kernel(M)if singular else[]
  cy=P.cycles(M);assert all(ok for v,ok in cy.values());cycle_products+=len(cy)
  if singular:
   assert len(ker)==1;obligations.append(i);cp=corner_patterns(row);assert cp==s['zero_touch_corner_patterns'];corners+=4
   assert s['orthogonal_spectator_rank_budget']=={str(d):d+2-n for d in(28,29)if d+2>=n}
  else:assert 'zero_touch_corner_patterns'not in s
  local=Counter();accept=[];dom=[];visits=0
  if n==31:
   assert s['trials']==s['accepted']==s['children']==[]
   status='DIMENSION29_SINGULAR_WEIGHT_OBLIGATION'if singular else'REJECT_RANK_EXCEEDS_DIMENSION29'
  else:
   dom,visits=D.domains(row);nodes+=visits;profiles+=len(dom);saved={tuple(t['contacts']):t for t in s['trials']}
   assert len(saved)==len(s['trials'])and set(saved)==set(dom),(i,len(saved),len(dom))
   for c in sorted(dom):
    ext=D.extension(row,c);t=saved[c]
    if D.h4_paid(ext):ts='PAID_ENTIRE_H4'
    elif(G:=finite_matrix(ext))is None:ts='REJECT_FINITE_CYCLE_FIELD'
    else:
     H=[z[:-1]for z in G[:-1]];K=ker if H==M else kernel(H)if singular else[]
     product=[P.dot(v,[G[j][-1]for j in range(n)])for v in K]
     if any(x!=Z for x in product):ts='REJECT_INHERITED_RELATION';relations+=1
     else:
      p=D.sparse_inertia(G,ext['white']);q=D.sparse_inertia(G,ext['white'],True);inertias+=2
      assert p==t['physical']and q==t['conjugate']
      if p[1]>1:ts='REJECT_TWO_NEGATIVE'
      elif q[1]:ts='REJECT_CONJUGATE'
      elif p[2]>1:ts='REJECT_MULTIPLE_LORENTZ_RELATIONS'
      else:
       assert p[1]==1 and p[2]==q[2]
       if singular:assert p[2]==1
       k=kernel(G)if p[2]else[]
       if k and(all(R.sign(x)>=0 for x in k[0])or all(R.sign(x)<=0 for x in k[0])):ts='REJECT_NONNEGATIVE_ROOT_DEPENDENCE'
       else:
        cc=P.cycles(G);assert all(ok for v,ok in cc.values());cycle_products+=len(cc);ts='ACCEPT_NECESSARY_UNPAID';accept.append(c)
    assert ts==t['status'],(i,c,ts,t['status']);local[ts]+=1;trials[ts]+=1
   assert set(accept)==set(map(tuple,s['accepted']))and len(accept)==len(s['accepted'])
   children={tuple(c['contacts']):c for c in s['children']};assert len(children)==len(s['children'])and set(children)==set(accept)
   for c in accept:
    child=children[c];comparisons+=check_map(D.extension(row,c),rows[child['index']]['row'],child['map']);maps+=1;reached.add(child['index'])
   status=('EXPANDED_SINGULAR'if singular else'EXPANDED_NONSINGULAR')if accept else('SEALED_SINGULAR_WEIGHT_OBLIGATION'if singular else'REJECT_SEALED_NONSINGULAR_RANK')
  assert status==s['status'];counts[status]+=1
  records.append(dict(index=i,nodes=n,singular=singular,status=status,domain_size=len(dom),independent_DFS_nodes=visits,domain_sha256=sha256(json.dumps(sorted(dom),separators=(',',':')).encode()).hexdigest(),trial_counts=dict(sorted(local.items())),accepted_count=len(accept)))
  if i%50==0:print('independent finite queue',i,'profiles',profiles,'singular obligations',len(obligations),flush=True)
 assert reached==set(range(N))and obligations==data['weight_obligations']and dict(counts)==data['counts']and dict(trials)==data['trial_counts']
 paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',input_path,f'scripts/{D.TAG}.py',f'results/{D.TAG}.json','results/fresh_sqrt5_sealed_subcore_closure_20261015.json','results/fresh_sqrt5_sealed_subcore_closure_20261015_review.json','results/fresh_sqrt5_triangle_weights_20261015_review.json','results/arithmetic_dimension_bounds_and_effective_high_dimension_search.md']
 if check_database:paths.append(data['database'])
 return dict(status='PASS_INDEPENDENT_COMPLETE_FINITE_ATTACHMENT_QUEUE_DIMENSION29'if complete else'PASS_INDEPENDENT_COMPLETED_FINITE_ATTACHMENT_BATCH_DIMENSION29',sources=P.manifest(paths),input_path=input_path,complete_queue=complete,database_checked=check_database,processed=processed,marked_states=N,pending=N-processed,literal_seeds=311,initial_marked_seeds=199,profiles=profiles,independent_DFS_nodes=nodes,exact_inertia_checks=inertias,cyclic_products=cycle_products,nonzero_inherited_relations=relations,seed_and_child_maps=maps,label_comparisons=comparisons,zero_touch_corner_patterns=corners,weight_obligations=len(obligations),counts=dict(counts),trial_counts=dict(trials),records=records,method='No producer imports. Reverse exact contact DFS with independent finite-forest leaf Schur tests; minimal radical gauge and sparse full congruence; fresh Gaussian nullspace; elementary cycles; complete literal seed/child maps and independently rebuilt entire finite corner types/orbits.',scope='Necessary finite-label attachment queue from311 literal nonsingular first extensions, restricted to dimension<=29. Every processed singular state, including nonsealed ones, exports a zero-contact obligation with its corner forests and spectator rank budget. Partial data do not imply exhaustion; actual geometry, payment, other seed families and a global dimension theorem remain separate.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--verify',action='store_true');ap.add_argument('--check-database',action='store_true');args=ap.parse_args();out=run(args.input,args.check_database);path=ROOT/f'results/{TAG}.json'
 if args.verify:assert out==P.read(str(path.relative_to(ROOT)))
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['processed'],out['marked_states'],out['counts'])
