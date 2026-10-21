#!/usr/bin/env python3
"""Atomic accelerated exact continuation of the immutable1433-core snapshot.
All new arithmetic uses the separate tuple-Fraction/minimal-gauge engine.
"""
import argparse,json,sqlite3,time
from collections import Counter,defaultdict
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_nonsingular_rank29_exhaustion_20261015_review as I
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J
D,P,R=I.D,I.P,I.R
Z=R.Z
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_nonsingular_rank29_fast_20261015'
SNAPSHOT='data/enumerations/fresh_sqrt5_nonsingular_rank29_exhaustion_20261015_batch1433.json'
DB=ROOT/f'data/enumerations/{TAG}_running.sqlite';OUT=ROOT/f'results/{TAG}.json'

def manifest():
 return P.manifest([f'scripts/{TAG}.py',f'scripts/{I.TAG}.py',f'results/{I.TAG}.json',SNAPSHOT,f'scripts/{J.TAG}.py',f'results/{J.TAG}.json',f'results/{D.TAG}.json','results/fresh_sqrt5_triangle_weights_20261015_review.json','results/arithmetic_dimension_bounds_and_effective_high_dimension_search.md'])

def analyze(M,white):
 # One field-valued elimination handles both embeddings. White pivots are
 # strictly positive in both embeddings. Kernel vectors are back-substituted.
 A=[r[:]for r in M];alive=set(range(len(A)));W=set(white);steps=[]
 while W:
  i=min(W,key=lambda x:sum(A[x][j]!=Z for j in W if j!=x));assert sum(A[i][j]!=Z for j in W if j!=i)<=1
  d=A[i][i];assert R.sign(d)>0 and R.sign(R.conj(d))>0
  neighbors=[j for j in alive if j!=i and A[i][j]!=Z];steps.append((i,d,[(j,A[i][j])for j in neighbors]))
  for k,j in enumerate(neighbors):
   v=R.div(A[i][j],d)
   for l in neighbors[k:]:A[j][l]=R.sub(A[j][l],R.mul(v,A[i][l]));A[l][j]=A[j][l]
  alive.remove(i);W.remove(i)
 ids=sorted(alive);S=[[A[i][j]for j in ids]for i in ids];pi=R.inertia(S);co=R.inertia([[R.conj(x)for x in r]for r in S]);pi[0]+=len(white);co[0]+=len(white)
 vectors=[]
 if pi[2]:
  for v in I.kernel(S):
   x=[Z]*len(M)
   for j,a in zip(ids,v):x[j]=a
   for j,d,terms in reversed(steps):x[j]=R.neg(R.div(P.sumfield(R.mul(a,x[k])for k,a in terms),d))
   assert all(P.dot(row,x)==Z for row in M);vectors.append(x)
 return pi,co,vectors

def h4_witness(row):
 result=[]
 for a,b in D.h4_paid(row):
  H=I.graph(row,row['white']+[a,b]);comp=next(c for c in I.components(H)if a in c);path=[a]
  while path[-1]!=b:path.append(next(j for j in H[path[-1]]if j not in path))
  result.append(dict(marked=[a,b],entire_H4=sorted(comp),path=path))
 return result

def evaluate(row):
 n=row['nodes'];assert n<=31 and row['white']==list(range(4,n));M=I.finite_matrix(row);assert M is not None
 pi,co,ker=analyze(M,row['white']);assert pi[1]==1 and pi[2]in(0,1)and co==[n-pi[2],0,pi[2]]
 singular=bool(pi[2]);out=dict(singular=singular,physical=pi,conjugate=co)
 if singular:
  out['zero_touch_corner_patterns']=I.corner_patterns(row)
  out['orthogonal_spectator_rank_budget']={str(d):d+2-n for d in(28,29)if d+2>=n}
 if n==31:return dict(**out,status='DIMENSION29_SINGULAR_WEIGHT_OBLIGATION'if singular else'REJECT_RANK_EXCEEDS_DIMENSION29',domain_counts={},trials=[],accepted=[])
 domain,nodes=D.domains(row);trials=[];accepted=[]
 for c in sorted(domain):
  ext=D.extension(row,c);t=dict(contacts=list(c))
  if(w:=h4_witness(ext)):t.update(status='PAID_ENTIRE_H4',witness=w)
  elif(G:=I.finite_matrix(ext))is None:t['status']='REJECT_FINITE_CYCLE_FIELD'
  else:
   if singular:
    H=[r[:-1]for r in G[:-1]];K=ker if H==M else I.kernel(H)
    products=[P.dot(v,[G[j][-1]for j in range(n)])for v in K]
   else:products=[]
   if any(x!=Z for x in products):t.update(status='REJECT_INHERITED_RELATION',products=[R.enc(x)for x in products])
   else:
    p,q,k=analyze(G,ext['white']);t.update(physical=p,conjugate=q)
    if p[1]>1:t['status']='REJECT_TWO_NEGATIVE'
    elif q[1]:t['status']='REJECT_CONJUGATE'
    elif p[2]>1:t['status']='REJECT_MULTIPLE_LORENTZ_RELATIONS'
    else:
     assert p[1]==1 and p[2]==q[2]
     if singular:assert p[2]==1
     if k and(all(R.sign(x)>=0 for x in k[0])or all(R.sign(x)<=0 for x in k[0])):t.update(status='REJECT_NONNEGATIVE_ROOT_DEPENDENCE',witness=[R.enc(x)for x in k[0]])
     else:t['status']='ACCEPT_NECESSARY_UNPAID';accepted.append(list(c))
  trials.append(t)
 status=('EXPANDED_SINGULAR'if singular else'EXPANDED_NONSINGULAR')if accepted else('SEALED_SINGULAR_WEIGHT_OBLIGATION'if singular else'REJECT_SEALED_NONSINGULAR_RANK')
 return dict(**out,status=status,domain_counts=dict(independent_reverse_DFS_nodes=nodes,profiles=len(domain)),trials=trials,accepted=accepted)

class Registry:
 def __init__(self,con):
  self.con=con;self.buckets=defaultdict(list);self.next_id=0
  for i,s in con.execute('SELECT id,row FROM cores ORDER BY id'):
   assert i==self.next_id;M=J.matrix(json.loads(s));self.buckets[J.invariant(M)].append((i,M));self.next_id+=1
 def add(self,row):
  M=J.matrix(row);key=J.invariant(M)
  for j,target in self.buckets[key]:
   mp=J.isomorphism(M,target)
   if mp is not None:J.check_map(M,target,mp);return j,mp
  j=self.next_id;self.next_id+=1;self.con.execute('INSERT INTO cores VALUES(?,?,NULL)',(j,json.dumps(row,separators=(',',':'))));self.buckets[key].append((j,M));return j,list(range(row['nodes']))

def export(con,state):
 records=[];counts=Counter();tc=Counter();obligations=[]
 for i,row,payload in con.execute('SELECT id,row,result FROM cores ORDER BY id'):
  result=json.loads(payload)if payload else None;rec=dict(index=i,row=json.loads(row),result=result);records.append(rec)
  assert(result is not None)==(i<state['next_core'])
  if result:
   counts[result['status']]+=1;tc.update(t['status']for t in result['trials'])
   if result['singular']:obligations.append(i)
  else:counts['PENDING']+=1
 complete=state['next_core']==len(records)
 data=dict(status='COMPLETE_ACCELERATED_FINITE_ATTACHMENT_QUEUE_DIMENSION29'if complete else'PARTIAL_ATOMIC_ACCELERATED_FINITE_ATTACHMENT_QUEUE_DIMENSION29',sources=state['sources'],state=state,dimension_bound=29,maximum_star_nodes=31,database=str(DB.relative_to(ROOT)),database_sha256=sha256(DB.read_bytes()).hexdigest(),marked_states=len(records),complete_queue=complete,counts=dict(counts),trial_counts=dict(tc),root_maps=[json.loads(r[0])for r in con.execute('SELECT payload FROM roots ORDER BY literal_index')],weight_obligations=obligations,records=records,scope='Exact dimension<=29 continuation from the independently reviewed immutable1433-state batch. Complete finite-corner/no-short contact domains, exact physical/conjugate inertia, root relations and full marked child maps. All singular stopping choices remain weight obligations, including expanded states and the31-node boundary. Orthogonal spectators obey the remaining rank budget. No all-rank, geometric-realization, weight-payment or dimension-improvement claim. The separate original producer remains unchanged; new sparse arithmetic requires independent review.')
 OUT.write_text(json.dumps(data,indent=2)+'\n');return data

def comparison(saved,new):
 # Independent gauges change relation vectors and domain counters, not decisions.
 for k in('status','singular','physical','conjugate','accepted'):assert saved[k]==new[k],k
 for k in('zero_touch_corner_patterns','orthogonal_spectator_rank_budget'):assert saved.get(k)==new.get(k),k
 a={tuple(t['contacts']):t for t in saved['trials']};b={tuple(t['contacts']):t for t in new['trials']};assert set(a)==set(b)
 for c in a:
  for k in('status','physical','conjugate'):assert a[c].get(k)==b[c].get(k),(c,k)


def verify():
 data=P.read(str(OUT.relative_to(ROOT)));assert data['sources']==manifest();I.verify_database(data)
 reached={m['target']for m in data['root_maps']};counts=Counter();tc=Counter();obligations=[]
 for rec in data['records']:
  i=rec['index'];assert i in reached
  if rec['result'] is None:counts['PENDING']+=1;continue
  new=evaluate(rec['row']);old=rec['result'];comparison(old,new);counts[new['status']]+=1;tc.update(t['status']for t in new['trials'])
  if new['singular']:obligations.append(i)
  children={tuple(c['contacts']):c for c in old['children']};assert set(children)==set(map(tuple,new['accepted']))and len(children)==len(old['children'])
  for c in new['accepted']:
   ch=children[tuple(c)];I.check_map(D.extension(rec['row'],c),data['records'][ch['index']]['row'],ch['map']);reached.add(ch['index'])
  if i%100==0:print('fast saved replay',i,flush=True)
 assert reached==set(range(len(data['records'])))and dict(counts)==data['counts']and dict(tc)==data['trial_counts']and obligations==data['weight_obligations']
 print('PASS_SAVED_ACCELERATED_FINITE_QUEUE',data['state']['next_core'],data['marked_states'])

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');ap.add_argument('--verify',action='store_true');ap.add_argument('--max-cores',type=int,default=10000);ap.add_argument('--max-seconds',type=float,default=600);args=ap.parse_args()
 if args.verify:verify();return
 sources=manifest();DB.parent.mkdir(parents=True,exist_ok=True)
 if not args.resume:assert not DB.exists(),'Use --resume; never overwrite an existing cursor.'
 con=sqlite3.connect(DB);con.execute('PRAGMA journal_mode=DELETE');con.execute('PRAGMA synchronous=FULL')
 if not args.resume:
  snap=P.read(SNAPSHOT);assert snap['state']['next_core']==1433 and len(snap['records'])==1649
  con.execute('BEGIN IMMEDIATE');con.execute('CREATE TABLE cores(id INTEGER PRIMARY KEY,row TEXT,result TEXT)');con.execute('CREATE TABLE roots(literal_index INTEGER PRIMARY KEY,payload TEXT)');con.execute('CREATE TABLE meta(id INTEGER PRIMARY KEY,value TEXT)')
  for r in snap['records']:con.execute('INSERT INTO cores VALUES(?,?,?)',(r['index'],json.dumps(r['row'],separators=(',',':')),json.dumps(r['result'],separators=(',',':'))if r['result']is not None else None))
  for r in snap['root_maps']:con.execute('INSERT INTO roots VALUES(?,?)',(r['literal_index'],json.dumps(r)))
  state=dict(sources=sources,next_core=1433,initial_literal_roots=311,initial_marked_roots=199,inherited_processed=1433,inherited_marked_states=1649,batches=[])
  con.execute('INSERT INTO meta VALUES(1,?)',(json.dumps(state),));con.commit()
 else:
  state=json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]);assert state['sources']==sources
 registry=Registry(con);start=state['next_core'];began=time.monotonic();stop='MAX_CORES'
 for k in range(args.max_cores):
  if k and time.monotonic()-began>=args.max_seconds:stop='TIME_GUARD_AFTER_COMPLETE_CORE';break
  i=state['next_core'];s=con.execute('SELECT row FROM cores WHERE id=?',(i,)).fetchone()
  if s is None:stop='QUEUE_EXHAUSTED';break
  row=json.loads(s[0]);z=evaluate(row);con.execute('BEGIN IMMEDIATE');children=[]
  for c in z['accepted']:
   j,mp=registry.add(D.extension(row,c));children.append(dict(index=j,contacts=c,map=mp))
  z['children']=children;con.execute('UPDATE cores SET result=? WHERE id=?',(json.dumps(z,separators=(',',':')),i));state['next_core']=i+1;con.execute('UPDATE meta SET value=? WHERE id=1',(json.dumps(state),));con.commit()
  if i%50==0:print('fast finite core',i,'nodes',row['nodes'],'states',registry.next_id,flush=True)
 state['batches'].append(dict(start=start,end=state['next_core'],stop_reason=stop,elapsed_seconds=round(time.monotonic()-began,3)))
 con.execute('UPDATE meta SET value=? WHERE id=1',(json.dumps(state),));con.commit();out=export(con,state);con.close()
 print(out['status'],state['next_core'],out['marked_states'],len(out['weight_obligations']))
if __name__=='__main__':main()
