#!/usr/bin/env python3
"""Separate integer Gale, full-face and link-factor replay of saved proposals."""
import argparse,json,os,time
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_rank5_gale_20261014';FULL=32767
BASES=list(combinations(range(15),5));CIRCUITS=list(combinations(range(15),6))
@lru_cache(None)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def checked(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d.get('sources',{}).items():assert sha(f)==h,f
 return d
def bits(m):
 while m:
  b=m&-m;yield b;m-=b
def det(A):
 A=[r[:]for r in A];n=len(A);s=1;last=1
 for k in range(n-1):
  pivot=next((j for j in range(k,n)if A[j][k]),None)
  if pivot is None:return 0
  if pivot!=k:A[k],A[pivot]=A[pivot],A[k];s=-s
  a=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    num=A[i][j]*a-A[i][k]*A[k][j];assert num%last==0;A[i][j]=num//last
   A[i][k]=0
  last=a
 return s*A[-1][-1]
def minimal(values):
 out=[]
 for m in sorted(set(values),key=lambda x:(x.bit_count(),x)):
  if not any(k&m==k for k in out):out.append(m)
 return sorted(out)
def replay(path,expected):
 assert sha(path)==expected;d=checked(path);r=d['geometry'];V=r['gale'];assert len(V)==15 and all(len(x)==5 for x in V)
 determinants={b:det([V[i]for i in b])for b in BASES};assert all(determinants.values())
 weights=[0]*15;facets=[]
 for C in CIRCUITS:
  cof=[(-1)**i*determinants[C[:i]+C[i+1:]]for i in range(6)]
  if all(x>0 for x in cof)or all(x<0 for x in cof):
   facets.append(FULL^sum(1<<i for i in C))
   for i,x in zip(C,cof):weights[i]+=abs(x)
 assert weights==r['positive_kernel']and all(x>0 for x in weights)
 assert all(sum(weights[i]*V[i][j]for i in range(15))==0 for j in range(5))
 assert sorted(facets)==r['facets']and all(x.bit_count()==9 for x in facets)
 assert min(map(abs,determinants.values()))==r['minimum_absolute_minor']and r['nonzero_five_minors']==3003
 assert all(any(f>>i&1 for f in facets)for i in range(15))
 faces={0}
 for f in facets:
  s=f
  while s:faces.add(s);s=(s-1)&f
 mm=[m for m in range(1,FULL+1)if m not in faces and all(m^b in faces for b in bits(m))]
 assert mm==r['minimal_nonfaces']and all(2<=m.bit_count()<=5 for m in mm)
 ridge=Counter(f^b for f in facets for b in bits(f));assert set(ridge.values())=={2}and len(ridge)==r['ridge_count']
 # Link(K,v) forbids precisely the minimal sets among M\{v}.
 for v in range(15):
  lm=minimal(m&~(1<<v)for m in mm);assert 0 not in lm
  ghosts=sum(m for m in lm if m.bit_count()==1)
  linkmm=[m for m in lm if m.bit_count()>1];assert all(not(m&ghosts)for m in linkmm)
  actual=0
  for f in facets:
   if f>>v&1:actual|=f^(1<<v)
  assert actual==FULL^ghosts^(1<<v)
  assert not any(all(m==k or not(m&k)for k in linkmm)for m in linkmm)
 assert r['simplex_boundary_link_factors']==[]
 order=r['canonical_permutation'];assert sorted(order)==list(range(15));mp={old:new for new,old in enumerate(order)}
 def move(m):return sum(1<<mp[i]for i in range(15)if m>>i&1)
 assert sorted(map(move,mm))==r['canonical_minimal_nonfaces']and sorted(map(move,facets))==r['canonical_vertices']
 gate=d['gate']
 if gate['gate']=='two_lanner_blocks':
  a,b=gate['blocks'];assert a in mm and b in mm and not(a&b)
  assert sorted((a.bit_count(),b.bit_count()))in([3,4],[3,5],[4,4],[4,5],[5,5])
  assert gate['support']==a|b and {m for m in mm if m&(a|b)==m}=={a,b}
 elif gate['gate']=='literal_induced_pattern':
  P=checked(gate['source'])
  if gate['pattern']=='seven_node_relation':p=P['seven_node_obstruction']
  else:p=next(x for x in P['patterns']if x['pattern']==gate['pattern'])
  assert p['verdict']=='EXACT_LOCAL_OBSTRUCTION'
  order=gate['permutation'];assert len(order)==p['order']and len(set(order))==len(order)
  u=sum(1<<i for i in order);assert u==gate['support'];mp={old:new for new,old in enumerate(order)}
  reduced=sorted(sum(1<<mp[i]for i in order if m>>i&1)for m in mm if m&u==m)
  assert reduced==p['minimal_nonfaces']
 else:
  assert gate['gate']=='complete_label_search'
  import fresh_h9_fifteen_facet_20261014_engine as E
  from fresh_h9_bugaenko_closure_20261014_large_first_audit import finite,determinant
  E.finite=finite;E.residue=determinant
  answer=json.loads(json.dumps(E.search(r,d['label_node_limit']if 'label_node_limit'in d else 3000)))
  assert answer==gate['search']
 return dict(case=d['case'],input=path,input_sha256=expected,status='PASS_SEPARATE_INTEGER_GALE_INCIDENCE_AND_LITERAL_GATE',
  vertices=len(facets),minimal_nonfaces=len(mm),five_minors=3003,six_circuits=5005,vertex_links=15,gate=gate['gate'],
  search_verdict=gate.get('search',{}).get('verdict'),search_states=gate.get('search',{}).get('nodes',0))
def atomic(p,d):
 tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(d,separators=(',',':'))+'\n');tmp.replace(p)
def run(jobs,follow):
 primary=f'results/{TAG}.json';folder=ROOT/f'results/{TAG}_audit_cases';folder.mkdir(exist_ok=True);records={}
 files=[f'scripts/{TAG}_independent.py','scripts/fresh_h9_fifteen_facet_20261014_engine.py',
  'scripts/fresh_h9_bugaenko_closure_20261014_large_first_audit.py','scripts/fresh_h9_bugaenko_surgery_20261014_rank_independent.py']
 sources={p:sha(p)for p in files};pending={};start=time.monotonic();initial=checked(primary)
 for p in sorted(folder.glob('case_*.json')):
  d=checked(str(p.relative_to(ROOT)));assert d['sources']==sources;records[d['result']['case']]=d['result']
 def save(final=False):
  data=json.loads((ROOT/primary).read_text());refs={r['case']:r for r in data['cases']}
  for i,r in records.items():assert r['input_sha256']==refs[i]['sha256']
  src=dict(sources)
  if final:src[primary]=sha(primary)
  out=dict(status='PASS_ALL_SAVED_DIRECT_RANK5_GALE_CASES'if final else'RUNNING_SEPARATE_GALE_AUDIT',sources=src,
   cases_checked=len(records),saved_primary_cases=len(refs),vertices_checked=sum(r['vertices']for r in records.values()),
   five_minors_checked=3003*len(records),six_circuits_checked=5005*len(records),vertex_links_checked=15*len(records),
   records=[records[i]for i in sorted(records)],scope='Separate integer Bareiss determinants, positive circuits, full downward face sets and minimalized link ideals. Every literal local/two-block application is rechecked. Their inherited complete obstruction tables remain dependencies. Global label cases, if any, reuse DFS with separate finite-tree and determinant arithmetic. No proposal-walk exhaustion or all-fifteen-facet classification claim.')
  atomic(ROOT/f'results/{TAG}_independent.json',out);return out
 print('PID',os.getpid(),'resumed',len(records),'jobs',jobs,flush=True)
 with ProcessPoolExecutor(max_workers=jobs)as pool:
  while True:
   data=json.loads((ROOT/primary).read_text());assert data['sources']==initial['sources']
   terminal=not data['status'].startswith('RUNNING');waiting=set(pending.values())
   todo=[r for r in data['cases']if r['case']not in records and r['case']not in waiting]
   for r in todo[:jobs-len(pending)]:pending[pool.submit(replay,r['path'],r['sha256'])]=r['case']
   if not pending:
    if terminal or not follow:break
    time.sleep(2);continue
   done,_=wait(pending,timeout=3,return_when=FIRST_COMPLETED)
   for f in done:
    i=pending.pop(f);r=f.result();records[i]=r;atomic(folder/f'case_{i:05d}.json',dict(sources=sources,result=r))
    if len(records)%25==0:print('checked',len(records),'seconds',round(time.monotonic()-start,1),flush=True)
   if done:save()
 return save(final=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=2);ap.add_argument('--follow',action='store_true');a=ap.parse_args()
 d=run(a.jobs,a.follow);print(d['status'],d['cases_checked'],flush=True)
