#!/usr/bin/env python3
"""Independent complete first-white domains and arithmetic/rank review.
Reversed variable/label order, separate finite-tree classifier and Fraction
field from the frozen independent contact checker; no producer imports.
"""
import argparse,json,sqlite3
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_quad_contacts_20261015_review as E
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_all_quad_contacts_20261015_review';PRIMARY='fresh_sqrt5_all_quad_contacts_20261015'
K=E.K;ALPHA=(10,6,5,4,3,2);TOKEN={3:(K(1),(0,0,0)),4:(K(1),(1,0,0)),5:(E.PHI,(0,0,0)),6:(K(1),(0,1,0)),10:(K(1),(0,0,1))}

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def pair(i,j):return tuple(sorted((i,j)))
@lru_cache(maxsize=100000)
def finite(n,word):return E.finite_word.__wrapped__(n,word)

def domains(row):
 n=row['nodes'];g=E.edges(row);c=[2]*n;out=[];nodes=0
 order=list(reversed(range(n)))
 corners=[tuple(sorted(row['white']+[n,a,b]))for a,b in product((0,2),(1,3))]
 def okay():return all(finite(len(ids),tuple(g.get(pair(i,j),2)for i,j in combinations(ids,2)))for ids in corners)
 def dfs(k):
  nonlocal nodes
  nodes+=1
  if k==n:
   if any(v!=2 for v in c):out.append(tuple(c))
   return
  i=order[k]
  for m in ALPHA:
   c[i]=m;g[i,n]=m
   if okay():dfs(k+1)
  del g[i,n];c[i]=2
 dfs(0);assert len(out)==len(set(out));return out,nodes

def field_pass(n,g):
 # Only finite-edge circuits are checked; dotted variables stay symbolic.
 adj={i:[]for i in range(n)}
 for(i,j),m in g.items():
  if m not in(0,2):adj[i].append((j,{3:0,4:1,5:0,6:2,10:4}[m]));adj[j].append((i,{3:0,4:1,5:0,6:2,10:4}[m]))
 colors={}
 for a in reversed(range(n)):
  if a in colors:continue
  colors[a]=0;todo=[a]
  for i in todo:
   for j,m in adj[i]:
    want=colors[i]^m
    if j in colors:
     if colors[j]!=want:return False
    else:colors[j]=want;todo.append(j)
 return True

def inv(M):
 n=len(M);B=[list(r)+[K(i==j)for j in range(n)]for i,r in enumerate(M)]
 for k in range(n):
  p=next(i for i in range(k,n)if B[i][k]);B[k],B[p]=B[p],B[k];v=B[k][k];B[k]=[z/v for z in B[k]]
  for i in range(n):
   if i==k:continue
   v=B[i][k];B[i]=[z-v*w for z,w in zip(B[i],B[k])]
 assert all(B[i][j]==(i==j)for i in range(n)for j in range(n))
 return[r[n:]for r in B]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),K())
def mv(M,v):return[dot(r,v)for r in M]
def scale_product(*terms):
 f=K(1);powers=[0,0,0]
 for a,e in terms:
  f=f*a;powers=[x+y for x,y in zip(powers,e)]
 assert all(p%2==0 for p in powers)
 for d,p in zip(E.RADS,powers):f=f*d**(p//2)
 return f

def prepare(row):
 n=row['nodes'];g=E.edges(row);adj=E.graph(range(n),g);assert len(E.components(adj))==1
 # Tree-product scaling from the largest-index root; no minimal-color scaling.
 paths={n-1:(K(1),(0,0,0))};todo=[n-1]
 while todo:
  i=todo.pop();f,e=paths[i]
  for j in sorted(adj[i],reverse=True):
   if j not in paths:
    h,k=TOKEN[g[pair(i,j)]];paths[j]=(f*h,tuple(x+y for x,y in zip(e,k)));todo.append(j)
 M=[[K()for _ in range(n)]for _ in range(n)]
 for i in range(n):M[i][i]=2*scale_product(paths[i],paths[i])
 for(i,j),m in g.items():
  if m!=2:M[i][j]=M[j][i]=-scale_product(paths[i],paths[j],TOKEN[m])
 assert M==E.tree_scaled(n,g)
 I=row['white']+[0,1];J=[2,3];P=E.restrict(M,I);PI=inv(P)
 assert E.inertia(P)==[len(I),0,0]and E.inertia([[z.conj()for z in r]for r in P])==[len(I),0,0]
 B=[[M[i][j]for i in I]for j in J];U=[mv(PI,b)for b in B]
 S=[[M[i][j]-dot(B[a],U[b])for b,j in enumerate(J)]for a,i in enumerate(J)]
 ii=E.inertia(S);jj=E.inertia([[z.conj()for z in r]for r in S]);assert ii[1]==1 and jj[1]==0
 return dict(n=n,paths=paths,I=I,J=J,PI=PI,B=B,U=U,S=S,core_inertia=[ii[0]+len(I),ii[1],ii[2]],core_conjugate=[jj[0]+len(I),jj[1],jj[2]])

def test(cache,contacts):
 n=cache['n'];a=next(i for i,m in enumerate(contacts)if m!=2);f,e=cache['paths'][a];h,k=TOKEN[contacts[a]];s=(f*h,tuple(x+y for x,y in zip(e,k)))
 diag=2*scale_product(s,s);v=[K()if m==2 else-scale_product(cache['paths'][i],s,TOKEN[m])for i,m in enumerate(contacts)]
 w=[v[i]for i in cache['I']];z=mv(cache['PI'],w);rho=diag-dot(w,z)
 assert rho.sign()>0 and rho.conj().sign()>0
 q=[v[j]-dot(u,w)for j,u in zip(cache['J'],cache['U'])]
 S=[cache['S'][i]+[q[i]]for i in range(2)]+[q+[rho]]
 physical=E.inertia(S);conjugate=E.inertia([[z.conj()for z in r]for r in S])
 assert physical[1]>=1
 p=[physical[0]+len(w),physical[1],physical[2]];c=[conjugate[0]+len(w),conjugate[1],conjugate[2]]
 return p,c

def run():
 data=read('results/'+PRIMARY+'.json');assert data['status']=='COMPLETE_ALL_QUAD_ONE_WHITE_DOMAINS'
 base=read('results/fresh_sqrt5_all_quad_unions_20261015.json')['labelled_union_supports'];classes=read('results/fresh_sqrt5_all_quad_unions_20261015_classes.json')['marked_isomorphism_classes']
 assert len(data['records'])==len(classes)==8298 and data['state']['next_class']==8298
 assert [r['class_index']for r in data['records']]==list(range(8298))
 # Verify the atomically committed producer records and cursor read-only.
 db=ROOT/data['database'];assert sha(data['database'])==data['database_sha256']
 con=sqlite3.connect('file:'+str(db)+'?mode=ro',uri=True);assert con.execute('PRAGMA integrity_check').fetchone()[0]=='ok'
 assert json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])==data['state']
 stored=[json.loads(r[0])for r in con.execute('SELECT value FROM records ORDER BY class_index')];con.close();assert stored==data['records']
 records=[];counts=Counter();trialcounts=Counter();domains_total=0;nodes_total=0;fulltests=0;affinechecks=0
 for rec in data['records']:
  ci=rec['class_index'];ui=classes[ci]['representative'];assert ui==rec['union_index'];row=base[ui];n=row['nodes'];assert n==rec['nodes']and(0 in row['labels'])==rec['dotted']
  status=rec['status'];counts[status]+=1;out=dict(class_index=ci,union_index=ui,status=status)
  if status=='REJECT_CORE_AFFINE_PRINCIPAL_SUBDIAGRAM':
   out['positive_kernel']=E.check_affine(row,rec['affine_witness']);affinechecks+=1
  else:
   dom,visits=domains(row);nodes_total+=visits
   saved={tuple(t['contacts']):t for t in rec['trials']};assert len(saved)==len(rec['trials'])and set(dom)==set(saved)
   assert len(dom)==rec['structural_profile_count'];domains_total+=len(dom)
   cache=None if rec['dotted']else prepare(row)
   if cache is not None:
    assert(cache['core_inertia'][2]>0)==rec['singular'];out.update(core_inertia=cache['core_inertia'],core_conjugate=cache['core_conjugate'])
   accepted=[];local=Counter();inertias=Counter()
   for contact in sorted(dom):
    trial=saved[contact];ts=trial['status'];trialcounts[ts]+=1;local[ts]+=1;ext=E.extension(row,contact);okay=field_pass(n+1,E.edges(ext))
    if not okay:assert ts=='REJECT_FIXED_FINITE_CYCLE_FIELD';continue
    assert ts!='REJECT_FIXED_FINITE_CYCLE_FIELD'
    if rec['dotted']:
     assert ts in('DOTTED_PARAMETER_UNFILTERED','REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM')
    else:
     p,c=test(cache,contact);fulltests+=1;inertias[str((tuple(p),tuple(c)))]+=1
     if p[1]>1:assert ts=='REJECT_TWO_NEGATIVE_DIRECTIONS';continue
     if c[1]>0:assert ts=='REJECT_CONJUGATE_NOT_PSD';continue
     assert ts in('PASS_EXTRA_WHITE_NECESSARY_TESTS','REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM')
    if ts=='REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM':E.check_affine(ext,trial['affine_witness']);affinechecks+=1
    else:accepted.append(contact)
   assert set(accepted)==set(map(tuple,rec['accepted']))and len(accepted)==len(rec['accepted'])
   expected='DOTTED_STRUCTURAL_ONE_WHITE_DOMAIN'if rec['dotted']else('FINITE_CORE_WITH_ONE_WHITE_EXTENSIONS'if accepted else'SEALED_SINGULAR_FINITE_CORE'if rec['singular']else'REJECT_SEALED_NONSINGULAR_FULL_STAR_RANK')
   assert status==expected
   out.update(domain_count=len(dom),domain_sha256=sha256(json.dumps(sorted(dom),separators=(',',':')).encode()).hexdigest(),independent_DFS_nodes=visits,accepted_count=len(accepted),accepted_sha256=sha256(json.dumps(sorted(accepted),separators=(',',':')).encode()).hexdigest(),trial_counts=dict(sorted(local.items())),full_inertia_counts=dict(sorted(inertias.items())))
  records.append(out)
  if len(records)%100==0:print('independent cores',len(records),'profiles',domains_total,'exactfinite',fulltests,flush=True)
 assert dict(counts)==data['counts']and dict(trialcounts)==data['profile_counts']and domains_total==81772
 paths=[f'scripts/{TAG}.py',f'scripts/{PRIMARY}.py',f'results/{PRIMARY}.json','scripts/fresh_sqrt5_quad_contacts_20261015_review.py','results/fresh_sqrt5_quad_contacts_20261015_review.json','results/fresh_sqrt5_all_quad_unions_20261015.json','results/fresh_sqrt5_all_quad_unions_20261015_classes.json','results/fresh_sqrt5_all_quad_unions_20261015_review_complete.json']
 for p in paths:
  if p.endswith('.json'):
   for q,h in read(p).get('sources',{}).items():assert sha(q)==h,(p,q)
 return dict(status='PASS_INDEPENDENT_ALL8298_QUAD_ONE_WHITE_DOMAINS_AND_EXACT_FILTERS',sources={p:sha(p)for p in paths},core_counts=dict(sorted(counts.items())),trial_counts=dict(sorted(trialcounts.items())),complete_cores=8298,structural_profiles=domains_total,independent_DFS_nodes=nodes_total,finite_exact_inertia_checks=fulltests,positive_affine_kernel_checks=affinechecks,records=records,method='Reverse-order complete domains from the separate finite-tree classifier; DFS tree-product scaling and Fraction3x3 congruence inertia after an independently inverted spherical corner. All positive affine kernels and producer atomic cursor/records verified. No producer imports.',scope='Complete necessary single-white domains only. Physical and conjugate tests are necessary, not sufficient for arbitrary simultaneous extensions. Finite empty domains give exact orthogonal sealing/full-star rank rejection; empty dotted domains give symbolic orthogonal sealing but leave distance/rank unknown. No new dimension bound.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['core_counts'],out['structural_profiles'])
