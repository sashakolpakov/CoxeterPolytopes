#!/usr/bin/env python3
"""Exact necessary-label DFS with modular rank obstructions, bounded nodes."""
import argparse,json,time,os
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
from fresh_h9_highbase_construction_20261014_cell import finite,codes

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_surgery_20261014'
PAIRS=list(combinations(range(14),2));PI={p:i for i,p in enumerate(PAIRS)}
P=241;SQ={d:next(x for x in range(P)if x*x%P==d)for d in[2,3,5]}
COS={2:0,3:pow(2,-1,P),4:SQ[2]*pow(2,-1,P)%P,5:(1+SQ[5])*pow(4,-1,P)%P,6:SQ[3]*pow(2,-1,P)%P}
def tup(m):return tuple(i for i in range(14)if m>>i&1)
def det(A):
 A=[a[:]for a in A];d=1
 for k in range(len(A)):
  pivot=next((j for j in range(k,len(A))if A[j][k]%P),None)
  if pivot is None:return 0
  if pivot!=k:A[k],A[pivot]=A[pivot],A[k];d=-d
  c=A[k][k]%P;d=d*c%P;ci=pow(c,-1,P)
  for j in range(k+1,len(A)):
   a=A[j][k]*ci%P
   for i in range(k+1,len(A)):A[j][i]=(A[j][i]-a*A[k][i])%P
 return d%P
def minor(labels,a):
 return det([[1 if i==j else -COS[labels[PI[tuple(sorted((i,j)))]]]%P for j in a]for i in a])

def search(row,limit,extension_limit):
 M=list(map(tup,row['canonical_minimal_nonfaces']));V=list(map(tup,row['canonical_vertices']))
 con=[(tuple(PI[p]for p in combinations(a,2)),codes(len(a)))for a in M if len(a)>=3]
 vp=[tuple(PI[p]for p in combinations(v,2))for v in V]
 affected=[[j for j,pv in enumerate(vp)if set(a)&set(pv)]for a,_ in con]
 dotted={PI[a]for a in M if len(a)==2};covered={p for a,_ in con for p in a}
 free=sorted(set(range(91))-dotted-covered);cliques=[a for a in combinations(range(14),11)if all(PI[p]not in dotted for p in combinations(a,2))]
 clique_pairs=[tuple(PI[p]for p in combinations(a,2))for a in cliques]
 labels=[2]*91;nodes=terminal_count=extension_count=rank_rejects=0;dead=set();candidate=None;limited=False;witnesses={};trace=[]
 def complete():
  nonlocal terminal_count,extension_count,rank_rejects,candidate,limited
  terminal_count+=1;options=[]
  for p in free:
   aa=[]
   for m in range(2,8):
    labels[p]=m
    if all(finite(9,tuple(labels[j]for j in vv))for vv in vp):aa.append(m)
   labels[p]=2
   if not aa:return False
   options.append(aa)
  possible=1
  for a in options:possible*=len(a)
  if possible>extension_limit:limited=True;trace.append([terminal_count,'EXTENSION_LIMIT',possible]);return None
  local=0
  for choice in product(*options):
   for p,m in zip(free,choice):labels[p]=m
   if not all(finite(9,tuple(labels[j]for j in vv))for vv in vp):continue
   extension_count+=1;local+=1;found=False
   for a,pp in zip(cliques,clique_pairs):
    if any(labels[p]>=7 for p in pp):continue
    residue=minor(labels,a)
    if residue:
     code=tuple(labels[p]for p in pp);key=(a,code,residue)
     witnesses[key]=witnesses.get(key,0)+1;rank_rejects+=1;found=True;break
   if not found:candidate=labels[:];break
  for p in free:labels[p]=2
  trace.append([terminal_count,local,candidate is not None])
  return candidate is not None
 def walk(left,bits):
  nonlocal nodes,limited
  nodes+=1
  if nodes>limit:limited=True;return None
  if not left:return complete()
  key=(left,tuple((i,labels[i])for i in range(91)if bits>>i&1))
  if key in dead:return False
  candidates=[k for k in range(len(con))if left>>k&1]
  k=max(candidates,key=lambda k:(sum(bits>>p&1 for p in con[k][0]),-sum(not(bits>>p&1)for p in con[k][0]),len(con[k][0]),-k))
  pp,oo=con[k];new=[p for p in pp if not(bits>>p&1)];nbits=bits|sum(1<<p for p in new);rest=left&~(1<<k);known=[(i,p)for i,p in enumerate(pp)if bits>>p&1]
  for code in oo:
   if any(code[i]!=labels[p]for i,p in known):continue
   for p,m in zip(pp,code):labels[p]=m
   valid=all(finite(9,tuple(labels[p]for p in vp[vi]))for vi in affected[k])
   if valid:
    for j in candidates:
     if j==k:continue
     jj,options=con[j];known_j=[(i,p)for i,p in enumerate(jj)if nbits>>p&1]
     if known_j and not any(all(c[i]==labels[p]for i,p in known_j)for c in options):valid=False;break
   if valid:
    result=walk(rest,nbits)
    if result is not False:
     for p in new:labels[p]=2
     return result
   for p in new:labels[p]=2
  dead.add(key);return False
 result=walk((1<<len(con))-1,0)
 verdict='RANK_SURVIVOR_METRIC_REQUIRED'if result is True else('LIMIT_UNRESOLVED'if result is None else'EXACT_LABEL_AND_RANK_OBSTRUCTION')
 out=dict(case=row['case'],verdict=verdict,nodes=nodes,complete_Lanner_labelings=terminal_count,finite_extensions=extension_count,rank_rejections=rank_rejects,trace=trace,
  rank_witnesses=[dict(support=list(a),labels=list(c),residue=r,multiplicity=n)for(a,c,r),n in witnesses.items()])
 if candidate is not None:out.update(label_vector=candidate,dotted_pairs=[list(PAIRS[p])for p in sorted(dotted)])
 return out

def run(ids,nodes,ext):
 inc=json.loads((ROOT/f'results/{TAG}_incidence.json').read_text());rows=[];t=time.time();print('PID',os.getpid(),flush=True)
 for i in ids:
  r=search(inc['cases'][i],nodes,ext);rows.append(r);print('case',i,r['verdict'],'nodes',r['nodes'],'labelings',r['complete_Lanner_labelings'],'extensions',r['finite_extensions'],'seconds',round(time.time()-t,2),flush=True)
 paths=[f'scripts/{TAG}_rank_search.py',f'results/{TAG}_incidence.json','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='EXACT_BOUNDED_LABEL_AND_RANK_SEARCH',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},case_ids=ids,node_limit=nodes,free_extension_limit=ext,prime=P,square_roots=SQ,cosine_residues=COS,cases=rows,scope='Complete finite-label categories in Lanner and vertex tests. Nonzero fixed 11-minor modulo241 proves rank>10 independently of all dotted distances and every variable >=7 edge. Only untruncated finished cases are excluded; a surviving category needs actual metric/large labels.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--ids',default='0,2,3,31,11,28,318,636,637');p.add_argument('--nodes',type=int,default=5000);p.add_argument('--extensions',type=int,default=200000);p.add_argument('--verify',action='store_true');a=p.parse_args();ids=list(map(int,a.ids.split(',')));o=run(ids,a.nodes,a.extensions);path=ROOT/f'results/{TAG}_rank_search.json'
 if a.verify:assert o==json.loads(path.read_text())
 else:path.write_text(json.dumps(o,indent=2)+'\n')
 print(Counter(r['verdict']for r in o['cases']))
