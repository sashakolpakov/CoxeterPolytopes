#!/usr/bin/env python3
"""Complete case-zero closure with literal symmetries and proof trees."""
import argparse,json,os
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
from fresh_h9_highbase_construction_20261014_cell import finite,codes

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_bugaenko_closure_20261014'
OLD='fresh_h9_bugaenko_surgery_20261014';PAIRS=list(combinations(range(14),2));PI={p:i for i,p in enumerate(PAIRS)}
PRIME=241;SQ={d:next(i for i in range(PRIME)if i*i%PRIME==d)for d in[2,3,5]}
COS={2:0,3:pow(2,-1,PRIME),4:SQ[2]*pow(2,-1,PRIME)%PRIME,5:(1+SQ[5])*pow(4,-1,PRIME)%PRIME,6:SQ[3]*pow(2,-1,PRIME)%PRIME}
def tup(m):return tuple(i for i in range(14)if m>>i&1)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def residue(code):
 A=[[int(i==j)for j in range(11)]for i in range(11)]
 for(i,j),m in zip(combinations(range(11),2),code):A[i][j]=A[j][i]=-COS[m]%PRIME
 out=1
 for i in range(11):
  j=next((j for j in range(i,11)if A[j][i]),None)
  if j is None:return 0
  if j!=i:A[i],A[j]=A[j],A[i];out=-out
  pivot=A[i][i];out=out*pivot%PRIME;u=pow(pivot,-1,PRIME)
  for j in range(i+1,11):
   c=A[j][i]*u%PRIME
   for k in range(i+1,11):A[j][k]=(A[j][k]-c*A[i][k])%PRIME
 return out%PRIME

class Gates:
 def __init__(self,row):
  self.M=list(map(tup,row['canonical_minimal_nonfaces']));self.V=list(map(tup,row['canonical_vertices']))
  self.vp=[tuple(PI[p]for p in combinations(v,2))for v in self.V]
  self.dotted={PI[m]for m in self.M if len(m)==2}
  self.cliques=[(a,tuple(PI[p]for p in combinations(a,2)))for a in combinations(range(14),11)if all(PI[p]not in self.dotted for p in combinations(a,2))]
  self.orthogonal=[(i,j,tuple(PI[tuple(sorted((a,b)))]for a in A for b in B))for i,A in enumerate(self.M)for j,B in enumerate(self.M)if i<j and not set(A)&set(B)]
  self.det_cache={}
 def bad_vertex(self,labels):
  return next((i for i,pp in enumerate(self.vp)if not finite(9,tuple(labels[p]for p in pp))),None)
 def proof(self,labels,known):
  bad=self.bad_vertex(labels)
  if bad is not None:return dict(kind='NONSPHERICAL_VERTEX',vertex=bad)
  for i,j,pp in self.orthogonal:
   if all(p in known and labels[p]==2 for p in pp):return dict(kind='ORTHOGONAL_LANNER',missing_faces=[i,j])
  for a,pp in self.cliques:
   if not all(p in known and labels[p]<=6 for p in pp):continue
   c=tuple(labels[p]for p in pp)
   if c not in self.det_cache:self.det_cache[c]=residue(c)
   if self.det_cache[c]:return dict(kind='RANK_AT_LEAST_ELEVEN',support=a,residue=self.det_cache[c])
  return None
 def complete(self,state):
  labels=[2]*91;known=set()
  for p,m in state:labels[PI[tuple(p)]]=m;known.add(PI[tuple(p)])
  free=sorted(set(range(91))-known-self.dotted);domains=[]
  for p in free:
   aa=[]
   for m in range(2,8):
    labels[p]=m
    if self.bad_vertex(labels)is None:aa.append(m)
   labels[p]=2;assert aa
   domains.append(dict(pair=p,labels=aa))
  for d in domains:
   if len(d['labels'])==1:known.add(d['pair']);labels[d['pair']]=d['labels'][0]
  unresolved={d['pair']:d['labels']for d in domains if len(d['labels'])>1}
  nodes=[];survivors=[];counts=Counter()
  def rec(left):
   idx=len(nodes);nodes.append(None);pr=self.proof(labels,known)
   if pr is not None:nodes[idx]=pr;counts[pr['kind']]+=1;return idx
   if not left:
    nodes[idx]=dict(kind='METRIC_REQUIRED',labels=labels[:]);survivors.append(labels[:]);return idx
   # First complete the most nearly fixed 11-minor; this is only an ordering.
   available=[(sum(p not in known for p in pp),pp)for a,pp in self.cliques if all(labels[p]<=6 for p in pp if p in known)]
   target=min(available,key=lambda q:q[0])[1]if available else []
   p=min(left,key=lambda p:(p not in target,len(unresolved[p]),p));known.add(p);children=[]
   for m in unresolved[p]:
    labels[p]=m;children.append(rec(left-{p}))
   labels[p]=2;known.remove(p);nodes[idx]=dict(kind='BRANCH',pair=p,labels=unresolved[p],children=children);return idx
  assert rec(set(unresolved))==0
  return dict(domains=domains,tree=nodes,terminal_counts=dict(counts),survivors=survivors)

def run():
 ip=f'results/{OLD}_incidence.json';lp=f'results/{OLD}_labels_0_25.json'
 inc=json.loads((ROOT/ip).read_text());ld=json.loads((ROOT/lp).read_text())
 for source in[inc,ld]:
  for p,h in source['sources'].items():assert sha(p)==h
 row=inc['cases'][0];table=ld['cases'][0];assert not table['state_limit_reached']and table['unprocessed_missing_faces']==0
 states=table['terminal_states'];assert len(states)==864
 M=row['canonical_minimal_nonfaces'];classes=defaultdict(list)
 for i in range(14):classes[tuple(j for j,m in enumerate(M)if m>>i&1)].append(i)
 blocks=list(classes.values());domain=[tuple(p)for p,m in states[0]];di={p:i for i,p in enumerate(domain)}
 group=[];actions=[];verts=set(row['canonical_vertices'])
 for images in product(*(list(permutations(b))for b in blocks)):
  perm=list(range(14))
  for b,im in zip(blocks,images):
   for a,z in zip(b,im):perm[a]=z
  move=lambda m:sum(1<<perm[i]for i in tup(m))
  assert {move(m)for m in M}==set(M)and{move(v)for v in verts}==verts
  group.append(perm);actions.append([di[tuple(sorted((perm[a],perm[b])))]for a,b in domain])
 assert len({tuple(p)for p in group})==288
 vectors=[tuple(m for p,m in s)for s in states];lookup={v:i for i,v in enumerate(vectors)};orbits=[];todo=set(lookup)
 while todo:
  s=min(todo);orb={tuple(s[j]for j in a)for a in actions};assert orb<=todo
  todo-=orb;orbits.append(dict(representative=lookup[s],members=sorted(lookup[v]for v in orb)))
 assert len(orbits)==6 and all(len(o['members'])==144 for o in orbits)
 gate=Gates(row);out=[]
 for i,o in enumerate(orbits):
  record=gate.complete(states[o['representative']]);out.append(dict(**o,**record));print('orbit',i,'tree',len(record['tree']),'terminal',record['terminal_counts'],'survivors',len(record['survivors']),flush=True)
 survivors=sum(len(o['survivors'])for o in out);files=[ip,lp,f'scripts/{TAG}.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 return dict(status='PASS_COMPLETE_CASE0_CLOSURE'if survivors==0 else'CASE0_METRIC_CANDIDATES',sources={p:sha(p)for p in files},case=0,prime=PRIME,square_roots=SQ,cosine_residues=COS,pair_order=PAIRS,twin_blocks=blocks,verified_literal_group=group,initial_partial_labelings=864,orbits=out,metric_candidate_count=survivors,scope='All864 complete-table necessary assignments, quotiented only by explicitly verified incidence automorphisms. Every remaining finite meeting-label category is branched or rejected by a necessary spherical, orthogonal-negative-index, or fixed-minor test. Dotted distances remain arbitrary. A zero result excludes this incidence; any survivor still needs a Lorentzian Gram metric.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();print('PID',os.getpid(),flush=True);d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
