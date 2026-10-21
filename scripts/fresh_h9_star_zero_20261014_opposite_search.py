#!/usr/bin/env python3
"""Resumable star0 exact necessary-label traversal from the empty assignment.

Explicit DFS frames retain all pending siblings. Every resumed batch validates
its stored paths and cursors. Full labels still need omitted-contact and metric completion.
"""
import argparse,json
from collections import Counter,defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h9_highbase_construction_20261014_cell import codes
from independent_h9_14_record30_a1_labels import finite_component

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_zero_20261014_opposite_search'
GEOM='results/fresh_h9_tree_cube_wedges_20261014.json'
GATE='results/fresh_h9_star_zero_20261014_opposite_gates.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read_checked(p,seen=None):
 if seen is None:seen=set()
 obj=json.loads((ROOT/p).read_text())
 if p in seen:return obj
 seen.add(p)
 for q,h in obj.get('sources',{}).items():
  assert sha(q)==h,q
  if q.endswith('.json'):read_checked(q,seen)
 return obj

def setup():
 g=read_checked(GEOM)['representatives'][0]['geometry']
 gate=read_checked(GATE)
 assert gate['case']==0 and gate['status']=='PASS_STAR_ZERO_WITH_OPPOSITE_SINGLE_FAMILIES'
 ms=[tuple(m) for m in g['minimal_nonfaces']];ts=[m for m in ms if len(m)==3]
 vs=[tuple(z['incident']) for z in g['vertices']]
 pairs=list(combinations(range(18),2));ix={p:i for i,p in enumerate(pairs)}
 T=[tuple(ix[p] for p in combinations(t,2)) for t in ts]
 V=[tuple(ix[p] for p in combinations(v,2)) for v in vs]
 pv=[tuple(j for j,v in enumerate(V) if i in v) for i in range(len(pairs))]
 clauses=[]
 for aa,bb in combinations(ms,2):
  A,B=set(aa),set(bb)
  if A&B or sum(set(c)<=A|B for c in ms)!=2:continue
  cross=tuple(sorted(ix[tuple(sorted((i,j)))] for i in aa for j in bb))
  assert all(pv[i] for i in cross)
  clauses.append(dict(blocks=[aa,bb],cross_pairs=cross))
 opts=sorted(codes(3),key=lambda c:(max(c),sum(c),c.count(2),c))
 domains=[list(opts) for _ in ts]
 for a,b in [b for b in g['blocks'] if len(b)==2]:
  pm=list(range(18));pm[a],pm[b]=b,a
  assert {tuple(sorted(pm[i] for i in v)) for v in vs}==set(vs)
  k=next(k for k,t in enumerate(ts) if a in t and b in t)
  leaf=next(i for i in ts[k] if i not in (a,b));cp=list(combinations(ts[k],2))
  u=cp.index(tuple(sorted((leaf,a))));v=cp.index(tuple(sorted((leaf,b))))
  domains[k]=[c for c in domains[k] if c[u]<=c[v]]
 @lru_cache(2000000)
 def finite(edges):
  labels=defaultdict(lambda:2);adj=defaultdict(set)
  for k,m in edges:
   a,b=pairs[k];labels[a,b]=m;adj[a].add(b);adj[b].add(a)
  left=set(adj)
  while left:
   start=min(left);co={start};q=[start];left.remove(start)
   for a in q:
    for b in adj[a]&left:left.remove(b);co.add(b);q.append(b)
   if len(co)>2 and not finite_component(co,labels):return False
  return True
 families=gate['families']
 old_family_count=0
 for f in families:
  f['conditions']=[(ix[a,b],'eq',m) for a,b,m in f['fixed']]+[(ix[a,b],'ge',m) for a,b,m in f['bounds']]
 keys={k for f in families for k,op,m in f['conditions']}
 incid={k:[(j,T[j].index(k)) for j in range(len(T)) if k in T[j]] for k in keys}
 assert all(incid.values())
 def family_hit(frame,new=True):
  arr,dom=frame['assignment'],frame['domains']
  fs=families if new else families[:old_family_count]
  # Prefer the new proof so its actual contribution is visible.
  if new:fs=families[old_family_count:]+families[:old_family_count]
  for f in fs:
   proof=[]
   for k,op,m in f['conditions']:
    if arr[k]:possible={arr[k]}
    else:
     possible=None
     for j,pos in incid[k]:
      val={row[pos] for row in dom[j]};possible=val if possible is None else possible&val
     if not possible:break
    if not all(x==m if op=='eq' else x>=m for x in possible):break
    proof.append(dict(pair=pairs[k],relation=op,bound=m,possible=sorted(possible)))
   else:return dict(family=f['id'],origin=f['origin'],conditions=proof)
  return None
 def bridge(frame):
  arr=frame['assignment'];tested={};vk=[tuple((i,arr[i]) for i in v if arr[i]>2) for v in V]
  def can(k):
   if arr[k]:return (arr[k]>2,-1)
   if k not in tested:
    tested[k]=(True,-1)
    for j in pv[k]:
     if not finite(tuple(sorted((*vk[j],(k,3))))):tested[k]=(False,j);break
   return tested[k]
  for j,c in enumerate(clauses):
   proof=[]
   for k in c['cross_pairs']:
    ok,vi=can(k)
    if ok:break
    proof.append((k,vi))
   else:return dict(clause=j,forced_cross_witnesses=proof)
  return None
 def select(frame):
  remaining=frame['remaining'];dom=frame['domains'];arr=frame['assignment']
  return min(remaining,key=lambda k:(len(dom[k]),-sum(arr[i]!=0 for i in T[k]),-sum(len(set(T[k])&set(T[j])) for j in remaining if j!=k),k))
 def child(frame,k,c):
  arr=frame['assignment'][:];dom=frame['domains'];remaining=frame['remaining']
  assert all(not arr[i] or arr[i]==m for i,m in zip(T[k],c))
  changed=[i for i in T[k] if not arr[i]]
  for i,m in zip(T[k],c):arr[i]=m
  affected=sorted({j for i in changed if arr[i]>2 for j in pv[i]})
  for j in affected:
   if not finite(tuple((i,arr[i]) for i in V[j] if arr[i]>2)):return None,'spherical_vertex'
  nxt=[j for j in remaining if j!=k];nd=dom[:]
  for j in nxt:
   nd[j]=[o for o in dom[j] if all(not arr[i] or arr[i]==m for i,m in zip(T[j],o))]
   if not nd[j]:return None,'empty_triangle_domain'
  return dict(assignment=arr,domains=nd,remaining=nxt,
      path=frame['path']+[dict(triangle=k,labels=list(c))],entered=False,selected=None,next_choice=0),None
 initial=dict(assignment=[0]*153,domains=domains,remaining=list(range(len(ts))),path=[],entered=False,selected=None,next_choice=0)
 return dict(g=g,ts=ts,pairs=pairs,clauses=clauses,families=families,initial=initial,
     finite=finite,family_hit=family_hit,bridge=bridge,select=select,child=child)

def validate_frames(env,stack):
 for frame in stack:
  base=env['initial']
  for step in frame['path']:
   k=env['select'](base);assert k==step['triangle']
   assert tuple(step['labels']) in [tuple(x) for x in base['domains'][k]]
   base,why=env['child'](base,k,tuple(step['labels']));assert base is not None and why is None
  for key in ('assignment','domains','remaining','path'):
   assert json.loads(json.dumps(base[key]))==json.loads(json.dumps(frame[key])),key
  if frame['selected'] is not None:
   assert frame['selected']==env['select'](frame)
   assert 0<=frame['next_choice']<=len(frame['domains'][frame['selected']])
 for parent,child in zip(stack,stack[1:]):
  assert child['path'][:-1]==parent['path']
  step=child['path'][-1];k=parent['selected'];assert k==step['triangle']
  assert tuple(parent['domains'][k][parent['next_choice']-1])==tuple(step['labels'])

def run(batch,budget,resume=None,verify=False):
 env=setup()
 if resume:
  prev=read_checked(resume);assert prev['status'] in ('PAUSED_EXACT_PENDING_DFS','RUNNING_RESUMABLE_DFS','FULL_MANDATORY_BRIDGE_SURVIVOR')
  stack=prev['pending_frames'];recovery=prev['recovery'];validate_frames(env,stack)
  if prev['status']=='FULL_MANDATORY_BRIDGE_SURVIVOR':
   assert env['family_hit'](stack[-1]) is not None,'A saved survivor cannot be skipped without a new proved gate'
   stack[-1]['entered']=False
  inherited_additional=prev['cumulative_additional_nodes']
 else:stack=[env['initial']];recovery={'starts_at_empty_assignment':True};inherited_additional=0
 sources=[GEOM,GATE,f'scripts/{TAG}.py','scripts/fresh_h9_highbase_construction_20261014_cell.py','scripts/independent_h9_14_record30_a1_labels.py']
 if resume:sources.append(resume)
 manifest={p:sha(p) for p in sources};prunes=Counter();nodes=0;solution=None;samples=[];trace=sha256();last_saved=-1
 def snapshot(status):
  return dict(status=status,sources=manifest,case=0,settings=dict(batch=batch,node_budget=budget,resume_input=resume),
      previously_completed_prefix_nodes=0,new_nodes=nodes,cumulative_additional_nodes=inherited_additional+nodes,
      prunes=dict(prunes),trace_sha256=trace.hexdigest(),prune_samples=samples,solution=solution,
      recovery=recovery,pending_frames=stack,proved_families=env['families'],triangle_supports=env['ts'],pairs=env['pairs'],
      scope='Exact necessary-label traversal from the empty state with explicit resumable pending frames. Category7 means every finite integer label>=7; dotted distances remain free. A survivor still requires omitted meeting contacts and Lorentzian metric completion. No old completed search prefix is assumed.')
 def save_running():
  path=ROOT/f'results/{TAG}_running_batch_{batch:03d}.json';tmp=path.with_suffix('.tmp')
  tmp.write_text(json.dumps(snapshot('RUNNING_RESUMABLE_DFS'),indent=2)+'\n');tmp.replace(path)
 def rejected(reason,frame,proof=None):
  prunes[reason]+=1;row=dict(path=frame['path'],reason=reason,proof=proof)
  trace.update(json.dumps(row,sort_keys=True,separators=(',',':')).encode())
  if len(samples)<40:samples.append(row)
 while stack:
  frame=stack[-1]
  if not frame['entered']:
   if nodes>=budget:break
   if nodes and nodes%1000==0 and last_saved!=nodes:
    if not verify:save_running()
    last_saved=nodes
    print('new_nodes',nodes,'depth',len(frame['path']),'pending_frames',len(stack),'prunes',dict(prunes),flush=True)
   nodes+=1
   hit=env['family_hit'](frame)
   if hit is not None:rejected('certified_'+hit['family'],frame,hit);stack.pop();continue
   hit=env['bridge'](frame)
   if hit is not None:rejected('orthogonal_negative_blocks',frame,hit);stack.pop();continue
   if not frame['remaining']:
    solution=dict(pair_labels=[list(p)+[m] for p,m in zip(env['pairs'],frame['assignment']) if m],triangle_path=frame['path'])
    break
   k=env['select'](frame)
   if frame['selected'] is None:frame['selected']=k;frame['next_choice']=0
   else:assert frame['selected']==k
   frame['entered']=True
  k=frame['selected']
  if frame['next_choice']>=len(frame['domains'][k]):stack.pop();continue
  c=tuple(frame['domains'][k][frame['next_choice']]);frame['next_choice']+=1
  nxt,why=env['child'](frame,k,c)
  if nxt is None:rejected(why,frame,dict(triangle=k,labels=list(c)));continue
  stack.append(nxt)
 status='FULL_MANDATORY_BRIDGE_SURVIVOR' if solution else 'PAUSED_EXACT_PENDING_DFS' if stack else 'EXACT_ALL_MANDATORY_LABELS_AND_BRIDGES_EMPTY'
 out=json.loads(json.dumps(snapshot(status)))
 if not verify:
  path=ROOT/f'results/{TAG}_batch_{batch:03d}.json';assert not path.exists(),path
  path.write_text(json.dumps(out,indent=2)+'\n')
  running=ROOT/f'results/{TAG}_running_batch_{batch:03d}.json'
  if running.exists():running.unlink()
 return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--batch',type=int,required=True);p.add_argument('--budget',type=int,default=50000);p.add_argument('--resume');p.add_argument('--verify',action='store_true');a=p.parse_args()
 if a.verify:
  saved=json.loads((ROOT/f'results/{TAG}_batch_{a.batch:03d}.json').read_text())
  a.budget=saved['settings']['node_budget'];a.resume=saved['settings']['resume_input']
 out=run(a.batch,a.budget,a.resume,a.verify)
 if a.verify:assert out==saved
 print(out['status'],'new_nodes',out['new_nodes'],'pending',len(out['pending_frames']),out['prunes'],flush=True)
