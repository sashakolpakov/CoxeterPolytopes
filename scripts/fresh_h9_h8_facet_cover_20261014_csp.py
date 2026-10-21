#!/usr/bin/env python3
"""Bounded numerical necessary vertex-to-ball assignment test, not placement."""
import argparse,json
from collections import deque
from hashlib import sha256
from pathlib import Path
from fresh_h9_h8_facet_cover_20261014_geometry import data,SOURCE,BASE,CAP
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_facet_cover_20261014_csp'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def solve(d,limit=100000):
 ds,bs,dc,rho=d['D'],d['bound'],d['DC'],d['rho'];n=41;b=12
 allowed={}
 for i in range(n):
  for j in range(n):
   if i==j:continue
   allowed[i,j]=[sum(1<<y for y in range(b)if ds[i,j]<=bs[x,y]+1e-30 and ds[i,j]+1e-30>=max(0,dc[x,y]-rho[x]-rho[y]))for x in range(b)]
 def propagate(dom,queue=None):
  q=deque(queue if queue is not None else[(i,j)for i in range(n)for j in range(n)if i!=j])
  while q:
   i,j=q.popleft();new=sum(1<<x for x in range(b)if(dom[i]>>x)&1 and allowed[i,j][x]&dom[j])
   if new!=dom[i]:
    if not new:return None
    dom[i]=new;q.extend((k,i)for k in range(n)if k!=i and k!=j)
  return dom
 nodes=0;cutoff=False
 initial=propagate([(1<<b)-1]*n)
 def go(dom):
  nonlocal nodes,cutoff
  nodes+=1
  if nodes>limit:cutoff=True;return None
  if all(x.bit_count()==1 for x in dom):return dom
  i=min((i for i in range(n)if dom[i].bit_count()>1),key=lambda i:(dom[i].bit_count(),-sum(float(ds[i,j])for j in range(n))))
  for x in range(b):
   if not(dom[i]>>x)&1:continue
   dd=dom.copy();dd[i]=1<<x;dd=propagate(dd,[(j,i)for j in range(n)if j!=i])
   if dd is not None:
    ret=go(dd)
    if ret:return ret
   if cutoff:return None
  return None
 r=go(initial)if initial else None
 assignment=None if r is None else[x.bit_length()-1 for x in r]
 if r:
  assert all(allowed[i,j][assignment[i]]&(1<<assignment[j])for i in range(n)for j in range(n)if i!=j)
 return dict(u=str(d['u']),t=str(d['t']),status='SAT_NECESSARY_RELAXATION'if r else'NODE_LIMIT'if cutoff else'NUMERICAL_UNSAT',nodes=nodes,node_limit=limit,
   initial_domain_sizes=None if initial is None else[x.bit_count()for x in initial],assignment=assignment,
   source_diameter=str(max(ds)),union_diameter=str(max(bs)),center_diameter=str(max(dc)),ball_radii=[str(x)for x in rho])
def run():
 rows=[]
 for u in[None]+[str(.7642+(.99613-.7642)*i/20)for i in range(21)]:
  d=data(u=u);r=solve(d);r['literal_cap']=u is None;rows.append(r);print(r['u'][:12],r['status'],r['nodes'],flush=True)
 paths=[SOURCE,BASE,CAP,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py','scripts/fresh_h9_h8_facet_cover_20261014_geometry.py',f'scripts/{TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_NECESSARY_ASSIGNMENT_TEST',sources={p:sha(p)for p in paths},precision_decimal=80,rows=rows,
  scope='Literal cap plus21 u-grid points. SAT only assigns vertices to balls under necessary pair-distance bounds; it is neither an isometric placement nor full-polytope coverage. Numerical UNSAT is not an exact certificate. No angle compatibility imposed.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print('PASS',flush=True)
