#!/usr/bin/env python3
"""Bounded fixed-fiber R1/R5 constructive profile CSP scout."""
from itertools import combinations
from collections import Counter,defaultdict
from pathlib import Path
import argparse,json,time
import z3

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_weighted_pentagonal_20261014_csp'
FULL=[sum(1<<j for j in c)for c in combinations(range(8),3)]
def sgn(a,b):
 x,y=2*a+b,b
 if not x:return(y>0)-(y<0)
 if not y:return(x>0)-(x<0)
 if x*y>0:return(x>0)-(x<0)
 n=x*x-5*y*y;assert n
 return((x>0)-(x<0))*((n>0)-(n<0))
def relation(A,B,S):
 u=(A^B).bit_count();r2=(A&B&S).bit_count();r1=((A^B)&S).bit_count();r0=S.bit_count()-r1-r2
 a=8*u+20*r0+4*r2-4*r1;b=-4*u-12*r0+4*r2+4*r1
 if(a,b)==(64,-32):return'2'
 if(a,b)==(96,-48):return'3'
 if(a,b)==(48,-16):return'5'
 t=sgn(a-128,b+64)
 return'dotted'if t>0 else'parallel'if not t else'forbidden'

def run(dim,seconds):
 stats={}
 for h in[2,4,6,8]:
  S=(1<<h)-1;stats[h]=dict(Counter(relation(a,b,S)for a in FULL for b in FULL))
 nodes=[p for p in range(1<<dim)if p.bit_count()%2]
 # Separate full-coordinate indicators permit direct exact relation cuts.
 solver=z3.Solver();solver.set(timeout=int(seconds*1000))
 choices=[z3.Int(f'c{i}')for i in range(len(nodes))]
 for c in choices:solver.add(c>=0,c<len(FULL))
 # All profile alternatives are exact; pair conflicts are added lazily.
 rounds=0;cuts=0;t0=time.time();result='unknown';saved=None
 while time.time()-t0<seconds:
  solver.set(timeout=max(1,int((seconds-(time.time()-t0))*1000)))
  sat=solver.check();rounds+=1
  if sat!=z3.sat:result=str(sat);break
  model=solver.model();pick=[model[c].as_long()for c in choices]
  failures=[]
  for i,a in enumerate(nodes):
   for j in range(i):
    if relation(FULL[pick[i]],FULL[pick[j]],a^nodes[j])=='forbidden':failures.append((i,j))
  if not failures:result='sat';saved=[dict(probe=p,full_mask=FULL[k])for p,k in zip(nodes,pick)];break
  for i,j in failures:solver.add(z3.Or(choices[i]!=pick[i],choices[j]!=pick[j]))
  cuts+=len(failures)
 return dict(status='BOUNDED_EXACT_RELATION_CSP_SCOUT',dimension=dim,profile_count=56,
  R5_nodes=len(nodes),R1_nodes=1<<(dim-1),stats=stats,result=result,rounds=rounds,
  generated_conflict_clauses=cuts,seconds=time.time()-t0,selection=saved,
  scope='Fixed unit/heavy signs, R1 at even golden parity and one k3 R5 cap at each odd probe. SAT is an exact finite cap-compatibility selection, not global spatial coverage; timeout is inconclusive.')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--dim',type=int,default=3);p.add_argument('--seconds',type=float,default=30);a=p.parse_args()
 data=run(a.dim,a.seconds);out=ROOT/f'results/{TAG}_d{a.dim}.json';out.write_text(json.dumps(data,indent=2)+'\n')
 print(json.dumps({k:v for k,v in data.items()if k!='selection'},indent=2))
 print(data['selection'])
