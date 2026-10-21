#!/usr/bin/env python3
"""Numerical discovery: two orthogonal five-wall sets with S5 symmetry.
Seed cross Gram is -(2-b)I-aJ. All ten normals are independent.
The quartic norm equation in b is solved for every 0/1/phi profile orbit.
"""
import argparse,itertools,json,math
from pathlib import Path
import numpy as np
from numpy.polynomial import Polynomial as Poly
from scripts.fresh_h9_multiblock_seed_20261013 import compat,geometry
ROOT=Path(__file__).resolve().parents[1];PHI=(1+5**.5)/2
VALUES=(0.,1.,PHI);STATES=list(itertools.product(VALUES,repeat=2))

def compositions(total,k):
 if k==1:yield (total,);return
 for q in range(total+1):
  for rest in compositions(total-q,k-1):yield (q,)+rest

def profiles():
 rows=[]
 for counts in compositions(5,9):
  ts=sum(n*(x+y) for n,(x,y) in zip(counts,STATES));tq=sum(n*(x+y)**2 for n,(x,y) in zip(counts,STATES));ds=sum(n*(x-y) for n,(x,y) in zip(counts,STATES));dq=sum(n*(x-y)**2 for n,(x,y) in zip(counts,STATES));multi=math.factorial(5)//math.prod(math.factorial(n) for n in counts)
  rows.append((counts,ts,tq,ds,dq,multi))
 return rows

def parameters(a,P):
 b=Poly([0,1]);p=Poly([5*a,-1]);q=Poly([4,-1]);r=Poly([4+5*a,-1]);pts={}
 for i,(counts,ts,tq,ds,dq,multi) in enumerate(P):
  poly=tq*p*q*r-a*ts*ts*q*r+dq*b*p*r-a*ds*ds*b*p-4*b*p*q*r
  for root in poly.roots():
   if abs(root.imag)>1e-7:continue
   bb=float(root.real)
   if not(1e-6<bb<min(4,a)-1e-6):continue
   key=round(bb,8)
   if key not in pts:pts[key]={'b':bb,'profile_indices':set(),'polynomial':poly.coef.tolist(),'polynomial_profile':i}
   pts[key]['profile_indices'].add(i)
 for point in pts.values():
  bb=point['b'];ii=[]
  for i,(counts,ts,tq,ds,dq,multi) in enumerate(P):
   norm=tq/(2*bb)-a*ts*ts/(2*bb*(5*a-bb))+dq/(2*(4-bb))-a*ds*ds/(2*(4-bb)*(4-bb+5*a))
   if abs(norm-2)<1e-6:ii.append(i)
  point['profile_indices']=ii;point['cap_count']=sum(P[i][-1] for i in ii)
 return sorted(pts.values(),key=lambda x:(-x['cap_count'],x['b']))

def bank(a,b,P,ii):
 cross=-(2-b)*np.eye(5)-a*np.ones((5,5));G=np.block([[2*np.eye(5),cross],[cross,2*np.eye(5)]]);I=np.linalg.inv(G);cov=[];states=[]
 for i in ii:
  labels=tuple(j for j,n in enumerate(P[i][0]) for _ in range(n))
  for row in sorted(set(itertools.permutations(labels))):
   x=[STATES[j][0] for j in row];y=[STATES[j][1] for j in row];cov.append(-np.array(x+y));states.append(row)
 cov=np.array(cov);R=cov@I
 assert np.max(abs(np.sum((R@G)*R,axis=1)-2))<1e-5
 return G,R,states

def run(out,limit,attempts):
 P=profiles();records=[];census=[];rng=np.random.default_rng(1013)
 for m in [3,4,5,6,7,8,9,10,12]:
  a=2*np.cos(np.pi/m);pts=parameters(a,P);census.append({'m':m,'a':a,'parameter_points':len(pts),'top_cap_count':pts[0]['cap_count'] if pts else 0});print(census[-1],flush=True)
  for point in pts[:limit]:
   b=point['b'];G,R,states=bank(a,b,P,point['profile_indices']);C=compat(R@G@R.T);selections=set()
   for trial in range(attempts):
    rem=np.ones(len(R),dtype=bool);sel=[]
    while np.any(rem):
     avail=np.flatnonzero(rem)
     j=avail[np.argmax(np.sum(C[np.ix_(avail,avail)],axis=1))] if trial==0 else rng.choice(avail)
     sel.append(int(j));rem&=C[j]
    selections.add(tuple(sorted(sel)))
   best=None
   for sel in selections:
    geo=geometry(G,R,list(sel));key=(geo['max_klein_square'],geo['outside'])
    if best is None or key<best[0]:best=(key,sel,geo)
   rec={'m':m,'a':a,'b':b,'cap_count':len(R),'edges':int(np.sum(C)//2),'selections':len(selections),'max_selected':max(map(len,selections)),'selected_indices':list(best[1]),'selected_states':[list(states[i]) for i in best[1]],'parameter_polynomial':point['polynomial'],'parameter_profile_counts':list(P[point['polynomial_profile']][0]),'profile_indices':point['profile_indices'],'geometry':best[2]};records.append(rec)
   print({k:rec[k] for k in ['m','b','cap_count','max_selected','geometry']},flush=True)
   out.write_text(json.dumps({'scope':'Numerical S5 exchange seed, a=2cos(pi/m), m=3,4,5,6,7,8,9,10,12; all1287 orbit equations for norm2 caps with seed pair magnitudes0,1,phi; b roots in(0,min(4,a)); top parameter points by cap count, randomized greedy selections. No exact acceptance/exhaustive completion claim.','limit_per_label':limit,'selection_attempts':attempts,'census':census,'records':records},indent=2)+'\n')
 print({'done':len(records),'best':min(r['geometry']['max_klein_square'] for r in records)},flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=16);p.add_argument('--attempts',type=int,default=30);p.add_argument('--output',default='data/enumerations/fresh_h9_exchange_seed_20261013_scout.json');args=p.parse_args();run(ROOT/args.output,args.limit,args.attempts)
