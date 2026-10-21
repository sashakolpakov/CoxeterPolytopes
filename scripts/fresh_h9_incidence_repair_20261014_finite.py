#!/usr/bin/env python3
"""Propose new nonright contacts on old right pairs, not old-edge relabeling."""
from pathlib import Path
from itertools import combinations,product
from fractions import Fraction as F
import json,hashlib,argparse,time
import numpy as np
from fresh_h9_incidence_repair_20261014 import source,ROOT
OUT=ROOT/'results/fresh_h9_incidence_repair_20261014_finite.json'

def covers(row,G):
 P=set(map(tuple,row['parallel_pairs']));C=[set(a['nodes'])for a in row['affine_circuits']if not any(i in a['nodes']and j in a['nodes']for i,j in P)]
 Z=[p for p in combinations(range(15),2)if not G[p[0]][p[1]]and any(set(p)<=c for c in C)]
 req=[set(p for p in Z if set(p)<=c)for c in C];forced=set(next(iter(z))for z in req if len(z)==1);remaining=[p for p in Z if p not in forced]
 for k in range(len(forced),len(forced)+8):
  ans=[]
  for comb in combinations(remaining,k-len(forced)):
   S=forced|set(comb)
   if all(z&S for z in req):ans.append(sorted(S))
  if ans:return ans,sorted(forced)
 raise AssertionError('cover bound')

def main():
 p=argparse.ArgumentParser();p.add_argument('--all-binary',action='store_true');a=p.parse_args();row,G,_=source();S,forced=covers(row,G);P=set(map(tuple,row['parallel_pairs']));Q=np.array(G,float);Q/=np.sqrt(np.diag(Q)[:,None]*np.diag(Q)[None,:]);Bs=[B for B in combinations(range(15),9)if not any(i in B and j in B for i,j in P)];B=np.array(Bs);rows=[]
 print('right-only covers',len(S),'size',len(S[0]),'forced',forced,'candidate9',len(Bs),flush=True)
 for si,edges in enumerate(S):
  assignments=list(product([3,4],repeat=len(edges)))if a.all_binary else[(3,)*len(edges),(4,)*len(edges),(5,)*len(edges)]
  for labels in assignments:
   H=Q.copy()
   for(i,j),m in zip(edges,labels):H[i,j]=H[j,i]=-np.cos(np.pi/m)
   vals=np.linalg.eigvalsh(H[B[:,:,None],B[:,None,:]]);good=np.where(vals[:,0]>1e-9)[0]
   vertices=[list(Bs[i])for i in good];covered=set(i for v in vertices for i in v)
   rows.append({'support_index':si,'labels':labels,'finite9_count':len(good),'covered_facets':sorted(covered),'vertices':vertices})
  print(si,'PD9max',max(r['finite9_count']for r in rows if r['support_index']==si),'allfacets',sum(len(r['covered_facets'])==15 for r in rows if r['support_index']==si),flush=True)
 d={'status':'NUMERICAL_NEW_FINITE_CONTACT_DIAGRAM_SCOUT','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[ROOT/'results/fresh_h9_incidence_surgery_20261013.json',Path(__file__).resolve()]},'case':row['id'],'right_only_minimum_covers':S,'forced_pairs':forced,'profiles':rows,'scope':'All minimum right-only affine covers; either uniform labels3/4/5 or binary3/4 profiles. Exact old labels retained, old parallels free. Additional right-to-dotted changes cannot create new spherical9 blocks; this scout tests necessary vertex availability only.'}
 OUT.write_text(json.dumps(d,indent=2)+'\n')
if __name__=='__main__':main()
