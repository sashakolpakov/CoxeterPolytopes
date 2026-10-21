#!/usr/bin/env python3
"""Joint incidence proposals built from three disjoint Lanner triangles.

Weakening old finite contacts opens vertices; new finite contacts then
create higher missing faces. Optional additional dotted pairs can only
delete the spherical9 blocks tested here. This is a bounded numerical scout.
"""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import argparse,json,hashlib
import numpy as np
from fresh_h9_incidence_repair_20261014 import source,family,ROOT
OUT=ROOT/'results/fresh_h9_incidence_repair_20261014_lanner.json'

def peel(vertices):
 V=set(map(tuple,vertices));layers=[]
 while V:
  C=Counter(B for v in V for B in combinations(v,8));bad={v for v in V if any(C[b]<2 for b in combinations(v,8))}
  if not bad:break
  layers.append(len(bad));V-=bad
 return sorted(V),layers

def proposals(G,H,F,maxcost=4):
 tri=[]
 for T in combinations(range(15),3):
  E=list(combinations(T,2))
  if any(p in F for p in E):continue
  missing=[p for p in E if H[p[0]][p[1]]==0]
  if not 1<=len(missing)<=2 or any(G[i][j]for i,j in missing):continue
  old=[]
  for i,j in E:
   if (i,j)not in missing:
    sq=H[i][j]*H[i][j]/(H[i][i]*H[j][j]);old.append(3 if float(sq)<.3 else 4)
  for ms in product([3,4],repeat=len(missing)):
   if sum(1/m for m in old+list(ms))>=1-1e-12:continue
   tri.append((T,tuple((i,j,m)for(i,j),m in zip(missing,ms))))
 ans={}
 for a,b,c in combinations(range(len(tri)),3):
  Ts=[tri[k]for k in[a,b,c]]
  if len(set(i for T,_ in Ts for i in T))<9:continue
  E=tuple(sorted(e for _,es in Ts for e in es))
  if len(E)>maxcost:continue
  ans.setdefault(E,[list(t[0])for t in Ts])
 return sorted(ans.items(),key=lambda z:(len(z[0]),z[0]))

def main():
 p=argparse.ArgumentParser();p.add_argument('--patterns',type=int,default=72);p.add_argument('--per-pattern',type=int,default=64);p.add_argument('--cost',type=int,default=4);a=p.parse_args();row,G,ps=source();rows=[]
 for bi,pr in enumerate(ps[:a.patterns]):
  H,F,_=family(row,G,pr);cand=proposals(G,H,F,a.cost);B=[b for b in combinations(range(15),9)if not any(i in b and j in b for i,j in F)];Ba=np.array(B);Q=np.array(H,float);Q/=np.sqrt(np.diag(Q)[:,None]*np.diag(Q)[None,:]);tests=[]
  for additions,triangles in cand[:a.per_pattern]:
   M=Q.copy()
   for i,j,m in additions:M[i,j]=M[j,i]=-np.cos(np.pi/m)
   ev=np.linalg.eigvalsh(M[Ba[:,:,None],Ba[:,None,:]]);vertices=[list(B[k])for k in np.where(ev[:,0]>1e-9)[0]];core,layers=peel(vertices)
   tests.append({'additions':additions,'lanner_triangles':triangles,'vertex_count':len(vertices),'peel_layers':layers,'core_count':len(core),'core_vertices':core})
  print('pattern',bi,'proposals',len(cand),'tested',len(tests),'coremax',max((t['core_count']for t in tests),default=0),flush=True)
  rows.append({'pattern_index':bi,'support_index':pr['support_index'],'new_dotted_mask':pr['new_dotted_mask'],'available_proposals':len(cand),'tests':tests})
  d={'status':'NUMERICAL_JOINT_FINITE_CONTACT_AND_DOTTED_INCIDENCE_SCOUT','sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in[Path(__file__).resolve(),ROOT/'scripts/fresh_h9_incidence_repair_20261014.py',ROOT/'results/fresh_h9_incidence_surgery_20261013.json']},'scope':vars(a),'cases':rows}
  OUT.write_text(json.dumps(d,indent=2)+'\n')
if __name__=='__main__':main()
