#!/usr/bin/env python3
"""Bounded coherent replacement scout for the final prism corner."""
import json,sys
from pathlib import Path
from itertools import product,combinations
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,z3
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S

def run():
 _,_,_,_,_,names,_,_,R,G,center=S.load();d=S.V.read(S.PREFIX+'_roof_walk_scout.json');pd=S.V.read(S.PREFIX+'_profile_scout.json');free=d['free'];old=d['all_new']+[pd['candidates'][i]for i in pd['selected']];u=R[names.index(free)]/np.sqrt(2);caps=[];stats=[];oldr=np.array([x['normal']for x in old]);seen=set()
 corners=[[n for j,n in enumerate(L.T)if j!=ci//2]+[L.BASES[ci%2]]+[n for n in S.COMMON if n!=free]for ci in range(8)]
 for ci in (6,7):
  edge=corners[ci];E=R[[names.index(n)for n in edge]];M=E@G@E.T;Mi=np.linalg.inv(M);A=np.vstack([E,R[names.index(free)]]);v=center-(center@G@A.T)@np.linalg.inv(A@G@A.T)@A;v/=np.sqrt(-v@G@v)
  labels=np.array(list(product((2,3,4,5),repeat=8)));bs=-2*np.cos(np.pi/labels);qs=np.einsum('bi,ij,bj->b',bs,Mi,bs);mask=qs<2-1e-7;profiles=[(ll.tolist(),b,q)for ll,b,q in zip(labels[mask],bs[mask],qs[mask])]
  for j in range(8):
   for m in(6,8,10,12,16,24,32):
    ll=[2]*8;ll[j]=m;b=-2*np.cos(np.pi/np.array(ll));q=b@Mi@b
    if q<2-1e-7:profiles.append((ll,b,q))
  count=0
  for ll,b,q in profiles:
   part=b@Mi@E;qr=2-q
   for aj,anchor in enumerate(old):
    ar=np.array(anchor['normal']);aa=-(v@G@ar);bb=-(u@G@ar);delta=aa*aa-bb*bb
    if delta<=1e-9:continue
    pp=part@G@ar
    for m in(2,3,4,5,6,8,10,12,16,24,32):
     z=(pp+2*np.cos(np.pi/m))/np.sqrt(qr)
     if z<-1e-8:continue
     s=(aa*z+bb*np.sqrt(delta+z*z))/delta;r=part+np.sqrt(qr)*(s*v-np.sqrt(1+s*s)*u)
     if abs(r@G@r-2)>1e-5 or not all(L.gate(p)for p in R@G@r):continue
     conflict=[i for i,p in enumerate(oldr@G@r)if not L.gate(p)]
     if len(conflict)>1:continue
     key=tuple(np.round(r,8))
     if key in seen:continue
     seen.add(key);caps.append({'corner':ci,'edge':edge,'profile':ll,'anchor_existing':aj,'anchor_label':m,'normal':r.tolist(),'conflicts_with_roof':conflict});count+=1
  stats.append({'corner':ci,'profiles':len(profiles),'accepted':count});print(stats[-1],flush=True)
 bank=old+caps;BR=np.array([x['normal']for x in bank]);P=BR@G@BR.T;covers=[]
 for rr in BR:
  cover=[]
  for ci,edge in enumerate(corners):
   E=R[[names.index(n)for n in edge]];p=E@G@rr;M=np.block([[E@G@E.T,p[:,None]],[p[None,:],np.array([[2.]])]])
   if np.linalg.eigvalsh(M)[0]>1e-7:cover.append(ci)
  covers.append(cover)
 conflicts=[(i,j)for i,j in combinations(range(len(bank)),2)if not L.gate(P[i,j])]
 sol=z3.Optimize();sol.set(timeout=60000);X=[z3.Bool('branch'+str(i))for i in range(len(bank))]
 for i,j in conflicts:sol.add(z3.Or(z3.Not(X[i]),z3.Not(X[j])))
 sol.add(X[0])
 for ci in range(8):sol.add(z3.Or([X[i]for i,c in enumerate(covers)if ci in c]))
 sol.maximize(z3.Sum([z3.If(X[i],1,0)for i in range(len(old))]));sol.minimize(z3.Sum([z3.If(x,1,0)for x in X]));status=sol.check();out={'status':'NUMERICAL_COHERENT_REPLACEMENT_SCOUT','stats':stats,'source':S.PREFIX+'_profile_scout.json','bank':bank,'covers':covers,'conflicts':conflicts,'solver':str(status)}
 if status==z3.sat:
  model=sol.model();ids=[i for i,x in enumerate(X)if z3.is_true(model.eval(x))];out['selected']=ids;out['geometry']=L.geom(np.vstack([R,BR[ids]]),G,center,[names.index(n)for n in S.COMMON if n!=free]);print('selected',ids,'geometry',out['geometry'],flush=True)
 print('result',status,'bank',len(bank),flush=True);(S.V.ROOT/(S.PREFIX+'_branch_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':run()
