#!/usr/bin/env python3
"""Spherical attachments on the two remaining prism-corner edges."""
import json,sys
from itertools import product,combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,z3
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S;LABELS=[2,3,4,5,6,8,10,12,16]
def run():
 _,_,_,_,_,names,_,_,R,G,center=S.load();d=S.V.read(S.PREFIX+'_roof_walk_scout.json');free=d['free'];old=d['all_new'];u=R[names.index(free)]/np.sqrt(2);roots=np.vstack([R,np.array([x['normal']for x in old])]);caps=[];stats=[]
 for ci in(6,7):
  sides=L.T[:3]+[L.BASES[ci%2]];edge=sides+[n for n in S.COMMON if n!=free];E=R[[names.index(n)for n in edge]];M=E@G@E.T;Mi=np.linalg.inv(M)
  A=np.vstack([E,R[names.index(free)]]);v=center-(center@G@A.T)@np.linalg.inv(A@G@A.T)@A;v/=np.sqrt(-v@G@v)
  labels=np.array(list(product((2,3,4,5),repeat=8)),dtype=int);bs=-2*np.cos(np.pi/labels);qs=np.einsum('bi,ij,bj->b',bs,Mi,bs);mask=qs<2-1e-7;profiles=[(ll.tolist(),b,q)for ll,b,q in zip(labels[mask],bs[mask],qs[mask])]
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
    for m in LABELS:
     want=-2*np.cos(np.pi/m);z=(pp-want)/np.sqrt(qr)
     if z<-1e-8:continue
     s=(aa*z+bb*np.sqrt(delta+z*z))/delta;c=np.sqrt(1+s*s);r=part+np.sqrt(qr)*(s*v-c*u);ps=roots@G@r
     if abs(r@G@r-2)>1e-5 or not all(L.gate(p)for p in ps):continue
     caps.append({'corner':ci,'edge':edge,'profile':ll,'anchor_existing':aj,'anchor_label':m,'normal':r.tolist()});count+=1
  stats.append({'corner':ci,'spherical_profiles':len(profiles),'accepted':count});print(stats[-1],flush=True)
 selected=[]
 if caps:
  A=np.array([c['normal']for c in caps]);P=A@G@A.T;opt=z3.Optimize();X=[z3.Bool('x'+str(i))for i in range(len(caps))]
  for i,j in combinations(range(len(caps)),2):
   if not L.gate(P[i,j]):opt.add(z3.Or(z3.Not(X[i]),z3.Not(X[j])))
  opt.maximize(z3.Sum([z3.If(z3.Or([X[j]for j,c in enumerate(caps)if c['corner']==i]),1,0)for i in(6,7)]));opt.maximize(z3.Sum([z3.If(x,1,0)for x in X]));assert opt.check()==z3.sat;model=opt.model();selected=[i for i,x in enumerate(X)if z3.is_true(model.eval(x))]
 fixed=[names.index(n)for n in S.COMMON if n!=free];geo=L.geom(np.vstack([roots,np.array([caps[i]['normal']for i in selected])])if selected else roots,G,center,fixed)
 out={'status':'NUMERICAL_SPHERICAL_PROFILE_ROOF_CONTINUATION','source':S.PREFIX+'_roof_walk_scout.json','stats':stats,'candidates':caps,'selected':selected,'H5_geometry':geo,'scope':'Bounded spherical profiles and finite angles to the five new roof anchors. No exhaustive field/label claim.'}
 (S.V.ROOT/(S.PREFIX+'_profile_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('selected',[(caps[i]['corner'],caps[i]['profile'],caps[i]['anchor_existing'],caps[i]['anchor_label'])for i in selected]);print(geo)
if __name__=='__main__':run()
