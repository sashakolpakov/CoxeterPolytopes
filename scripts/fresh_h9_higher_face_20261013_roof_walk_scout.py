#!/usr/bin/env python3
"""Continue the coupled roof along the remaining original prism vertices."""
import json,sys
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,z3
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S;LABELS=[2,3,4,5,6,8,10,12,16,24,32]
def run():
 _,_,_,_,_,names,_,_,R,G,center=S.load();d=S.V.read(S.PREFIX+'_coupled_scout.json')['cases'][0];free=d['free'];u=R[names.index(free)]/np.sqrt(2);points=[];sides=[]
 for k in range(4):
  for base in L.BASES:
   ss=[n for j,n in enumerate(L.T)if j!=k]+[base];sides.append(ss);A=R[[names.index(n)for n in ss+S.COMMON]];v=center-(center@G@A.T)@np.linalg.inv(A@G@A.T)@A;v/=np.sqrt(-v@G@v);points.append(v)
 D=-np.array(points)@G@np.array(points).T;old=[{'corner':1,'sinh':d['initial_sinh'],'normal':d['initial_normal'],'anchor':'wall17','label':16}]+[d['coupled_candidates'][i]for i in d['selected']];steps=[]
 for iteration in range(3):
  roots=np.vstack([R,np.array([x['normal']for x in old])]);covered=[]
  for i,ss in enumerate(sides):
   E=R[[names.index(n)for n in ss+[n for n in S.COMMON if n!=free]]]
   for row in old:
    r=np.array(row['normal']);M=np.block([[E@G@E.T,(E@G@r)[:,None]],[(r@G@E.T)[None,:],np.array([[2.]])]])
    if np.linalg.eigvalsh(M)[0]>1e-7:covered.append(i);break
  missing=[i for i in range(8)if i not in covered];caps=[];seen=set()
  for i in missing:
   for j,anchor in enumerate(old):
    sj=anchor['sinh'];A=np.sqrt(1+sj*sj);B=D[i,anchor['corner']]*sj;delta=B*B-A*A
    if delta<=1e-9:continue
    for m in LABELS:
     c=np.cos(np.pi/m);s=(B*c+A*np.sqrt(delta+c*c))/delta;cc=np.sqrt(1+s*s);r=np.sqrt(2)*(s*points[i]-cc*u)
     if not all(L.gate(p)for p in roots@G@r):continue
     key=tuple(np.round(r,8))
     if key in seen:continue
     seen.add(key);caps.append({'corner':i,'sinh':float(s),'normal':r.tolist(),'anchor_existing':j,'anchor_label':m})
  row={'iteration':iteration,'covered':covered,'missing':missing,'candidates':caps}
  if not caps:steps.append(row);break
  A=np.array([r['normal']for r in caps]);P=A@G@A.T;opt=z3.Optimize();X=[z3.Bool(f'x{iteration}_{i}')for i in range(len(caps))]
  for i,j in combinations(range(len(caps)),2):
   if not L.gate(P[i,j]):opt.add(z3.Or(z3.Not(X[i]),z3.Not(X[j])))
  opt.maximize(z3.Sum([z3.If(z3.Or([X[j]for j,c in enumerate(caps)if c['corner']==i]),1,0)for i in missing]));assert opt.check()==z3.sat;model=opt.model();chosen=[i for i,x in enumerate(X)if z3.is_true(model.eval(x))]
  row['selected']=chosen;old.extend(caps[i]for i in chosen);steps.append(row)
 fixed=[names.index(n)for n in S.COMMON if n!=free];geo=L.geom(np.vstack([R,np.array([r['normal']for r in old])]),G,center,fixed)
 out={'status':'NUMERICAL_COUPLED_ROOF_WALK','source':S.PREFIX+'_coupled_scout.json','free':free,'steps':steps,'all_new':old,'final_H5_geometry':geo,'scope':'Finite anchor labels on the prescribed original-vertex families, retaining the selected coherent roof. Numerical only.'}
 (S.V.ROOT/(S.PREFIX+'_roof_walk_scout.json')).write_text(json.dumps(out,indent=2)+'\n')
 print([(s['iteration'],s['covered'],len(s['candidates']),s.get('selected'))for s in steps]);print(geo)
if __name__=='__main__':run()
