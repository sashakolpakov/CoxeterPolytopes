#!/usr/bin/env python3
"""Couple free dotted-side endpoint families to a new pi/16 prism cap."""
import json,sys
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import z3
from scripts import fresh_h9_higher_face_20261013_linked_scout as L
S=L.S;LABELS=[2,3,4,5,6,8,10,12,16]
def run():
 _,_,_,_,_,names,_,_,R,G,center=S.load();corners=[];points=[]
 for k in range(4):
  for base in L.BASES:
   sides=[n for j,n in enumerate(L.T)if j!=k]+[base];corners.append(sides)
   A=R[[names.index(n)for n in sides+S.COMMON]];M=A@G@A.T;p=center-(center@G@A.T)@np.linalg.inv(M)@A;p/=np.sqrt(-p@G@p);points.append(p)
 D=-np.array(points)@G@np.array(points).T;out=[]
 for free in S.COMMON:
  if free=='compact2':continue
  u=R[names.index(free)]/np.sqrt(2);b=R[names.index('wall17')];v0=points[1]
  s0=2*np.cos(np.pi/16)/(-np.sqrt(2)*(v0@G@b));c0=np.sqrt(1+s0*s0);r0=np.sqrt(2)*(s0*v0-c0*u)
  assert all(L.gate(p)for p in R@G@r0)
  caps=[]
  for i,v in enumerate(points):
   if i==1:continue
   A=c0;B=D[1,i]*s0;delta=B*B-A*A
   if delta<=0:continue
   for m in LABELS:
    c=np.cos(np.pi/m);s=(B*c+A*np.sqrt(delta+c*c))/delta;cc=np.sqrt(1+s*s);r=np.sqrt(2)*(s*v-cc*u)
    assert abs(r@G@r-2)<1e-6 and abs(r@G@r0+2*c)<1e-6
    ps=R@G@r
    if r@G@center>=0 or not all(L.gate(p)for p in ps):continue
    caps.append({'corner':i,'sides':corners[i],'anchor_label':m,'sinh':float(s),'normal':r.tolist(),'finite_old':[names[j]for j,p in enumerate(ps)if p>-2+1e-8]})
  matrix=np.array([c['normal']for c in caps]);pair=matrix@G@matrix.T;conflict=[(i,j)for i,j in combinations(range(len(caps)),2)if not L.gate(pair[i,j])]
  opt=z3.Optimize();X=[z3.Bool(f'{free}_{i}')for i in range(len(caps))]
  for i,j in conflict:opt.add(z3.Or(z3.Not(X[i]),z3.Not(X[j])))
  opt.maximize(z3.Sum([z3.If(z3.Or([X[j]for j,c in enumerate(caps)if c['corner']==i]),1,0)for i in range(8)]));opt.maximize(z3.Sum([z3.If(x,1,0)for x in X]));assert opt.check()==z3.sat
  model=opt.model();chosen=[i for i,x in enumerate(X)if z3.is_true(model.eval(x))];fixed=[names.index(n)for n in S.COMMON if n!=free]
  roots=np.vstack([R,r0,matrix[chosen]]);geo=L.geom(roots,G,center,fixed)
  row={'free':free,'initial_corner':1,'initial_anchor':'wall17','initial_label':16,'initial_sinh':float(s0),'initial_normal':r0.tolist(),'coupled_candidates':caps,'conflicts':conflict,'selected':chosen,'H5_geometry':geo}
  print(free,'candidates',len(caps),'corners',sorted({c['corner']for c in caps}),'selected',[(caps[i]['corner'],caps[i]['anchor_label'])for i in chosen],geo,flush=True);out.append(row)
 (S.V.ROOT/(S.PREFIX+'_coupled_scout.json')).write_text(json.dumps({'status':'NUMERICAL_NEW_ANCHOR_COUPLED_PRISM_CAPS','labels':LABELS,'cases':out,'scope':'Selected continuous old-dotted families solved at prescribed pair angles to a newpi16cap. Numerical proposals only.'},indent=2)+'\n')
if __name__=='__main__':run()
