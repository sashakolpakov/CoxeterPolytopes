#!/usr/bin/env python3
"""All existing spherical8 faces, genuinely trivalent new finite-node profiles."""
import json,itertools,time
from pathlib import Path
from fractions import Fraction
import numpy as np
ROOT=Path(__file__).resolve().parents[1];PRE='data/enumerations/fresh_h9_component_bridge_20261013'

def rootsof(A,B,C):
 if abs(A)<1e-12:return [-C/B]if abs(B)>1e-12 else[]
 d=B*B-4*A*C
 if d< -1e-8:return[]
 d=np.sqrt(max(d,0));return[(-B-d)/(2*A),(-B+d)/(2*A)]

def run():
 d=json.loads((ROOT/(PRE+'_faces_scout.json')).read_text());G=np.array(d['gram']);R=np.array(d['normals']);names=d['names'];h=np.array(d['h']);center=np.array(d['center']);targets=[np.array(d['target'])]
 od=json.loads((ROOT/'data/enumerations/fresh_h9_prism_completion_20261013_orbit_obstruction_point.json').read_text());targets.append(np.array([float(Fraction(x['a']))+float(Fraction(x['b']))*(1+np.sqrt(5))/2 for x in od['point']]))
 labels=(2,3,4,5,6,8,10,12);values=-2*np.cos(np.pi/np.array(labels));cov=R@G;target=targets[1];w=-np.linalg.solve(G,target)
 rows=[];stats={'faces':0,'trivalent_profiles':0,'line_solutions':0,'target_cutters':0,'seed_compatible':0,'at_most_two_optional_conflicts':0};seen=set();start=time.time()
 for fi,face in enumerate(d['spherical8_faces']):
  ids=np.array(face);A=R[ids];M=A@G@A.T;Mi=np.linalg.inv(M)
  proj=lambda v:v-(A@G@v)@Mi@A
  t=proj(center);qt=t@G@t
  if qt>=-1e-8:continue
  t/=np.sqrt(-qt);e=proj(w);e+=(e@G@t)*t;qe=e@G@e
  if qe<=1e-12:continue
  e/=np.sqrt(qe);stats['faces']+=1
  exterior=[j for j in range(len(R))if j not in face]
  for inds in itertools.combinations(range(8),3):
   p=np.zeros(8);p[list(inds)]=-1.;N=2-p@Mi@p
   if N<=1e-7:continue
   base=p@Mi@A;rho=np.sqrt(N);stats['trivalent_profiles']+=1
   for branch in(-1,1):
    aa=rho*(t+branch*e)/2;bb=rho*(-t+branch*e)/2
    pa,pb,pc=cov@aa,cov@base,cov@bb
    for j in exterior:
     for m,val in zip(labels,values):
      for z in rootsof(pa[j],pb[j]-val,pc[j]):
       if z<=0 or not np.isfinite(z):continue
       r=base+z*aa+bb/z;stats['line_solutions']+=1
       if r@h<=1e-7 or min(r@y for y in targets)>=-1e-7:continue
       stats['target_cutters']+=1
       pp=cov@r;ok=(pp<-2-1e-7)|(np.min(abs(pp[:,None]-values),axis=1)<1e-6)
       if not all(ok[:16]):continue
       stats['seed_compatible']+=1
       bad=np.flatnonzero(~ok).tolist()
       if len(bad)>2:continue
       stats['at_most_two_optional_conflicts']+=1
       key=tuple(np.round(r,7))
       if key in seen:continue
       seen.add(key)
       rows.append({'face_index':fi,'edge':face,'nonright_positions':list(inds),'anchor':j,'anchor_label':int(m),'branch':branch,'parameter':float(z),'normal':r.tolist(),'norm':float(r@G@r),'target_cuts':[-float(r@y)for y in targets],'conflicts':bad,'finite_old':[names[i]for i in np.flatnonzero(pp>-2+1e-6)]})
  if fi%25==0:print(fi,stats,'accepted',len(rows),'seconds',round(time.time()-start,1),flush=True)
 out={'status':'BOUNDED_NUMERICAL_TRIVALENT_ENDPOINT_SCOUT','source':PRE+'_faces_scout.json','stats':stats,'labels':list(labels),'candidates':rows,'scope':'All numerically found old102 spherical8 subsets. New node has exactly3 nonright edge pairings, allpi3, and one old finite anchor. Numeric proposals only; at most2 optional old walls may be dropped. Both known escapes tested.'}
 (ROOT/(PRE+'_trivalent_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',stats,'accepted',len(rows),flush=True)
if __name__=='__main__':run()
