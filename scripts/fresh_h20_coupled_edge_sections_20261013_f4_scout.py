#!/usr/bin/env python3
"""Rational candidate-center census plus numerical bad-vertex proposals."""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations
import json,hashlib,math
import numpy as np
import fresh_h20_cusp_sections_20261013 as P
ROOT=P.ROOT;TAG='fresh_h20_coupled_edge_sections_20261013'
INPUT=f'results/{TAG}_scout.json';OUT=ROOT/f'results/{TAG}_f4_scout.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 raw=json.loads((ROOT/P.RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in r]for r in raw['basis_gram']]
 W=[list(map(int,r))for r in raw['wall_coordinates_as_columns_transposed']];C=[P.mv(Q,w)for w in W];G=[[P.dot(w,c)for c in C]for w in W]
 prior=json.loads((ROOT/INPUT).read_text());walk=json.loads((ROOT/P.WALK).read_text());h=json.loads(P.OUT.read_text())['strict_center']
 Qn=np.array(Q,float);Wn=np.array(W,float);Cn=np.array(C,float);Gn=np.array(G,float)
 rows=[];counts=[]
 for si,seed in enumerate(prior['selected_edges']):
  if seed['template']!='A1_A3_to_F4':continue
  E=seed['edge_walls'];cc=seed['component'];I=P.inv([[G[i][j]for j in E]for i in E]);centers={};candidates=[]
  for l in range(210):
   if l in E or any(G[l][i]for i in cc):continue
   v=P.mv(I,[G[i][l]for i in E]);p=[F(W[l][d])-sum(v[k]*W[j][d]for k,j in enumerate(E))for d in range(22)]
   q=P.pair(p,p,Q);assert any(p)
   record={'wall':l,'projection':list(map(str,p)),'norm':str(q)}
   if q>=0:record['status']='NOT_TIMELIKE'
   else:
    p=P.scale(-1/P.pair(p,h,Q),p);ps=[P.dot(p,c)for c in C]
    if any(ps[i]>=0 for i in range(210)if i not in E):record['status']='OUTSIDE_OPEN_EDGE'
    else:record['status']='INTERIOR_CENTER';centers.setdefault(tuple(p),[]).append(l)
   candidates.append(record)
  c=walk['corner_orbits'][seed['corner']];edge=c['edge_orbits'][seed['edge_index']]
  right=np.array(edge['endpoint_primitive_ray'],float);right/=-(right@Qn@np.array(h,float))
  In=np.array(I,float);a=np.zeros(20)
  for i,coef,rad in seed['pairings']:a[E.index(i)]=float(F(coef))*math.sqrt(float(F(rad)))
  v=In@a;w=v@Wn[E];gamma=1-a@v;assert abs(gamma-.05)<1e-8
  current=[]
  for point,ls in centers.items():
   pn=np.array(point,float);qp=pn@Qn@pn;z0=right-(right@Qn@pn/qp)*pn;qz=z0@Qn@z0;pp=Cn@pn
   for sg in[-1,1]:
    z=sg*math.sqrt(gamma/qz)*z0;e=w+z;ep=Cn@e;duals=In@Wn[E]-np.outer(v/gamma,z);bads=[]
    for k,release in enumerate(E):
     dp=Cn@duals[k];bounds=[(pp[j]/dp[j],j)for j in range(210)if j not in E and dp[j]<-1e-9]
     if not bounds:continue
     tt,l=min(bounds);x=pn-tt*duals[k];xp=Cn@x;active=np.flatnonzero(abs(xp)<1e-7).tolist()
     if x@Qn@x>=-1e-8 or len(active)!=20 or max(xp)>1e-7:continue
     B=Gn[np.ix_(active,active)]-np.outer(ep[active],ep[active])
     if np.linalg.eigvalsh(B)[0]<1e-7:continue
     for ii,jj in combinations(range(20),2):
      p=B[ii,jj];r=p*p/(B[ii,ii]*B[jj,jj]);kind=None;score=0
      if p>1e-7:kind='obtuse';score=p
      elif p<-1e-7 and 1e-5<r<.25-1e-5:kind='between_0_and_1/4';score=min(r,.25-r)
      elif p<-1e-7 and .25+1e-5<r<.5-1e-5:kind='between_1/4_and_1/2';score=min(r-.25,.5-r)
      if kind:bads.append({'release':release,'new_wall':l,'active':active,'bad_pair':[active[ii],active[jj]],'kind':kind,'score':float(score)})
    record={'seed':si,'point':list(map(str,point)),'possible_extra_F4_face_walls':ls,'radical_sign':sg}
    if bads:record['witness']=max(bads,key=lambda v:v['score'])
    else:record['status']='UNRESOLVED'
    current.append(record)
  rows+=current;counts.append({'seed':si,'old_wall_projection_census':candidates,'distinct_centers':len(centers),'branches':len(current)})
  print(si,len(centers),len(current),sum('witness'not in r for r in current),flush=True)
 out={'status':'EXACT_RATIONAL_CENTER_ENUMERATION_WITH_NUMERICAL_VERTEX_WITNESSES',
  'sources':{p:sha(p)for p in[INPUT,P.RAW,P.WALK,'results/'+P.TAG+'.json','scripts/'+Path(__file__).name]},
  'families':counts,'cases':rows}
 OUT.write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':run()
