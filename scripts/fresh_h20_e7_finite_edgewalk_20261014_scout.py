#!/usr/bin/env python3
"""Complete finite attachment templates; bounded numerical endpoint scout."""
import json,math,gzip
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import numpy as np
import fresh_h20_e7_residual_gate_20261014_exact_data as X

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_finite_edgewalk_20261014_scout'

def q5sign(a,b):
 if not b:return (a>0)-(a<0)
 if not a or (a>0)==(b>0):return (b>0)-(b<0)
 n=a*a-5*b*b;assert n
 return ((a>0)-(a<0))*((n>0)-(n<0))

def run():
 data=X.load();parent=data['parent'];S=data['S'];face=parent['face_facet_walls'];K=X.K
 N=np.array([[X.value(c)for c in v]for v in data['N']]);Q=np.array(data['Q'],float)
 H=np.array(list(map(X.value,data['H'])));E=np.array(list(map(X.value,data['E'])))
 W=np.array([[X.value(c)for c in v]for v in data['W']]);norms=np.array(list(map(float,data['base']['norms'])))
 G7=K.inverse(data['I']);FG=[[F(x)for x in row]for row in parent['face_gram']]
 scale2={i:F(2)if i in face[11:13]else F(1)for i in S+face}
 fullids=S+face;sp={i:j for j,i in enumerate(S)};fp={i:j for j,i in enumerate(face)}
 def gram(i,j):
  if i in sp and j in sp:return G7[sp[i]][sp[j]]
  if i in fp and j in fp:return FG[fp[i]][fp[j]]
  return F(0)
 mag={3:1.,4:math.sqrt(2),5:(1+math.sqrt(5))/2,6:math.sqrt(3)}
 squares={3:(F(1),F(0)),4:(F(2),F(0)),5:(F(3,2),F(1,2)),6:(F(3),F(0))}
 def components(ids):
  rem=set(ids);out=[]
  while rem:
   seen={min(rem)};todo=list(seen)
   for i in todo:
    for j in sorted(rem-seen):
     if gram(i,j):seen.add(j);todo.append(j)
   out.append(sorted(seen));rem-=seen
  return out
 edges=[];hits=[];stats=Counter();templates_total=0;tails_total=0;trials=0
 for vi,row in enumerate(parent['face_finite_vertices']):
  active=row['ambient_incident'];B=[[gram(i,j)for j in active]for i in active];BI=K.inverse(B)
  V=np.array([X.value(c)for c in K.vector(row['ambient_point'])]);q=-float(F(row['norm']));v=V/math.sqrt(q)
  for jj,j in enumerate(S):
   pos=active.index(j);ids=[i for i in active if i!=j];ii=[k for k in range(20)if k!=pos]
   inv=[[BI[a][b]-BI[a][pos]*BI[pos][b]/BI[pos][pos]for b in ii]for a in ii]
   diag=[scale2[i]*inv[k][k]for k,i in enumerate(ids)]
   assert min(diag)>=F(1,2)
   scales=np.sqrt([float(scale2[i])for i in ids]);invnp=np.array(inv,float)*np.outer(scales,scales)
   comps=components(ids);which={i:k for k,c in enumerate(comps)for i in c};tails=[i for c in comps if len(c)==1 for i in c]
   profiles=[dict(contacts=[],cost=['0','0'],residual=['2','0'])]
   for sz in [1,2,3]:
    for positions in combinations(range(19),sz):
     if len({which[ids[i]]for i in positions})!=sz:continue
     if sum(diag[i]for i in positions)>=2:continue
     for labels in product([3,4,5,6],repeat=sz)if sz<3 else[(3,3,3)]:
      ca=sum(diag[i]*squares[m][0]for i,m in zip(positions,labels));cb=sum(diag[i]*squares[m][1]for i,m in zip(positions,labels))
      if q5sign(2-ca,-cb)<=0:continue
      profiles.append(dict(contacts=[[ids[i],m]for i,m in zip(positions,labels)],cost=list(map(str,[ca,cb])),residual=list(map(str,[2-ca,-cb]))))
   degree=float(data['I'][jj][jj]);e=W[jj]/math.sqrt(degree)
   null=v-e;nullpairs=N@Q@null;nonshared=[i for i in range(87)if i not in ids]
   edgekind='open_numeric'if max(nullpairs[nonshared])<=1e-8 else'old_wall_closes_before_null_numeric'
   # The m>=7 tails remain symbolic. Only these labels are sampled here.
   sampled=profiles+[dict(contacts=[[i,m]],tail=True)for i in tails for m in [7,8,10,12,16,24,30,60]]
   erec=dict(vertex=vi,released=j,support=ids,components=comps,finite_profiles=profiles,
             I2_tails=[dict(wall=i,labels='all integers m>=7')for i in tails],edge_status=edgekind)
   edges.append(erec);templates_total+=len(profiles);tails_total+=len(tails);stats[edgekind]+=1
   for pi,pf in enumerate(sampled):
    b=np.zeros(19)
    for i,m in pf['contacts']:b[ids.index(i)]=-2*math.cos(math.pi/m)
    coeff=invnp@b;p=coeff@N[ids];cost=b@coeff;res=2-cost
    if res<=1e-10:continue
    sh=math.sqrt(res);pj=float(p@Q@N[j])
    for anchor in [2,3,4,5,6,7,8,10,12,16,24,30,60,'parallel']:
     cc=2. if anchor=='parallel'else 2*math.cos(math.pi/anchor)
     ch=(pj+cc)*math.sqrt(degree)/sh
     if ch<=1+1e-10:continue
     sn=math.sqrt(ch*ch-1);r=p+sh*(sn*v-ch*e);trials+=1
     pair=N@Q@r
     if max(pair)>2e-7 or r@Q@H>=-1e-8:continue
     co=-pair/np.sqrt(2*norms);labels=[];fail=False
     for z in co:
      if abs(z)<1e-7:labels.append(2)
      elif z>=1-1e-8:labels.append('parallel_or_dotted')
      elif z<.5-1e-7:fail=True;break
      else:
       m=math.pi/math.acos(min(z,1));mi=round(m)
       if abs(m-mi)>2e-6:fail=True;break
       labels.append(mi)
     if fail:continue
     point=v-(sn/ch)*e;xp=N@Q@point
     if max(xp)>2e-8 or any(xp[i]>=-1e-9 for i in nonshared):continue
     ans=dict(vertex=vi,released=j,support=ids,profile=pf,anchor=anchor,normal=r.tolist(),point=point.tolist(),
              normal_norm=float(r@Q@r),point_norm=float(point@Q@point),escape_pairing=float(r@Q@E),all_labels=labels)
     hits.append(ans);print('NUMERIC_ENDPOINT',vi,j,pf['contacts'],anchor,'cuts',ans['escape_pairing'],flush=True)
 print('counts',len(edges),templates_total,tails_total,trials,len(hits),dict(stats),flush=True)
 return dict(status='EXACT_TEMPLATE_SCHUR_TESTS_WITH_BOUNDED_NUMERICAL_ENDPOINT_SCOUT',
             sources={p:X.digest(p)for p in [X.BASE,X.PARENT,X.FACE,f'scripts/{TAG}.py','scripts/fresh_h20_e7_residual_gate_20261014_exact_data.py']},
             counts=dict(edges=len(edges),finite_pool_profiles=templates_total,symbolic_I2_tails=tails_total,numerical_trials=trials,numerical_endpoints=len(hits),edge_status=dict(stats)),
             edge_templates=edges,candidates=hits,
             scope='Finite spherical contact templates are enumerated with exact rational inverse diagonals and Qsqrt5 Schur signs; symbolic I2 tails are retained. Numerical conic scout fixes the omitted-wall contact to a bounded label list. It is not a full all-parameter or all-label endpoint exclusion. No candidate exactified by this script.')

if __name__=='__main__':
 d=run()
 with (ROOT/f'results/{TAG}.json.gz').open('wb') as out:
  with gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as packed:
   packed.write((json.dumps(d,indent=2)+'\n').encode())
