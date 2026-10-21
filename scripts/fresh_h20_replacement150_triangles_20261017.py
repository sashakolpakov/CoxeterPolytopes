#!/usr/bin/env python3
"""Finite triangular replacements of150; numerical ambient gates."""
import argparse,hashlib,json,math,collections
import numpy as np
from scipy.linalg import null_space
import fresh_h20_replacement150_boundary_20261017 as A
O,L,ROOT=A.O,A.L,A.ROOT
TAG='fresh_h20_replacement150_triangles_20261017'
def run():
 D,F,Y,CT,J,q,m,T,ret=A.geometry();ret=[j for j in ret if j!=150];G=dict(T=T,m=m);whole=Y@m@Y.T+CT@J@CT.T
 cs=[]
 for edge in [[142,145],[142,152],[145,152]]:
  W=T+edge;c=O.census(W,whole[np.ix_(W,W)]);c['edge']=edge;cs.append(c);print('CENSUS',edge,len(c['profiles']),c['parametric_singletons'],flush=True)
 rows=O.P.profiles(dict(endpoint_censuses=cs));assert len(rows)==63;print('INTERSECTION',len(rows),rows,flush=True)
 common=set.intersection(*(set(c['parametric_singletons'])for c in cs));assert common=={151};rows += [dict(id=f'single{j}_{lab}',contacts=[[j,lab]])for j in sorted(common)for lab in list(range(3,25))+[32,64,128,256,1024]]
 U=Y[[142,145,152]];gram=U@m@U.T;assert abs(np.linalg.det(gram))>1;inv=np.linalg.inv(gram);ker=null_space(U@m)[:,0];kn=ker@m@ker
 B=json.load(open('results/fresh_h20_replacement150_boundary_20261017.json'));cont=B['candidates'];out=[];stats=collections.Counter();V=np.array([[1]+z['point']for z in B['base_boundary']['poly_vertices']if z['kind']=='finite'and 150 not in z['active']]);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
 for row in rows:
  dt=np.zeros(len(T));f=np.zeros(3)
  for j,lab in row['contacts']:
   mag=2*math.cos(math.pi/lab)
   if j in T:dt[T.index(j)]=-mag
   else:f[[142,145,152].index(j)]=-mag*math.sqrt(q[j]/2)
  offset=CT@J@dt;f-=offset[[142,145,152]];p=(inv@f)@U;h=2-dt@J@dt;rad=(h-p@m@p)/kn
  if rad<0:stats['negative_rad']+=1;continue
  for sign in [-1,1]:
   r=p+sign*math.sqrt(rad)*ker;co=-(r@m@Y.T+offset)/np.sqrt(2*q);labels=[O.R.classify(float(x))for x in co];bad=[j for j in ret if labels[j]in('obtuse','angle_gap')];stats['real']+=1
   if bad:stats['native_fail']+=1;continue
   assert max(r@m@V.T)<-1e-8
   assert abs(r@m@r+dt@J@dt-2)<1e-7
   YY=np.r_[Y,[r]];CC=np.r_[CT,[dt]];qq=np.r_[q,2];base=L.boundary(YY,m,ret+[175],T)
   rec=dict(profile=row,sign=sign,preserved_finite_vertex_products=(r@m@V.T).tolist(),residual4=r.tolist(),T_pairings=dt.tolist(),base=base,compatible=[])
   print('REPLACEMENT',row,sign,'base',base['poly_vertex_counts'],len(base['open_edges']),flush=True)
   for z in cont:
    rr=np.array(z['recipe']['residual4']);dd=np.array(z['recipe']['T_pairings']);cross=float(-(r@m@rr+dt@J@dd)/2);lab=O.R.classify(cross)
    if lab in('obtuse','angle_gap'):continue
    bd=L.boundary(np.r_[YY,[rr]],m,ret+[175,176],T);rec['compatible'].append(dict(index=z['index'],mutual=lab,cosine=cross,boundary=bd));print('PAIR',z['index'],lab,bd['poly_vertex_counts'],len(bd['open_edges']),bd['maximum_squared_radius'],flush=True)
   out.append(rec)
 print('DONE',dict(stats),len(out),flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',A.PAIR,A.BELT,f'scripts/{O.TAG}.py',f'scripts/{O.P.TAG}.py',f'scripts/{L.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_TRIANGULAR_REPLACEMENTS_OF_OLD150',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},censuses=cs,intersection=rows,frame_gram=gram.tolist(),kernel_square=float(kn),stats=dict(stats),candidates=out,scope='Complete fixed finite triangular-replacement profile intersection; singleton151 labels3..24,32,64,128,256,1024 are only a bounded sample. Both residual norm branches and all172 retained native gates tested numerically, followed by all18 recorded continuation caps and full target-H3 boundary recount. No exact promotion or all-label negative conclusion.')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');z=a.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if z.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(p.read_text()))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
