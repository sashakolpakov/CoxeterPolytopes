#!/usr/bin/env python3
"""Bounded finite-contact continuation and cross-bridge checks after replacing150."""
import argparse,hashlib,json,math,collections
import numpy as np
import fresh_h20_replacement150_boundary_20261017 as A
O,L,ROOT=A.O,A.L,A.ROOT
TAG='fresh_h20_replacement150_continuation_20261017'
TRI='results/fresh_h20_replacement150_triangles_20261017.json'
CASES=[([152,172],176,152),([173,174],176,174),([152,176],174,176),([174,176],172,None),([171,172],173,171)]
def vertices(b,m):
 V=np.array([[1]+x['point']for x in b['poly_vertices']if x['kind']=='finite']);return V/np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
def catalogue(Y,CT,J,m,T,edge):
 whole=Y@m@Y.T+CT@J@CT.T;W=T+edge;c=O.census(W,whole[np.ix_(W,W)]);c['edge']=edge
 rows=list(c['profiles'])+[dict(id=f'single{j}_{lab}',contacts=[[j,lab]])for j in c['parametric_singletons']for lab in range(3,13)]
 return c,rows

def run():
 D,F,Y,CT,J,q,m,T,ret=A.geometry();ret=[j for j in ret if j!=150];R=O.P.P.load(TRI)['candidates'][20];C=O.P.P.load(f'results/{A.TAG}.json')['candidates'][14]['recipe'];assert R['profile']['contacts']==[]and C['contacts']==[]and C['label']==3
 for z in[R,C]:Y=np.r_[Y,[z['residual4']]];CT=np.r_[CT,[z['T_pairings']]];q=np.r_[q,2]
 ret += [175,176];G=dict(T=T,m=m);base=next(x['boundary']for x in R['compatible']if x['index']==14);V=vertices(base,m);st=collections.Counter();out=[];near=[];cens=[]
 for ci,(edge,anchor,shared)in enumerate(CASES):
  c,rows=catalogue(Y,CT,J,m,T,edge);cens.append(c)
  for row in rows:
   for sg,p,a,b,dt,off in O.conic(c,row,G,Y,CT,J,q):
    aa=a@m@Y[anchor];bb=b@m@Y[anchor];cc=p@m@Y[anchor]+off[anchor]
    for lab in range(2,7):
     for k in O.B.roots(aa,cc+math.sqrt(2*q[anchor])*(0 if lab==2 else math.cos(math.pi/lab)),bb):
      st['real']+=1;r=p+k*a+b/k
      if max(r@m@V.T)>=-1e-8:st['oldV']+=1;continue
      YY=np.r_[Y,[r]];CC=np.r_[CT,[dt]];full=YY@m@YY.T+CC@J@CC.T
      if shared is not None and min(np.linalg.eigvalsh(full[np.ix_(T+[shared,anchor,177],T+[shared,anchor,177])]))<=1e-8:st['joint_not_PD']+=1;continue
      co=-(r@m@Y.T+off)/np.sqrt(2*q);labs=[O.R.classify(float(x))for x in co];bad=[j for j in ret if labs[j]in('obtuse','angle_gap')]
      rec=dict(case=ci,contacts=row['contacts'],branch=sg,label=lab,k=float(k),residual4=r.tolist(),T_pairings=dt.tolist(),bad=[dict(wall=j,cosine=float(co[j]),kind=labs[j])for j in bad])
      if bad:
       st['gate']+=1
       if len(bad)<=2:near.append(rec)
       continue
      bnd=L.boundary(YY,m,ret+[177],T);rec['boundary']=bnd;out.append(rec);st['ok']+=1
  print('CASE',ci,dict(st),flush=True)
 # The next wall must also meet the existing opposite edge173/174. This
 # tests finite cross-bridges to each right-T accepted continuation177.
 c,rows=catalogue(Y,CT,J,m,T,[173,174]);forms=[(row,form)for row in rows for form in O.conic(c,row,G,Y,CT,J,q)]
 select=[(i,z)for i,z in enumerate(out)if not any(j in T for j,lab in z['contacts'])];bs=collections.Counter();bridges=[];bridge_trials=[]
 for ix,(si,seed)in enumerate(select):
  r177=np.array(seed['residual4']);dt177=np.array(seed['T_pairings']);YY=np.r_[Y,[r177]];CC=np.r_[CT,[dt177]];qq=np.r_[q,2];rr=ret+[177];VV=vertices(seed['boundary'],m)
  for row,(sg,p,a,b,dt,off)in forms:
   aa=a@m@r177;bb=b@m@r177;cc=p@m@r177+dt@J@dt177
   for lab in range(2,7):
    for k in O.B.roots(aa,cc+2*(0 if lab==2 else math.cos(math.pi/lab)),bb):
     bs['real']+=1;r=p+k*a+b/k
     if max(r@m@VV.T)>=-1e-8:bs['oldV']+=1;continue
     co=-(r@m@YY.T+CC@J@dt)/np.sqrt(2*qq);ls=[O.R.classify(float(x))for x in co];bad=[j for j in rr if ls[j]in('obtuse','angle_gap')]
     bridge_trials.append(dict(seed=si,profile=row['id'],branch=sg,k=float(k),label=lab,bad=[dict(wall=j,cosine=float(co[j]),kind=ls[j])for j in bad]))
     if bad:bs['gate']+=1;continue
     bnd=L.boundary(np.r_[YY,[r]],m,rr+[178],T);bridges.append(dict(seed=si,contacts=row['contacts'],branch=sg,k=float(k),label=lab,residual4=r.tolist(),T_pairings=dt.tolist(),boundary=bnd));bs['ok']+=1
  if ix%10==0:print('BRIDGE',ix,dict(bs),flush=True)
 # Cost of deleting the next sole blocker171. This is a diagnostic only.
 release=[j for j in ret if j!=171];delete=L.boundary(Y,m,release,T);rel=[]
 for i,z in enumerate(near):
  if [u['wall']for u in z['bad']]!=[171]:continue
  r=np.array(z['residual4']);dt=np.array(z['T_pairings']);co=-(r@m@Y.T+CT@J@dt)/np.sqrt(2*q);ls=[O.R.classify(float(x))for x in co];assert all(ls[j]not in('obtuse','angle_gap')for j in release)
  rel.append(dict(near_index=i,boundary=L.boundary(np.r_[Y,[r]],m,release+[177],T)))
 paths=[f'scripts/{TAG}.py',TRI,f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',A.PAIR,A.BELT,f'scripts/{O.TAG}.py',f'scripts/{L.TAG}.py']
 return dict(status='BOUNDED_NUMERICAL_REPAIRED150_BELT_AND_CROSS_BRIDGE_DIAGNOSTIC',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},retained_ids=ret,base_replacement_index=20,base_continuation_index=14,base_boundary=base,cases=[list(z)for z in CASES],censuses=cens,statistics=dict(st),candidates=out,at_most_two_bad_gates=near,cross_bridge=dict(census=c,seed_indices=[i for i,z in select],statistics=dict(bs),trials_passing_old_vertex_gate=bridge_trials,candidates=bridges),release171=dict(retained_ids=release,deletion_only=delete,candidates=rel),scope='Numerical finite-contact families only: full fixed endpoint catalogues plus singleton labels3..12 and anchor labels2..6. The cross-bridge tests all67 right-T survivors. Full native offset gates and old-finite-vertex strictness are necessary filters; all accepted continuations receive a complete numerical H3 boundary recount. Releasing171 is separately labelled as a replacement diagnostic. No exact promotion, all-label negative claim, or finite-volume conclusion.')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');z=a.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if z.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(p.read_text()))
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],d['cross_bridge']['statistics'],flush=True)
