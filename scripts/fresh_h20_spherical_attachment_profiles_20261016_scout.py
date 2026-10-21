#!/usr/bin/env python3
"""Bounded full-ambient attachment scout; all-parameter old-edge necessary tests."""
import argparse,hashlib,itertools,json
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_spherical_attachment_profiles_20261016_scout';INPUT='results/fresh_h20_spherical_attachment_profiles_20261016.json'
def load(p):return json.loads((ROOT/p).read_text())
def geometry():
 D=load(INPUT);g=D['geometry'];X=np.array(g['residual_normals']);M=np.array(g['residual_metric']);q=np.array(g['native_norms']);CS=np.array(g['old_support_pairings']);J=np.array(g['support_inverse']);S=g['support'];T=S+[151]
 H=np.array(g['weak_residual_center']);w=X[151];h=H-(H@M@w)/2*w;h/=np.sqrt(-h@M@h)
 Z=null_space(np.array([w@M,h@M])).T;ev,V=np.linalg.eigh(Z@M@Z.T);assert min(ev)>0
 F=np.r_[[h],(V/np.sqrt(ev)[None,:]).T@Z];m=np.diag([-1.,1.,1.,1.]);assert np.max(abs(F@M@F.T-m))<1e-7
 Y=(X@M@F.T)@m
 # Norms are evaluated from native norms minus positive-space projections.
 c151=X@M@w;qt=q-np.einsum('ij,jk,ik->i',CS,J,CS)-c151*c151/2
 # Recompute the seven small positive projections without cancellation from
 # large old roots, using independent old actual-facet list.
 sides=[95,142,145,150,152,156,157]
 for j in sides:qt[j]=Y[j]@m@Y[j]
 a=Y[:,1:];b=Y[:,0];live=[j for j in range(158)if j not in T and np.linalg.norm(a[j])>1e-9]
 aa=a[live];bb=b[live];lens=np.linalg.norm(aa,axis=1)
 lp=linprog([0,0,0,-1],A_ub=np.c_[aa/lens[:,None],np.ones(len(live))],b_ub=bb/lens,bounds=[(None,None)]*3+[(0,None)],method='highs');assert lp.success and lp.x[3]>0
 verts=HalfspaceIntersection(np.c_[aa,-bb],lp.x[:3]).intersections
 return D,dict(X=X,M=M,q=q,CS=CS,J=J,S=S,T=T,F=F,m=m,Y=Y,qt=qt,c151=c151,sides=sides,A=aa,b=bb,ids=live,vertices=verts,center=lp.x[:3])
def poly(G,y):
 A=np.r_[G['A'],[y[1:]]];b=np.r_[G['b'],y[0]];lens=np.linalg.norm(A,axis=1)
 lp=linprog([0,0,0,-1],A_ub=np.c_[A/lens[:,None],np.ones(len(A))],b_ub=b/lens,bounds=[(None,None)]*3+[(0,None)],method='highs')
 if not lp.success or lp.x[3]<1e-10:return dict(empty=True)
 vv=HalfspaceIntersection(np.c_[A,-b],lp.x[:3]).intersections;r2=np.sum(vv*vv,axis=1)
 return dict(vertices=len(vv),finite=int(np.sum(r2<1-1e-6)),exterior=int(np.sum(r2>1+1e-6)),near_sphere=int(np.sum(abs(r2-1)<=1e-6)),maximum_squared_radius=float(max(r2)))
def offsets(G,row):
 ds=np.zeros(16);d151=0
 for j,m in row['contacts']:
  d=2*np.cos(np.pi/m)
  if j==151:d151=d
  else:ds[G['S'].index(j)]=d
 return -G['CS']@G['J']@ds-d151/2*G['c151']
def gates(G,ys,offset):
 pairs=ys@G['m']@G['Y'].T+offset
 co=-pairs/np.sqrt(2*G['q'])
 co[:,G['T']]=0 # prescribed finite profile is checked separately
 ok=(abs(co)<2e-6)|(co>=1-2e-6)
 for m in [3,4,5,6,7,8,9,10,12,16,20,24]:ok|=abs(co-np.cos(np.pi/m))<2e-6
 return np.all(ok,axis=1),co

def run():
 D,G=geometry();rows=[];candidates=[];stats=Counter();labelset=[2,3,4,5,6,8,10]
 exterior=G['vertices'][np.sum(G['vertices']**2,axis=1)>1+1e-6]
 for row in D['census']['profiles']:
  h=row['numeric_residual_norm'];offset=offsets(G,row);domains={};ratios={}
  for j in G['sides']:
   ratios[str(j)]=float(offset[j]**2/(h*G['qt'][j]))
   dom=[]
   for m in labelset:
    target=-np.sqrt(2*G['q'][j])*(0 if m==2 else np.cos(np.pi/m));p=target-offset[j]
    if p<1e-8 and p*p<h*G['qt'][j]-1e-10:dom.append(m)
   if dom:domains[j]=dom
  rows.append(dict(profile=row['id'],contacts=row['contacts'],minimum_projection_ratio=ratios,finite_old_edge_label_domains={str(k):v for k,v in domains.items()},eligible_old_neighbors=len(domains)))
  for ids in itertools.combinations(domains,3):
   U=G['Y'][list(ids)];gram=U@G['m']@U.T
   if abs(np.linalg.det(gram))<1e-12:stats['singular_frames']+=1;continue
   inv=np.linalg.inv(gram);ker=null_space(U@G['m'])[:,0];kn=ker@G['m']@ker
   for labels in itertools.product(*(domains[j]for j in ids)):
    f=np.array([-np.sqrt(2*G['q'][j])*(0 if m==2 else np.cos(np.pi/m))-offset[j]for j,m in zip(ids,labels)])
    p=(inv@f)@U;rad=(h-p@G['m']@p)/kn;stats['patterns']+=1
    if rad < -1e-9:continue
    ys=np.array([p+sgn*np.sqrt(max(rad,0))*ker for sgn in [-1,1]])
    mask,co=gates(G,ys,offset)
    for i in np.where(mask&(ys[:,0]>0))[0]:
     y=ys[i];cuts=exterior@y[1:]-y[0];stats['accepted']+=1
     candidates.append(dict(profile=row['id'],contacts=row['contacts'],extra_old_sides=list(ids),extra_labels=list(labels),branch=[-1,1][i],residual4=y.tolist(),offset=offset.tolist(),old_exterior_cuts=np.where(cuts>1e-7)[0].tolist(),face=poly(G,y)))
  print('PROFILE',row['id'],'oldneighbors',len(domains),'accepted',len(candidates),flush=True)
 # Unrestricted all-other-dotted rays: use ideal directions inside the old
 # H3 cone and solve the exact-in-k quadratic inequalities numerically.
 directions=[]
 for z in exterior:
  v=z/np.linalg.norm(z)
  if not any(np.linalg.norm(v-w)<1e-7 for w in directions):directions.append(v)
 all_dotted=[];notT=[j for j in range(158)if j not in G['T']]
 for row in [D['census']['profiles'][i]for i in D['census']['E8_profiles']]:
  h=row['numeric_residual_norm'];offset=offsets(G,row)
  for vi,v in enumerate(directions):
   aa=np.sqrt(h)/2*(np.r_[1,v]@G['m']@G['Y'].T);bb=np.sqrt(h)/2*(np.r_[-1,v]@G['m']@G['Y'].T);cc=offset+np.sqrt(2*G['q'])
   if max(aa[notT])>=-1e-10:continue
   k=1.
   for j in notT:
    disc=cc[j]**2-4*aa[j]*bb[j]
    if disc>0:
     roots=np.roots([aa[j],cc[j],bb[j]]);k=max(k,float(max(roots)))
   k*=1.000001
   y=np.sqrt(h)/2*np.r_[k-1/k,(k+1/k)*v];mask,co=gates(G,y[None,:],offset)
   assert mask[0],(row['id'],vi,min(co[0,notT]))
   all_dotted.append(dict(profile=row['id'],direction=vi,k=float(k),residual4=y.tolist(),minimum_nonsupport_cosine=float(min(co[0,notT])),face=poly(G,y)))
 src=[INPUT,f'scripts/{TAG}.py']
 return dict(status='BOUNDED_FULL_AMBIENT_SPHERICAL_PROFILE_SCOUT',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},old_face=dict(vertices=G['vertices'].tolist(),maximum_squared_radius=float(max(np.sum(G['vertices']**2,axis=1)))),profile_old_edge_domains=rows,label_domain=labelset,statistics=dict(stats),candidates=candidates,E8_all_other_dotted_caps=all_dotted,geometry=dict(residual4_basis=G['F'].tolist(),residual4_metric=G['m'].tolist()),scope='Numerical discovery only.57finiteprofiles,7bounded additional old finite labels, three actual old finite-neighbor contacts. E8all-dotted rays use finitely many old-exterior directions. Oldedge minimum projection test is an all-angle necessary inequality, pending exact reconstruction; no global closure exclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  old=load(str(out.relative_to(ROOT)))
  from fresh_h20_release144_cap_20261016_residual_scout import close
  assert close(d,old)
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],'E8rays',len(d['E8_all_other_dotted_caps']),flush=True)
