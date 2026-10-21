#!/usr/bin/env python3
"""Quartic coordinated E8-partner rotation / old-edge endpoint experiment.

No sampled rotation grid: both new-partner contacts determine a rotated
2-vector, and equality of its length gives a polynomial of degree at most4.
All geometry remains numerical, conditional on the frozen exact-promotable
162-wall recipe. The old-edge labels2..5 are a stated bounded domain.
"""
import argparse,hashlib,itertools,json,math
import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog
from scipy.spatial import HalfspaceIntersection
import fresh_h20_e8_partner_endpoints_20261017 as P
ROOT=P.ROOT;TAG=P.TAG+'_coordinated';INPUT='results/'+P.TAG+'.json'
OLD_EDGES=[(95,142),(95,157),(142,152),(145,152),(145,157)]
def poly(Y,T):
 ids=[i for i in range(len(Y))if i not in T and np.linalg.norm(Y[i,1:])>1e-9]
 A=Y[ids,1:];b=Y[ids,0];norm=np.linalg.norm(A,axis=1)
 lp=linprog([0,0,0,-1],A_ub=np.c_[A/norm[:,None],np.ones(len(ids))],b_ub=b/norm,bounds=[(None,None)]*3+[(0,None)],method='highs')
 if not lp.success or lp.x[3]<=1e-9:return None
 v=HalfspaceIntersection(np.c_[A,-b],lp.x[:3]).intersections;rr=np.sum(v*v,axis=1)
 return dict(vertices=len(v),finite=int(sum(rr<1-1e-6)),exterior=int(sum(rr>1+1e-6)),near_sphere=int(sum(abs(rr-1)<=1e-6)),maximum_squared_radius=float(max(rr)),center=lp.x[:3].tolist())
def run():
 D=P.load(INPUT);G,Y,CT,JT,q,H,live,verts=P.geometry();m=G['m'];r0=Y[158];W=Y[159:162];M=W@m@W.T
 a=-np.linalg.solve(M,np.ones(3));hh=a@W;hh/=math.sqrt(-(hh@m@hh));t=-W[0]@m@hh
 assert np.max(abs(W@m@hh+t))<1e-9
 E=null_space(np.array([r0@m,hh@m])).T;ee,V=np.linalg.eigh(E@m@E.T);E=(V/np.sqrt(ee)[None,:]).T@E
 spatial=W@m@E.T
 assert np.max(abs(t*hh+spatial@E-W))<1e-9
 off=CT[159:162]@JT@CT[:158].T
 stats={'old_edge_label_patterns':0,'quartics':0,'positive_real_roots':0,'coherent_Coxeter':0};records=[];near=[]
 notT=[j for j in range(158)if j not in G['T']]
 for oldedge in OLD_EDGES:
  U=Y[list(oldedge)];gi=U@m@U.T;inv=np.linalg.inv(gi)
  h0=np.array([1.,0,0,0]);h0-=((h0@m@U.T)@inv)@U;h0/=math.sqrt(-(h0@m@h0))
  e0=null_space(np.r_[U@m,[h0@m]])[:,0];e0/=math.sqrt(e0@m@e0)
  for labels in itertools.product(range(2,6),repeat=2):
   f=np.array([-math.sqrt(2*q[j])*(0 if lab==2 else math.cos(math.pi/lab))for j,lab in zip(oldedge,labels)])
   p=(inv@f)@U;h=2-p@m@p
   if h<=1e-9:continue
   stats['old_edge_label_patterns']+=1
   for sign in [-1,1]:
    va=math.sqrt(h)/2*(h0+sign*e0);vb=math.sqrt(h)/2*(-h0+sign*e0)
    # Laurent coefficient order is -1,0,1. Multiplication by k^2 after
    # squaring therefore turns convolution directly into degree0..4.
    rh=np.array([vb@m@hh,p@m@hh,va@m@hh]);rv=np.array([[vb@m@e,p@m@e,va@m@e]for e in E])
    for ci,cen in enumerate(D['endpoint_censuses']):
     pair=[j-159 for j in cen['edge']];SS=spatial[pair];si=np.linalg.inv(SS)
     for row in cen['profiles']:
      if any(j in G['T']for j,lab in row['contacts']):continue
      target=np.zeros(2)
      for j,lab in row['contacts']:target[cen['edge'].index(j)]=-2*math.cos(math.pi/lab)
      cc=-t*np.tile(rh,(2,1));cc[:,1]+=target
      rz=si@cc;pol=sum((np.convolve(v,v)for v in rz),np.zeros(5))-sum((np.convolve(v,v)for v in rv),np.zeros(5))
      scale=max(abs(pol));assert scale>1e-10;pol/=scale
      degree=4
      while degree and abs(pol[degree])<1e-11:degree-=1
      rr=np.roots(pol[:degree+1][::-1]);stats['quartics']+=1
      for root in rr:
       if abs(root.imag)>1e-7*(1+abs(root.real))or root.real<=0:continue
       k=float(root.real);r=p+k*va+vb/k;v=r@m@E.T;z=si@(target-t*(r@m@hh))
       if np.linalg.norm(v)<1e-9:continue
       angle=math.atan2(v[0]*z[1]-v[1]*z[0],v@z);rot=np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
       WW=t*hh+spatial@rot@E;stats['positive_real_roots']+=1
       assert max(abs(WW[pair]@m@r-target))<2e-5
       assert abs(r@m@r-2)<2e-5
       cosp=-(WW@m@Y[:158].T+off)/np.sqrt(2*q[:158])[None,:]
       cosr=-(r@m@Y[:159].T)/np.sqrt(2*q[:159])
       others=[j for j in range(159)if j not in G['T'] and j not in oldedge]
       remaining=next(j for j in range(3)if j not in pair);cross=-WW[remaining]@m@r/2
       # Entire unprescribed old contacts are required dotted in this lane.
       margin=min(float(min(cosp[:,notT].ravel())-1),float(min(cosr[others])-1),float(cross-1))
       why=dict(partner_old=float(min(cosp[:,notT].ravel())),new_old=float(min(cosr[others])),unpaired_partner=float(cross))
       rec=dict(old_edge=list(oldedge),old_labels=list(labels),branch=sign,new_edge=cen['edge'],endpoint_profile=row['id'],endpoint_contacts=row['contacts'],polynomial_degree=degree,polynomial_coefficients=pol.tolist(),positive_parameter=k,rotation_angle=angle,rotated_partners=WW.tolist(),new_residual4=r.tolist(),minimum_dotted_margin=margin,gate_minima=why)
       near.append(rec)
       if margin>=-1e-6:
        face=poly(np.r_[Y[:159],WW,[r]],G['T'])
        if face is not None:rec['face']=face;records.append(rec);stats['coherent_Coxeter']+=1
 print('QUARTIC',stats,flush=True)
 near.sort(key=lambda r:-r['minimum_dotted_margin'])
 src=[f'scripts/{TAG}.py',INPUT,P.RECIPE,f'scripts/{P.TAG}.py']
 return dict(status='BOUNDED_COORDINATED_ENDPOINT_QUARTIC_EXPERIMENT',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},fixed_E8_root=158,rotated_partners=[159,160,161],old_edges=[list(x)for x in OLD_EDGES],label_domain=list(range(2,6)),statistics=stats,candidates=records,all_solved_real_branches=near,scope='Every positive real solution of the numerical quartics in the stated finite old-label domain and all zero-T endpoint profiles is retained. This is not an exact root-isolation certificate; roots with tiny imaginary parts use tolerance. Both labels are necessarily2 on oldlabel6 edges by spherical3-node classification. Other profiles, unbounded oldI2labels, and more than one added cap remain outside scope. The E8 root is fixed but partners rotate, so this is an alternative to the frozen162 recipe, not a retaining extension.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,P.load(str(out.relative_to(ROOT))))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
