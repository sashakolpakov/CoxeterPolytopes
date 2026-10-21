#!/usr/bin/env python3
"""Free dotted self-image search for the only unreduced neighbor sphere.

All17 other finite-neighbor relations and the facet175 angle are copied from
old9. Every additional relation is required dotted, with at most one finite
anchor. Circle feasibility uses complete angular arcs, not point sampling.
This is a numerical proposal search, not a universal exact obstruction.
"""
import json,numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog,minimize
from fresh_h20_facet_bridge_20261013_scout import setup,ROOT
from fresh_h20_coupled_completion_20261013 import read,num
from fresh_h20_coupled_completion_20261013_second_exact import OUT as SOURCE
from fresh_dimension_bridge_two_angle_20261013_independent import B

OUT='data/enumerations/fresh_h20_changed_assembly_20261013_free_scout.json'

def arcs(a,b,c):
 intervals=[(0.,2*np.pi)]
 for aa,bb,cc in zip(a,b,c,strict=True):
  amp=np.hypot(aa,bb)
  if amp<1e-10:
   if cc<-1e-8:return []
   continue
  if cc>=amp-1e-10:continue
  if cc<-amp-1e-10:return []
  phi=np.arctan2(bb,aa);theta=np.arccos(np.clip(cc/amp,-1,1))
  # The forbidden arc is (phi-theta,phi+theta), periodically translated.
  for k in range(-2,3):
   lo,hi=phi-theta+2*k*np.pi,phi+theta+2*k*np.pi;new=[]
   for left,right in intervals:
    if hi<=left or lo>=right:new.append((left,right));continue
    if lo>left:new.append((left,min(lo,right)))
    if hi<right:new.append((max(hi,left),right))
   intervals=new
   if not intervals:return []
 return [(l,r)for l,r in intervals if r-l>1e-9]

Q,R,h,neighbors,*_=setup();cov=R@Q;qn=np.einsum('ij,jk,ik->i',R,Q,R);j=9;I=[i for i in neighbors if i!=j]+[175]
N=null_space(cov[I]);ev,E=np.linalg.eigh(N.T@Q@N);assert min(ev)>0;N=N@E@np.diag(1/np.sqrt(ev))
c=R[j]-N@(N.T@Q@R[j]);radius2=qn[j]-c@Q@c;L=N*np.sqrt(radius2)
J=[i for i in range(176)if i not in I];A=cov[J]@L;b=-np.sqrt(qn[j]*qn[J])-cov[J]@c
base=linprog(np.zeros(3),A_ub=A,b_ub=b,bounds=[(-1,1)]*3,method='highs')
print('all-dotted cube feasible',base.success,flush=True)
target=np.array([num(B,B.parse(x))for x in read(SOURCE)['strict_escape']['vector']]);tc=target@Q@c;tv=target@Q@L
labels=list(range(2,31));records=[];cutters=[]
for ii,k in enumerate(J):
 a=A[ii];aa=a@a
 if aa<1e-12:continue
 Z=null_space(a[None,:]);assert Z.shape==(3,2)
 for m in labels:
  beta=-np.sqrt(qn[j]*qn[k])*np.cos(np.pi/m)-cov[k]@c;u0=beta*a/aa;rho2=1-u0@u0
  if rho2<=1e-9:continue
  T=np.sqrt(rho2)*Z;others=np.arange(len(J))!=ii;AB=A[others]@T;bb=b[others]-A[others]@u0
  intervals=arcs(AB[:,0],AB[:,1],bb)
  if not intervals:continue
  cutAB=np.vstack([AB,-tv@T]);cutb=np.r_[bb,tc+tv@u0-1e-6]
  cut_intervals=arcs(cutAB[:,0],cutAB[:,1],cutb)
  theta=sum((cut_intervals or intervals)[0])/2;u=u0+T@np.array([np.cos(theta),np.sin(theta)]);r=c+L@u
  record={'anchor':k,'label':m,'arcs':intervals,'cutting_arcs':cut_intervals,'normal':r.tolist(),'norm':float(r@Q@r),
          'target_pairing':float(r@Q@target),'max_nondotted_residual':float(max(A[others]@u-b[others]))}
  records.append(record)
  if cut_intervals:cutters.append(record)
  print('arc',k,m,'cuts',bool(cut_intervals),'pair',record['target_pairing'],flush=True)
out={'status':'NUMERICAL_FREE_DOTTED_NEIGHBOR9_SPHERE_SEARCH','source':SOURCE,'fixed_rows':I,'copied_row':9,'radius_squared':float(radius2),
     'all_dotted_cube_feasible':bool(base.success),'anchor_labels':labels,'continuous_families':records,'target_cutting_families':cutters,
     'scope':'All17 other finite-neighbor pairings and175 relation copied; at most one additional finite anchor from labels2..30, all remaining158 relations dotted. Floating arc calculations only.'}
(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
