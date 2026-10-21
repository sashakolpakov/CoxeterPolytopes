#!/usr/bin/env python3
"""Small exact facet/edge diagnostic for the two sqrt22 caps.

The input gate certifies every wall outside the22 listed neighbors is dotted
to the cap. Such a wall cannot cut its intrinsic hyperbolic facet. This checks
all21 outgoing edges at its saved finite vertex, not a global chamber census.
"""
import argparse,json,hashlib
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from fresh_h22_incidence_mixed_height_20261013 import inertia
from fresh_h22_incidence_mixed_height_20261013_gap import roots
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_incidence_mixed_height_20261013'
BASE='data/enumerations/'+PRE+'.json';GAP='data/enumerations/'+PRE+'_gap.json';RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';OUT='data/enumerations/'+PRE+'_frontier.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),F(0))
def mv(A,v):return[dot(r,v)for r in A]
def pair(A,a,b):return dot(a,mv(A,b))
def inv(A):
 n=len(A);B=[[F(x)for x in row]+[F(i==j)for j in range(n)]for i,row in enumerate(A)]
 for j in range(n):
  k=next(i for i in range(j,n)if B[i][j]);B[j],B[k]=B[k],B[j];p=B[j][j];B[j]=[x/p for x in B[j]]
  for i in range(n):
   if i!=j:
    p=B[i][j];B[i]=[x-p*y for x,y in zip(B[i],B[j])]
 return [r[n:]for r in B]
def run():
 b=read(BASE);c=read(GAP)
 for d in(b,c):
  for p,s in d['sources'].items():assert sha(p)==s
 raw=read(RAW)['coxeter_input'];Q=[[-F(x)for x in r]for r in raw['basis_gram']];W=[[F(x)for x in r]for r in raw['wall_coordinates_as_columns_transposed']]
 u=list(map(F,c['new_old_projection']));h=list(map(F,b['center']));D=[[F(x)for x in r]for r in b['null_projections']];X=[[F(x)/2 for x in r]for r in b['sqrt3_twice_projections']]
 I=c['actual_finite_vertices']['old_incident'];G=[[pair(Q,W[i],W[j])for j in I]for i in I];J=inv(G)
 assert all(x>=0 for r in J for x in r)and inertia(G)==[21,0,0]
 assert [i for i,w in enumerate(W)if pair(Q,u,w)==-1]==[14]
 bvec=[pair(Q,W[i],W[14])for i in I];aa=mv(J,bvec)
 assert all(a<=0 for a in aa)
 edges=[];finite=[];duals=[]
 for n,i in enumerate(I):
  d=[sum(J[n][j]*W[k][t]for j,k in enumerate(I))for t in range(22)];duals.append(d)
  assert pair(Q,d,u)==0 and pair(Q,d,d)==J[n][n]>0
  assert [pair(Q,d,W[k])for k in I]==[F(i==k)for k in I]
  assert pair(Q,d,W[14])==aa[n]
  gap=11*J[n][n]-20*aa[n]*aa[n]
  kind='finite'if aa[n]<0 and gap<0 else'ideal'if aa[n]<0 and gap==0 else'open'
  row={'released_old_wall':i,'dual_square':str(J[n][n]),'wall14_dual_pair':str(aa[n]),'comparison_11d2_minus20a2':str(gap),'status':kind}
  if kind=='finite':
   t=-1/aa[n];x=[a-t*z for a,z in zip(u,d)];q=pair(Q,x,x)+F(200,11);assert q<0 and pair(Q,x,h)<0
   oldp=[pair(Q,x,w)for w in W];zero=[j for j,p in enumerate(oldp)if p==0]
   expect=sorted([k for k in I if k!=i]+[14]);assert zero==expect and max(oldp)<=0
   dp=max(pair(Q,x,w)for w in D);xp=max(pair(Q,x,w)for w in X)
   assert dp<0 and xp<0
   a_lower=(roots(6)[0]+roots(3)[0])/2
   assert a_lower*dp+20/roots(11)[0]<0 and xp+F(10,11)*roots(66)[1]<0
   B=[[pair(Q,W[j],W[k])for k in expect]+[pair(Q,u,W[j])]for j in expect]
   B.append([pair(Q,u,W[k])for k in expect]+[F(2)])
   assert inertia(B)==[22,0,0]
   row['point_old_coordinates']=list(map(str,x));row['point_height']='+/-20/sqrt22';row['point_square']=str(q);row['old_incident']=expect
   finite.append(i)
  edges.append(row)
 assert finite==[12,13,15] and sum(r['status']=='open'for r in edges)==18
 assert all(r['status']!='ideal'for r in edges)
 n=I.index(1);d=duals[n];assert J[n][n]==8 and aa[n]==0
 # eta=(u-sqrt(5/22)*d, +/-20/sqrt22) is future null.
 assert F(-20)+F(5,22)*8+F(200,11)==0
 # Old retained20 plus the newcap vanish; released1 isstrict, wall14=-1.
 assert pair(Q,d,W[14])==0
 return {'status':'PASS_EXACT_SQRT22_CAP_FACET_THREE_FINITE_EIGHTEEN_OPEN_EDGES','sources':{p:sha(p)for p in(BASE,GAP,RAW)},'scope':'Identical intrinsic-facet calculation for both newcaps and allthree1336 branches. It covers the21 edges from the saved capvertex, not the whole chamber frontier.','finite_neighbors':I+[14],'start_point':'(u,+/-20/sqrt22)','start_square':'-20/11','edge_parameterization':'v(t)=(u-t*d_i,+/-20/sqrt22), t>=0; d_i is the old21 Gram dual; wall14 pairing=-1-t*a_i.','edges':edges,'finite_endpoint_count_per_cap':3,'open_edge_count_per_cap':18,'ideal_endpoint_count_per_cap':0,'explicit_open_edge':{'released_old_wall':1,'dual_old_coordinates':list(map(str,d)),'dual_square':'8','null_endpoint':'(u-sqrt(5/22)*d,+/-20/sqrt22)','null_endpoint_incident':{'old':[i for i in I if i!=1],'new':'the same-sign sqrt22 cap'},'null_endpoint_wall14_pair':'-1','proof_other_walls_strict':'Every other wall is strictly dotted to this cap, by the frozen gap certificate. Its projection to cap-perp is timelike and pairs negatively with the start vertex, so it is strictly negative at every future null point of the facet.','interpretation':'The incident21 Gram is positive definite, not an affine cusp diagram. This is an open ideal end of the edge, not a bounded cusp vertex.'},'necessary_next_step':'Any finite acute-compatible completion retaining the newcap needs at least one added wall with a finite-angle or parallel relation to it, or must replace it. Adding only dotted walls leaves all18 openedges unchanged.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
